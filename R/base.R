#### Hidden functions ####

.extractparms <- function(kmobj, allvars=FALSE){
  outlist = list()
  matlist = c("beta", "lambda", "r")
  veclist = c("sigsq.eps")
  if (kmobj$est.h) matlist = c("h.hat", matlist)
  if (kmobj$varsel & allvars) matlist = c(matlist, "delta") # rhat useful?
  if(allvars) veclist = c(veclist)

  for (matparm in matlist){
    width = ncol(kmobj[[matparm]])
    if(!is.null(width)){
      outlist[[matparm]] = as.data.frame(as.matrix(kmobj[[matparm]]))
      names(outlist[[matparm]]) = paste0(matparm, "_", 1:width)
    }
  }
  for (vecparm in veclist){
    outlist[[vecparm]] = as.data.frame(as.matrix(kmobj[[vecparm]]))
    names(outlist[[vecparm]]) = paste0(vecparm)
  }
  outdf = as.data.frame(do.call("cbind", outlist))
  outdf
}

.diag_par <- function(kmobj.list){
  #' @importFrom rstan Rhat ess_bulk ess_tail
  #'
  getparmvec <- function(lst, parm){
    lst[[parm]]
  }
  odf = lapply(kmobj.list, .extractparms)
  nms = names(odf[[1]])
  esstabs_b = sapply(odf, function(x) sapply(x, rstan::ess_bulk))
  esstabs_t = sapply(odf, function(x) sapply(x, rstan::ess_tail))
  rhatvec = list()
  for(nm in nms){
    nmmat = do.call("cbind", lapply(odf, function(x) getparmvec(x, nm)))
    rhatvec[[nm]] = rstan::Rhat(nmmat)
  }
  data.frame(ess_bulk = apply(esstabs_b, 1, sum), ess_tail = apply(esstabs_t, 1, sum), rhat = as.numeric(rhatvec))
}

.diag <- function(kmobj){
  #' @importFrom rstan ess_bulk ess_tail
  odf = .extractparms(kmobj)
  ess_b = sapply(odf, rstan::ess_bulk)
  ess_t = sapply(odf, rstan::ess_tail)
  data.frame(ess_bulk=ess_b, ess_tail=ess_t)
}

.predictivemean <- function(object, ptype=c("mean", "sd.fit"), ...){
  #' @importFrom bkmr SamplePred
  #' @importFrom stats sd
  sf = bkmr::SamplePred(object, ...)
  if (ptype[1] == "mean"){
    meansample = as.numeric(drop(apply(sf, 1, mean)))
    return(meansample)
  }
  if (ptype[1] == "sd.fit"){
    sdsample = as.numeric(drop(apply(sf, 1, sd)))
    return(sdsample)
  }
}


#### User visible functions ####

kmbayes_diag <- function(kmobj,...) {
  #' kmbayes_diag: rstan diagnostics
  #'
  #' @description Give MCMC diagnostistics from the \code{rstan} package
  #' using the \code{\link[rstan]{Rhat}}, \code{\link[rstan]{ess_bulk}},
  #' and \code{\link[rstan]{ess_tail}} functions. Note that r-hat is only
  #' reported for \code{bkmrfit.list} objects from \code{\link[bkmrhat]{kmbayes_parallel}}
  #'
  #' @param kmobj Either an object from \code{\link[bkmr]{kmbayes}} or
  #' from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param ... unused
  #'
  #' @export
  #'
  #' @examples
  #' \donttest{
  #' set.seed(111)
  #' dat <- bkmr::SimData(n = 50, M = 4)
  #' y <- dat$y
  #' Z <- dat$Z
  #' X <- dat$X
  #' set.seed(111)
  #' Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
  #' future::plan(strategy = future::multiprocess)
  #' fitkm.list <- kmbayes_parallel(nchains=4, y = y, Z = Z, X = X, iter = 5000,
  #'   verbose = FALSE, varsel = TRUE)
  #' kmbayes_diag(fitkm.list)
  #' kmbayes_diag(fitkm.list[[1]]) # just the first chain
  #' }
  #'
  if (inherits(kmobj, "bkmrfit.list")){
    message("Parallel chains\n")
    res = .diag_par(kmobj)
  }
  if (inherits(kmobj, "bkmrfit")){
    message("Single chain\n")
    res = .diag(kmobj)
  }
  res
}


kmbayes_parallel <- function(nchains=4, ...) {
  #' kmbayes_parallel: parallel bkmr
  #'
  #' @description Fit parallel chains from the \code{\link[bkmr]{kmbayes}} function.
  #' These chains leverage parallel processing from the \code{future} package, which
  #' can speed fitting and enable diagnostics that rely on multiple Markov
  #' chains from dispersed initial values.
  #'
  #' @param nchains number of parallel chains
  #' @param ... arguments to kmbayes
  #'
  #' @return a "bkmrfit.list" object, which is just an R list object in which each entry is a "bkmrfit" object \code{\link[bkmr]{kmbayes}}
  #' @importFrom rstan Rhat ess_bulk ess_tail
  #' @import future bkmr
  #' @export
  #'
  #' @examples
  #' \donttest{
  #' set.seed(111)
  #' dat <- bkmr::SimData(n = 50, M = 4)
  #' y <- dat$y
  #' Z <- dat$Z
  #' X <- dat$X
  #' set.seed(111)
  #' Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
  #' future::plan(strategy = future::multiprocess)
  #' fitkm.list <- kmbayes_parallel(nchains=4, y = y, Z = Z, X = X, iter = 5000,
  #'   verbose = FALSE, varsel = TRUE)
  #' }
  ff <- list()
  for (ii in 1:nchains) {
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      bkmr::kmbayes(...)
    })
  }
  #res <- future_lapply(ff, FUN = value)
  res <- values(ff)
  class(res) = c("bkmrfit.list", class(res))
  res
}

comb_bkmrfits <- function(fitkm.list, burnin=0, reorder=TRUE) {
  #' comb_bkmrfits: combine fits across multiple chains
  #'
  #' @description Combine multiple chains comprising BKMR fits at different starting
  #' values.
  #'
  #' @details Chains are not combined fully sequentially
  #'
  #' @param fitkm.list output from kmbayes_par
  #' @param burnin add in custom burnin (number of burnin iterations per chain)
  #' @param reorder ensures that the first half of the combined chain contains
  #'  only the first half of each individual chain - this allows unaltered use
  #'  of standard functions from bkmr package, which automatically trims the first
  #'  half of the iterations. This can be used for posterior summaries, but
  #'  certain diagnostics may not work well (autocorrelation, effective sample size)
  #'  so the diagnostics should be done on the individual chains
  #'  #' @param ... arguments to \code{\link[bkmrhat]{as.mcmc.bkmrfit}}

  #' @return a \code{bkmrplusfit} object, which inherits from \code{bkmrfit}
  #' (from the \code{\link[bkmr]{kmbayes}} function)
  #' with multiple chains combined into a single object and additional parameters
  #' given by \code{chain} and \code{iters}, which index the specific chains and
  #' iterations for each posterior sample in the \code{bkmrplusfit} object
  #' @export
  #'
  #' @examples
  #' \donttest{
  #' # following example from https://jenfb.github.io/bkmr/overview.html
  #' set.seed(111)
  #' library(bkmr)
  #' dat <- bkmr::SimData(n = 50, M = 4)
  #' y <- dat$y
  #' Z <- dat$Z
  #' X <- dat$X
  #' set.seed(111)
  #' Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
  #' future::plan(strategy = future::multiprocess)
  #' # run 4 parallel Markov chains
  #' fitkm.list <- kmbayes_parallel(nchains=4, y = y, Z = Z, X = X, iter = 5000,
  #'   verbose = FALSE, varsel = TRUE)
  #' bigkm = comb_bkmrfits(fitkm.list)
  #' ests = ExtractEsts(bigkm)
  #' ExtractPIPs(bigkm)
  #' pred.resp.univar <- PredictorResponseUnivar(fit = bigkm)
  #' risks.overall <- OverallRiskSummaries(fit = bigkm, y = y, Z = Z, X = X,
  #'   qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, method = "exact")
  #'
  #' # additional objects that are not in a standard bkmrfit object:
  #' summary(bigkm$iters)
  #' table(bigkm$chain)
  #' }
  #'
  #'
  #'
  kmoverall = fitkm.list[[1]]
  nchains = length(fitkm.list)
  fnm = names(kmoverall)
  kmN = length(kmoverall$y)
  kmIter = length(kmoverall$sigsq.eps)
  #c("h.hat", "beta", "lambda", "sigsq.eps", "r", "acc.r", "acc.lambda", "delta",
  # "acc.rdelta", "move.type", "est.h", "time1", "time2", "iter", "family",
  # "starting.values", "control.params", "X", "Z", "y", "ztest", "data.comps", "varsel")
  # names(kmoverall)
  chainspecific = c("starting.values")
  kmoverall$chain = rep(1:nchains, each=kmIter)
  kmoverall$iters = rep(1:kmIter, times=nchains)
  autoburn = which(kmoverall$iters <= ceiling(kmIter/2))
  autonotburn = which(kmoverall$iters > ceiling(kmIter/2))
  getparm <- function(lst, parm){
    lst[[parm]]
  }
  getparmmat <- function(lst, parm){
    lst[[parm]][(burnin+1):kmIter,,drop=FALSE]
  }
  getparmvec <- function(lst, parm){
    lst[[parm]][(burnin+1):kmIter]
  }
  for (matparm in c("h.hat", "beta", "lambda", "r", "acc.r", "acc.lambda", "delta")){
    tmp = do.call("rbind", lapply(fitkm.list, FUN=getparmmat, parm=matparm))
    kmoverall[[matparm]] = rbind(tmp[autoburn,,drop=FALSE], tmp[autonotburn,,drop=FALSE])
  }
  for (vecparm in c("sigsq.eps", "acc.rdelta", "move.type", "iters")){
    temp = do.call("c", lapply(fitkm.list, FUN=getparmvec, parm=vecparm))
    kmoverall[[vecparm]] = c(tmp[autoburn], tmp[autonotburn])
  }
  for (sumparm in c("iter")){
    kmoverall[[sumparm]] = do.call("sum", lapply(fitkm.list, FUN=getparm, parm=sumparm))
  }
  class(kmoverall) = c("bkmrplusfit", class(kmoverall))
  kmoverall
}


as.mcmc.bkmrfit <- function(kmobj, iterstart=1, thin=1){
  #' as.mcmc.bkmrfit: convert bkmrfit to mcmc object
  #'
  #' @description Converts a \code{kmrfit} (from the bkmr package) into
  #' an \code{\link[coda]{mcmc}} object from the \code{coda} package.
  #'
  #' @param kmobj object of type kmrfit (from bkmr package)
  #' @param iterstart first iteration to use (e.g. for implementing burnin)
  #' @param thin keep 1/thin % of the total iterations (at regular intervals)
  #'
  #' @return An \code{\link[coda]{mcmc}} object
  #' @importFrom coda mcmc as.mcmc
  #' @export
  #'
  #' @examples
  #'  #' # following example from https://jenfb.github.io/bkmr/overview.html
  #'  \donttest{
  #' set.seed(111)
  #' library(coda)
  #' library(bkmr)
  #' dat <- bkmr::SimData(n = 50, M = 4)
  #' y <- dat$y
  #' Z <- dat$Z
  #' X <- dat$X
  #' set.seed(111)
  #' fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 5000, verbose = FALSE, varsel = FALSE)
  #' mcmcobj = as.mcmc(fitkm, iterstart=2501)
  #' summary(mcmcobj) # posterior summaries of model parameters
  #' # compare with default from bkmr package, which omits first 1/2 of chain
  #' summary(fitkm)
  #' # note this only works on multiple chains (see kmbayes_parallel)
  #' # gelman.diag(mcmcobj)
  #' # lots of functions in the coda package to use
  #' traceplot(mcmcobj)
  #' # will also fail with delta functions (when using variable selection)
  #' try(geweke.plot(mcmcobj))
  #' }
  df = .extractparms(kmobj, allvars = TRUE)
  nr = nrow(df)
  mcmc(df[seq(iterstart, nr, by=thin),], start=iterstart, thin=thin)
}


as.mcmc.list.bkmrfit.list <- function(kmobj, ...){
  #' as.mcmc.list.bkmrfit.list: convert multi-chain bkmrfit to mcmc.list object
  #'
  #' @description Converts a \code{kmrfit.list} (from the bkmrhat package) into
  #' an \code{\link[coda]{mcmc.list}} object from the \code{coda} package.
  #'
  #' @param kmobj object of type kmrfit.list (from bkmrhat package)
  #' @param ... arguments to \code{\link[bkmrhat]{as.mcmc.bkmrfit}}
  #'
  #' @return An \code{\link[coda]{mcmc.list}} object
  #' @importFrom coda mcmc.list as.mcmc.list
  #' @export
  #'
  #' @examples
  #'  #' # following example from https://jenfb.github.io/bkmr/overview.html
  #'  \donttest{
  #' set.seed(111)
  #' library(coda)
  #' dat <- bkmr::SimData(n = 50, M = 4)
  #' y <- dat$y
  #' Z <- dat$Z
  #' X <- dat$X
  #' set.seed(111)
  #' Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
  #' future::plan(strategy = future::multiprocess)
  #' # run 4 parallel Markov chains
  #' fitkm.list <- kmbayes_parallel(nchains=4, y = y, Z = Z, X = X, iter = 5000,
  #'   verbose = FALSE, varsel = FALSE)
  #' mcmcobj = as.mcmc.list(fitkm.list)
  #' summary(mcmcobj)
  #' # Gelman/Rubin diagnostics won't work on certain objects,
  #' # like delta parameters (when using variable selection),
  #' # so the rstan version of this will work better (does not give errors)
  #'  try(gelman.diag(mcmcobj))
  #' # lots of functions in the coda package to use
  #' plot(mcmcobj)
  #' # both of these will also fail with delta functions (when using variable selection)
  #' try(gelman.plot(mcmcobj))
  #' try(geweke.plot(mcmcobj))
  #' }
  res = lapply(kmobj, as.mcmc.bkmrfit, ...)
  mcmc.list(res)
}


predict.bkmrfit <- function(object, ptype=c("mean", "sd.fit"), ...){
  #' predict.bkmrfit: posterior mean/sd predictions
  #'
  #' @description Provides observation level predictions based on the posterior mean, or, alternatively, yields the posterior standard deviations of predictions for an observation. This function is useful for interfacing with ensemble machine learning packages such as `SuperLearner`, which utilize only point estimates.
  #'
  #' @param object fitted object of class inheriting from "bkmrfit".
  #' @param ptype "mean" or "sd.fit", where "mean" yields posterior
  #'  mean prediction for every observation in the data, and "sd.fit"
  #'  yields the posterior standard deviation for every observation in
  #'  the data.
  #' @param ... arguments to bkmr::SamplePred
  #' @importFrom stats predict
  #' @importFrom bkmr SamplePred
  #' @importFrom stats sd
  #' @return vector of predictions the same length as the outcome in the bkmrfit object
  #' @export
  #'
  #' @examples
  #' # following example from https://jenfb.github.io/bkmr/overview.html
  #' \donttest{
  #' library(bkmr)
  #' set.seed(111)
  #' dat <- bkmr::SimData(n = 50, M = 4)
  #' y <- dat$y
  #' Z <- dat$Z
  #' X <- dat$X
  #' set.seed(111)
  #' fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 5000, verbose = FALSE,
  #'   varsel = TRUE)
  #' postmean = predict(fitkm)
  #' postmean2 = predict(fitkm, Znew=Z/2)
  #' # mean difference in posterior means
  #' mean(postmean-postmean2)
  #' }
  sf = bkmr::SamplePred(object, ...)
  if (ptype[1] == "mean"){
    postmean = as.numeric(drop(apply(sf, 2, mean)))
    return(postmean)
  }
  if (ptype[1] == "sd.fit"){
    postsd = as.numeric(drop(apply(sf, 2, sd)))
    return(postsd)
  }
}

