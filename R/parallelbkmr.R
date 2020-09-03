#### Parallel processing ####


kmbayes_parallel <- function(nchains=4, ...) {
  #' Run multiple BKMR chains in parallel
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
  #' \dontrun{
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
  res <- values(ff)
  class(res) <- c("bkmrfit.list", class(res))
  res
}

comb_bkmrfits <- function(fitkm.list, burnin=0, reorder=TRUE) {
  #' Combine multiple BKMR chains
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
  #' \dontrun{
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
  kmoverall <- fitkm.list[[1]]
  nchains <- length(fitkm.list)
  kmIter <- length(kmoverall$sigsq.eps)
  #c("h.hat", "beta", "lambda", "sigsq.eps", "r", "acc.r", "acc.lambda", "delta",
  # "acc.rdelta", "move.type", "est.h", "time1", "time2", "iter", "family",
  # "starting.values", "control.params", "X", "Z", "y", "ztest", "data.comps", "varsel")
  kmoverall$chain <- rep(1:nchains, each=kmIter)
  kmoverall$iters <- rep(1:kmIter, times=nchains)
  autoburn <- which(kmoverall$iters <= ceiling(kmIter/2))
  autonotburn <- which(kmoverall$iters > ceiling(kmIter/2))
  getparm <- function(lst, parm) {
    lst[[parm]]
  }
  getparmmat <- function(lst, parm) {
    lst[[parm]][(burnin+1):kmIter, , drop=FALSE]
  }
  getparmvec <- function(lst, parm) {
    lst[[parm]][(burnin+1):kmIter]
  }
  for (matparm in c("h.hat", "beta", "lambda", "r", "acc.r", "acc.lambda", "delta", "ystar")) {
    tmp <- do.call("rbind", lapply(fitkm.list, FUN=getparmmat, parm=matparm))
    kmoverall[[matparm]] <- rbind(tmp[autoburn, , drop=FALSE],
                                  tmp[autonotburn, , drop=FALSE])
  }
  for (vecparm in c("sigsq.eps", "acc.rdelta", "move.type", "iters")) {
    tmp <- do.call("c", lapply(fitkm.list, FUN=getparmvec, parm=vecparm))
    kmoverall[[vecparm]] <- c(tmp[autoburn], tmp[autonotburn])
  }
  for (sumparm in c("iter")) {
    kmoverall[[sumparm]] <- do.call("sum", lapply(fitkm.list, FUN=getparm, parm=sumparm))
  }
  class(kmoverall) <- c("bkmrplusfit", class(kmoverall))
  kmoverall
}



OverallRiskSummaries_parallel <- function(x, ...){
  #' Overall summary by chain
  #' @param x bkmrfit.list object from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param ... arguments to \code{\link[bkmr]{OverallRiskSummaries}}
  #'
  #' @return data.frame with all chains together
  #' @importFrom bkmr OverallRiskSummaries
  #' @export
  #'
  ff <- list()
  nchains = length(x)
  for (ii in 1:nchains) {
    xii = x[[ii]]
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      df = suppressWarnings(bkmr::OverallRiskSummaries(xii, ...))
      df$chain=ii
      df
    })
  }
  res <- values(ff)
  as.data.frame(do.call("rbind", res))
}

PredictorResponseUnivar_parallel <- function(x, ...){
  #' Univariate predictor response summary by chain
  #' @param x bkmrfit.list object from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param ... arguments to \code{\link[bkmr]{PredictorResponseUnivar}}
  #'
  #' @return data.frame with all chains together
  #' @importFrom bkmr PredictorResponseUnivar
  #' @export
  #'
  ff <- list()
  nchains = length(x)
  for (ii in 1:nchains) {
    xii = x[[ii]]
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      df = suppressWarnings(bkmr::PredictorResponseUnivar(xii, ...))
      df$chain=ii
      df
    })
  }
  res <- values(ff)
  as.data.frame(do.call("rbind", res))
}


PredictorResponseBivar_parallel <- function(x, ...){
  #' Bivariate predictor response by chain
  #' @param x bkmrfit.list object from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param ... arguments to \code{\link[bkmr]{PredictorResponseBivar}}
  #'
  #' @return data.frame with all chains together
  #' @importFrom bkmr PredictorResponseBivar
  #' @export
  #'
  ff <- list()
  nchains = length(x)
  for (ii in 1:nchains) {
    xii = x[[ii]]
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      df = suppressWarnings(bkmr::PredictorResponseBivar(xii, ...))
      df$chain=ii
      df
    })
  }
  res <- values(ff)
  as.data.frame(do.call("rbind", res))
}


SingVarRiskSummaries_parallel <- function(x, ...){
  #' Single variable summary by chain
  #' @param x bkmrfit.list object from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param ... arguments to \code{\link[bkmr]{SingVarRiskSummaries}}
  #'
  #' @return data.frame with all chains together
  #' @importFrom bkmr SingVarRiskSummaries
  #' @export
  #'
  ff <- list()
  nchains = length(x)
  for (ii in 1:nchains) {
    xii = x[[ii]]
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      df = suppressWarnings(bkmr::SingVarRiskSummaries(xii, ...))
      df$chain=ii
      df
    })
  }
  res <- values(ff)
  as.data.frame(do.call("rbind", res))
}


#ExtractSamps_parallel <- function(x, ...){
#  #' Extract posterior samples by chain
#  #' @param x bkmrfit.list object from \code{\link[bkmrhat]{kmbayes_parallel}}
#  #' @param ... arguments to \code{\link[bkmr]{ExtractSamps}}
#  #'
#  #' @return data.frame with all chains together
#  #' @importFrom bkmr ExtractSamps
#  #' @export
#  #'
#  ff <- list()
#  nchains = length(x)
#  for (ii in 1:nchains) {
#    xii = x[[ii]]
#    ff[[ii]] <- future({
#      cat(paste("Chain", ii, "\n"))
#      df = suppressWarnings(bkmr::ExtractSamps(xii, ...))
#      df$chain=ii
#      df
#    })
#  }
#  res <- values(ff)
#  as.data.frame(do.call("rbind", res))
#}

ExtractPIPs_parallel <- function(x, ...){
  #' Posterior inclusion probabilities by chain
  #' @param x bkmrfit.list object from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param ... arguments to \code{\link[bkmr]{CalcPIPs}}
  #'
  #' @return data.frame with all chains together
  #' @importFrom bkmr ExtractPIPs
  #' @export
  #'
  ff <- list()
  nchains = length(x)
  for (ii in 1:nchains) {
    xii = x[[ii]]
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      df = suppressWarnings(data.frame(bkmr::ExtractPIPs(xii)))
      df$chain=ii
      df
    })
  }
  res <- values(ff)
  as.data.frame(do.call("rbind", res))
}


SamplePred_parallel <- function(x, ...){
  #' Posterior inclusion probabilities by chain
  #' @param x bkmrfit.list object from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param ... arguments to \code{\link[bkmr]{CalcPIPs}}
  #'
  #' @return data.frame with all chains together
  #' @importFrom bkmr SamplePred
  #' @export
  #'
  ff <- list()
  nchains = length(x)
  for (ii in 1:nchains) {
    xii = x[[ii]]
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      df = suppressWarnings(as.data.frame(bkmr::SamplePred(xii, ...)))
      df$chain=ii
      df
    })
  }
  res <- values(ff)
  as.data.frame(do.call("rbind", res))
}

# may also want parallel implementations of:
# CalcWithinGroupPIPs


# this should go into another file
predict.bkmrfit <- function(object, ptype=c("mean", "sd.fit"), ...) {
  #' Posterior mean/sd predictions
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
  #' \dontrun{
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
  sf <- bkmr::SamplePred(object, ...)
  if (ptype[1] == "mean") {
    postmean <- as.numeric(drop(apply(sf, 2, mean)))
    return(postmean)
  }
  if (ptype[1] == "sd.fit") {
    postsd <- as.numeric(drop(apply(sf, 2, sd)))
    return(postsd)
  }
}
