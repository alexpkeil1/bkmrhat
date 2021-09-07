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
  #' @importFrom stats runif
  # #' @importFrom future value
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
  #'
  #' future::plan(strategy = future::multisession, workers=2)
  #' # only 50 iterations fit to save installation time
  #' fitkm.list <- kmbayes_parallel(nchains=2, y = y, Z = Z, X = X, iter = 50,
  #'   verbose = FALSE, varsel = TRUE)
  #' closeAllConnections()
  #' }
  ff <- list()
  ss = round(runif(nchains) * .Machine$integer.max)
  for (ii in 1:nchains) {
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      bkmr::kmbayes(...)
    }, seed=ss[ii])
  }
  res <- future::value(ff)
  class(res) <- c("bkmrfit.list", class(res))
  res
}

kmbayes_combine <- function(fitkm.list, burnin=NULL, excludeburnin=FALSE, reorder=TRUE) {
  #' Combine multiple BKMR chains
  #'
  #' @description Combine multiple chains comprising BKMR fits at different starting
  #' values.
  #'
  #' @details Chains are not combined fully sequentially
  #'
  #' @param fitkm.list output from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param burnin (numeric, or default=NULL) add in custom burnin (number of burnin iterations per chain).
  #' If NULL, then default to half of the chain
  #' @param excludeburnin (logical, default=FALSE) should burnin iterations be excluded from the final chains?
  #' Note that all bkmr package functions automatically exclude burnin from calculations.
  #' @param reorder (logical, default=TRUE) ensures that the first half of the combined chain contains
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
  #' @name kmbayes_combine
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
  #'
  #' future::plan(strategy = future::multisession, workers=2)
  #' # run 4 parallel Markov chains (low iterations used for illustration)
  #' fitkm.list <- kmbayes_parallel(nchains=2, y = y, Z = Z, X = X, iter = 500,
  #'   verbose = FALSE, varsel = TRUE)
  #' # use bkmr defaults for burnin, but keep them
  #' bigkm = kmbayes_combine(fitkm.list, excludeburnin=FALSE)
  #' ests = ExtractEsts(bigkm) # defaults to keeping second half of samples
  #' ExtractPIPs(bigkm)
  #' pred.resp.univar <- PredictorResponseUnivar(fit = bigkm)
  #' risks.overall <- OverallRiskSummaries(fit = bigkm, y = y, Z = Z, X = X,
  #'   qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, method = "exact")
  #'
  #' # additional objects that are not in a standard bkmrfit object:
  #' summary(bigkm$iters) # note that this reflects how fits are re-ordered to reflect burnin
  #' table(bigkm$chain)
  #' }
  #'
  #' closeAllConnections()
  #'
  kmoverall <- fitkm.list[[1]]
  nchains <- length(fitkm.list)
  kmIter <- length(kmoverall$sigsq.eps)
  #c("h.hat", "beta", "lambda", "sigsq.eps", "r", "acc.r", "acc.lambda", "delta",
  # "acc.rdelta", "move.type", "est.h", "time1", "time2", "iter", "family",
  # "starting.values", "control.params", "X", "Z", "y", "ztest", "data.comps", "varsel")
  kmoverall$chain <- rep(1:nchains, each=kmIter)
  kmoverall$iters <- rep(1:kmIter, times=nchains)
  for(kmfitidx in 1:nchains){
    fitkm.list[[kmfitidx]]$iters <- 1:kmIter
  }
  burnend <- ifelse(is.null(burnin), ceiling(kmIter/2), burnin)
  autoburn <- which(kmoverall$iters <= burnend)
  if(excludeburnin) autoburn = integer(0L)
  autonotburn <- which(kmoverall$iters > burnend)
  getparm <- function(lst, parm) {
    lst[[parm]]
  }
  getparmmat <- function(lst, parm) {
    lst[[parm]]#[(burnin+1):kmIter, , drop=FALSE]
  }
  getparmvec <- function(lst, parm) {
    lst[[parm]]#[(burnin+1):kmIter]
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

#' @rdname kmbayes_combine
#' @export
comb_bkmrfits <- kmbayes_combine

kmbayes_combine_lowmem <- function(fitkm.list, burnin=NULL, excludeburnin=FALSE, reorder=TRUE) {
  #' Combine multiple BKMR chains in lower memory settings
  #'
  #' @description Combine multiple chains comprising BKMR fits at different starting
  #' values. This function writes some results
  #' to disk, rather than trying to process fully within memory which, in some cases,
  #' will result in avoiding "out of memory" errors that can happen with kmbayes_combine.
  #'
  #' @details Chains are not combined fully sequentially (see "reorder")
  #'
  #' @param fitkm.list output from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param burnin (numeric, or default=NULL) add in custom burnin (number of burnin iterations per chain).
  #' If NULL, then default to half of the chain
  #' @param excludeburnin (logical, default=FALSE) should burnin iterations be excluded from the final chains?
  #' Note that all bkmr package functions automatically exclude burnin from calculations.
  #' @param reorder (logical, default=TRUE) ensures that the first half of the combined chain contains
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
  #' @importFrom data.table fwrite fread as.data.table
  #' @export
  #' @name kmbayes_combine_lowmem
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
  #'
  #' future::plan(strategy = future::multisession, workers=2)
  #' # run 4 parallel Markov chains (low iterations used for illustration)
  #' fitkm.list <- kmbayes_parallel(nchains=2, y = y, Z = Z, X = X, iter = 500,
  #'   verbose = FALSE, varsel = TRUE)
  #' # use bkmr defaults for burnin, but keep them
  #' bigkm = kmbayes_combine_lowmem(fitkm.list, excludeburnin=FALSE)
  #' ests = ExtractEsts(bigkm) # defaults to keeping second half of samples
  #' ExtractPIPs(bigkm)
  #' pred.resp.univar <- PredictorResponseUnivar(fit = bigkm)
  #' risks.overall <- OverallRiskSummaries(fit = bigkm, y = y, Z = Z, X = X,
  #'   qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, method = "exact")
  #'
  #' # additional objects that are not in a standard bkmrfit object:
  #' summary(bigkm$iters) # note that this reflects how fits are re-ordered to reflect burnin
  #' table(bigkm$chain)
  #' }
  #'
  #' closeAllConnections()
  #'
  kmoverall <- fitkm.list[[1]]
  nchains <- length(fitkm.list)
  kmIter <- length(kmoverall$sigsq.eps)
  #c("h.hat", "beta", "lambda", "sigsq.eps", "r", "acc.r", "acc.lambda", "delta",
  # "acc.rdelta", "move.type", "est.h", "time1", "time2", "iter", "family",
  # "starting.values", "control.params", "X", "Z", "y", "ztest", "data.comps", "varsel")
  kmoverall$chain <- rep(1:nchains, each=kmIter)
  kmoverall$iters <- rep(1:kmIter, times=nchains)
  for(kmfitidx in 1:nchains){
    fitkm.list[[kmfitidx]]$iters <- 1:kmIter
  }
  burnend <- ifelse(is.null(burnin), ceiling(kmIter/2), burnin)
  autoburn <- which(kmoverall$iters <= burnend)
  if(excludeburnin) autoburn = integer(0L)
  autonotburn <- which(kmoverall$iters > burnend)
  getparm <- function(lst, parm) {
    lst[[parm]]
  }
  getparmmat <- function(lst, parm) {
    lst[[parm]]
  }
  getparmvec <- function(lst, parm) {
    lst[[parm]]
  }
  for (matparm in c("h.hat", "ystar")) {
    tfa = tempfile()
    tfb = tempfile()
    for(k in 1:nchains){
      tmp = getparmmat(fitkm.list[[k]], parm=matparm)
      if(!is.null(tmp)){
        tmp = data.table::as.data.table(tmp)
        if(k == 1){
          data.table::fwrite(tmp[autoburn, , drop=FALSE], file = tfa, row.names=FALSE, verbose=FALSE)
          data.table::fwrite(tmp[autonotburn, , drop=FALSE], file = tfb, row.names=FALSE, verbose=FALSE)
        }
        if(k > 1){
          data.table::fwrite(tmp[autoburn, , drop=FALSE], file = tfa, append=TRUE, verbose=FALSE)
          data.table::fwrite(tmp[autonotburn, , drop=FALSE], file = tfb, append=TRUE, verbose=FALSE)
        }
      }
    }
    if(!is.null(tmp)){
      kmoverall[[matparm]] <- rbind(
        as.matrix(fread(tfa, header=FALSE)),
        as.matrix(fread(tfb, header=FALSE))
      )
    }else{
      kmoverall[[matparm]] <- tmp
    }
    rm("tmp")
  }
  for (matparm in c("beta", "lambda", "r", "acc.r", "acc.lambda", "delta")) {
    tmp <- do.call("rbind", lapply(fitkm.list, FUN=getparmmat, parm=matparm))
    kmoverall[[matparm]] <- rbind(tmp[autoburn, , drop=FALSE],
                                  tmp[autonotburn, , drop=FALSE])
    rm("tmp")
  }
  for (vecparm in c("sigsq.eps", "acc.rdelta", "move.type", "iters")) {
    tmp <- do.call("c", lapply(fitkm.list, FUN=getparmvec, parm=vecparm))
    kmoverall[[vecparm]] <- c(tmp[autoburn], tmp[autonotburn])
    rm("tmp")
  }
  for (sumparm in c("iter")) {
    kmoverall[[sumparm]] <- do.call("sum", lapply(fitkm.list, FUN=getparm, parm=sumparm))
  }
  class(kmoverall) <- c("bkmrplusfit", class(kmoverall))
  kmoverall
}

#' @rdname kmbayes_combine_lowmem
#' @export
comb_bkmrfits_lowmem <- kmbayes_combine_lowmem
