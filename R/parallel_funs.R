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
  #' \donttest{
  #' set.seed(111)
  #' dat <- bkmr::SimData(n = 50, M = 4)
  #' y <- dat$y
  #' Z <- dat$Z
  #' X <- dat$X
  #' set.seed(111)
  #' Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
  #' future::plan(strategy = future::multiprocess, workers=2)
  #' # only 50 iterations fit to save installation time
  #' fitkm.list <- kmbayes_parallel(nchains=2, y = y, Z = Z, X = X, iter = 50,
  #'   verbose = FALSE, varsel = TRUE)
  #' closeAllConnections()
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

kmbayes_combine <- function(fitkm.list, burnin=0, reorder=TRUE) {
  #' Combine multiple BKMR chains
  #'
  #' @description Combine multiple chains comprising BKMR fits at different starting
  #' values.
  #'
  #' @details Chains are not combined fully sequentially
  #'
  #' @param fitkm.list output from \code{\link[bkmrhat]{kmbayes_parallel}}
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
  #' Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
  #' future::plan(strategy = future::multiprocess, workers=2)
  #' # run 4 parallel Markov chains (low iterations used for illustration)
  #' fitkm.list <- kmbayes_parallel(nchains=2, y = y, Z = Z, X = X, iter = 500,
  #'   verbose = FALSE, varsel = TRUE)
  #' bigkm = kmbayes_combine(fitkm.list)
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

#' @rdname kmbayes_combine
#' @export
comb_bkmrfits <- kmbayes_combine

