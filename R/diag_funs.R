#### MCMC diagnostics ####


#' MCMC diagnostics using rstan
#'
#' @description Give MCMC diagnostistics from the \code{rstan} package
#' using the \code{\link[rstan]{Rhat}}, \code{\link[rstan]{ess_bulk}},
#' and \code{\link[rstan]{ess_tail}} functions. Note that r-hat is only
#' reported for \code{bkmrfit.list} objects from \code{\link[bkmrhat]{kmbayes_parallel}}
#'
#' @param kmobj Either an object from \code{\link[bkmr]{kmbayes}} or
#' from \code{\link[bkmrhat]{kmbayes_parallel}}
#' @param ... arguments to \code{\link[rstan]{monitor}}
#'
#' @export
#'
#' @name kmbayes_diagnose
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
#' fitkm.list <- kmbayes_parallel(nchains=2, y = y, Z = Z, X = X, iter = 1000,
#'   verbose = FALSE, varsel = TRUE)
#' kmbayes_diag(fitkm.list)
#' kmbayes_diag(fitkm.list[[1]]) # just the first chain
#'
#' closeAllConnections()
#' }
kmbayes_diagnose <- function(kmobj, ...) {
  if (inherits(kmobj, "bkmrfit.list")) {
    message("Parallel chains\n")
    res <- .diag_par(kmobj, ...)
  }
  if (inherits(kmobj, "bkmrfit")) {
    message("Single chain\n")
    res <- .diag(kmobj, ...)
  }
  res
}


#' @rdname kmbayes_diagnose
#' @export
kmbayes_diag <- kmbayes_diagnose
