#### Interfaces with coda package ####

as.mcmc.bkmrfit <- function(x, iterstart=1, thin=1, ...) {
  #' Convert bkmrfit to mcmc object for coda MCMC diagnostics
  #'
  #' @description Converts a \code{kmrfit} (from the bkmr package) into
  #' an \code{\link[coda]{mcmc}} object from the \code{coda} package. The
  #' \code{coda} package enables many different types of single chain MCMC
  #' diagnostics, including \code{\link[coda]{geweke.diag}}, \code{\link[coda]{traceplot}} and
  #' \code{\link[coda]{effectiveSize}}. Posterior summarization is also available,
  #' such as \code{\link[coda]{HPDinterval}} and \code{\link[coda]{summary.mcmc}}.
  #'
  #' @param x object of type kmrfit (from bkmr package)
  #' @param iterstart first iteration to use (e.g. for implementing burnin)
  #' @param thin keep 1/thin % of the total iterations (at regular intervals)
  #' @param ... unused
  #'
  #' @return An \code{\link[coda]{mcmc}} object
  #' @importFrom coda mcmc as.mcmc
  #' @export
  #'
  #' @examples
  #'
  #' # following example from https://jenfb.github.io/bkmr/overview.html
  #'  \donttest{
  #' set.seed(111)
  #' library(coda)
  #' library(bkmr)
  #' dat <- bkmr::SimData(n = 50, M = 4)
  #' y <- dat$y
  #' Z <- dat$Z
  #' X <- dat$X
  #' set.seed(111)
  #' fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 500, verbose = FALSE,
  #'   varsel = FALSE)
  #' mcmcobj <- as.mcmc(fitkm, iterstart=251)
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
  requireNamespace("coda")
  df <- .extractparms(x, allvars = TRUE)
  nr <- nrow(df)
  mcmc(df[seq(iterstart, nr, by=thin), ], start=iterstart, thin=thin)
}


as.mcmc.list.bkmrfit.list <- function(x, ...) {
  #' Convert multi-chain bkmrfit to mcmc.list for coda MCMC diagnostics
  #'
  #' @description Converts a \code{kmrfit.list} (from the bkmrhat package) into
  #' an \code{\link[coda]{mcmc.list}} object from the \code{coda} package.The
  #' \code{coda} package enables many different types of MCMC diagnostics,
  #' including \code{\link[coda]{geweke.diag}}, \code{\link[coda]{traceplot}} and
  #' \code{\link[coda]{effectiveSize}}. Posterior summarization is also available,
  #' such as \code{\link[coda]{HPDinterval}} and \code{\link[coda]{summary.mcmc}}.
  #' Using multiple chains is necessary for certain MCMC diagnostics, such as
  #' \code{\link[coda]{gelman.diag}}, and \code{\link[coda]{gelman.plot}}.
  #'
  #' @param x object of type kmrfit.list (from bkmrhat package)
  #' @param ... arguments to \code{\link[bkmrhat]{as.mcmc.bkmrfit}}
  #'
  #' @return An \code{\link[coda]{mcmc.list}} object
  #' @importFrom coda mcmc.list as.mcmc.list
  #' @export
  #'
  #' @examples
  #' # following example from https://jenfb.github.io/bkmr/overview.html
  #'  \donttest{
  #' set.seed(111)
  #' library(coda)
  #' dat <- bkmr::SimData(n = 50, M = 4)
  #' y <- dat$y
  #' Z <- dat$Z
  #' X <- dat$X
  #' set.seed(111)
  #' Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
  #' future::plan(strategy = future::multiprocess, workers=2)
  #' # run 2 parallel Markov chains (more usually better)
  #' fitkm.list <- kmbayes_parallel(nchains=2, y = y, Z = Z, X = X, iter = 1000,
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
  #'
  #' closeAllConnections()
  #' }
  requireNamespace("coda")
  res <- lapply(x, as.mcmc.bkmrfit, ...)
  mcmc.list(res)
}

