.checkver <- function(print=FALSE){
  #' @importFrom utils packageVersion
  ver = packageVersion('bkmr')[[1]]
  if(print) message(paste('bkmr package version', ver))
  minver = list(c(0L,2L,9000L))
  class(minver) = c("package_version", "numeric_version")
  ver >= minver
}

.add_bkmrfits <- function(fitkm.list, trim=TRUE) {
  # combine two bkmr fits of possibly different lengths
  burnin=0
  reorder=FALSE
  kmoverall <- fitkm.list[[1]]
  nchains <- length(fitkm.list)
  iters = sapply(fitkm.list, function(x) length(x$sigsq.eps))
  kmIter <- sum(iters)
  #c("h.hat", "beta", "lambda", "sigsq.eps", "r", "acc.r", "acc.lambda", "delta",
  # "acc.rdelta", "move.type", "est.h", "time1", "time2", "iter", "family",
  # "starting.values", "control.params", "X", "Z", "y", "ztest", "data.comps", "varsel")
  getparm <- function(lst, parm) {
    lst[[parm]]
  }
  getparmmat <- function(lst, parm) {
    lst[[parm]]
  }
  getparmvec <- function(lst, parm) {
    lst[[parm]]
  }
  for (matparm in c("h.hat", "beta", "lambda", "r", "acc.r", "acc.lambda", "delta", "ystar")) {
    tmp <- do.call("rbind", lapply(fitkm.list, FUN=getparmmat, parm=matparm))
    # cut first iteration from second chain, which is just the starting values
    if(trim) kmoverall[[matparm]] <- tmp[-(iters[[1]] + 1),,drop=FALSE]
    if(!trim) kmoverall[[matparm]] <- tmp[,,drop=FALSE]
  }
  for (vecparm in c("sigsq.eps", "acc.rdelta", "move.type", "iters")) {
    tmp <- do.call("c", lapply(fitkm.list, FUN=getparmvec, parm=vecparm))
    if(trim) kmoverall[[vecparm]] <- tmp[-(iters[[1]] + 1)]
    if(!trim) kmoverall[[vecparm]] <- tmp
  }
  for (sumparm in c("iter")) {
    if(trim) kmoverall[[sumparm]] <- -1 + do.call("sum", lapply(fitkm.list, FUN=getparm, parm=sumparm))
    if(!trim) kmoverall[[sumparm]] <- do.call("sum", lapply(fitkm.list, FUN=getparm, parm=sumparm))
  }
  class(kmoverall) <- c("bkmrfit.continued", class(kmoverall))
  kmoverall
}

.ensuremat <- function(x){
  if(!is.null(dim(x)[1])) return(x)
  cbind(x,x)[,1,drop=FALSE]
}

#' Continue sampling from existing bkmr fit
#'
#' Use this when you've used MCMC sampling with the \code{\link[bkmr]{kmbayes}} function, but you did not take enough samples and do not want to start over.
#'
#' Note this does not fully start from the prior values of the MCMC chains. The \code{\link[bkmr]{kmbayes}} function does not allow full specification of the kernel function parameters, so this will restart the chain at the last values of all fixed effect parameters, and start the kernel `r` parmeters at the arithmetic mean of all `r` parameters from the last step in the previous chain.
#' @param fit output from \code{\link[bkmr]{kmbayes}}
#' @param ... arguments to \code{\link[bkmrhat]{kmbayes_continue}}
#'
#' @return a `bkmrfit.continued` object, which inherits from `bkmrfit` objects similar to \code{\link[bkmr]{kmbayes}} output, and which can be used to make inference using functions from the `bkmr` package.
#' @export
#' @seealso \code{\link[bkmrhat]{kmbayes_parallel}}
#' @import bkmr
#' @importFrom utils modifyList
#'
#' @examples
#' set.seed(111)
#' dat <- bkmr::SimData(n = 50, M = 4)
#' y <- dat$y
#' Z <- dat$Z
#' X <- dat$X
#' \dontrun{
#' fitty1 = bkmr::kmbayes(y=y,Z=Z,X=X, est.h=TRUE, tier=100)
#' # do some diagnostics here to see if 100 iterations (default) is enough
#' # add 100 additional iterations (for illustration - still will not be enough)
#' fitty2 = kmbayes_continue(fitty1, iter=100)
#' cobj = as.mcmc(fitty2)
#' varnames(cobj)
#'
#' }
#'
kmbayes_continue <- function(fit, ...){
  bkmrvernew = .checkver()
  eps = 1e-9
  last.iter = fit$iter
  ending.values = sapply(names(fit$starting.values), function(x) .ensuremat(fit[[x]])[last.iter,])
  if(!fit$est.h) ending.values$h.hat = ending.values$h.hat + 1
  if(!bkmrvernew){
    # using old version of bkmr
    if(fit$est.h){
      cm = pmax(eps, colMeans(fit$h.hat))
      gt0 = ending.values$h.hat>0
      if(sum(1-gt0)>0){
        message("Modifying h.hat starting values to meet kmbayes initial value constraints (this isn't a perfect continuation)")
      }
      ending.values$h.hat = ending.values$h.hat*gt0 + cm*(1-gt0)
    }
    message("Modifying r starting values to meet kmbayes initial value constraints (this isn't a perfect continuation)")
    message("This issue can be fixed by updating bkmr to the development version via: install.packages('devtools'); devtools::install_github('jenfb/bkmr')")
    if(sum(ending.values$delta)>0) ending.values$r = mean(ending.values$r[which(ending.values$delta==1)])
    if(sum(ending.values$delta)==0) ending.values$r = eps
  }
  newstart = list(starting.values = ending.values)
  unfit = unclass(fit)
  keepargs = names(as.list(args(kmbayes)))
  oldargs = unfit[names(unfit)[which(names(unfit) %in% keepargs)]]
  newargs = list(...)
  newargs = modifyList(oldargs, newargs)
  newargs = modifyList(newargs, newstart)
  newargs$iter = newargs$iter + 1 # account for restart
  fit2 = do.call(kmbayes, newargs)
  res = list(fit, fit2)
  .add_bkmrfits(res, trim=TRUE)
}

#' Continue sampling from existing bkmr_parallel fit
#'
#' Use this when you've used MCMC sampling with the \code{\link[bkmrhat]{kmbayes_parallel}} function, but you did not take enough samples and do not want to start over.
#' @param fitkm.list output from \code{\link[bkmrhat]{kmbayes_parallel}}
#' @param ... arguments to \code{\link[bkmrhat]{kmbayes_continue}}
#'
#' @return a `bkmrfit.list` object, which is just a list of `bkmrfit` objects similar to \code{\link[bkmrhat]{kmbayes_parallel}}
#' @export
#' @seealso \code{\link[bkmrhat]{kmbayes_parallel}}
#' @import future
#'
#' @examples
#' set.seed(111)
#' dat <- bkmr::SimData(n = 50, M = 4)
#' y <- dat$y
#' Z <- dat$Z
#' X <- dat$X
#' \dontrun{
#' Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
#' future::plan(strategy = future::multiprocess, workers=2)
#' fitty1p = kmbayes_parallel(nchains=2, y=y,Z=Z,X=X)
#'
#' fitty2p = kmbayes_parallel_continue(fitty1p, iter=3000)
#' cobj = as.mcmc.list(fitty2p)
#' plot(cobj)
#' }
kmbayes_parallel_continue <- function(fitkm.list, ...) {
  ff <- list()
  nchains = length(fitkm.list)
  for (ii in 1:nchains) {
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      kmbayes_continue(fitkm.list[[ii]], ...)
    })
  }
  res <- values(ff)
  class(res) <- c("bkmrfit.list", class(res))
  res
}
