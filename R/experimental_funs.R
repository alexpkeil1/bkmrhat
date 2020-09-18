.add_bkmrfits <- function(fitkm.list) {
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
    lst[[parm]][,,drop=FALSE]
  }
  getparmvec <- function(lst, parm) {
    lst[[parm]]
  }
  for (matparm in c("h.hat", "beta", "lambda", "r", "acc.r", "acc.lambda", "delta", "ystar")) {
    tmp <- do.call("rbind", lapply(fitkm.list, FUN=getparmmat, parm=matparm))
    kmoverall[[matparm]] <- tmp
  }
  for (vecparm in c("sigsq.eps", "acc.rdelta", "move.type", "iters")) {
    tmp <- do.call("c", lapply(fitkm.list, FUN=getparmvec, parm=vecparm))
    kmoverall[[vecparm]] <- tmp
  }
  for (sumparm in c("iter")) {
    kmoverall[[sumparm]] <- do.call("sum", lapply(fitkm.list, FUN=getparm, parm=sumparm))
  }
  class(kmoverall) <- c("bkmrcontfit", class(kmoverall))
  kmoverall
}

# experimental functions

.ensuremat <- function(x){
  if(!is.null(dim(x)[1])) return(x)
  cbind(x,x)[,1,drop=FALSE]
}





#' Continue sampling from existing bkmr fit
#'
#' Use this when you've used MCMC sampling with the \code{\link[bkmr]{kmbayes}} function, but you did not take enough samples and do not want to start over.
#'
#' Note this does not fully start from the prior values of the MCMC chains. The \code{\link[bkmr]{kmbayes}} does not allow full specification of the kernel function parameters, so
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
#' fitty1 = bkmr::kmbayes(y=y,Z=Z,X=X)
#'
#' fitty2 = kmbayes_continue(fitty1, iter=3000)
#' cobj = as.mcmc(fitty2)
#' plot(cobj)
#' }
#'
kmbayes_continue <- function(fit, ...){
  eps = 1e-9
  message("kmbayes cannot be continued perfectly at the last point\n
          but this function will get close")
  last.iter = fit$iter
  ending.values = sapply(names(fit$starting.values), function(x) .ensuremat(fit[[x]])[last.iter,])
  if(!fit$est.h) ending.values$h.hat = ending.values$h.hat + 1
  ending.values$r = mean(ending.values$r + eps)
  newstart = list(starting.values = ending.values)
  unfit = unclass(fit)
  keepargs = names(as.list(args(kmbayes)))
  oldargs = unfit[names(unfit)[which(names(unfit) %in% keepargs)]]
  newargs = list(...)
  newargs = modifyList(oldargs, newargs)
  newargs = modifyList(newargs, newstart)
  fit2 = do.call(kmbayes, newargs)
  res = list(fit, fit2)
  class(res) <- c("bkmrfit.continued", class(res))
  res
  .add_bkmrfits(res)
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
#' fitty1p = kmbayes_parallel(y=y,Z=Z,X=X)
#'
#' fitty2p = kmbayes_parallel_continue(fitty1p, iter=3000)
#' cobj = as.mcmc(fitty2)
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
