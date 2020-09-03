#### Hidden functions ####

.extractparms <- function(kmobj, allvars=FALSE) {
  outlist <- list()
  matlist <- c("lambda", "r")
  veclist <- NULL
  if (!all(kmobj$X == 0)) matlist <- c("beta", matlist)
  if (kmobj$family!="binomial") veclist <- c("sigsq.eps")
  if (kmobj$family=="binomial") matlist <- c("ystar", matlist)
  if (kmobj$est.h) matlist <- c("h.hat", matlist)
  if (kmobj$varsel & allvars) matlist <- c(matlist, "delta") # rhat useful?

  for (matparm in matlist) {
    width <- ncol(kmobj[[matparm]])
    if (!is.null(width)) {
      outlist[[matparm]] <- as.data.frame(as.matrix(kmobj[[matparm]]))
      if(width>1){
        names(outlist[[matparm]]) <- paste0(matparm, 1:width)
      } else{
        names(outlist[[matparm]]) <- paste0(matparm)
      }
    }
  }
  for (vecparm in veclist) {
    outlist[[vecparm]] <- as.data.frame(as.matrix(kmobj[[vecparm]]))
    names(outlist[[vecparm]]) <- paste0(vecparm)
  }
  names(outlist) <- NULL
  outdf <- as.data.frame(do.call("cbind", outlist))
  outdf
}

.diag_par <- function(kmobj.list, ...) {
  #' @importFrom rstan Rhat ess_bulk ess_tail
  #'
  getparmvec <- function(lst, parm) {
    lst[[parm]]
  }
  odf <- lapply(kmobj.list, .extractparms)
  nms <- names(odf[[1]])
  arr <- array(NA, c(nrow(odf[[1]]), length(odf), ncol(odf[[1]])))
  for (ll in seq_len(length(odf))) {
    arr[, ll, ] <- as.matrix(odf[[ll]])
  }
  dimnames(arr) <- list(NULL, seq_len(length(odf)), nms)
  rstan::monitor(arr, ...)
}

.diag <- function(kmobj, ...) {
  #' @importFrom rstan ess_bulk ess_tail
  odf <- .extractparms(kmobj)
  arr <- array(NA, c(nrow(odf), 1, ncol(odf)))
  arr[, 1, ] <- as.matrix(odf)
  dimnames(arr) <- list(NULL, 1, names(odf))
  rstan::monitor(arr, ...)
}

.predictivemean <- function(object, ptype=c("mean", "sd.fit"), ...) {
  #' @importFrom bkmr SamplePred
  #' @importFrom stats sd
  sf <- bkmr::SamplePred(object, ...)
  if (ptype[1] == "mean") {
    meansample <- as.numeric(drop(apply(sf, 1, mean)))
    return(meansample)
  }
  if (ptype[1] == "sd.fit") {
    sdsample <- as.numeric(drop(apply(sf, 1, sd)))
    return(sdsample)
  }
}
