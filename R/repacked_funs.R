#### repackaged functions from bkmr package ####
# todo:
# may also want parallel implementations of:
# CalcWithinGroupPIPs

OverallRiskSummaries_parallel <- function(x, ...){
  #' Overall summary by chain
  #' @param x bkmrfit.list object from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param ... arguments to \code{\link[bkmr]{OverallRiskSummaries}}
  #'
  #' @return data.frame with all chains together
  #' @importFrom bkmr OverallRiskSummaries
  #' @importFrom future value
  #' @importFrom stats runif
  #' @export
  #'
  ff <- list()
  nchains = length(x)
  ss = round(runif(nchains) * .Machine$integer.max)
  for (ii in 1:nchains) {
    xii = x[[ii]]
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      df = suppressWarnings(bkmr::OverallRiskSummaries(xii, ...))
      df$chain=ii
      df
    }, seed=ss[ii])
  }
  res <- value(ff)
  as.data.frame(do.call("rbind", res))
}

PredictorResponseUnivar_parallel <- function(x, ...){
  #' Univariate predictor response summary by chain
  #' @param x bkmrfit.list object from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param ... arguments to \code{\link[bkmr]{PredictorResponseUnivar}}
  #'
  #' @return data.frame with all chains together
  #' @importFrom bkmr PredictorResponseUnivar
  #' @importFrom stats runif
  #' @importFrom future value
  #' @export
  #'
  ff <- list()
  nchains = length(x)
  ss = round(runif(nchains) * .Machine$integer.max)
  for (ii in 1:nchains) {
    xii = x[[ii]]
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      df = suppressWarnings(bkmr::PredictorResponseUnivar(xii, ...))
      df$chain=ii
      df
    }, seed=ss[ii])
  }
  res <- value(ff)
  as.data.frame(do.call("rbind", res))
}


PredictorResponseBivar_parallel <- function(x, ...){
  #' Bivariate predictor response by chain
  #' @param x bkmrfit.list object from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param ... arguments to \code{\link[bkmr]{PredictorResponseBivar}}
  #'
  #' @return data.frame with all chains together
  #' @importFrom bkmr PredictorResponseBivar
  #' @importFrom stats runif
  #' @importFrom future value
  #' @export
  #'
  ff <- list()
  nchains = length(x)
  ss = round(runif(nchains) * .Machine$integer.max)
  for (ii in 1:nchains) {
    xii = x[[ii]]
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      df = suppressWarnings(bkmr::PredictorResponseBivar(xii, ...))
      df$chain=ii
      df
    }, seed=ss[ii])
  }
  res <- value(ff)
  as.data.frame(do.call("rbind", res))
}


SingVarRiskSummaries_parallel <- function(x, ...){
  #' Single variable summary by chain
  #' @param x bkmrfit.list object from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param ... arguments to \code{\link[bkmr]{SingVarRiskSummaries}}
  #'
  #' @return data.frame with all chains together
  #' @importFrom bkmr SingVarRiskSummaries
  #' @importFrom future value
  #' @importFrom stats runif
  #' @export
  #'
  ff <- list()
  nchains = length(x)
  ss = round(runif(nchains) * .Machine$integer.max)
  for (ii in 1:nchains) {
    xii = x[[ii]]
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      df = suppressWarnings(bkmr::SingVarRiskSummaries(xii, ...))
      df$chain=ii
      df
    }, seed=ss[ii])
  }
  res <- value(ff)
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
#    }, seed=TRUE)
#  }
#  res <- value(ff)
#  as.data.frame(do.call("rbind", res))
#}

ExtractPIPs_parallel <- function(x, ...){
  #' Posterior inclusion probabilities by chain
  #' @param x bkmrfit.list object from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param ... arguments to \code{\link[bkmr]{CalcPIPs}}
  #'
  #' @return data.frame with all chains together
  #' @importFrom bkmr ExtractPIPs
  #' @importFrom future value
  #' @importFrom stats runif
  #' @export
  #'
  ff <- list()
  nchains = length(x)
  ss = round(runif(nchains) * .Machine$integer.max)
  for (ii in 1:nchains) {
    xii = x[[ii]]
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      df = suppressWarnings(data.frame(bkmr::ExtractPIPs(xii)))
      df$chain=ii
      df
    }, seed=ss[ii])
  }
  res <- value(ff)
  as.data.frame(do.call("rbind", res))
}


SamplePred_parallel <- function(x, ...){
  #' Posterior samples of E(Y|h(Z),X,beta) by chain
  #' @param x bkmrfit.list object from \code{\link[bkmrhat]{kmbayes_parallel}}
  #' @param ... arguments to \code{\link[bkmr]{CalcPIPs}}
  #'
  #' @return data.frame with all chains together
  #' @importFrom bkmr SamplePred
  #' @importFrom future value
  #' @importFrom stats runif
  #' @export
  #'
  ff <- list()
  nchains = length(x)
  ss = round(runif(nchains) * .Machine$integer.max)
  for (ii in 1:nchains) {
    xii = x[[ii]]
    ff[[ii]] <- future({
      cat(paste("Chain", ii, "\n"))
      df = suppressWarnings(as.data.frame(bkmr::SamplePred(xii, ...)))
      df$chain=ii
      df
    }, seed=ss[ii])
  }
  res <- value(ff)
  as.data.frame(do.call("rbind", res))
}


