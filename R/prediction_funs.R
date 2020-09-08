### functions for prediction


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
  requireNamespace("stats")
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

