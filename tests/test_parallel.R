cat("Testing parallel")
library(bkmrhat)
set.seed(111)
dat <- bkmr::SimData(n = 50, M = 4)
y <- dat$y
Z <- dat$Z
X <- dat$X
set.seed(111)
Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
future::plan(strategy = future::sequential)
fitkm.list <- kmbayes_parallel(nchains=2, y = y, Z = Z, X = X, iter = 10,
                               verbose = FALSE, varsel = TRUE)
PredictorResponseBivar_parallel(fitkm.list)
PredictorResponseUnivar_parallel(fitkm.list)
SingVarRiskSummaries_parallel(fitkm.list)
ExtractPIPs_parallel(fitkm.list)
OverallRiskSummaries_parallel(fitkm.list)
SamplePred_parallel(fitkm.list)
comb_bkmrfits(fitkm.list)
kmbayes_diag(fitkm.list)
closeAllConnections()

future::plan(strategy = future::multiprocess, workers=2)
fitkm.list <- kmbayes_parallel(nchains=2, y = y, Z = Z, X = X, iter = 10,
                               verbose = FALSE, varsel = TRUE)
