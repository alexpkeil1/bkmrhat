cat("Testing parallel")
library(bkmrhat)
set.seed(111)
dat <- bkmr::SimData(n = 50, M = 4)
y <- dat$y
Z <- dat$Z
X <- dat$X
set.seed(111)

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

future::plan(strategy = future::multisession, workers=2)
fitkm.list <- kmbayes_parallel(nchains=2, y = y, Z = Z, X = X, iter = 10,
                               verbose = FALSE, varsel = TRUE)

kmbayes_combine(fitkm.list, burnin=0)
kmbayes_combine(fitkm.list, burnin=8, reorder = TRUE, excludeburnin=FALSE)
kmbayes_combine(fitkm.list, burnin=8, reorder = TRUE, excludeburnin=TRUE)
kmbayes_combine(fitkm.list, burnin=8, reorder = FALSE, excludeburnin=TRUE)
