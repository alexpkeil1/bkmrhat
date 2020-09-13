cat("Testing coda")
library(bkmrhat)
set.seed(111)
dat <- bkmr::SimData(n = 50, M = 4)
y <- dat$y
Z <- dat$Z
X <- dat$X
set.seed(111)
Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
future::plan(strategy = future::multiprocess)
fitkm.list <- kmbayes_parallel(nchains=2, y = y, Z = Z, X = X, iter = 10,
                               verbose = FALSE, varsel = TRUE)

as.mcmc.list(fitkm.list)

as.mcmc(fitkm.list[[1]])
as.mcmc(comb_bkmrfits(fitkm.list))
