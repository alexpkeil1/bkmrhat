cat("Testing parallel")
library(bkmrhat)
set.seed(111)
dat <- bkmr::SimData(n = 100, M = 4)
expit <- function(mu) 1/(1+exp(-mu))
y <- rbinom(100, 1, expit(dat$y))
Z <- dat$Z
X <- dat$X
set.seed(111)
Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
future::plan(strategy = future::multiprocess)
fitkm.list <- kmbayes_parallel(nchains=2, y = y, Z = Z, X = X, iter = 10,
                               verbose = FALSE, varsel = TRUE, family="binomial")
sinkit = kmbayes_diag(fitkm.list)
xx = comb_bkmrfits(fitkm.list)
sinkit = suppressWarnings(predict(xx))
sinkit = suppressWarnings(predict(xx, ptype="sd.fit"))

closeAllConnections()
