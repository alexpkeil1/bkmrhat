cat("Testing diag")
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

kmbayes_diagnose(fitkm.list)
combobj = comb_bkmrfits(fitkm.list)
stopifnot(inherits(combobj, "bkmrfit"))
kmbayes_diagnose(combobj)
closeAllConnections()
