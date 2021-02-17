cat("Testing for non duplicated results in seeded simulations")

library(bkmrhat)
set.seed(111)
dat <- bkmr::SimData(n = 50, M = 4)
y <- dat$y
Z <- dat$Z
X <- dat$X

future::plan(strategy = future::sequential)
fitkm.list <- kmbayes_parallel(nchains=2, y = y, Z = Z, X = X, iter = 10,
                               verbose = FALSE, varsel = TRUE)

rr = as.mcmc.list(fitkm.list)
tt = all.equal(as.data.frame(rr[[1]]), as.data.frame(rr[[2]]))
stopifnot(mode(tt)=="character")


fitkm.list[[1]] = fitkm.list[[2]]
badpred = SamplePred_parallel(fitkm.list, method="exact", Znew = dat$Z[1,], Xnew = dat$X[1,])
ff = tapply(badpred$znew1, badpred$chain, mean)
tt = all.equal(ff[1], ff[2])
stopifnot(mode(tt)=="character")
