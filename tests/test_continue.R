library(bkmrhat)

set.seed(111)
dat <- bkmr::SimData(n = 50, M = 5, ind=1:3, Zgen="realistic")
y <- dat$y
Z <- dat$Z
X <- cbind(dat$X, rnorm(50))

# run 10 initial iterations for a model with only 2 exposures
Z2 = Z
kmfitbma.start <- suppressWarnings(bkmr::kmbayes(y = y, Z = Z2, X = X, iter = 10, verbose = FALSE, varsel = TRUE, est.h = TRUE))

# run 20 additional iterations
moreiterations = suppressWarnings(kmbayes_continue(kmfitbma.start, iter=20))
res = kmbayes_diag(moreiterations)

#bkmr::TracePlot(moreiterations, par="r", comp=5)
#bkmr::TracePlot(moreiterations, par="beta", comp=1)
#bkmr::TracePlot(moreiterations, par="h", comp=50)


stopifnot(kmfitbma.start$iter<moreiterations$iter)
stopifnot(all(kmfitbma.start$sigsq.eps %in% moreiterations$sigsq.eps))
stopifnot(all(kmfitbma.start$beta[,1] %in% moreiterations$beta[,1]))
stopifnot(all(kmfitbma.start$r[,1] %in% moreiterations$r[,1]))
stopifnot(all(kmfitbma.start$h.hat[,1] %in% moreiterations$h.hat[,1]))
stopifnot(ncol(kmfitbma.start$beta) == ncol(moreiterations$beta))
stopifnot(ncol(kmfitbma.start$r) == ncol(moreiterations$r))
stopifnot(ncol(kmfitbma.start$h.hat) == ncol(moreiterations$h.hat))


# now in paralelel
kmfitbma.start2 <- suppressWarnings(kmbayes_parallel(nchains=2,y = y, Z = Z2, X = X, iter = 10, verbose = FALSE, varsel = TRUE, est.h = FALSE))

# run 20 additional iterations
moreiterations2 = suppressWarnings(kmbayes_parallel_continue(kmfitbma.start2, iter=20))
res2 = kmbayes_diag(moreiterations2)


stopifnot(kmfitbma.start2[[1]]$iter < moreiterations2[[1]]$iter)
stopifnot(all(kmfitbma.start2[[1]]$sigsq.eps %in% moreiterations2[[1]]$sigsq.eps))
stopifnot(all(kmfitbma.start2[[1]]$beta[,1] %in% moreiterations2[[1]]$beta[,1]))
stopifnot(all(kmfitbma.start2[[1]]$r[,1] %in% moreiterations2[[1]]$r[,1]))
stopifnot(all(kmfitbma.start2[[1]]$h.hat[,1] %in% moreiterations2[[1]]$h.hat[,1]))
stopifnot(ncol(kmfitbma.start2[[1]]$beta) == ncol(moreiterations2[[1]]$beta))
stopifnot(ncol(kmfitbma.start2[[1]]$r) == ncol(moreiterations2[[1]]$r))
stopifnot(ncol(kmfitbma.start2[[1]]$h.hat) == ncol(moreiterations2[[1]]$h.hat))


# just see if it will work with probit model
y <- 1.0*(dat$y>median(dat$y))
fitty1 = suppressWarnings(bkmr::kmbayes(y=y,Z=Z,X=X, est.h=TRUE, iter=5, family="binomial"))
# do some diagnostics here to see if 1000 iterations (default) is enough
# add 3000 additional iterations
fitty2 = suppressWarnings(kmbayes_continue(fitty1, iter=5))
stopifnot(ncol(fitty1$ystar[,1]) %in% ncol(fitty2$ystar[,1]))



# force old version
kmfitbma.start2 = kmfitbma.start
kmfitbma.start2$delta = kmfitbma.start2$delta*0
moreiterations = suppressWarnings(kmbayes_continue(kmfitbma.start2, iter=20))
res = kmbayes_diag(moreiterations)
