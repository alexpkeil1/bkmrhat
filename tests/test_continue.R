library(bkmrhat)

set.seed(111)
dat <- bkmr::SimData(n = 50, M = 5, ind=1:3, Zgen="realistic")
y <- dat$y
Z <- dat$Z
X <- cbind(dat$X, rnorm(50))

# run 10 initial iterations for a model with only 2 exposures
Z2 = Z[,1:2]
kmfitbma.start <- suppressWarnings(bkmr::kmbayes(y = y, Z = Z2, X = X, iter = 10, verbose = FALSE, varsel = TRUE, est.h = TRUE))

# run 20 additional iterations
moreiterations = suppressWarnings(kmbayes_continue(kmfitbma.start, iter=20))

stopifnot(kmfitbma.start$iter<moreiterations$iter)
stopifnot(all(kmfitbma.start$sigsq.eps %in% moreiterations$sigsq.eps))
stopifnot(all(kmfitbma.start$beta[,1] %in% moreiterations$beta[,1]))
stopifnot(all(kmfitbma.start$r[,1] %in% moreiterations$r[,1]))
stopifnot(all(kmfitbma.start$h.hat[,1] %in% moreiterations$h.hat[,1]))
stopifnot(ncol(kmfitbma.start$beta) == ncol(moreiterations$beta))
stopifnot(ncol(kmfitbma.start$r) == ncol(moreiterations$r))
stopifnot(ncol(kmfitbma.start$h.hat) == ncol(moreiterations$h.hat))
