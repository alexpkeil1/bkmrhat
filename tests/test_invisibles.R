#### testing invisible functions ####
library(bkmrhat)

set.seed(111)
dat <- bkmr::SimData(n = 50, M = 5, ind=1:3, Zgen="realistic")
y <- dat$y
Z <- dat$Z
X <- cbind(dat$X, rnorm(50))

# run 10 initial iterations for a model with only 2 exposures
Z2 = Z[,1:2]
kmfitbma.start <- suppressWarnings(bkmr::kmbayes(y = y, Z = Z2, X = X, iter = 10, verbose = FALSE, varsel = TRUE, est.h = TRUE))

#.extractparms
suppressWarnings(bkmrhat:::.extractparms(kmfitbma.start))

# .diag_par
suppressWarnings(bkmrhat:::.diag_par(list(kmfitbma.start, kmfitbma.start)))

# .diag
suppressWarnings(bkmrhat:::.diag(kmfitbma.start))

# .predictivemean
suppressWarnings( bkmrhat:::.predictivemean(kmfitbma.start))
suppressWarnings(bkmrhat:::.predictivemean(kmfitbma.start, ptype='sd.fit'))


#.checkver
bkmrhat:::.checkver(print=TRUE)


#.add_bkmrfits
bkmrhat:::.add_bkmrfits(list(kmfitbma.start,kmfitbma.start), trim=FALSE)

