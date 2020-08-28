`bkmrhat`: Diagnostics and multi-chain tools with Bayesian kernel machine regression (bkmr)

Bayesian kernel machine regression (BKMR) is a semi-parametric approach to Bayesian GLMs. The `bkmr` package implements BKMR under identity and probit links, but does not implement standard Bayesian diagnostics or interface with R packages that can do these diagnostics, nor does it allow easy implementation of parallel chains, which is one of the easiest ways to speed up Markov chain Monte Carlo on modern computers. This package fixes those problems.


# Quick start
    devtools::install_github("alexpkeil1/bkmrhat")

	set.seed(111)
	library(coda)
	library(bkmr)
	library(bkmrhat)
	dat <- bkmr::SimData(n = 50, M = 4)
	y <- dat$y
	Z <- dat$Z
	X <- dat$X
	set.seed(111)
	Sys.setenv(R_FUTURE_SUPPORTSMULTICORE_UNSTABLE="quiet")
	future::plan(strategy = future::multiprocess)
	# run 4 parallel Markov chains
	fitkm.list <- kmbayes_parallel(nchains=4, y = y, Z = Z, X = X, iter = 5000,
	  verbose = FALSE, varsel = FALSE)
	mcmcobj = as.mcmc.list(fitkm.list)
	summary(mcmcobj)
	# Gelman/Rubin diagnostics won't work on certain objects,
	# like delta parameters (when using variable selection),
	# so the rstan version of this will work better (does not give errors)
	 try(gelman.diag(mcmcobj))
	# lots of functions in the coda package to use
	plot(mcmcobj)
	# both of these will also fail with delta functions (when using variable selection)
	try(gelman.plot(mcmcobj))
	try(geweke.plot(mcmcobj))



	fitkm.list <- kmbayes_parallel(nchains=4, y = y, Z = Z, X = X, iter = 5000,
								   verbose = FALSE, varsel = TRUE)
	mcmcobj = as.mcmc.list(fitkm.list)
	summary(mcmcobj)
	# Gelman/Rubin diagnostics won't work on certain objects,
	# like delta parameters (when using variable selection),
	# so the rstan version of this will work better (does not give errors)
	try(gelman.diag(mcmcobj))
	# lots of functions in the coda package to use
	plot(mcmcobj)
	# both of these will also fail with delta functions (when using variable selection)
	try(gelman.plot(mcmcobj))
	try(geweke.plot(mcmcobj))
