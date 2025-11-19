`bkmrhat` v1.1.5




Diagnostics and multi-chain tools with Bayesian kernel machine regression (bkmr)

Bayesian kernel machine regression (BKMR) is a semi-parametric approach to Bayesian GLMs. The `bkmr` package implements BKMR under identity and probit links, but does not implement standard Bayesian diagnostics or interface with R packages that can do these diagnostics, nor does it allow easy implementation of parallel chains, which is one of the easiest ways to speed up Markov chain Monte Carlo on modern computers. This package fixes those problems.

Note that this package functions best using the development verision of the bkmr package (instructions [here](https://github.com/jenfb/bkmr/blob/master/README.md))


# Quick start

#### install
    install.packages("devtools")
    devtools::install_github("alexpkeil1/bkmrhat", build_vignettes = TRUE)

#### simulate data from a `bkmr` package function
	set.seed(111)
	library(coda)
	library(bkmr)
	library(bkmrhat)
	dat <- bkmr::SimData(n = 50, M = 4)
	y <- dat$y
	Z <- dat$Z
	X <- dat$X
	set.seed(111)
	

#### fit BKMR model via four parallel MCMC chains (10,000 iterations total)
	future::plan(strategy = future::multisession)
	# run 4 parallel Markov chains
	fitkm.list <- kmbayes_parallel(nchains=4, y = y, Z = Z, X = X, iter = 2500,
	  verbose = FALSE, varsel = TRUE)

#### run diagnostics from `rstan` package
    # rstan has a few excellent diagnostics (modern r-hat and effective sample size)
    diagres = kmbayes_diag(fitkm.list)
    # estimate of standard error (gives margin of error for reporting estimates)
    diagres[,"se_mean"]

#### convert to `coda` package object for other diagnostics
	mcmcobj = as.mcmc.list(fitkm.list)
	# lots of functions in the coda package to use
	# get info on rejection probabilities (from h() function in bkmr - won't be correct for discrete parameters like delta parameters)
	  rejectionRate(mcmcobj)
	# check trace plots for obvious issues between/within chains
	  plot(mcmcobj)
	# autocorrelation for efficiency issues
	  acfplot(mcmcobj)
	# posterior correlation for
	  crosscorr(mcmcobj)
	# Gelman/Rubin diagnostics for convergence (multivariate may fail when including delta functions from variable selection)
	  try(gelman.diag(mcmcobj, multivariate = FALSE))  
    # batch standard error: a different approach to estimating standard error
      batchSE(kmres)
    # effective size, another way to assess whether enough samples have been made (I prefer the `rstan` implementations)
      effectiveSize(mcmcobj)
      try(geweke.diag(mcmcobj))
    # posterior summary
      summary(mcmcobj)
    # HPD intervals by chain
      HPDinterval(mcmcobj)
	  

#### combine chains and use `bkmr` native posterior functions
    combobj = comb_bkmrfits(fitkm.list)
    combobj$iter
    summary(combobj) # rejection rates will be slightly off for multi-chain objects
    # mean difference function from `bkmr` package (default burnin of half of total number of iterations)
      mdiff = OverallRiskSummaries(fit = combobj,
                                      qs = seq(.05, 0.95, by = 0.1),
                                      q.fixed = 0.5, method = "exact")
    mdiff
