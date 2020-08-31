#### actions on package load ####

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Diagnostics and parallel chain functioning for Bayesian kernel machine regression")
  #require("coda")
}
