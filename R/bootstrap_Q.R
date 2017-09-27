
bootstrap_Q <- function(p, c0, c1, n_boots, n_cores) {

  # define mixture sampling function
  n <- length(p)
  f <- function() {
    samp <- sample_beta_mix(n, c0)
    ll0 <- loglike(samp, c0)
    ll1 <- loglike(samp, c1)
    qw <- 2 * (ll1 - ll0)
    return(qw)
  }

  # sample and calculate Qs in parallel
  doMC::registerDoMC(cores = n_cores)
  Qw <- foreach::foreach(i = 1:n_boots, .combine=c) %dopar% (f())
  Qo <- 2 * (loglike(p, c1) - loglike(p, c0))
  Q <- list(Qo = Qo, Qw = Qw)
  return(Q)
}
