
.bootstrap_Q <- function(x, c0, c1, n_boots, n_cores) {
  n <- length(x)
  f <- function() {
    samp <- sample_beta_mix(n, c0)
    ll0 <- loglike(samp, c0)
    ll1 <- loglike(samp, c1)
    qw <- 2 * (ll1 - ll0)
    return(qw)
  }
  doMC::registerDoMc(cores = n_cores)
  Qw <- foreach(i = 1:n_boots, .combine=c) %dopar% (f())
  Qo <- 2 * (loglike(x, c1) - loglike(x, c0))
  Q <- list(Qo = Qo, Qw = Qw)
  return(Q)
}
