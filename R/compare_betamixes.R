#' @title Compute statistics for model selection
#' @importFrom foreach "%dopar%"
#' @description Compute the p-value for the improvement of one uniform-beta mixture over another
#' with respect to a given set of p-values
#' @param p A numeric vector of p-values
#' @param mixm0 A \code{betamix} object that defines a uniform-beta mixture to be tested against
#' @param mixm1 A \code{betamix} object that defines a uniform-beta mixture to be tested
#' @param n_boots Number of bootstraps to perform
#' @param n_cores Number of cores to use to perform the bootstrapping
#' @export
compare_betamixes <- function(p, mixm0, mixm1, n_boots, n_cores) {

  # define mixture sampling function
  n <- length(p)
  f <- function() {
    samp <- sample_betamix(mixm0, n)
    ll0 <- loglike_ubeta(samp, mixm0)
    ll1 <- loglike_ubeta(samp, mixm1)
    qw <- 2 * (ll1 - ll0)
    return(qw)
  }

  # sample and calculate Qs in parallel
  doMC::registerDoMC(cores = n_cores)
  Qw <- foreach::foreach(i = 1:n_boots, .combine=c) %dopar% (f())
  Qo <- 2 * (loglike_ubeta(p, mixm1) - loglike_ubeta(p, mixm0))

  # calculate p-value
  sig <- sum(Qw > Qo) / n_boots

  return(sig)
}
