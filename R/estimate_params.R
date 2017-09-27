
estimate_params <- function(p, n_boots = 1000, alpha = 0.01,
                             n_cores = 1, subsample = 1, max_comp = 5) {

  # subsample p-values if necessary
  if (!((subsample <= 1) & (subsample > 0))) {
    stop("error: subsample must be on the interval (0, 1]")
  }
  if (subsample < 1) {
    s <- round(length(p) * subsample)
    p <- sample(p, s)
  }

  # two possible models
  coefs0 <- c("l0" = 1)
  coefs1 <- add_beta(p, coefs0)

  q <- bootstrap_Q(p, coefs0, coefs1, n_boots, n_cores)
  sig <- sum(q$Qw > q$Qo) / n_boots

  # add components until the stopping critera are met
  n_comp <- 1
  while ((sig < alpha) & (n_comp <= max_comp)) {

    # return a uniform distribution if the mixture is uniform
    if (coefs1["l0"] == 1) {
      return(coefs0)
    }
    if ((coefs1[paste0("r", n_comp)] == 1) & (coefs1[paste0("s", n_comp)] == 1)) {
      return(coefs0)
    }
    # return the previous model if lambda = 0 for the new component
    if (coefs1[paste0("l", n_comp)] == 0) {
      return(coefs0)
    }

    # update coefficients
    coefs0 <- coefs1
    print(coefs0)
    coefs1 <- add_beta(p, coefs0)
    q <- bootstrap_Q(p, coefs0, coefs1, n_boots, n_cores)
    sig <- sum(q$Qw > q$Qo) / n_boots
    n_comp <- n_comp + 1
  }
  return(coefs0)
}
