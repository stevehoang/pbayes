
.estimate_params <- function(x, n_boots = 1000, alpha = 0.01,
                             n_cores = 1, subsample = 1, max_comp = Inf) {

  # subsample p-values if necessary
  if (!((subsample <= 1) & (subsample > 0))) {
    stop("error: subsample must be on the interval (0, 1]")
  }
  if (subsample < 1) {
    s <- round(length(x) * subsample)
    x <- sample(x, s)
  }

  coefs0 <- c("l0" = 1)
  coefs1 <- add_beta(x, coefs0)
  q <- .bootstrap_Q(x, coefs0, coefs1, n_boots, n_cores)
  sig <- sum(q$Qw > q$Qo) / n_boots
  n_comp <- 1

  while ((sig > alpha) & (n_comp <= max_comp)) {
    coefs0 <- coefs1
    print(coefs0)
    coefs1 <- .add_beta(x, coefs0)
    q <- .bootstrap_Q(x, coefs0, coefs1, n_boots, n_cores)
    sig <- sum(q$Qw > q$Qo) / n_boots
    n_comp <- n_comp + 1
  }

  return(coefs0)

}
