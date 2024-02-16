#' @title Fit a uniform-beta mixture to a distribution of p-values.
#' @description Estimates the parameters of a uniform-beta mixture distribution,
#' given some data. The mixture permits multiple betas and a single uniform distribution.
#' @param p A numeric vector of p-values
#' @param n_boots A number providing the number of bootstraps used
#' to calculate the convergence statistic.
#' @param alpha A number representing the convergence criterion for formulating the mixture.
#' @param n_cores A number representing the number of cores to use
#' for the bootstrap calculation.
#' @param subsample The proportion of p to subsample to calculate the beta mixture.
#' @param max_comp A number representing the maximum number of non-uniform
#' components to include in the mixture distribution.
#' @param min_null The lower bound on the null distribution mixing fraction.
#' @param opt_method Optimization method. Options are "L-BFGS-B" or "SANN".
#' @param monotone Logical. Only use monotonically-decreasing beta components in the
#' mixture model.
#' @param ... Additional parameters to be passed to \code{bbmle::mle2}.
#' export
fit_betamix <- function(p, n_boots = 1000, alpha = 0.01,
                         n_cores = 1, subsample = 1,
                         max_comp = 4, min_null = 0.4,
                         opt_method = "L-BFGS-B",
                         monotone = TRUE, ...) {

  # check for valid p-values
  if (!is.numeric(p)) {
    stop("p is not numeric.")
  }
  if (!(all(p >= 0) & all(p <= 1))) {
    stop(paste0("The input does not look like a vector of p-values. ",
                "Some values are not in [0,1]."))
  }

  # subsample p-values
  if (!((subsample <= 1) & (subsample > 0))) {
    stop("subsample must be on the interval (0, 1]")
  }
  if (subsample < 1) {
    s <- round(length(p) * subsample)
    p <- sample(p, s)
  }

  # two possible models m0 and m1
  mixm0 <- list("l0" = 1)
  class(mixm0) <- "betamix"
  mixm1 <- add_beta(p, mixm0, opt_method = opt_method, monotone = monotone,
                    min_null = min_null, ...)

  # compare models
  sig <- compare_betamixes(p, mixm0, mixm1, n_boots, n_cores)

  # add components until the stopping criteria are met
  n_comp <- 1
  while ((sig < alpha) & (n_comp <= max_comp)) {

    # return a uniform distribution if the mixture is uniform
    if (mixm1["l0"] == 1) {
      return(mixm0)
    }
    if ((mixm1[paste0("r", n_comp)] == 1) & (mixm1[paste0("s", n_comp)] == 1)) {
      return(mixm0)
    }

    # return the previous model if lambda = 0 for the new component
    if (mixm1[paste0("l", n_comp)] == 0) {
      return(mixm0)
    }

    # update parameters
    mixm0 <- mixm1
    print(mixm0)
    mixm1 <- add_beta(p, mixm0, opt_method = opt_method, monotone = monotone,
                      min_null = min_null, ...)
    sig <- compare_betamixes(p, mixm0, mixm1, n_boots, n_cores)
    n_comp <- n_comp + 1
  }

  return(mixm0)
}
