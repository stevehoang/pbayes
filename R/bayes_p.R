#' @title Map p-values to posterior probabilities
#' @description Given a vector of p-values from independent hypothesis
#' tests, calculate the Bayesian posterior probability of the alternative
#' hypothesis being true. The method implemeted is approximately that
#' which is described by Erikson et al., 2006 and Allison et al., 2002.
#' The only difference being that the mixture of beta distribtutions used
#' to model the distribution of p-values only permits monotonically decreasing
#' betas.
#' @param p A numeric vector of p-values'
#' @param n_boots A number providing the number of bootstraps used
#' to calculate the convergence statistic.
#' @param alpha A number representing the convergence statistic.
#' @param n_cores A number representing the number of cores to use
#' for the bootstrap calculation.
#' @param subsample The proportion of p to subsample to calculate
#' the beta mixture
#' @param level_p Logical. Apply \code{Polyfit::levelPvalues} to \code{p}
#' before the fitting procedure.
#' @param max_comp A number representing the maximum number of non-uniform
#' components to include in the mixture distribution.
#' @export
bayes_p <- function(p, n_boots = 1000, alpha = 0.01, n_cores = 1,
                    subsample = 1, level_p = FALSE, max_comp = 5) {

  # level p-values if required
  if (level_p) {
    p <- Polyfit::levelPValues(p)$pValueCorr
  }

  # check for evidence of true positives
  h <- hist(p)$counts
  if (max(h) != h[1]) {
    print("Warning: Small p-values are not strongly overreprensented.
          Inspection of p-value histogram is recommended")
  }

  # fit a beta mixture to the data
  model <- estimate_params(p, n_boots = n_boots, alpha = alpha,
                           n_cores = n_cores, subsample = subsample,
                           max_comp = max_comp)

  # calculate the posterior probabilities of the alternative truth
  pp <- purrr::map_dbl(p, ~ alt_posterior(., model=model))

  res <- list(p_vals = p,
              post_prob = pp,
              mixture_model = model)
}
