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
