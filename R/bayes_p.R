bayes_p <- function(p, n_boots = 1000, alpha = 0.01, n_cores = 1,
                    subsample = 1, level_p = FALSE) {

  # make zeros very small positive reals
  p[p == 0] <- .Machine[["double.xmin"]]

  # level p-values if required
  if (level_p) {
    p <- Polyfit::levelPValues(p)$pValueCorr
  }

  # check for evidence of null hypothesis rejections
  h <- hist(p)$counts
  if (max(h) != h[1]) {
    print("Warning: Small p-values are not strongly overreprensented.
          Inspection of p-value histogram is recommended")
  }

  model <- .estimate_params(p, n_boots = n_boots, alpha = alpha,
                            n_cores = n_cores, subsample = subsample)

  pp <- purrr::map_dbl(p, ~ alt_posterior(., model=model))

  non_unif <- sum(grepl("^l", names(model))) - 1
  null_frac - model["l0"]
  if (null_frac != 1) {
    for (i in 1:non_unif) {
      l <- paste0("l", i)
      r <- paste0("r", i)
      s <- paste0("s", i)
      if (model[r] > model[s]) {
        null_frac <- null_frac + model[l]
      }
    }
  }
  de_frac <- 1 - null_frac
  res <- list(p_vals = p,
              post_prob = pp,
              mixture_model = model)
}
