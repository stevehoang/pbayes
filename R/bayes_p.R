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
#' @param opt_method Optimization method (see \code{?optim} for options).
#' @param ... Additional parameters to be passed to \code{bbmle::mle2}.
#' @export
bayes_p <- function(p, n_boots = 1000, alpha = 0.01, n_cores = 1,
                    subsample = 1, level_p = FALSE, max_comp = 5,
                    opt_method = "L-BFGS-B", ...) {

  # level p-values if required
  if (level_p) {
    p <- Polyfit::levelPValues(p)$pValueCorr
  }

  # check for evidence of true positives
  h <- hist(p, plot = FALSE)$counts
  if (max(h) != h[1]) {
    print("Warning: Small p-values are not strongly overreprensented.
          Inspection of p-value histogram is recommended")
  }

  # fit a beta mixture to the data
  model <- estimate_params(p, n_boots = n_boots, alpha = alpha,
                           n_cores = n_cores, subsample = subsample,
                           max_comp = max_comp, opt_method = opt_method, ...)

  # calculate the posterior probabilities of the alternative truth
  pp <- purrr::map_dbl(p, ~ alt_posterior(., model=model))

  res <- list(p_vals = p,
              post_prob = pp,
              mixture_model = model)

  class(res) <- "bayes_p"

  return(res)
}


plot.bayes_p <- function(bp, ggtheme = ggplot2::theme_bw()) {

  tp_1 <- data.frame(p_vals = bp$p_vals, post_prob = bp$post_prob)

  p1 <- ggplot2::ggplot(tp_1, ggplot2::aes(x = p_vals, y = post_prob)) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "original p-values", y = "posterior probabilities",
                  title = "Probability mapping") +
    ggtheme
  p2 <- ggplot2::ggplot(tp_1, ggplot2::aes(x = p_vals)) +
    ggplot2::geom_histogram(bins = 50) +
    ggplot2::labs(x = "original p-values", title = "Original p-values") +
    ggtheme
  p3 <- ggplot2::ggplot(tp_1, ggplot2::aes(x = post_prob)) +
    ggplot2::geom_histogram(bins = 50) +
    ggplot2::labs(x = "posterior probabilities",
                  title = "Posterior probabilities") +
    ggtheme

  # sample from each component separately
  coefs <- bp$mixture_model
  fracs <- coefs[grepl("^l", names(bp$mixture_model))]
  samps <- quantile(c(0, 1), probs = seq(0, 1, 0.01))
  res <- list()
  for (i in 1:length(fracs)) {
    if (i == 1) {
      res[[paste("component", i - 1)]] <- dunif(samps) #* fracs[i]
    } else {
      res[[paste("component", i - 1)]] <- dbeta(samps, coefs[[paste0("r", i - 1)]],
                                                coefs[[paste0("s", i - 1)]]) #* fracs[i]
    }
  }
  res <- lapply(names(res), function(x) data.frame(val = res[[x]], component = x,
                                                   p_vals = samps,
                                                   stringsAsFactors = FALSE))
  res <- do.call(dplyr::bind_rows, res)

  # res %<>% dplyr::rowwise() %>%
    # dplyr::mutate(val = val * length(bp$p_vals))

  label <- lapply(names(fracs), function(x) paste(sub("l", "component ", x), "=", signif(fracs[x], 3)))
  label <- paste(label, collapse = "\n")
  label <- paste0("mixing fractions:\n", label)

  p4 <- ggplot2::ggplot(res, ggplot2::aes(x = p_vals, y = val)) +
    ggplot2::geom_line(data = res, ggplot2::aes(color = component, y = val),
                       size = 1) +
    ggplot2::geom_label(data = NULL, ggplot2::aes(x = Inf, y = Inf, label = label),
                        hjust = 1, vjust = 1) +
    ggplot2::labs(x = "original p-values", y = "density", color = "mixture\ncomponents",
                  title = "Mixture model") +
    ggsci::scale_color_d3() +
    ggtheme

  egg::ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)

}
