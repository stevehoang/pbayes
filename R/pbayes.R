#' @title Map p-values to posterior probabilities
#' @description Given a vector of p-values from independent hypothesis
#' tests, calculate the Bayesian posterior probability of the alternative
#' hypothesis being true. The method implemeted is similar to the
#' method described by Erikson et al., 2010 and Allison et al., 2002.
#' @param p A numeric vector of p-values
#' @param n_boots A number providing the number of bootstraps used
#' to calculate the convergence statistic.
#' @param alpha A number representing the convergence statistic for fitting
#' the uniform-beta mixture model.
#' @param n_cores A number representing the number of cores to use
#' for the bootstrap calculation.
#' @param subsample The proportion of p to subsample to calculate
#' the beta mixture
#' @param max_comp A number representing the maximum number of non-uniform
#' components to include in the mixture distribution.
#' @param min_null The lower bound on the null distribution mixing fraction.
#' @param opt_method Optimization method. Options are "L-BFGS-B" or "SANN".
#' @param monotone Logical. Only use monotonically-decreasing beta components in
#' the mixture model
#' @param mask_flagpole Logical. During the fitting procedure, mask p-values
#' that are equal to 1. These values can be overrepresented due to p-value
#' calulations on discrete distributions.
#' @param level_pvals Logical. Apply a cubic polynomial fit to a portion of the
#' p-value histogram to level non-uniform trends.
#' @param ... Additional parameters to be passed to \code{bbmle::mle2}.
#' @references
#' Allison, D. B., et al. (2002). A mixture model approach for the analysis of
#' microarray gene expression data. Computational Statistics & Data Analysis,
#' 39(1), 1-20. https://doi.org/10.1016/S0167-9473(01)00046-9
#'
#' Erikson S., et al. (2010). Composite hypothesis testing: an approach built
#' on intersection-union tests and Bayesian posterior probabilities. In
#' Guerra, R., and Goldstein, D. R., (Ed.), Meta-analysis and Combining
#' Information in Genetics and Genomics. (pp. 83-93). Chapman & Hall/CRC.
#' @return A list of 1) the original p-values, 2) the posterior probabilities
#' corresponding to each p-value, and 3) parameters of the fitted mixture model.
#' @export
pbayes <- function(p, n_boots = 1000, alpha = 0.01, n_cores = 1, subsample = 1,
                   max_comp = 4, min_null = 0.4, opt_method = "L-BFGS-B",
                   monotone = TRUE, mask_flagpole = TRUE,
                   level_pvals = FALSE, ...) {

  # check for valid p-values
  if (!is.numeric(p)) {
    stop("p is not numeric.")
  }
  if (!(all(p >= 0) & all(p <= 1))) {
    stop("The input does not look like a vector of p-values.
         Some values are not in [0,1].")
  }

  # check for valid opt_method
  if (opt_method != "L-BFGS-B" & opt_method != "SANN") {
    stop("Invalid opt_method")
  }

  # handle flagpole
  if (mask_flagpole) {
    p_orig <- p
    p <- p[p < 0.999]
  } else {
    p_orig <- p
  }

  # make sure there are enough values
  if (length(p) < 100) {
    warning(paste0("The length of vector p is less than 100. This may not ",
                   "provide sufficient data for accurate analysis."))
  }

  # level p-values if necessary
  if (level_pvals) {
    p = level_p(p)
  }

  # QC check: Does it look like there are true positives?
  h <- hist(p, plot = FALSE)$counts
  if (max(h) != h[1]) {
    warning("Small p-values are not strongly overreprensented.
            Inspection of p-value histogram is recommended")
  }

  # fit a monotonic beta mixture to the data
  mixm <- fit_betamix(p, n_boots = n_boots, alpha = alpha,
                       n_cores = n_cores, subsample = subsample,
                       max_comp = max_comp, min_null = min_null,
                       opt_method = opt_method,
                       monotone = TRUE, ...)

  # calculate posterior probabilities
  pp <- sapply(p_orig, function(x) alt_posterior(x, mixm))

  # generate output
  res <- list(p_value = p_orig,
              posterior_prob = pp,
              mixture_model = mixm)

  class(res) <- "pbayes"

  return(res)

}


#' Plot Bayesian Probability Mapping Results
#'
#' Visualizes the results of Bayesian probability mapping for p-values
#' by plotting the fitted uniform-beta mixture distribution alongside
#' the histogram of the original p-values. This function provides a
#' graphical representation of how well the mixture model fits the
#' data and allows for the inspection of posterior probabilities associated
#' with each p-value.
#'
#' @param pb A `pbayes` object containing the results from Bayesian probability
#'   mapping, including the fitted mixture model and original p-values. Must be
#'   of class 'pbayes'.
#' @param hist_breaks Integer, the number of breaks (bins) to use for the histogram
#'   of p-values. Defaults to 70.
#' @param linewidth Numeric, specifying the line width for the mixture model components.
#'   Defaults to 1.
#' @param brewer_pal Character string, specifying the color palette from the RColorBrewer
#'   package to use for the different components of the mixture model. Defaults to "Set2".
#'
#' @return A ggplot object representing the fitted mixture model and the distribution
#'   of p-values. This object can be further modified using ggplot2 functions or printed
#'   to display the plot.
#'
#' @examples
#' \dontrun{
#'   # Assuming you have a pbayes object `pb` from previous analysis
#'   plot.pbayes(pb)
#'
#'   # Customizing the plot
#'   plot.pbayes(pb, hist_breaks = 100, linewidth = 2, brewer_pal = "Dark2")
#' }
#'
#' @export
plot.pbayes <- function(pb, hist_breaks = 70,
                        linewidth = 1, brewer_pal = "Set2") {

  # Ensure the input is of class pbayes
  if (!inherits(pb, "pbayes")) {
    stop("Input must be of class 'pbayes'.")
  }

  # get the number of components
  n_comp <- sum(grepl("^l", names(pb$mixture_model)))

  # loop over the components to make the curves
  sample_length <- 200
  xtemp <- seq(0, 1, length.out = sample_length) # x values for each loop
  xvals <- rep(xtemp, n_comp) # x values for plotting
  components <- character(sample_length * n_comp)
  yvals <- numeric(sample_length * n_comp)
  for (i in 1:n_comp) {
    comp_index <- i - 1
    segment_start <- (comp_index * sample_length) + 1
    segment_end <- i * sample_length
    if (comp_index == 0) { # if the uniform component
      yvals[segment_start:segment_end] <-
        dbeta(xtemp, 1, 1) * pb$mixture_model["l0"]
      components[segment_start:segment_end] <- rep("uniform", length(xtemp))
    }
    else { # if another component
      yvals[segment_start:segment_end] <-
        dbeta(xtemp, pb$mixture_model[paste0("r", comp_index)],
              pb$mixture_model[paste0("s", comp_index)]) *
        pb$mixture_model[paste0("l", comp_index)]
      components[segment_start:segment_end] <- rep(paste0("beta ", comp_index),
                                                   length(xtemp))
    }

  }

  # assemble for plotting
  curves <- data.frame(xvals, yvals, components)
  curves$components <-  factor(curves$components)
  curves$components <- relevel(curves$components, "uniform")

  # make plot
  p <- ggplot2::ggplot(data.frame(xvals = pb$p_value), ggplot2::aes(x = xvals)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            breaks = seq(0, 1, length.out = hist_breaks)) +
    ggplot2::geom_line(data = curves,
                       ggplot2::aes(color = components, y = yvals, x = xvals),
                       linewidth = linewidth) +
    ggplot2::scale_color_brewer(palette = brewer_pal) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "p-value", y = "density")

  return(p)

}

