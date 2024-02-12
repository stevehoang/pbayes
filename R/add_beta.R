#' @title Add a new component to a mixture distribution.
#' @description Adds a monotonically decreasing beta distribution to an
#' existing distribution.
#' @param p A vector of p-values
#' @param mixm A \code{betamix} object representing the beta mixture model.
#' This is typically created by the fit_beta_mix function.
#' @param opt_method Optimization method. Options are "L-BFGS-B" or "SANN".
#' @param monotone Only use monotonically-decreasing beta components in the
#' mixture model
#' @param min_null The lower bound on the null distribution mixing fraction.
#' @param ... Additional parameters to be passed to \code{bbmle::mle2}.
#' @export
add_beta <- function(p, mixm, opt_method = "L-BFGS-B", monotone = TRUE,
                     min_null = 0.4, suppress_warn = TRUE, ...) {

  # check for valid p-values
  if (!is.numeric(p)) {
    stop("p is not numeric.")
  }
  if (!(all(p >= 0) & all(p <= 1))) {
    stop("The input does not look like a vector of p-values.
         Some values are not in [0,1].")
  }

  # get number of components
  n_comp <- sum(grepl("^l", names(mixm)))

  # create a set of possible starting points to try for optimization
  n_points <- 4
  n_try <- n_points^2
  r <- seq(0.01, 0.99, length.out = n_points)
  s <- seq(1.01, 10, length.out = n_points)
  rs <- expand.grid(r, s)
  rs <- rs[order(rev(rs[,2]), rs[,1]), ]

  # set starting parameters
  frac <- 1 - ((2 * sum(p > 0.5)) / length(p)) # pick a reasonable mixing frac
  init <- mixm
  init[[paste0("l", n_comp)]] <- frac
  init[[paste0("r", n_comp)]] <- rs[1, 1] # 0.5
  init[[paste0("s", n_comp)]] <- rs[1, 2] # 5

  # normalize the mixing fractions
  norm <- init[grepl("^l", names(init))]
  norm <- sum(unlist(norm))
  init[grepl("^l", names(init))] <-
    lapply(init[grepl("^l", names(init))], function(x) x / norm)

  # set optimization boundaries for mixing fractions and
  # beta distribution parameters
  upper_bounds <- rep(1, length(init))
  names(upper_bounds) <- names(init)
  upper_bounds[grepl("^s", names(upper_bounds))] <- 100

  lower_bounds <- rep(.Machine$double.xmin, length(init))
  names(lower_bounds) <- names(init)
  lower_bounds[grepl("^l", names(lower_bounds))] <- 0
  # assume uniform component is no less than min_null of the total
  lower_bounds[grepl("l0", names(lower_bounds))] <- min_null

  if (monotone) {
    lower_bounds[grepl("^s", names(lower_bounds))] <- 1
  } else {
    upper_bounds[grepl("^r", names(upper_bounds))] <- 100
  }

  # construc negative log-likelihood function
  nll <- .construct_nll(names(init))

  # calculate new parameters
  data <- list(p = p)
  if (opt_method == "SANN") { # simulated annealing
    est <- bbmle::mle2(minuslogl = nll, start = init, data = data,
                       method = "SANN", gr=.sann_generate, ...)
  } else {
    # try different starting points for optimization if it fails to converge
    for (i in 1:n_try) {
      init[[paste0("r", n_comp)]] <- rs[i, 1]
      init[[paste0("s", n_comp)]] <- rs[i, 2]
      if (suppress_warn) {
        est <- try(suppressWarnings(bbmle::mle2(minuslogl = nll, start = init,
                                                data = data, method = opt_method,
                                                lower = lower_bounds,
                                                upper = upper_bounds, ...)),
                   silent = TRUE)

      } else {
        est <- try(bbmle::mle2(minuslogl = nll, start = init, data = data,
                               method = opt_method, lower = lower_bounds,
                               upper = upper_bounds, ...), silent = TRUE)
      }
      if (class(est) != "try-error") {
        break
      }
      print(paste0("Attempting optimization with new initial values. Attempt ",
                   i))
    }
    # If L-BFGS-B fails ...
    if (class(est) == "try-error") {
      stop('MLE did not converge. Try running with opt_method = "SANN"')
    }
  }

  # get new coefficients
  new_coefs <- est@fullcoef

  # normalize the mixing fractions and return the result
  norm_fac <- sum(new_coefs[grepl("^l", names(new_coefs))])
  new_coefs[grepl("^l", names(new_coefs))] <-
    new_coefs[grepl("^l", names(new_coefs))] / norm_fac
  new_coefs <- new_coefs[order(names(new_coefs))]

  class(new_coefs) <- "betamix"

  return(new_coefs)
}
