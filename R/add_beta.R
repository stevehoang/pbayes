#' @title Add a new component to a mixture distribution.
#' @description Adds a monotonically decreasing beta distribution to an
#' existing distribution.
#' @param p A numeric vector of p-values
#' @param coefs A named numeric vector of coefficients representing the
#' initial model.
#' @param sann Logical. Use simulated annealing (defaults to FALSE).
#' @param ... Additional parameters to be passed to \code{bbmle::mle2}.
add_beta <- function(p, coefs, sann = FALSE, ...) {

  # get number of components
  n_comp <- sum(grepl("^l", names(coefs)))
  data <- list(p = p)

  # set starting parameters
  init <- as.list(coefs)
  init[[paste0("l", n_comp)]] <- 0.1
  init[[paste0("r", n_comp)]] <- 1
  init[[paste0("s", n_comp)]] <- 1

  # normalize the mixing fractions
  norm <- init[grepl("^l", names(init))] %>%
    unlist() %>%
    sum()
  init[grepl("^l", names(init))] <-
    lapply(init[grepl("^l", names(init))], function(x) x / norm)

  # set boundaries for mixing fractions and
  # beta distribution parameters
  upper_bounds <- rep(1, length(init))
  names(upper_bounds) <- names(init)
  upper_bounds[grepl("^s", names(upper_bounds))] <- 100

  lower_bounds <- rep(.Machine$double.xmin, length(init))
  names(lower_bounds) <- names(init)
  lower_bounds[grepl("^l", names(lower_bounds))] <- 0
  lower_bounds[grepl("^s", names(lower_bounds))] <- 1

  # construc negative log-likelihood function
  nll <- construct_nll(names(init))

  # calculate new coefficients
  if (sann == TRUE) { # simulated annealing
    est <- bbmle::mle2(minuslogl = nll, start = init, data = data,
                       method = "SANN", gr=sann_generate, ...)
  } else {
    est <- try(bbmle::mle2(minuslogl = nll, start = init, data = data, method = "L-BFGS-B",
         lower = lower_bounds, upper = upper_bounds, ...))
    if (class(est) == "try-error") {
      stop("Error in MLE. Try running with sann = TRUE.")
    }
  }

  # get new coefficients
  new_coefs <- est@fullcoef
  new_coefs[grepl("^l", names(new_coefs))] %<>% `/`(., sum(.))
  new_coefs <- new_coefs[order(names(new_coefs))]

  return(new_coefs)
}
