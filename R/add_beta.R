.add_beta <- function(ps, coefs) {

  n_comp <- sum(grepl("^l", names(coefs)))
  data <- list(x = ps)

  # set boundaries for mixing fractions and
  # beta distribution parameters
  upper_bounds <- c(rep(1, n_comp + 1),
                    rep(1, n_comp),
                    rep(Inf, ncomp))
  names(upper_bounds) <- names(data)

  lower_bounds <- c(rep(0, n_comp * 2 + 1),
                    rep(1, n_comp))
  names(upper_bounds) <- names(data)

  # set starting parameters
  start[[paste0("l", n_comp)]] <- 0.5
  start[[paste0("r", n_comp)]] <- 1
  start[[paste0("s", n_comp)]] <- 1

  # calculate new coefficients
  est <- mle2(minuslogl=negLogLike, start=start, data=data, method = "L-BFGS-B",
       lower = lower, upper = upper) # this should be modified to use SANN or some such

  new_coefs <- est@fullcoef
  new_coefs[grepl("^l", names(new_coefs))] %<>% `/`(., sum(.))
  new_coefs <- new_coefs[order(names(new_coefs))]

  retrun(new_coefs)
}
