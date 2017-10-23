#' @title Sample a mixture of beta distributions
#' @description Draw a random sample from a beta mixture model
#' @param n Number of observations
#' @param coefs A named list providing the parameters of the beta mixture
sample_beta_mix <- function(n, coefs) {

  # get mixing fractions of components
  fracs <- coefs[grepl("^l", names(coefs))]
  fracs %<>% sort()
  n_comp <- length(fracs)

  # calculate intervals proprotional to the mixing fractions
  intervals <- c()
  for (i in 1:(n_comp - 1)) {
    if (i == 1) {
      intervals <- c(intervals, fracs[i])
    } else {
      intervals <- c(intervals, fracs[i] + intervals[i - 1])
    }
  }

  # sample the mixture
  u <- runif(n) # random uniform draws on [0,1]
  samps <- c()
  for (i in u) {
    ind <- sum(intervals < i) + 1 # ID the interval
    comp <- names(fracs)[ind] # get the component associated with the interval
    # sample from the appropriate component
    if (comp == "l0") {
      s <- rbeta(1, 1, 1)
    } else {
      s <- rbeta(1, coefs[sub("l", "r", comp)], coefs[sub("l", "s", comp)])
    }
    samps <- c(samps, s)
  }
  return(samps)
}
