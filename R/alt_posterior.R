#' @title Compute the posterior probability of the alternative hypothesis
#' @description Maps a vector of p-values to a vector of posterior probabilities
#' given a beta mixture model.
#' @param p A numeric vector of p-values
#' @param mixm A \code{betamix} object representing the beta mixture model
alt_posterior <- function(p, mixm) {

  # check for valid p-values
  if (!is.numeric(p)) {
    stop("p is not numeric.")
  }
  if (!(all(p >= 0) & all(p <= 1))) {
    stop("The input does not look like a vector of p-values.
         Some values are not in [0,1].")
  }

  # check for valid mixture model object
  if (class(mixm) != "betamix") {
    stop("Mixture model must be of class betamix")
  }

  # return 0 if there is only a uniform component
  if (length(mixm) == 1) {
    return(0)
  }

  # set inital values for posterior prob calculation
  denom <- mixm[["l0"]] * dbeta(p, 1, 1) * p
  numer <- 0

  # iterate through model components and grab the measure
  non_unif <- sum(grepl("^l", names(mixm))) - 1
  for (i in 1:non_unif) {
    l <- paste0("l", i)
    r <- paste0("r", i)
    s <- paste0("s", i)
    p_comp <- unname(mixm[[l]] * dbeta(p, mixm[[r]], mixm[[s]]) * p)
    numer <- numer + p_comp
  }

  # calculate the odds, probability
  odds <- numer / denom
  p_alt <- odds / (1 + odds)
  return(p_alt)
}
