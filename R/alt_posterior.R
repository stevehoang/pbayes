
alt_posterior <- function(p, model) {

  # return 0 if there is only a uniform component
  if (length(model) == 1) {
    return(0)
  }

  # set inital values for posterior prob calculation
  denom <- model["l0"] * dbeta(p, 1, 1) * p
  numer <- 0

  # iterate through model componetnts
  non_unif <- sum(grepl("^l", names(model))) - 1
  for (i in 1:non_unif) {
    l <- paste0("l", i)
    r <- paste0("r", i)
    s <- paste0("s", i)
    p_comp <- unname(model[l] * dbeta(p, model[r], model[s]) * p)
    if (model[r] > model[s]) {
      denom <- denom + p_comp
    } else {
      numer <- numer + p_comp
    }
  }

  # calculate the odds, probability
  odds <- numer / denom
  p_alt <- odds / (1 + odds)
  return(p_alt)
}
