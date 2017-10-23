#' @title Calculate the value of the log likelihood function
#' @description Given a distribution of p-values and mixture model coefficients
#' calculate the value of the log likelihood function
#' @param p A numeric vector of p-values
#' @param model The beta mixture model represented as a named list of distribution parameters.
loglike <- function(p, coefs) {

  # get model parameters
  lambda <- coefs[grepl("^l", names(coefs))]
  r <- coefs[grepl("^r", names(coefs))]
  s <- coefs[grepl("^s", names(coefs))]

  # change zeros in lambda and normalize
  lnorm <- lambda / sum(lambda)

  # calculate the log-likelihood
  n_nonunif <- length(lambda) - 1
  if (n_nonunif > 0) {
    ll <- rowSums(sapply(1:n_nonunif, function (i) {
      lnorm[i + 1] * dbeta(p, r[i], s[i])
    }))
    ll[is.infinite(ll)] <- .Machine$double.xmax
    ll <- sum(log(lnorm[1] + ll))
  } else {
    ll <- log(lnorm[1]) * length(p)
  }
  return(ll)
}
