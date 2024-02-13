#' @title Calculate the value of the log likelihood function.
#' @description Given a distribution of p-values and parameters for a beta-uniform mixture,
#' calculate the value of the log likelihood function.
#' @param p A numeric vector of p-values
#' @param mixm An object of class betamix representing the uniform-beta mixture model.
#' @export
loglike_ubeta <- function(p, mixm) {

  # get model parameters
  lambda <- unlist(mixm[grepl("^l", names(mixm))])
  r <- unlist(mixm[grepl("^r", names(mixm))])
  s <- unlist(mixm[grepl("^s", names(mixm))])

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
