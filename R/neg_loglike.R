.neg_loglike <- function(x, ...) {

  # get model parameters
  args <- list(...)
  args <- args[sort(names(args))]
  lambda <- args[grepl("^l", names(args))] %>%
    do.call(c, .)
  r <- args[grepl("^r", names(args))] %>%
    do.call(c, .)
  s <- args[grepl("^s", names(args))] %>%
    do.call(c, .)

  lnorm <- lambda / sum(lambda)

  n_nonunif <- length(lambda) - 1
  if (n_nonunif > 0) {
    ll <- rowSums(sapply(1:n_nonunif, function (i) {
      lnorm[i + 1] * dbeta(x, r[i], s[i])
    }))
    ll <- ll + (x * lnorm[1])
    ll[is.infinite(ll)] <- .Machine$double.xmax
    ll <- sum(log(ll))
    ll <- ll * -1
  } else {
    ll <- sum(log(lnorm[1]) * x) * -1
  }
  return(ll)
}
