construct_nll <- function(arg_names) {

  # define negative log-likelihood function
  # formals will be reassigned after the definition
  neg_loglike <- function() {

    args <- as.list(match.call())
    p <- args[grepl("^p", names(args))] %>%
      do.call(c, .)
    lambda <- args[grepl("^l", names(args))] %>%
      do.call(c, .)
    r <- args[grepl("^r", names(args))] %>%
      do.call(c, .)
    s <- args[grepl("^s", names(args))] %>%
      do.call(c, .)

    # change zeros in lambda and normalize
    lnorm <- lambda / sum(lambda)

    n_nonunif <- length(lambda) - 1
    if (n_nonunif > 0) {
      ll <- rowSums(sapply(1:n_nonunif, function (i) {
        lnorm[i + 1] * dbeta(p, r[i], s[i])
      }))
      ll[is.infinite(ll)] <- .Machine$double.xmax
      ll <- sum(log(lnorm[1] + ll))
      nll <- ll * -1
    } else {
      nll <- log(lnorm[1]) * length(p) * -1
    }
    return(nll)
  }

  # reassign formals
  f <- neg_loglike
  new_forms <- paste(arg_names, collapse = " = ,")
  new_forms <- paste0("alist(p = , ", new_forms, " = )")
  formals(f) <- eval(parse(text = new_forms))
  return(f)
}
