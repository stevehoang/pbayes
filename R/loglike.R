.loglike <- function(x, coefs) {
  ind <- grepl("^l", names(coefs))
  coefs[ind] <- coefs[ind] / sum(coefs[ind])
  ll_call <- paste(names(coefs), coefs, sep="=")
  ll_call <- paste(ll_call, collapse = ",")
  ll <- eval(parse(text = paste0("neg_loglike(x,", ll_call, ")"))) * -1
  return(ll)
}
