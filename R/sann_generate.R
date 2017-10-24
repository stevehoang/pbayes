
sann_generate <- function(...) {
  params <- list(...)
  params <- params[-length(params)]
  n_comp <- ( length(params) - 1 ) / 3
  names(params) <- c("l0", paste0("l", 1:n_comp),
                     paste0("r", 1:n_comp), paste0("s", 1:n_comp))
  indicator <- sample(c(0,1), 1)
  frac <- runif(1)
  ind <- sample(1:length(params), 1)
  if (grepl("^s", names(params)[ind])) {
    val <- 1 + 99 * frac
  } else {
    val <- frac
  }
  # print(p aste(names(params), val))
  print(val)
  print (unlist(params))
  params[[ind]] <- val
  # print (unlist(params))
  norm <- params[grepl("^l", names(params))] %>%
    unlist() %>%
    sum()
  params[grepl("^l", names(params))] <-
    lapply(params[grepl("^l", names(params))], function(x) x / norm)
  res <- unlist(params)
  print(res)
  return(res)
}
