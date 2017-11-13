#' @title Generate a new point for simulated annealing
#' @param ... The arguments to the negative log likelihood function
sann_generate <- function(...) {
  params <- list(...)
  params <- params[-length(params)]
  n_comp <- ( length(params) - 1 ) / 3
  names(params) <- c("l0", paste0("l", 1:n_comp),
                     paste0("r", 1:n_comp), paste0("s", 1:n_comp))
  ind <- sample(1:length(params), 1) # select a parameter to change
  frac <- runif(1) # the change factor
  if (grepl("^s", names(params)[ind])) {
    val <- 1 + 99 * frac
  } else {
    val <- frac
  }
  params[[ind]] <- val
  norm <- params[grepl("^l", names(params))] %>%
    unlist() %>%
    sum()
  params[grepl("^l", names(params))] <-
    lapply(params[grepl("^l", names(params))], function(x) x / norm)
  res <- unlist(params)
  return(res)
}
