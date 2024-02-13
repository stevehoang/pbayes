#' @title Sample from a uniform-beta mixture distribution
#' @description Random draws from a uniform-beta mixture model
#' @param mixm A betamix object defining the beta mixture
#' @param n Number of draws
#' @export
sample_betamix <- function(mixm, n) {

  # check for valid mixture model object
  if (class(mixm) != "betamix") {
    stop("Mixture model must be of class betamix")
  }

  # get mixing fractions of components
  fracs <- sort(unlist(mixm[grepl("^l", names(mixm))]))
  intervals <- c(0, cumsum(fracs))

  # use the intervals to segment a random uniform vector
  u <- runif(n)
  cuts <- cut(u, intervals)
  intervals2 <- levels(cuts)
  names(intervals2) <- names(intervals)[-1]

  # loop over the intervals and sample the appropriate component
  samps <- c()
  for(i in 1:length(intervals2)) {
    nsamp <- sum(cuts == intervals2[i])
    comp <- names(intervals2)[i]
    if (comp == "l0") {
      s <- rbeta(nsamp, 1, 1)
    } else {
      s <- rbeta(nsamp, mixm[[sub("l", "r", comp)]],
                 mixm[[sub("l", "s", comp)]])
    }
    samps <- c(samps, s)
  }

  return(sample(samps))
}
