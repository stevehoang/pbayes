#' @title Levels the uniform component of a p-value distribution.
#' @description The uniform component of a p-distribution is leveled using a cubic
#' polynomial. The method is similar to that of \code{Polyfit::levelPValues}.
#' @param p A numeric vector of p-values
#' @export
level_p <- function(p) {

  # check for valid p-values
  if (!is.numeric(p)) {
    stop("p is not numeric.")
  }
  if (!(all(p >= 0) & all(p <= 1))) {
    stop("The input does not look like a vector of p-values. Some values are not in [0,1].")
  }

  # get center section from p-value distribution
  h <- hist(p, breaks = seq(0, 1, length.out = 100), plot = FALSE)
  ind <- (h$mids > 0.2) & (h$mids < 0.95) # this is the section we want to fit to
  x <- h$mids[ind]
  y <- h$counts[ind]

  # fit linear model
  m <- lm(y ~ poly(x, degree = 3, raw = TRUE))
  fitted <- coef(m)

  # auc function
  auc <- function(x_val) {
    area <- fitted[1] * x_val +
      (fitted[2] / 2) * x_val^2 +
      (fitted[3] / 3) * x_val^3 +
      (fitted[4] / 4) * x_val^4
    return(area)
  }

  # calculate leveled p-values
  total_auc <- auc(1)
  p_corrected <- sapply(p, function(x) auc(x) / total_auc)

  return(p_corrected)

}
