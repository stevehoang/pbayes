library(testthat)
library(pbayes)

# Test that level_p returns numeric vector
test_that("level_p returns a numeric vector", {
  p <- runif(100, 0, 1) # Generate 100 random p-values
  expect_true(is.numeric(level_p(p)))
})

# Test input validation for non-numeric inputs
test_that("level_p stops with non-numeric input", {
  expect_error(level_p("not numeric"), "p is not numeric")
  expect_error(level_p(c(TRUE, FALSE)), "p is not numeric")
})

# Test the effect of leveling on a skewed distribution
test_that("level_p modifies skewed distribution of p-values", {
  p <- c(rep(0.01, 25), runif(75, 0.2, 1)) # Skewed distribution
  original_mean <- mean(p)
  corrected_p <- level_p(p)
  corrected_mean <- mean(corrected_p)
  expect_true(corrected_mean != original_mean)
})
