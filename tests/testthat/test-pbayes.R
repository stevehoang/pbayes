library(testthat)
library(pbayes)

# Test for correct input
test_that("pbayes returns correct structure with valid input", {
  p <- runif(100, 0, 1) # Generate 100 random p-values
  result <- pbayes(p)
  expect_true(is.list(result))
  expect_true("pbayes" %in% class(result))
  expect_equal(length(result$p_value), 100)
  expect_equal(length(result$posterior_prob), 100)
  expect_true(is.numeric(result$p_value))
  expect_true(is.numeric(result$posterior_prob))
})

# Test for input validation
test_that("pbayes stops with non-numeric input", {
  expect_error(pbayes("not numeric"))
  expect_error(pbayes(c(NA, 0.1, 0.5)))
})

# Test for small sample size warning
test_that("pbayes warns on small sample size", {
  p <- runif(50, 0, 1) # Generate 50 random p-values, less than recommended 100
  expect_warning(pbayes(p))
})
