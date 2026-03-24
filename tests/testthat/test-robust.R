
test_that("robust estimation works", {
  
  # Generate test data
  set.seed(123)
  data <- simulate_meta_data(k = 30)
  
  # Fit robust model
  result <- meta_robust(data$yi, data$vi)
  
  # Basic checks
  expect_true(is.list(result))
  expect_true(!is.na(result$estimate))
  expect_true(result$se > 0)
  expect_true(result$converged)
})

test_that("contamination models work", {
  
  # Generate contaminated data
  set.seed(123)
  data <- simulate_meta_data(k = 30, contamination = 0.2)
  
  # Test different contamination models
  result_t <- meta_robust(data$yi, data$vi, contamination = "t-mixture")
  result_none <- meta_robust(data$yi, data$vi, contamination = "none")
  
  # Robust estimate should exist
  expect_true(!is.na(result_t$estimate))
  expect_true(!is.na(result_none$estimate))
})

