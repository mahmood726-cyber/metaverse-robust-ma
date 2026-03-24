
context("Complete metaverse package tests")

library(testthat)
library(metaverse)

# Test robust estimation
test_that("All robust estimators work correctly", {
  
  set.seed(123)
  data <- simulate_meta_data(k = 30, tau2 = 0.1)
  
  # Test all methods
  methods <- c("MM", "M", "tau", "S", "ML", "REML")
  
  for(method in methods) {
    result <- meta_robust(data$yi, data$vi, method = method)
    
    expect_s3_class(result, "meta_robust")
    expect_true(!is.na(result$estimate))
    expect_true(result$se > 0)
    expect_true(!is.null(result$converged))
    expect_true(result$converged)
  }
})

# Test contamination models
test_that("Contamination models handle outliers", {
  
  set.seed(123)
  data <- simulate_meta_data(k = 30, contamination = 0.2)
  
  contam_types <- c("t-mixture", "slash", "box-cox")
  
  for(type in contam_types) {
    result <- meta_robust(data$yi, data$vi, contamination = type)
    
    expect_true(!is.null(result$contamination_fit))
    expect_true(any(result$outlier_prob > 0.5))
  }
})

# Test variable selection methods
test_that("All selection methods work", {
  
  set.seed(123)
  data <- simulate_meta_data(k = 50, moderators = 5)
  
  X <- as.matrix(data[, grep("mod", names(data))])
  
  methods <- c("knockoff", "spike_slab", "lasso", "elastic_net")
  
  for(method in methods) {
    result <- select_moderators(X, data$yi, data$vi, method = method)
    
    expect_s3_class(result, "selection_result")
    expect_true(is.integer(result$selected))
    expect_true(!is.null(result$method))
  }
})

# Test inference methods
test_that("Inference methods provide valid intervals", {

  set.seed(123)
  data <- simulate_meta_data(k = 30)

  # Create model object with explicit method
  rob <- meta_robust(data$yi, data$vi, method = "MM")
  model <- new("metaverse",
               data = data,
               model = list(method = "MM"),
               robust = unclass(rob))

  # Test bootstrap (reduced R for speed)
  boot_result <- quantify_uncertainty(model, method = "bootstrap", R = 50)

  expect_true(!is.null(boot_result$intervals))
  expect_true(all(!is.na(boot_result$se)))
})

# Test sensitivity analyses
test_that("Sensitivity analyses run without errors", {

  set.seed(123)
  data <- simulate_meta_data(k = 20)

  rob <- meta_robust(data$yi, data$vi, method = "MM")

  model <- new("metaverse",
               data = data,
               model = list(method = "MM"),
               robust = rob)

  # Test outlier analysis (lightweight)
  sens_result <- assess_sensitivity(model, analyses = c("outliers"))

  expect_s3_class(sens_result, "sensitivity_result")
  expect_true("outliers" %in% names(sens_result))
})

# Performance tests
test_that("Package handles large datasets efficiently", {
  
  set.seed(123)
  data_large <- simulate_meta_data(k = 1000)
  
  time_taken <- system.time({
    result <- meta_robust(data_large$yi, data_large$vi)
  })
  
  expect_true(time_taken["elapsed"] < 5)  # Should complete in < 5 seconds
})

# Edge cases
test_that("Package handles edge cases gracefully", {

  # Small k -- should still produce a result
  data_small <- simulate_meta_data(k = 5, seed = 42)
  result_small <- meta_robust(data_small$yi, data_small$vi)
  expect_true(!is.na(result_small$estimate))

  # Perfect homogeneity -- tau estimator should give tau2 near 0
  data_homog <- data.frame(yi = rep(0.5, 20), vi = rep(0.01, 20))
  result_homog <- meta_robust(data_homog$yi, data_homog$vi, method = "tau")
  expect_true(!is.null(result_homog$tau2))
  expect_true(result_homog$tau2 < 0.01)

  # Mismatched lengths should error
  expect_error(meta_robust(c(1, 2, 3), c(0.1, 0.2)))

  # Negative variance should error
  expect_error(meta_robust(c(1, 2, 3), c(0.1, -0.1, 0.1)))
})

