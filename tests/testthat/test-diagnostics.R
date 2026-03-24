# tests/testthat/test-diagnostics.R
# Tests for diagnose() function

context("diagnose: model diagnostics")

make_diag_model <- function(k = 20, seed = 42) {
  data <- simulate_meta_data(k = k, theta = 0.4, tau2 = 0.1, seed = seed)
  rob <- meta_robust(data$yi, data$vi, method = "M")
  rob$data <- data  # attach data for list-based access too
  model <- new("metaverse",
               data = data,
               model = list(method = "M"),
               robust = unclass(rob))
  model
}

test_that("diagnose returns all expected components", {
  model <- make_diag_model()
  d <- diagnose(model)

  expect_s3_class(d, "meta_diagnostics")
  expect_equal(d$n_studies, 20)
  expect_true(!is.na(d$estimate))
  expect_true(d$se > 0)

  # Normality test
  expect_true(!is.null(d$normality_test))
  expect_true(!is.null(d$normality_test$p_value))

  # Heterogeneity
  expect_true(!is.null(d$heterogeneity$Q))
  expect_true(d$heterogeneity$Q >= 0)
  expect_true(!is.null(d$heterogeneity$I2))

  # Outliers
  expect_true(!is.null(d$outliers))
})

test_that("diagnose print method works", {
  model <- make_diag_model()
  d <- diagnose(model)
  expect_output(print(d), "Diagnostics")
  expect_output(print(d), "Heterogeneity")
})

test_that("diagnose works on list-based meta_robust result", {
  data <- simulate_meta_data(k = 10, seed = 99)
  rob <- meta_robust(data$yi, data$vi, method = "M")
  rob$data <- data
  d <- diagnose(rob)
  expect_s3_class(d, "meta_diagnostics")
  expect_equal(d$n_studies, 10)
})

test_that("diagnose handles k = 3 (Shapiro-Wilk minimum)", {
  data <- simulate_meta_data(k = 3, seed = 77)
  rob <- meta_robust(data$yi, data$vi, method = "M")
  model <- new("metaverse",
               data = data,
               model = list(method = "M"),
               robust = unclass(rob))
  d <- diagnose(model)
  expect_true(!is.null(d$normality_test))
})
