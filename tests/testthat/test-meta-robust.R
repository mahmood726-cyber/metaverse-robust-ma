# tests/testthat/test-meta-robust.R
# Comprehensive tests for meta_robust() and estimation methods

context("meta_robust: estimation methods")

# ============================================================
# Helper: generate a clean dataset for reuse
# ============================================================
make_data <- function(k = 30, theta = 0.5, tau2 = 0.1, seed = 42, ...) {
  simulate_meta_data(k = k, theta = theta, tau2 = tau2, seed = seed, ...)
}

# ============================================================
# 1. Basic functionality for every method
# ============================================================
test_that("meta_robust runs for all 6 methods without error", {
  data <- make_data()
  methods <- c("MM", "M", "tau", "S", "ML", "REML")
  for (m in methods) {
    res <- meta_robust(data$yi, data$vi, method = m)
    expect_s3_class(res, "meta_robust", info = paste("method =", m))
    expect_true(is.numeric(res$estimate), info = m)
    expect_true(!is.na(res$estimate), info = m)
    expect_true(res$se > 0, info = m)
    expect_true(!is.null(res$ci.lower), info = m)
    expect_true(!is.null(res$ci.upper), info = m)
    expect_true(res$ci.lower < res$ci.upper, info = m)
    expect_true(!is.null(res$p_value), info = m)
    expect_true(res$p_value >= 0 && res$p_value <= 1, info = m)
  }
})

test_that("M-estimator converges on clean data", {
  data <- make_data(k = 50)
  res <- meta_robust(data$yi, data$vi, method = "M")
  expect_true(res$converged)
  # Should be reasonably close to true theta = 0.5
  expect_true(abs(res$estimate - 0.5) < 0.5)
})

test_that("MM-estimator produces both S and M stage results", {
  data <- make_data(k = 30)
  res <- meta_robust(data$yi, data$vi, method = "MM")
  expect_true(!is.null(res$s_estimate))
  expect_true(!is.null(res$s_scale))
  expect_true(res$converged)
})

test_that("S-estimator returns a scale estimate", {
  data <- make_data(k = 20)
  res <- s_estimator(data$yi, data$vi, 1/data$vi)
  expect_true(!is.na(res$scale))
  expect_true(res$scale > 0)
})

test_that("tau estimator returns heterogeneity stats", {
  data <- make_data(k = 30, tau2 = 0.15)
  res <- meta_robust(data$yi, data$vi, method = "tau")
  expect_true(!is.null(res$tau2))
  expect_true(res$tau2 >= 0)
  expect_true(!is.null(res$I2))
  expect_true(res$I2 >= 0 && res$I2 <= 100)
  expect_true(!is.null(res$H2))
  expect_true(!is.null(res$prediction_interval))
  expect_length(res$prediction_interval, 2)
  expect_true(res$prediction_interval[1] < res$prediction_interval[2])
})

test_that("REML and ML estimators return results", {
  data <- make_data(k = 25)
  for (m in c("ML", "REML")) {
    res <- meta_robust(data$yi, data$vi, method = m)
    expect_true(!is.na(res$estimate), info = m)
    expect_true(res$se > 0, info = m)
  }
})

# ============================================================
# 2. Contamination models
# ============================================================
test_that("t-mixture contamination model works", {
  data <- make_data(k = 30, contamination = 0.2)
  res <- meta_robust(data$yi, data$vi, contamination = "t-mixture")
  expect_true(!is.null(res$contamination_fit))
  expect_true(!is.null(res$outlier_prob))
  expect_length(res$outlier_prob, nrow(data))
  expect_true(all(res$outlier_prob >= 0 & res$outlier_prob <= 1))
})

test_that("slash contamination model works", {
  data <- make_data(k = 30, contamination = 0.1)
  res <- meta_robust(data$yi, data$vi, contamination = "slash")
  expect_true(!is.null(res$contamination_fit))
  expect_true(!is.null(res$outlier_prob))
})

test_that("box-cox contamination model works", {
  # box-cox requires positive data
  data <- make_data(k = 30, theta = 2.0, tau2 = 0.05)
  res <- meta_robust(data$yi, data$vi, contamination = "box-cox")
  expect_true(!is.null(res$contamination_fit))
  expect_true(!is.null(res$contamination_fit$lambda))
})

test_that("contamination = 'none' does not alter estimates", {
  data <- make_data(k = 25)
  res_none <- meta_robust(data$yi, data$vi, contamination = "none")
  res_def  <- meta_robust(data$yi, data$vi)
  expect_equal(res_none$estimate, res_def$estimate)
  expect_equal(res_none$se, res_def$se)
})

# ============================================================
# 3. Moderators (meta-regression)
# ============================================================
test_that("meta_robust works with single moderator", {
  data <- make_data(k = 40, moderators = 1)
  mods <- as.matrix(data[, grep("mod", names(data))])
  res <- meta_robust(data$yi, data$vi, mods = mods, method = "M")
  expect_true(!is.na(res$estimate))
  expect_length(res$coefficients, 2)  # intercept + 1 moderator
})

test_that("meta_robust works with multiple moderators", {
  data <- make_data(k = 50, moderators = 3)
  mods <- as.matrix(data[, grep("mod", names(data))])
  res <- meta_robust(data$yi, data$vi, mods = mods, method = "M")
  expect_true(!is.na(res$estimate))
  expect_length(res$coefficients, 4)  # intercept + 3 moderators
})

# ============================================================
# 4. Custom weights and control parameters
# ============================================================
test_that("custom weights are accepted", {
  data <- make_data(k = 20)
  custom_w <- rep(1, nrow(data))
  res <- meta_robust(data$yi, data$vi, weights = custom_w)
  expect_true(!is.na(res$estimate))
})

test_that("control parameters override defaults", {
  data <- make_data(k = 20)
  res <- meta_robust(data$yi, data$vi, method = "M",
                     control = list(tol = 1e-3, maxiter = 5))
  # Should still produce a result (may not converge in 5 iter)
  expect_true(!is.na(res$estimate))
})

# ============================================================
# 5. Input validation
# ============================================================
test_that("mismatched lengths error", {
  expect_error(meta_robust(c(1, 2, 3), c(0.1, 0.2)),
               "same length")
})

test_that("negative variances error", {
  expect_error(meta_robust(c(1, 2, 3), c(0.1, -0.1, 0.1)),
               "positive")
})

test_that("zero variances error", {
  expect_error(meta_robust(c(1, 2, 3), c(0.1, 0, 0.1)),
               "positive")
})

test_that("mismatched moderator rows error", {
  expect_error(
    meta_robust(c(1, 2, 3), c(0.1, 0.2, 0.3), mods = matrix(1:4, 2, 2)),
    "must equal"
  )
})

test_that("negative weights error", {
  expect_error(
    meta_robust(c(1, 2, 3), c(0.1, 0.2, 0.3), weights = c(1, -1, 1)),
    "non-negative"
  )
})

# ============================================================
# 6. Edge cases
# ============================================================
test_that("k = 2 works without error", {
  yi <- c(0.3, 0.7)
  vi <- c(0.05, 0.08)
  res <- meta_robust(yi, vi, method = "M")
  expect_true(!is.na(res$estimate))
  expect_true(res$se > 0)
})

test_that("k = 3 works for MM estimator", {
  yi <- c(0.2, 0.4, 0.6)
  vi <- c(0.01, 0.02, 0.03)
  res <- meta_robust(yi, vi, method = "MM")
  expect_true(!is.na(res$estimate))
})

test_that("homogeneous data (all same yi) is handled", {
  yi <- rep(0.5, 10)
  vi <- rep(0.01, 10)
  res <- meta_robust(yi, vi, method = "M")
  expect_true(abs(res$estimate - 0.5) < 0.01)
})

test_that("highly heterogeneous data does not crash", {
  set.seed(99)
  yi <- rnorm(20, 0, 5)
  vi <- runif(20, 0.01, 0.5)
  res <- meta_robust(yi, vi, method = "MM")
  expect_true(!is.na(res$estimate))
})

test_that("large dataset (k=500) completes quickly", {
  data <- make_data(k = 500, seed = 7)
  timing <- system.time(res <- meta_robust(data$yi, data$vi, method = "M"))
  expect_true(timing["elapsed"] < 10)
  expect_true(!is.na(res$estimate))
})

# ============================================================
# 7. Print method
# ============================================================
test_that("print.meta_robust runs without error", {
  data <- make_data(k = 15)
  res <- meta_robust(data$yi, data$vi, method = "tau")
  expect_output(print(res), "Robust Meta-Analysis")
  expect_output(print(res), "Estimate:")
})

# ============================================================
# 8. Robustness to outliers
# ============================================================
test_that("MM-estimator resists a single extreme outlier", {
  set.seed(11)
  yi <- c(rnorm(19, 0.5, 0.1), 50)  # one extreme outlier

  vi <- rep(0.02, 20)
  res_mm <- meta_robust(yi, vi, method = "MM")
  # Should be close to 0.5, not pulled toward 50
  expect_true(abs(res_mm$estimate - 0.5) < 1.0)
  # Compare with naive weighted mean
  naive <- sum((1/vi) * yi) / sum(1/vi)
  expect_true(abs(res_mm$estimate - 0.5) < abs(naive - 0.5))
})
