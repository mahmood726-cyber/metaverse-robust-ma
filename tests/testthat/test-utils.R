# tests/testthat/test-utils.R
# Tests for utility functions: escalc_robust, create_network_data,
# meta_regression_power, extract_2x2_table

context("Utility functions")

# ============================================================
# 1. escalc_robust: OR
# ============================================================
test_that("escalc_robust computes log-OR correctly", {
  res <- escalc_robust(measure = "OR", ai = 10, bi = 40, ci = 5, di = 45)
  expected_logOR <- log((10 * 45) / (40 * 5))
  expected_vi <- 1/10 + 1/40 + 1/5 + 1/45
  expect_equal(res$yi, expected_logOR, tolerance = 1e-10)
  expect_equal(res$vi, expected_vi, tolerance = 1e-10)
})

test_that("escalc_robust applies continuity correction for zero cells", {
  res <- escalc_robust(measure = "OR", ai = 0, bi = 50, ci = 5, di = 45)
  # With add = 0.5 default, ai becomes 0.5
  expected_logOR <- log((0.5 * 45) / (50 * 5))
  expect_true(is.finite(res$yi))
  expect_true(is.finite(res$vi))
})

# ============================================================
# 2. escalc_robust: RR
# ============================================================
test_that("escalc_robust computes log-RR correctly", {
  res <- escalc_robust(measure = "RR", ai = 15, ci = 10, n1i = 100, n2i = 100)
  expected <- log((15/100) / (10/100))
  expect_equal(res$yi, expected, tolerance = 1e-10)
})

# ============================================================
# 3. escalc_robust: ZCOR
# ============================================================
test_that("escalc_robust computes Fisher z correctly", {
  res <- escalc_robust(measure = "ZCOR", ri = 0.5, ni = 50)
  expected_z <- 0.5 * log((1 + 0.5) / (1 - 0.5))
  expected_vi <- 1 / (50 - 3)
  expect_equal(res$yi, expected_z, tolerance = 1e-10)
  expect_equal(res$vi, expected_vi, tolerance = 1e-10)
})

# ============================================================
# 4. extract_2x2_table
# ============================================================
test_that("extract_2x2_table computes marginals and rates", {
  d <- data.frame(a = c(10, 20), b = c(40, 30), c = c(5, 15), d = c(45, 35))
  res <- extract_2x2_table(d, a, b, c, d)
  expect_equal(res$n1i, c(50, 50))
  expect_equal(res$n2i, c(50, 50))
  expect_equal(res$pi1, c(10/50, 20/50))
  expect_equal(res$pi2, c(5/50, 15/50))
})

# ============================================================
# 5. create_network_data
# ============================================================
test_that("create_network_data produces valid edge list", {
  d <- data.frame(
    study = c("S1", "S1", "S2", "S2", "S2"),
    treat = c("A", "B", "A", "B", "C"),
    outcome = c(0.3, 0.5, 0.4, 0.6, 0.2),
    se = c(0.1, 0.1, 0.15, 0.15, 0.12)
  )
  net <- create_network_data(d, "study", "treat", "outcome", "se")
  expect_true(is.data.frame(net))
  expect_true(nrow(net) > 0)
  expect_true(all(c("treat1", "treat2", "y_diff") %in% names(net)))
  # S1 has 1 comparison (A-B), S2 has 3 (A-B, A-C, B-C) = 4 total
  expect_equal(nrow(net), 4)
})

# ============================================================
# 6. meta_regression_power
# ============================================================
test_that("meta_regression_power computes power and required k", {
  res <- meta_regression_power(k = 50, R2 = 0.1, n_covariates = 3)
  expect_true(!is.null(res$power))
  expect_true(res$power > 0 && res$power <= 1)
  expect_true(!is.null(res$k_required))
})

test_that("meta_regression_power errors with too few studies", {
  expect_error(
    meta_regression_power(k = 3, R2 = 0.1, n_covariates = 5),
    "Not enough"
  )
})

# ============================================================
# 7. working_correlation_matrix
# ============================================================
test_that("working_correlation_matrix CS is symmetric and valid", {
  d <- data.frame(a = 1:5)
  R <- working_correlation_matrix(d, type = "CS", rho = 0.3)
  expect_equal(nrow(R), 5)
  expect_equal(ncol(R), 5)
  expect_true(isSymmetric(R))
  expect_equal(diag(R), rep(1, 5))
  expect_equal(R[1, 2], 0.3)
})

test_that("working_correlation_matrix AR1 decays correctly", {
  d <- data.frame(a = 1:4)
  R <- working_correlation_matrix(d, type = "AR1", rho = 0.5)
  expect_equal(R[1, 1], 1)
  expect_equal(R[1, 2], 0.5)
  expect_equal(R[1, 3], 0.25, tolerance = 1e-10)
  expect_equal(R[1, 4], 0.125, tolerance = 1e-10)
})
