# tests/testthat/test-causal.R
# Tests for causal_forest_meta() and transport()

context("Causal inference helpers")

test_that("causal_forest_meta returns CATE, importance, and subgroups", {
  set.seed(42)
  k <- 40
  X <- matrix(rnorm(k * 3), k, 3, dimnames = list(NULL, c("age", "dose", "year")))
  yi <- 0.3 + 0.5 * X[, 1] - 0.2 * X[, 2] + rnorm(k, 0, 0.2)
  vi <- rep(0.04, k)

  res <- causal_forest_meta(yi, vi, X)
  expect_length(res$cate, k)
  expect_length(res$importance, 3)
  expect_true(all(res$importance >= 0))
  expect_true(abs(sum(res$importance) - 1) < 0.01)  # should sum to 1
  expect_true(length(res$subgroups) > 0)
})

test_that("causal_forest_meta subgroups have correct structure", {
  set.seed(42)
  k <- 30
  X <- matrix(rnorm(k * 2), k, 2, dimnames = list(NULL, c("mod1", "mod2")))
  yi <- 0.5 + 0.3 * X[, 1] + rnorm(k, 0, 0.15)
  vi <- rep(0.03, k)

  res <- causal_forest_meta(yi, vi, X, min_leaf = 3)
  for (sg in res$subgroups) {
    expect_true(is.data.frame(sg))
    expect_true(all(c("group", "estimate", "se", "k") %in% names(sg)))
    expect_equal(nrow(sg), 2)  # low and high
  }
})

test_that("transport produces valid predictions for multiple targets", {
  set.seed(42)
  k <- 30
  X <- matrix(rnorm(k * 2), k, 2)
  yi <- 0.5 + 0.3 * X[, 1] + rnorm(k, 0, 0.2)
  vi <- rep(0.04, k)
  target <- matrix(c(1, 0, -1, 1), nrow = 2, byrow = TRUE)

  res <- transport(yi, vi, X, target)
  expect_length(res$estimate, 2)
  expect_length(res$se, 2)
  expect_true(all(res$se > 0))
  expect_true(all(res$ci_lower < res$ci_upper))
})
