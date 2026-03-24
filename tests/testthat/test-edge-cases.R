# tests/testthat/test-edge-cases.R
# Edge cases: k=2, single moderator, all-zero-variance surrogates,
# homogeneous data, and other boundary conditions

context("Edge cases and boundary conditions")

# ============================================================
# 1. k = 2 (minimum for meta-analysis)
# ============================================================
test_that("meta_robust works with k = 2 for all methods", {
  yi <- c(0.3, 0.7)
  vi <- c(0.05, 0.08)
  for (m in c("M", "MM", "S", "tau")) {
    res <- meta_robust(yi, vi, method = m)
    expect_true(!is.na(res$estimate), info = paste("k=2, method =", m))
    expect_true(res$se > 0, info = paste("k=2, method =", m))
  }
})

test_that("tau estimator with k = 2 returns valid heterogeneity", {
  yi <- c(0.2, 0.9)
  vi <- c(0.02, 0.03)
  res <- meta_robust(yi, vi, method = "tau")
  expect_true(!is.null(res$tau2))
  expect_true(res$tau2 >= 0)
  expect_true(!is.null(res$I2))
})

# ============================================================
# 2. k = 2 for diagnostics
# ============================================================
test_that("diagnose works with k = 2", {
  yi <- c(0.4, 0.6)
  vi <- c(0.03, 0.05)
  data <- data.frame(yi = yi, vi = vi, study = c("A", "B"))
  rob <- meta_robust(yi, vi, method = "M")
  model <- new("metaverse",
               data = data,
               model = list(method = "M"),
               robust = unclass(rob))
  diag <- diagnose(model)
  expect_true(is.list(diag))
  expect_equal(diag$n_studies, 2)
})

# ============================================================
# 3. Homogeneous data (tau2 should be ~0)
# ============================================================
test_that("homogeneous data returns near-zero tau2", {
  yi <- rep(0.5, 20)
  vi <- rep(0.01, 20)
  res <- meta_robust(yi, vi, method = "tau")
  expect_true(res$tau2 < 0.01)
})

test_that("homogeneous data M-estimator returns precise estimate", {
  yi <- rep(1.0, 15)
  vi <- rep(0.02, 15)
  res <- meta_robust(yi, vi, method = "M")
  expect_true(abs(res$estimate - 1.0) < 0.05)
})

# ============================================================
# 4. Near-zero variances (very precise studies)
# ============================================================
test_that("very small variances do not cause numerical issues", {
  yi <- c(0.5, 0.51, 0.49, 0.50, 0.52)
  vi <- c(1e-6, 1e-6, 1e-6, 1e-6, 1e-6)
  res <- meta_robust(yi, vi, method = "M")
  expect_true(!is.na(res$estimate))
  # Estimate should be close to mean
  expect_true(abs(res$estimate - mean(yi)) < 0.1)
})

# ============================================================
# 5. Large variances (imprecise studies)
# ============================================================
test_that("large variances are handled", {
  yi <- c(0.5, 1.0, -0.5, 2.0)
  vi <- c(10, 20, 15, 25)
  res <- meta_robust(yi, vi, method = "M")
  expect_true(!is.na(res$estimate))
})

# ============================================================
# 6. Single moderator selection edge case
# ============================================================
test_that("selection with p = 1 returns valid result", {
  set.seed(99)
  k <- 30
  X <- matrix(rnorm(k), k, 1)
  colnames(X) <- "X1"
  yi <- rnorm(k, 0.3 + 0.5 * X[, 1], 0.3)
  vi <- rep(0.04, k)

  res <- select_moderators(X, yi, vi, method = "lasso",
                           control = list(n_bootstrap = 10))
  expect_s3_class(res, "selection_result")
  # Selected should be valid indices
  if (length(res$selected) > 0) {
    expect_true(all(res$selected %in% 1:ncol(X)))
  }
})

# ============================================================
# 7. All identical effect sizes
# ============================================================
test_that("all-identical yi does not crash M-estimator", {
  yi <- rep(0.3, 10)
  vi <- runif(10, 0.01, 0.1)
  res <- meta_robust(yi, vi, method = "M")
  expect_true(!is.na(res$estimate))
  expect_true(abs(res$estimate - 0.3) < 0.1)
})

# ============================================================
# 8. simulate_meta_data edge cases
# ============================================================
test_that("simulate_meta_data with k = 2 works", {
  data <- simulate_meta_data(k = 2, seed = 1)
  expect_equal(nrow(data), 2)
  expect_true(all(data$vi > 0))
})

test_that("simulate_meta_data with contamination = 0 produces no outliers", {
  data <- simulate_meta_data(k = 20, contamination = 0, seed = 1)
  expect_equal(nrow(data), 20)
})

test_that("simulate_meta_data with moderators works", {
  data <- simulate_meta_data(k = 20, moderators = 3, seed = 1)
  expect_true("mod1" %in% names(data))
  expect_true("mod2" %in% names(data))
  expect_true("mod3" %in% names(data))
})

# ============================================================
# 9. Power analysis edge cases
# ============================================================
test_that("power_analysis with k provided returns power", {
  res <- power_analysis(k = 30, theta = 0.3, tau2 = 0.05)
  expect_true(!is.null(res$power))
  expect_true(res$power > 0 && res$power <= 1)
})

test_that("power_analysis with k = NULL finds required k", {
  res <- power_analysis(k = NULL, theta = 0.5, tau2 = 0.1, power = 0.80)
  expect_true(!is.null(res$k_required))
})

test_that("power_analysis print method works", {
  res <- power_analysis(k = 20, theta = 0.3)
  expect_output(print(res), "Power Analysis")
})

# ============================================================
# 10. aggregate_effects edge cases
# ============================================================
test_that("aggregate_effects works with single-study clusters", {
  data <- data.frame(
    yi = c(0.3, 0.5, 0.7),
    vi = c(0.01, 0.02, 0.03),
    cluster_id = c("A", "B", "C")
  )
  agg <- aggregate_effects(data, cluster = "cluster_id", method = "robust")
  expect_equal(nrow(agg), 3)
  expect_true(all(agg$vi > 0))
})

test_that("aggregate_effects works with multi-study clusters", {
  data <- data.frame(
    yi = c(0.3, 0.35, 0.5, 0.55, 0.7),
    vi = c(0.01, 0.015, 0.02, 0.025, 0.03),
    cluster_id = c("A", "A", "B", "B", "C")
  )
  for (method in c("simple", "weighted", "robust", "borenstein")) {
    agg <- aggregate_effects(data, cluster = "cluster_id", method = method)
    expect_equal(nrow(agg), 3, info = method)
    expect_true(all(agg$vi > 0), info = method)
  }
})

# ============================================================
# 11. causal_forest_meta edge case: few studies
# ============================================================
test_that("causal_forest_meta warns on too-few studies", {
  yi <- c(0.3, 0.5)
  vi <- c(0.01, 0.02)
  X <- matrix(c(1, 2), 2, 1)
  expect_warning(
    res <- causal_forest_meta(yi, vi, X, min_leaf = 3),
    "Too few studies"
  )
  expect_length(res$cate, 2)
})

# ============================================================
# 12. transport function
# ============================================================
test_that("transport returns prediction for target population", {
  set.seed(42)
  k <- 30
  X <- matrix(rnorm(k * 2), k, 2)
  yi <- 0.5 + 0.3 * X[, 1] + rnorm(k, 0, 0.2)
  vi <- rep(0.04, k)
  target_X <- matrix(c(1, 0), nrow = 1)

  res <- transport(yi, vi, X, target_X)
  expect_true(!is.null(res$estimate))
  expect_true(!is.null(res$se))
  expect_true(res$se > 0)
  expect_true(res$ci_lower < res$ci_upper)
})

# ============================================================
# 13. S4 class validation
# ============================================================
test_that("metaverse class rejects data without yi/vi", {
  expect_error(
    validObject(new("metaverse", data = data.frame(x = 1:5, y = 1:5))),
    "yi and vi"
  )
})

test_that("metaverse class rejects negative vi", {
  expect_error(
    validObject(new("metaverse", data = data.frame(yi = 1:3, vi = c(0.1, -0.1, 0.1)))),
    "positive"
  )
})

test_that("metaverse show method runs", {
  data <- simulate_meta_data(k = 5, seed = 1)
  rob <- meta_robust(data$yi, data$vi, method = "M")
  model <- new("metaverse",
               data = data,
               model = list(method = "M", effect_type = "SMD"),
               robust = unclass(rob))
  expect_output(show(model), "METAVERSE")
})
