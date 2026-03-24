# tests/testthat/test-sensitivity.R
# Tests for assess_sensitivity() and all sub-analyses

context("assess_sensitivity: sensitivity analysis suite")

# ============================================================
# Helper
# ============================================================
make_model <- function(k = 20, seed = 42) {
  data <- simulate_meta_data(k = k, theta = 0.5, tau2 = 0.1, seed = seed)
  rob <- meta_robust(data$yi, data$vi, method = "MM")
  new("metaverse",
      data = data,
      model = list(method = "MM"),
      robust = unclass(rob))
}

# ============================================================
# 1. Outlier analysis
# ============================================================
test_that("outlier analysis returns standardized residuals", {
  model <- make_model(k = 20)
  res <- assess_sensitivity(model, analyses = "outliers")
  expect_true("outliers" %in% names(res))
  out <- res$outliers
  expect_true(!is.null(out$standardized_residuals))
  expect_length(out$standardized_residuals, 20)
  expect_true(!is.null(out$max_z))
  expect_true(out$max_z > 0)
})

# ============================================================
# 2. Publication bias
# ============================================================
test_that("publication bias analysis runs all tests", {
  model <- make_model(k = 25)
  res <- assess_sensitivity(model, analyses = "publication_bias")
  expect_true("publication_bias" %in% names(res))
  pb <- res$publication_bias

  # Egger tests
  expect_true(!is.null(pb$egger_standard))
  expect_true(!is.null(pb$egger_robust))
  expect_true(!is.null(pb$egger_standard$p_value))

  # Begg
  expect_true(!is.null(pb$begg))
  expect_true(!is.null(pb$begg$p_value))

  # Trim and fill
  expect_true(!is.null(pb$trim_fill_l0))
  expect_true(!is.null(pb$trim_fill_l0$k0))
  expect_true(pb$trim_fill_l0$k0 >= 0)

  # P-curve
  expect_true(!is.null(pb$p_curve))

  # Contour data
  expect_true(!is.null(pb$contour_data))
})

# ============================================================
# 3. Influence analysis
# ============================================================
test_that("influence analysis computes leave-one-out and Cook's D", {
  model <- make_model(k = 15)
  res <- assess_sensitivity(model, analyses = "influence")
  expect_true("influence" %in% names(res))
  inf <- res$influence

  expect_true(!is.null(inf$loo))
  expect_equal(nrow(inf$loo), 15)
  expect_equal(ncol(inf$loo), 5)

  expect_true(!is.null(inf$cooks_distance))
  expect_length(inf$cooks_distance, 15)

  expect_true(!is.null(inf$dffits))
  expect_length(inf$dffits, 15)

  expect_true(!is.null(inf$covratio))
  expect_length(inf$covratio, 15)
})

# ============================================================
# 4. Fragility
# ============================================================
test_that("fragility analysis returns index and quotient", {
  model <- make_model(k = 15)
  res <- assess_sensitivity(model, analyses = "fragility")
  expect_true("fragility" %in% names(res))
  frag <- res$fragility

  expect_true(!is.null(frag$fragility_index))
  expect_true(!is.null(frag$fragility_quotient))
  expect_true(!is.null(frag$interpretation))
  expect_true(nchar(frag$interpretation) > 0)
})

# ============================================================
# 5. E-values
# ============================================================
test_that("e-values analysis returns all components", {
  model <- make_model(k = 20)
  res <- assess_sensitivity(model, analyses = "evalues")
  expect_true("evalues" %in% names(res))
  ev <- res$evalues

  expect_true(!is.null(ev$evalue_estimate))
  expect_true(ev$evalue_estimate >= 1)
  expect_true(!is.null(ev$evalue_ci))
  expect_true(!is.null(ev$interpretation))
  expect_true(!is.null(ev$bias_scenarios))
})

# ============================================================
# 6. Multiple analyses at once
# ============================================================
test_that("running multiple analyses works", {
  model <- make_model(k = 15)
  res <- assess_sensitivity(model, analyses = c("outliers", "evalues"))
  expect_true("outliers" %in% names(res))
  expect_true("evalues" %in% names(res))
})

# ============================================================
# 7. "all" keyword
# ============================================================
test_that("analyses = 'all' runs cumulative too", {
  model <- make_model(k = 10)
  res <- assess_sensitivity(model, analyses = "all")
  expect_true("cumulative" %in% names(res))
  cum <- res$cumulative
  expect_true(is.data.frame(cum))
  expect_equal(nrow(cum), 10)
  expect_true(all(c("estimate", "se", "ci_lower", "ci_upper") %in% names(cum)))
})

# ============================================================
# 8. sensitivity_to_heterogeneity standalone function
# ============================================================
test_that("sensitivity_to_heterogeneity returns a data frame and critical tau2", {
  set.seed(42)
  d <- simulate_meta_data(k = 20, theta = 0.3, tau2 = 0.05, seed = 42)
  res <- sensitivity_to_heterogeneity(d$yi, d$vi, tau2_range = seq(0, 0.3, 0.05))
  expect_true(is.data.frame(res$results))
  expect_equal(nrow(res$results), length(seq(0, 0.3, 0.05)))
})

# ============================================================
# 9. Edge case: k = 5 for sensitivity
# ============================================================
test_that("sensitivity runs on small k = 5", {
  model <- make_model(k = 5)
  # Should not error on small samples
  res <- assess_sensitivity(model, analyses = "outliers")
  expect_true(!is.null(res$outliers))
})
