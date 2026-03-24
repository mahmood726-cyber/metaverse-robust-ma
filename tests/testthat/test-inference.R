# tests/testthat/test-inference.R
# Tests for quantify_uncertainty() and inference methods

context("quantify_uncertainty: inference methods")

make_inference_model <- function(k = 20, seed = 42) {
  data <- simulate_meta_data(k = k, theta = 0.4, tau2 = 0.08, seed = seed)
  rob <- meta_robust(data$yi, data$vi, method = "MM")
  new("metaverse",
      data = data,
      model = list(method = "MM"),
      robust = unclass(rob))
}

# ============================================================
# 1. Conformal inference
# ============================================================
test_that("conformal inference returns intervals with coverage", {
  model <- make_inference_model(k = 30)
  res <- quantify_uncertainty(model, method = "conformal", conf.level = 0.90)
  expect_true(!is.null(res$q_hat))
  expect_true(!is.null(res$intervals))
  expect_true(is.data.frame(res$intervals))
  expect_true(!is.null(res$coverage))
  expect_equal(res$conf.level, 0.90)
})

# ============================================================
# 2. Permutation inference
# ============================================================
test_that("permutation inference returns p-value and CI", {
  model <- make_inference_model(k = 15)
  res <- quantify_uncertainty(model, method = "permutation", R = 50)
  expect_true(!is.null(res$p_value))
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  expect_true(!is.null(res$ci))
  expect_length(res$ci, 2)
  expect_true(res$ci[1] < res$ci[2])
})

# ============================================================
# 3. Bayesian inference
# ============================================================
test_that("bayesian inference returns posterior samples", {
  model <- make_inference_model(k = 15)
  res <- quantify_uncertainty(model, method = "bayesian",
                              control = list(n_iter = 500, burn_in = 100))
  expect_true(!is.null(res$mu_posterior))
  expect_true(length(res$mu_posterior) > 0)
  expect_true(!is.null(res$tau2_posterior))
  expect_true(!is.null(res$estimate))
  expect_true(!is.null(res$hpd_mu))
  expect_length(res$hpd_mu, 2)
})

# ============================================================
# 4. Selective inference (no selection done)
# ============================================================
test_that("selective inference handles missing selection gracefully", {
  model <- make_inference_model(k = 20)
  res <- quantify_uncertainty(model, method = "selective")
  expect_true(!is.null(res$message) || !is.null(res$adjusted_p))
})

# ============================================================
# 5. Method label is propagated
# ============================================================
test_that("inference result carries method name", {
  model <- make_inference_model(k = 15)
  res <- quantify_uncertainty(model, method = "permutation", R = 20)
  expect_equal(res$method, "permutation")
})
