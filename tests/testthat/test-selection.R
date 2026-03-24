# tests/testthat/test-selection.R
# Comprehensive tests for select_moderators() and related selection methods

context("select_moderators: variable selection methods")

# ============================================================
# Helper data
# ============================================================
make_selection_data <- function(k = 60, p = 5, seed = 42) {
  set.seed(seed)
  X <- matrix(rnorm(k * p), k, p)
  colnames(X) <- paste0("X", seq_len(p))
  # True effect on first variable only
  beta <- c(0.5, rep(0, p - 1))
  theta <- 0.3 + X %*% beta
  tau2 <- 0.05
  vi <- runif(k, 0.01, 0.08)
  yi <- rnorm(k, theta, sqrt(vi + tau2))
  list(X = X, yi = yi, vi = vi, beta_true = beta)
}

# ============================================================
# 1. Knockoff filter
# ============================================================
test_that("knockoff selection runs without error", {
  d <- make_selection_data(k = 80, p = 5)
  res <- select_moderators(d$X, d$yi, d$vi, method = "knockoff",
                           fdr = 0.2, control = list(randomize = FALSE))
  expect_s3_class(res, "selection_result")
  expect_true(is.integer(res$selected) || is.numeric(res$selected))
  expect_equal(res$method, "knockoff")
  # importance should have length p
  expect_length(res$importance, ncol(d$X))
})

test_that("knockoff with different knockoff types works", {
  d <- make_selection_data(k = 60, p = 4)
  for (ktype in c("sdp", "equi", "sequential")) {
    res <- select_moderators(d$X, d$yi, d$vi, method = "knockoff",
                             fdr = 0.2,
                             control = list(knockoff_type = ktype,
                                            randomize = FALSE))
    expect_s3_class(res, "selection_result", info = ktype)
  }
})

# ============================================================
# 2. Spike-and-slab MCMC
# ============================================================
test_that("spike_slab selection runs and returns MCMC output", {
  d <- make_selection_data(k = 50, p = 4)
  res <- select_moderators(d$X, d$yi, d$vi, method = "spike_slab",
                           control = list(n_iter = 500, burn_in = 100, thin = 1))
  expect_s3_class(res, "selection_result")
  expect_equal(res$method, "spike_slab")
  # PIP should be between 0 and 1
  expect_true(all(res$pip >= 0 & res$pip <= 1))
  # Should have coefficient means
  expect_length(res$coefficients, ncol(d$X))
  # Gamma samples should exist
  expect_true(!is.null(res$gamma_samples))
})

test_that("spike_slab with single moderator works", {
  d <- make_selection_data(k = 40, p = 1)
  res <- select_moderators(d$X, d$yi, d$vi, method = "spike_slab",
                           control = list(n_iter = 300, burn_in = 50, thin = 1))
  expect_s3_class(res, "selection_result")
  expect_length(res$pip, 1)
})

# ============================================================
# 3. LASSO with stability selection
# ============================================================
test_that("lasso selection works", {
  d <- make_selection_data(k = 60, p = 5)
  res <- select_moderators(d$X, d$yi, d$vi, method = "lasso",
                           control = list(n_bootstrap = 20))
  expect_s3_class(res, "selection_result")
  expect_equal(res$method, "lasso")
  expect_true(!is.null(res$selection_prob))
  expect_length(res$selection_prob, ncol(d$X))
})

# ============================================================
# 4. Elastic net
# ============================================================
test_that("elastic_net selection works", {
  d <- make_selection_data(k = 60, p = 5)
  res <- select_moderators(d$X, d$yi, d$vi, method = "elastic_net")
  expect_s3_class(res, "selection_result")
  expect_equal(res$method, "elastic_net")
  expect_true(!is.null(res$alpha))  # best alpha found
})

# ============================================================
# 5. SCAD and MCP
# ============================================================
test_that("scad selection works", {
  d <- make_selection_data(k = 50, p = 5)
  res <- select_moderators(d$X, d$yi, d$vi, method = "scad")
  expect_s3_class(res, "selection_result")
  expect_equal(res$method, "scad")
})

test_that("mcp selection works", {
  d <- make_selection_data(k = 50, p = 5)
  res <- select_moderators(d$X, d$yi, d$vi, method = "mcp")
  expect_s3_class(res, "selection_result")
  expect_equal(res$method, "mcp")
})

# ============================================================
# 6. E-values based selection
# ============================================================
test_that("evalues selection works", {
  d <- make_selection_data(k = 50, p = 4)
  res <- select_moderators(d$X, d$yi, d$vi, method = "evalues")
  expect_s3_class(res, "selection_result")
  expect_equal(res$method, "evalues")
  expect_true(!is.null(res$e_values))
  expect_true(!is.null(res$p_values))
  expect_length(res$e_values, ncol(d$X))
})

# ============================================================
# 7. Edge cases
# ============================================================
test_that("selection with single moderator does not crash", {
  set.seed(77)
  k <- 40
  X <- matrix(rnorm(k), k, 1)
  colnames(X) <- "X1"
  yi <- rnorm(k, 0.3 + 0.5 * X[, 1], 0.2)
  vi <- rep(0.04, k)

  for (m in c("lasso", "elastic_net", "evalues")) {
    res <- select_moderators(X, yi, vi, method = m)
    expect_s3_class(res, "selection_result", info = m)
  }
})

test_that("selection correctly maps variable names", {
  d <- make_selection_data(k = 60, p = 5)
  res <- select_moderators(d$X, d$yi, d$vi, method = "lasso",
                           control = list(n_bootstrap = 10, threshold = 0.01))
  if (length(res$selected) > 0) {
    expect_true(!is.null(res$variable_names))
    expect_true(all(res$variable_names %in% colnames(d$X)))
  }
})

# ============================================================
# 8. Print method
# ============================================================
test_that("print.selection_result runs without error", {
  d <- make_selection_data(k = 50, p = 3)
  res <- select_moderators(d$X, d$yi, d$vi, method = "lasso",
                           control = list(n_bootstrap = 10))
  expect_output(print(res), "Variable Selection Results")
  expect_output(print(res), "lasso")
})
