# tests/testthat/test-convert-effects.R
# Roundtrip and correctness tests for convert_effect_sizes()

context("convert_effect_sizes: conversions and roundtrips")

# ============================================================
# 1. Identity conversion
# ============================================================
test_that("identity conversion returns input unchanged", {
  for (type in c("SMD", "OR", "logOR", "r", "z", "f")) {
    res <- convert_effect_sizes(0.5, vi = 0.04, from = type, to = type)
    expect_equal(res$es, 0.5, info = type)
    expect_equal(res$var, 0.04, info = type)
  }
})

# ============================================================
# 2. SMD <-> logOR roundtrip
# ============================================================
test_that("SMD -> logOR -> SMD roundtrip is approximately correct", {
  d <- 0.6
  v <- 0.05
  step1 <- convert_effect_sizes(d, vi = v, from = "SMD", to = "logOR")
  step2 <- convert_effect_sizes(step1$es, vi = step1$var, from = "logOR", to = "SMD")
  expect_equal(step2$es, d, tolerance = 1e-10)
  expect_equal(step2$var, v, tolerance = 1e-10)
})

# ============================================================
# 3. SMD -> r and r -> SMD roundtrip
# ============================================================
test_that("SMD -> r -> SMD roundtrip works (equal groups)", {
  d <- 0.4
  v <- 0.03
  step1 <- convert_effect_sizes(d, vi = v, from = "SMD", to = "r")
  # r should be between -1 and 1

  expect_true(abs(step1$es) < 1)
  step2 <- convert_effect_sizes(step1$es, vi = step1$var, from = "r", to = "SMD")
  expect_equal(step2$es, d, tolerance = 1e-6)
})

# ============================================================
# 4. r <-> z (Fisher) roundtrip
# ============================================================
test_that("r -> z -> r roundtrip is exact", {
  r <- 0.35
  v_r <- 0.01
  step1 <- convert_effect_sizes(r, vi = v_r, from = "r", to = "z", n1 = 50)
  # Fisher's z should be atanh(r)
  expect_equal(step1$es, atanh(r), tolerance = 1e-10)
  step2 <- convert_effect_sizes(step1$es, vi = step1$var, from = "z", to = "r")
  expect_equal(step2$es, r, tolerance = 1e-10)
})

# ============================================================
# 5. SMD -> Hedges' g
# ============================================================
test_that("SMD -> g applies J correction factor", {
  d <- 0.8
  v <- 0.05
  n1 <- 15
  n2 <- 15
  df <- n1 + n2 - 2
  J <- 1 - (3 / (4 * df - 1))
  res <- convert_effect_sizes(d, vi = v, from = "SMD", to = "g", n1 = n1, n2 = n2)
  expect_equal(res$es, d * J, tolerance = 1e-10)
  expect_equal(res$var, v * J^2, tolerance = 1e-10)
})

# ============================================================
# 6. SMD -> f and f -> SMD roundtrip
# ============================================================
test_that("SMD -> f -> SMD roundtrip is exact", {
  d <- 0.5
  v <- 0.04
  step1 <- convert_effect_sizes(d, vi = v, from = "SMD", to = "f")
  expect_equal(step1$es, d / 2, tolerance = 1e-10)
  step2 <- convert_effect_sizes(step1$es, vi = step1$var, from = "f", to = "SMD")
  expect_equal(step2$es, d, tolerance = 1e-10)
})

# ============================================================
# 7. f -> eta2
# ============================================================
test_that("f -> eta2 is correct", {
  f <- 0.25
  res <- convert_effect_sizes(f, from = "f", to = "eta2")
  expected <- f^2 / (1 + f^2)
  expect_equal(res$es, expected, tolerance = 1e-10)
})

# ============================================================
# 8. f -> f2
# ============================================================
test_that("f -> f2 is correct", {
  f <- 0.3
  v <- 0.01
  res <- convert_effect_sizes(f, vi = v, from = "f", to = "f2")
  expect_equal(res$es, f^2, tolerance = 1e-10)
})

# ============================================================
# 9. OR -> SMD -> OR roundtrip
# ============================================================
test_that("OR -> SMD -> OR roundtrip via logOR is consistent", {
  lor <- 0.7
  v_lor <- 0.06
  step1 <- convert_effect_sizes(lor, vi = v_lor, from = "logOR", to = "SMD")
  step2 <- convert_effect_sizes(step1$es, vi = step1$var, from = "SMD", to = "logOR")
  expect_equal(step2$es, lor, tolerance = 1e-10)
})

# ============================================================
# 10. z -> SMD via r
# ============================================================
test_that("z -> SMD works", {
  z <- 0.4
  res <- convert_effect_sizes(z, vi = 0.01, from = "z", to = "SMD")
  expect_true(is.numeric(res$es))
  expect_true(!is.na(res$es))
})

# ============================================================
# 11. Invalid types error
# ============================================================
test_that("invalid type raises error", {
  expect_error(convert_effect_sizes(0.5, from = "XYZ", to = "SMD"), "must be one of")
  expect_error(convert_effect_sizes(0.5, from = "SMD", to = "XYZ"), "must be one of")
})

# ============================================================
# 12. NULL vi is handled
# ============================================================
test_that("conversion with vi = NULL returns only es", {
  res <- convert_effect_sizes(0.5, vi = NULL, from = "SMD", to = "logOR")
  expect_true(!is.null(res$es))
  expect_true(is.null(res$var))
})

# ============================================================
# 13. Vector inputs
# ============================================================
test_that("vectorized SMD -> logOR works", {
  d <- c(0.2, 0.5, 0.8)
  v <- c(0.01, 0.03, 0.05)
  res <- convert_effect_sizes(d, vi = v, from = "SMD", to = "logOR")
  expect_length(res$es, 3)
  expect_length(res$var, 3)
  # logOR = d * pi / sqrt(3)
  expect_equal(res$es, d * pi / sqrt(3), tolerance = 1e-10)
})
