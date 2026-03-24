# tests/testthat/test-visualization.R
# Tests for visualize() and individual plot functions
# These verify the functions run without error and return expected types.

context("visualize: plotting functions")

# ============================================================
# Helper
# ============================================================
make_plot_model <- function(k = 15, seed = 42) {
  data <- simulate_meta_data(k = k, theta = 0.4, tau2 = 0.08, seed = seed)
  rob <- meta_robust(data$yi, data$vi, method = "MM")
  new("metaverse",
      data = data,
      model = list(method = "MM"),
      robust = unclass(rob))
}

# ============================================================
# 1. Forest plot
# ============================================================
test_that("forest_plot_complete returns a ggplot", {
  model <- make_plot_model()
  p <- forest_plot_complete(model)
  expect_s3_class(p, "ggplot")
})

test_that("forest_plot_complete with sort_by = 'effect' works", {
  model <- make_plot_model()
  p <- forest_plot_complete(model, sort_by = "effect")
  expect_s3_class(p, "ggplot")
})

# ============================================================
# 2. Funnel plot
# ============================================================
test_that("funnel_plot_complete returns a ggplot", {
  model <- make_plot_model()
  p <- funnel_plot_complete(model)
  expect_s3_class(p, "ggplot")
})

test_that("funnel_plot_complete with contours disabled works", {
  model <- make_plot_model()
  p <- funnel_plot_complete(model, show_contours = FALSE)
  expect_s3_class(p, "ggplot")
})

# ============================================================
# 3. Diagnostic plots (base R)
# ============================================================
test_that("diagnostic_plots_complete runs without error", {
  model <- make_plot_model()
  # Redirect to null device to avoid opening plot windows in tests
  pdf(file = NULL)
  expect_silent(diagnostic_plots_complete(model))
  dev.off()
})

# ============================================================
# 4. Influence plot
# ============================================================
test_that("influence_plot returns a ggplot", {
  model <- make_plot_model(k = 10)
  p <- influence_plot(model)
  expect_s3_class(p, "ggplot")
})

# ============================================================
# 5. Cumulative plot
# ============================================================
test_that("cumulative_plot returns a ggplot", {
  model <- make_plot_model()
  p <- cumulative_plot(model)
  expect_s3_class(p, "ggplot")
})

# ============================================================
# 6. Radial (Galbraith) plot
# ============================================================
test_that("radial_plot returns a ggplot", {
  model <- make_plot_model()
  p <- radial_plot(model)
  expect_s3_class(p, "ggplot")
})

# ============================================================
# 7. L'Abbe plot
# ============================================================
test_that("labbe_plot returns a ggplot (approximate rates)", {
  model <- make_plot_model()
  # labbe_plot generates approximate rates when pi1/pi2 are missing
  p <- labbe_plot(model)
  expect_s3_class(p, "ggplot")
})

# ============================================================
# 8. Baujat plot
# ============================================================
test_that("baujat_plot returns a ggplot", {
  model <- make_plot_model(k = 10)
  p <- baujat_plot(model)
  expect_s3_class(p, "ggplot")
})

# ============================================================
# 9. visualize() dispatcher
# ============================================================
test_that("visualize dispatcher routes to correct plot type", {
  model <- make_plot_model(k = 10)

  p_forest <- visualize(model, type = "forest")
  expect_s3_class(p_forest, "ggplot")

  p_funnel <- visualize(model, type = "funnel")
  expect_s3_class(p_funnel, "ggplot")

  p_radial <- visualize(model, type = "radial")
  expect_s3_class(p_radial, "ggplot")
})

# ============================================================
# 10. Edge case: k = 3 plots
# ============================================================
test_that("plots work with k = 3", {
  data <- simulate_meta_data(k = 3, theta = 0.3, tau2 = 0.01, seed = 99)
  rob <- meta_robust(data$yi, data$vi, method = "M")
  model <- new("metaverse",
               data = data,
               model = list(method = "M"),
               robust = unclass(rob))

  p <- forest_plot_complete(model)
  expect_s3_class(p, "ggplot")

  p2 <- funnel_plot_complete(model)
  expect_s3_class(p2, "ggplot")
})

# ============================================================
# 11. meta_robust list input to .extract_model
# ============================================================
test_that("plot functions accept meta_robust list with $data", {
  data <- simulate_meta_data(k = 10, seed = 55)
  rob <- meta_robust(data$yi, data$vi, method = "M")
  rob$data <- data  # attach data

  p <- forest_plot_complete(rob)
  expect_s3_class(p, "ggplot")
})
