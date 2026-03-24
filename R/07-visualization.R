
#' @title Visualization Suite for Robust Meta-Analysis
#' @description Publication-quality plots using ggplot2 and base graphics
#' @import ggplot2
#' @importFrom graphics abline hist par qqline qqnorm
#' @importFrom grDevices rgb
#' @importFrom stats pnorm qnorm

#' Visualize Meta-Analysis Results
#'
#' Dispatcher function that creates publication-quality plots for meta-analysis
#' results, including forest, funnel, diagnostic, influence, cumulative,
#' radial, L'Abbe, and Baujat plots.
#'
#' @param model A metaverse S4 object or meta_robust list with \code{$data}.
#' @param type Character: plot type. One of "forest", "funnel", "diagnostic",
#'   "influence", "cumulative", "radial", "labbe", "baujat".
#' @param ... Additional arguments passed to the specific plot function.
#' @return A ggplot2 object (for most types) or NULL (for base-R diagnostic panel).
#' @export
visualize <- function(model,
                     type = c("forest", "funnel", "diagnostic", "influence",
                            "cumulative", "radial", "labbe", "baujat"),
                     ...) {

  type <- match.arg(type)

  plot_obj <- switch(type,
    "forest"     = forest_plot_complete(model, ...),
    "funnel"     = funnel_plot_complete(model, ...),
    "diagnostic" = diagnostic_plots_complete(model, ...),
    "influence"  = influence_plot(model, ...),
    "cumulative" = cumulative_plot(model, ...),
    "radial"     = radial_plot(model, ...),
    "labbe"      = labbe_plot(model, ...),
    "baujat"     = baujat_plot(model, ...)
  )

  plot_obj
}

# Helper to extract data from either metaverse S4 or meta_robust list
.extract_model <- function(model) {
  if (methods::is(model, "metaverse")) {
    list(
      data = model@data,
      est  = model@robust$estimate,
      se   = model@robust$se,
      ci_lo = model@robust$ci.lower,
      ci_hi = model@robust$ci.upper,
      weights = model@robust$weights
    )
  } else if (is.list(model) && !is.null(model$estimate)) {
    list(
      data = model$data,
      est  = model$estimate,
      se   = model$se,
      ci_lo = model$ci.lower,
      ci_hi = model$ci.upper,
      weights = model$weights
    )
  } else {
    stop("model must be a metaverse object or a meta_robust result")
  }
}

#' Forest plot
#' @export
forest_plot_complete <- function(model,
                                show_weights = TRUE,
                                show_prediction = TRUE,
                                sort_by = NULL,
                                annotate = TRUE,
                                theme = "classic") {

  m   <- .extract_model(model)
  data <- m$data
  est  <- m$est

  if (is.null(data) || nrow(data) == 0) stop("No data found in model object")

  yi  <- data$yi
  vi  <- data$vi
  sei <- sqrt(vi)
  n   <- length(yi)

  z_crit <- qnorm(0.975)
  lo <- yi - z_crit * sei
  hi <- yi + z_crit * sei

  ord <- if (!is.null(sort_by) && sort_by == "effect") order(yi) else seq_len(n)
  labels <- if ("study" %in% names(data)) as.character(data$study) else paste0("Study ", seq_len(n))

  df <- data.frame(
    study  = factor(labels[ord], levels = labels[ord]),
    yi     = yi[ord],
    lo     = lo[ord],
    hi     = hi[ord],
    weight = (1 / vi[ord]) / sum(1 / vi)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$yi, y = .data$study)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$lo, xmax = .data$hi), height = 0.2) +
    ggplot2::geom_vline(xintercept = est, linetype = "dashed", colour = "blue") +
    ggplot2::labs(x = "Effect Size", y = "", title = "Forest Plot") +
    ggplot2::theme_minimal()

  if (show_weights) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(x = max(df$hi) + 0.1,
                   label = paste0(round(.data$weight * 100, 1), "%")),
      hjust = 0, size = 3
    )
  }

  p
}

#' Funnel plot with contour enhancement
#' @export
funnel_plot_complete <- function(model,
                                show_contours = TRUE,
                                add_trim_fill = FALSE,
                                highlight_outliers = TRUE) {

  m    <- .extract_model(model)
  data <- m$data
  est  <- m$est

  if (is.null(data) || nrow(data) == 0) stop("No data found in model object")

  yi  <- data$yi
  sei <- sqrt(data$vi)

  df <- data.frame(yi = yi, sei = sei)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$yi, y = .data$sei)) +
    ggplot2::scale_y_reverse() +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(xintercept = est, linetype = "dashed") +
    ggplot2::labs(x = "Effect Size", y = "Standard Error", title = "Funnel Plot") +
    ggplot2::theme_minimal()

  if (show_contours) {
    se_range <- seq(0.001, max(sei) * 1.2, length.out = 200)
    for (level in c(0.10, 0.05, 0.01)) {
      z <- qnorm(1 - level / 2)
      cdf <- data.frame(se = se_range,
                         lo = est - z * se_range,
                         hi = est + z * se_range)
      p <- p +
        ggplot2::geom_line(data = cdf,
                           ggplot2::aes(x = .data$lo, y = .data$se),
                           linetype = "dotted", colour = "grey60") +
        ggplot2::geom_line(data = cdf,
                           ggplot2::aes(x = .data$hi, y = .data$se),
                           linetype = "dotted", colour = "grey60")
    }
  }

  if (highlight_outliers) {
    z95 <- qnorm(0.975)
    outlier_idx <- which(abs(yi - est) > z95 * sei)
    if (length(outlier_idx) > 0) {
      df_out <- data.frame(yi = yi[outlier_idx], sei = sei[outlier_idx])
      p <- p + ggplot2::geom_point(data = df_out,
                                    ggplot2::aes(x = .data$yi, y = .data$sei),
                                    colour = "red", size = 3, shape = 1)
    }
  }

  p
}

#' Diagnostic panel: QQ plot + residual plot + weight histogram
#' @export
diagnostic_plots_complete <- function(model, ...) {

  m   <- .extract_model(model)
  data <- m$data
  est  <- m$est
  wts  <- m$weights

  yi  <- data$yi
  vi  <- data$vi
  sei <- sqrt(vi)
  zi  <- (yi - est) / sei

  oldpar <- par(mfrow = c(1, 3))
  on.exit(par(oldpar))

  qqnorm(zi, main = "Normal Q-Q of Residuals", pch = 16, cex = 0.8)
  qqline(zi, col = "blue")

  plot(1 / sei, zi, pch = 16, cex = 0.8,
       xlab = "Precision (1/SE)", ylab = "Standardized Residual",
       main = "Residuals vs Precision")
  abline(h = c(-2, 0, 2), lty = c(2, 1, 2), col = c("red", "black", "red"))

  if (!is.null(wts) && length(wts) > 0) {
    hist(wts, breaks = 20, main = "Robust Weight Distribution",
         xlab = "Weight", col = "steelblue", border = "white")
  } else {
    wi <- 1 / vi
    hist(wi / sum(wi), breaks = 20, main = "Inverse-Variance Weights",
         xlab = "Relative Weight", col = "steelblue", border = "white")
  }

  invisible(NULL)
}

#' Influence (leave-one-out) plot
#' @export
influence_plot <- function(model, ...) {

  m   <- .extract_model(model)
  data <- m$data
  est  <- m$est

  yi <- data$yi
  vi <- data$vi
  n  <- length(yi)

  loo_est <- numeric(n)
  for (i in seq_len(n)) {
    fit_i <- meta_robust(yi[-i], vi[-i])
    loo_est[i] <- fit_i$estimate
  }

  labels <- if ("study" %in% names(data)) as.character(data$study) else paste0("Study ", seq_len(n))

  df <- data.frame(
    study = factor(labels, levels = rev(labels)),
    est   = loo_est
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$est, y = .data$study)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(xintercept = est, linetype = "dashed", colour = "blue") +
    ggplot2::labs(x = "Pooled Estimate (leave-one-out)", y = "",
                  title = "Influence Analysis") +
    ggplot2::theme_minimal()
}

#' Cumulative meta-analysis plot
#' @export
cumulative_plot <- function(model, ...) {

  m   <- .extract_model(model)
  data <- m$data

  yi <- data$yi
  vi <- data$vi
  n  <- length(yi)

  cum_est <- cum_lo <- cum_hi <- numeric(n)
  z_crit <- qnorm(0.975)

  for (i in seq_len(n)) {
    wi_i <- 1 / vi[1:i]
    mu_i <- sum(wi_i * yi[1:i]) / sum(wi_i)
    se_i <- sqrt(1 / sum(wi_i))
    cum_est[i] <- mu_i
    cum_lo[i]  <- mu_i - z_crit * se_i
    cum_hi[i]  <- mu_i + z_crit * se_i
  }

  labels <- if ("study" %in% names(data)) as.character(data$study) else paste0("Study ", seq_len(n))

  df <- data.frame(
    study = factor(labels, levels = labels),
    est   = cum_est,
    lo    = cum_lo,
    hi    = cum_hi
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$est, y = .data$study)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$lo, xmax = .data$hi), height = 0.2) +
    ggplot2::labs(x = "Cumulative Pooled Estimate", y = "",
                  title = "Cumulative Meta-Analysis") +
    ggplot2::theme_minimal()
}

#' Radial (Galbraith) plot
#' @export
radial_plot <- function(model, ...) {

  m   <- .extract_model(model)
  data <- m$data
  est  <- m$est

  yi  <- data$yi
  sei <- sqrt(data$vi)
  x   <- 1 / sei
  y   <- yi / sei

  df <- data.frame(precision = x, z_score = y)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$precision, y = .data$z_score)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_abline(intercept = 0, slope = est, colour = "blue") +
    ggplot2::geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", colour = "grey50") +
    ggplot2::labs(x = "1 / SE (Precision)", y = "yi / SE (Standardized Effect)",
                  title = "Radial (Galbraith) Plot") +
    ggplot2::theme_minimal()
}

#' L'Abbe plot (for binary outcomes)
#' @export
labbe_plot <- function(model, ...) {

  m   <- .extract_model(model)
  data <- m$data

  if (!all(c("pi1", "pi2") %in% names(data))) {
    message("L'Abbe plot requires pi1 (treatment rate) and pi2 (control rate).")
    message("Generating approximate rates from effect sizes.")
    data$pi2 <- pnorm(0)
    data$pi1 <- pnorm(data$yi)
  }

  wi    <- 1 / data$vi
  sizes <- 0.5 + 3 * wi / max(wi)

  df <- data.frame(pi2 = data$pi2, pi1 = data$pi1, w = sizes)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$pi2, y = .data$pi1)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$w)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::scale_size_identity() +
    ggplot2::labs(x = "Control Event Rate", y = "Treatment Event Rate",
                  title = "L'Abbe Plot") +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::theme_minimal()
}

#' Baujat plot (contribution vs influence)
#' @export
baujat_plot <- function(model, ...) {

  m   <- .extract_model(model)
  data <- m$data
  est  <- m$est

  yi <- data$yi
  vi <- data$vi
  n  <- length(yi)
  wi <- 1 / vi

  Q_i <- wi * (yi - est)^2

  influence_i <- numeric(n)
  for (i in seq_len(n)) {
    fit_i <- meta_robust(yi[-i], vi[-i])
    influence_i[i] <- abs(est - fit_i$estimate)
  }

  labels <- if ("study" %in% names(data)) as.character(data$study) else paste0("Study ", seq_len(n))

  df <- data.frame(Q_contrib = Q_i, influence = influence_i, label = labels)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$Q_contrib, y = .data$influence)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label),
                       hjust = -0.1, vjust = -0.2, size = 2.5) +
    ggplot2::labs(x = "Contribution to Q", y = "Influence on Pooled Estimate",
                  title = "Baujat Plot") +
    ggplot2::theme_minimal()
}
