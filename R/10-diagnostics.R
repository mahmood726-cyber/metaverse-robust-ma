
#' @title Model Diagnostics for Robust Meta-Analysis
#' @description Comprehensive diagnostic checks for fitted models
#' @importFrom stats shapiro.test ks.test

#' Run diagnostic suite on a meta_robust result or metaverse object
#'
#' Checks residual normality, heterogeneity, and convergence.
#'
#' @param model A meta_robust result or metaverse S4 object.
#' @param ... Additional arguments (unused).
#' @return A list of diagnostic results.
#' @export
diagnose <- function(model, ...) {

  # Extract data depending on object type

  if (methods::is(model, "metaverse")) {
    data <- model@data
    est  <- model@robust$estimate
    se   <- model@robust$se
    wts  <- model@robust$weights
    converged <- model@robust$converged
  } else if (is.list(model) && !is.null(model$estimate)) {
    data <- model$data
    est  <- model$estimate
    se   <- model$se
    wts  <- model$weights
    converged <- model$converged
  } else {
    stop("model must be a metaverse object or a meta_robust result")
  }

  yi  <- data$yi
  vi  <- data$vi
  sei <- sqrt(vi)
  n   <- length(yi)

  # 1. Standardized residuals
  zi <- (yi - est) / sei

  # 2. Normality test on residuals
  normality <- if (n >= 3 && n <= 5000) {
    shapiro.test(zi)
  } else {
    ks.test(zi, "pnorm")
  }

  # 3. Heterogeneity
  wi <- 1 / vi
  Q  <- sum(wi * (yi - est)^2)
  Q_df <- n - 1
  Q_p  <- pchisq(Q, df = Q_df, lower.tail = FALSE)
  I2   <- max(0, (Q - Q_df) / Q * 100)

  # 4. Number of outlier residuals
  n_outliers_2sd <- sum(abs(zi) > 2)
  n_outliers_3sd <- sum(abs(zi) > 3)

  # 5. Weight diagnostics
  weight_diag <- if (!is.null(wts) && length(wts) > 0) {
    list(
      min_weight = min(wts),
      max_weight = max(wts),
      median_weight = median(wts),
      n_downweighted = sum(wts < 0.5),
      prop_downweighted = mean(wts < 0.5)
    )
  } else {
    NULL
  }

  # 6. Convergence
  convergence_ok <- if (!is.null(converged)) converged else NA

  diagnostics <- list(
    n_studies = n,
    estimate  = est,
    se        = se,
    residuals_z = zi,
    normality_test = list(
      statistic = normality$statistic,
      p_value   = normality$p.value,
      method    = normality$method,
      normal    = normality$p.value > 0.05
    ),
    heterogeneity = list(
      Q     = Q,
      Q_df  = Q_df,
      Q_p   = Q_p,
      I2    = I2,
      significant = Q_p < 0.05
    ),
    outliers = list(
      n_2sd = n_outliers_2sd,
      n_3sd = n_outliers_3sd,
      indices_2sd = which(abs(zi) > 2),
      indices_3sd = which(abs(zi) > 3)
    ),
    weights = weight_diag,
    converged = convergence_ok
  )

  class(diagnostics) <- c("meta_diagnostics", "list")
  diagnostics
}

#' @export
print.meta_diagnostics <- function(x, ...) {
  cat("\nMeta-Analysis Diagnostics\n")
  cat("========================\n\n")

  cat("Studies:", x$n_studies, "\n")
  cat("Estimate:", round(x$estimate, 4), "(SE =", round(x$se, 4), ")\n\n")

  cat("Residual normality (", x$normality_test$method, "):\n", sep = "")
  cat("  p-value =", format.pval(x$normality_test$p_value, digits = 3))
  cat(if (x$normality_test$normal) " [OK]\n" else " [WARNING: non-normal]\n")

  cat("\nHeterogeneity:\n")
  cat("  Q =", round(x$heterogeneity$Q, 2),
      "  df =", x$heterogeneity$Q_df,
      "  p =", format.pval(x$heterogeneity$Q_p, digits = 3), "\n")
  cat("  I^2 =", round(x$heterogeneity$I2, 1), "%\n")

  cat("\nOutliers:\n")
  cat("  |z| > 2:", x$outliers$n_2sd, "\n")
  cat("  |z| > 3:", x$outliers$n_3sd, "\n")

  if (!is.null(x$weights)) {
    cat("\nRobust weights:\n")
    cat("  Down-weighted (w < 0.5):", x$weights$n_downweighted,
        "(", round(x$weights$prop_downweighted * 100, 1), "%)\n")
  }

  if (!is.na(x$converged)) {
    cat("\nConvergence:", if (x$converged) "YES" else "NO", "\n")
  }

  invisible(x)
}
