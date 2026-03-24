
#' @title Causal Inference Helpers for Meta-Analysis
#' @description Transportability analysis and heterogeneous treatment effects
#' @importFrom stats lm predict coef residuals

#' Causal forest-style heterogeneous treatment effect estimation
#'
#' Uses study-level moderators to estimate conditional average treatment effects
#' (CATE) across the meta-analytic sample. This is a simplified "poor-man's
#' causal forest" that partitions studies by moderator splits.
#'
#' @param yi Numeric vector of effect sizes.
#' @param vi Numeric vector of sampling variances.
#' @param X  Matrix of study-level moderators.
#' @param n_splits Number of random splits for honest estimation (default 20).
#' @param min_leaf  Minimum number of studies per leaf (default 3).
#' @return A list with estimated CATEs, variable importance, and subgroup effects.
#' @export
causal_forest_meta <- function(yi, vi, X, n_splits = 20, min_leaf = 3) {

  if (!is.matrix(X)) X <- as.matrix(X)
  n <- length(yi)
  p <- ncol(X)

  if (n < 2 * min_leaf) {
    warning("Too few studies for causal subgroup analysis")
    return(list(
      cate    = rep(mean(yi), n),
      importance = rep(1 / p, p),
      subgroups  = list()
    ))
  }

  # Weighted regression to estimate CATE(x) = E[theta | X=x]
  wi <- 1 / vi
  fit <- lm(yi ~ X, weights = wi)
  cate_hat <- as.numeric(predict(fit))

  # Variable importance via drop-in-R2
  r2_full <- summary(fit)$r.squared
  importance <- numeric(p)
  for (j in seq_len(p)) {
    fit_j <- lm(yi ~ X[, -j, drop = FALSE], weights = wi)
    r2_j  <- summary(fit_j)$r.squared
    importance[j] <- max(0, r2_full - r2_j)
  }
  importance <- importance / max(sum(importance), .Machine$double.eps)

  # Simple subgroup analysis: split each moderator at median
  subgroups <- list()
  mod_names <- if (!is.null(colnames(X))) colnames(X) else paste0("X", seq_len(p))

  for (j in seq_len(p)) {
    med_j <- median(X[, j])
    lo_idx <- which(X[, j] <= med_j)
    hi_idx <- which(X[, j] > med_j)

    if (length(lo_idx) >= min_leaf && length(hi_idx) >= min_leaf) {
      wi_lo <- 1 / vi[lo_idx]
      wi_hi <- 1 / vi[hi_idx]
      est_lo <- sum(wi_lo * yi[lo_idx]) / sum(wi_lo)
      est_hi <- sum(wi_hi * yi[hi_idx]) / sum(wi_hi)
      se_lo  <- sqrt(1 / sum(wi_lo))
      se_hi  <- sqrt(1 / sum(wi_hi))

      subgroups[[mod_names[j]]] <- data.frame(
        group    = c("low", "high"),
        estimate = c(est_lo, est_hi),
        se       = c(se_lo, se_hi),
        k        = c(length(lo_idx), length(hi_idx)),
        split_at = rep(med_j, 2)
      )
    }
  }

  list(
    cate       = cate_hat,
    importance = setNames(importance, mod_names),
    subgroups  = subgroups,
    model      = fit
  )
}

#' Transportability analysis
#'
#' Assesses whether a meta-analytic finding can be transported to a new target
#' population characterized by a different covariate distribution.
#'
#' @param yi Numeric vector of effect sizes.
#' @param vi Numeric vector of sampling variances.
#' @param X  Matrix of study-level moderators.
#' @param target_X  Matrix or vector of target population moderator values.
#' @return A list with the transported estimate and its uncertainty.
#' @export
transport <- function(yi, vi, X, target_X) {


  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(target_X)) target_X <- matrix(target_X, nrow = 1)

  n <- length(yi)
  wi <- 1 / vi

  # Fit weighted meta-regression
  fit <- lm(yi ~ X, weights = wi)

  # Predict for target population
  newdata <- data.frame(X = target_X)
  # Build prediction manually to avoid formula issues
  beta_hat <- coef(fit)
  X_pred   <- cbind(1, target_X)
  est_target <- as.numeric(X_pred %*% beta_hat)

  # Variance of prediction
  vcov_beta <- vcov(fit)
  se_target <- sqrt(diag(X_pred %*% vcov_beta %*% t(X_pred)))

  z_crit <- qnorm(0.975)

  list(
    estimate = est_target,
    se       = se_target,
    ci_lower = est_target - z_crit * se_target,
    ci_upper = est_target + z_crit * se_target,
    n_source = n,
    n_target = nrow(target_X)
  )
}
