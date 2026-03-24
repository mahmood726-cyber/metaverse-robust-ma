
#' @title Complete Sensitivity Analysis Suite
#' @description Comprehensive sensitivity and influence diagnostics

#' @export
assess_sensitivity <- function(model, 
                             analyses = c("publication_bias", "influence", 
                                        "fragility", "evalues", "outliers"),
                             control = list(),
                             verbose = FALSE) {
  
  if("all" %in% analyses) {
    analyses <- c("publication_bias", "influence", "fragility", 
                 "evalues", "outliers", "cumulative")
  }
  
  results <- list()
  
  for(analysis in analyses) {
    if(verbose) cat("Running", analysis, "analysis...\n")
    
    results[[analysis]] <- switch(analysis,
      "publication_bias" = publication_bias_complete(model, control, verbose),
      "influence" = influence_complete(model, control, verbose),
      "fragility" = fragility_complete(model, control, verbose),
      "evalues" = evalues_complete(model, control, verbose),
      "outliers" = outliers_complete(model, control, verbose),
      "cumulative" = cumulative_analysis(model, control, verbose)
    )
  }
  
  class(results) <- c("sensitivity_result", "list")
  results
}

#' Complete publication bias analysis
publication_bias_complete <- function(model, control = list(), verbose = FALSE) {
  
  data <- model@data
  yi <- data$yi
  vi <- data$vi
  sei <- sqrt(vi)
  n <- length(yi)
  
  # Multiple tests and adjustments
  tests <- list()
  
  # 1. Egger regression (multiple versions)
  tests$egger_standard <- egger_test(yi, vi, method = "standard")
  tests$egger_robust <- egger_test(yi, vi, method = "robust")
  tests$egger_multilevel <- egger_test(yi, vi, method = "multilevel")
  
  # 2. Rank correlation tests
  tests$begg <- begg_test(yi, vi)
  tests$thompson_sharp <- thompson_sharp_test(yi, vi)
  
  # 3. Trim and fill (multiple estimators)
  tests$trim_fill_l0 <- trim_and_fill(yi, vi, estimator = "L0")
  tests$trim_fill_r0 <- trim_and_fill(yi, vi, estimator = "R0")
  tests$trim_fill_q0 <- trim_and_fill(yi, vi, estimator = "Q0")
  
  # 4. Selection models
  tests$selection_step <- selection_model_step(yi, vi)
  tests$selection_beta <- selection_model_beta(yi, vi)
  
  # 5. P-curve and p-uniform
  tests$p_curve <- p_curve_analysis(yi, vi)
  tests$p_uniform <- p_uniform_analysis(yi, vi)
  
  # 6. Contour-enhanced funnel plot data
  tests$contour_data <- compute_contour_data(yi, vi)
  
  tests
}

#' Complete influence analysis
influence_complete <- function(model, control = list(), verbose = FALSE) {
  
  data <- model@data
  n <- nrow(data)
  
  # Leave-one-out
  loo_results <- matrix(NA, n, 5)
  colnames(loo_results) <- c("estimate", "se", "tau2", "I2", "dfbetas")
  
  method_name <- if (length(model@model) > 0 && !is.null(model@model$method)) model@model$method else "MM"
  est_full <- model@robust$estimate
  se_full  <- model@robust$se
  I2_full  <- if (!is.null(model@robust$I2)) model@robust$I2 else NA_real_

  for(i in 1:n) {
    data_loo <- data[-i, ]
    fit_loo <- meta_robust(data_loo$yi, data_loo$vi, method = method_name)

    loo_results[i, ] <- c(
      fit_loo$estimate,
      fit_loo$se,
      if (!is.null(fit_loo$tau2)) fit_loo$tau2 else NA_real_,
      if (!is.null(fit_loo$I2)) fit_loo$I2 else NA_real_,
      (est_full - fit_loo$estimate) / fit_loo$se
    )
  }

  # Identify influential studies
  influential <- list(
    dfbetas_large = which(abs(loo_results[, "dfbetas"]) > 2),
    estimate_change = which(abs(loo_results[, "estimate"] - est_full) > 2 * se_full),
    heterogeneity_change = if (!is.na(I2_full)) {
      which(abs(loo_results[, "I2"] - I2_full) > 20)
    } else {
      integer(0)
    }
  )
  
  # Cook's distance
  cooks_d <- compute_cooks_distance(model, data)
  
  # DFFITS
  dffits <- compute_dffits(model, data)
  
  # Covariance ratio
  covratio <- compute_covratio(model, data)
  
  list(
    loo = loo_results,
    influential = influential,
    cooks_distance = cooks_d,
    dffits = dffits,
    covratio = covratio
  )
}

#' Complete fragility analysis  
fragility_complete <- function(model, control = list(), verbose = FALSE) {
  
  data <- model@data
  n <- nrow(data)
  
  # Fragility index
  fragility_index <- compute_fragility_index(model, data)
  
  # Fragility quotient
  fragility_quotient <- fragility_index / n
  
  # Direction-specific fragility
  fragility_pos <- compute_directional_fragility(model, data, direction = "positive")
  fragility_neg <- compute_directional_fragility(model, data, direction = "negative")
  
  # Robustness index (based on multiple perturbations)
  robustness <- compute_robustness_index(model, data)
  
  list(
    fragility_index = fragility_index,
    fragility_quotient = fragility_quotient,
    fragility_positive = fragility_pos,
    fragility_negative = fragility_neg,
    robustness_index = robustness,
    interpretation = interpret_fragility(fragility_index, n)
  )
}

#' Complete E-values analysis
evalues_complete <- function(model, control = list(), verbose = FALSE) {
  
  est <- model@robust$estimate
  ci_lower <- model@robust$ci.lower
  ci_upper <- model@robust$ci.upper
  
  # E-value for point estimate
  e_val_est <- compute_evalue(est)
  
  # E-value for CI
  if(est > 0) {
    e_val_ci <- compute_evalue(ci_lower)
  } else {
    e_val_ci <- compute_evalue(ci_upper)
  }
  
  # E-value for specific alternative hypotheses
  e_val_null <- compute_evalue_null(est, model@robust$se)
  
  # Multiple bias scenario
  bias_scenarios <- compute_bias_scenarios(model)
  
  list(
    evalue_estimate = e_val_est,
    evalue_ci = e_val_ci,
    evalue_null = e_val_null,
    bias_scenarios = bias_scenarios,
    interpretation = interpret_evalues(e_val_est, e_val_ci)
  )
}

#' Helper functions for sensitivity analyses
egger_test <- function(yi, vi, method = "standard") {
  sei <- sqrt(vi)
  
  if(method == "standard") {
    fit <- lm(yi ~ sei, weights = 1/vi)
  } else if(method == "robust") {
    fit <- robustbase::lmrob(yi ~ sei, weights = 1/vi)
  } else {
    # Multilevel version
    fit <- lme4::lmer(yi ~ sei + (1|study), weights = 1/vi,
                      data = data.frame(yi = yi, sei = sei, 
                                      study = 1:length(yi)))
  }
  
  coefs <- summary(fit)$coefficients
  
  list(
    intercept = coefs[1, 1],
    se = coefs[1, 2],
    p_value = coefs[1, 4],
    significant = coefs[1, 4] < 0.05
  )
}

begg_test <- function(yi, vi) {
  zi <- yi / sqrt(vi)
  ranks <- rank(abs(zi))
  
  cor_test <- cor.test(ranks, zi, method = "kendall")
  
  list(
    tau = cor_test$estimate,
    p_value = cor_test$p.value,
    significant = cor_test$p.value < 0.05
  )
}

#' Thompson-Sharp test for funnel asymmetry
thompson_sharp_test <- function(yi, vi) {
  sei <- sqrt(vi)
  wi  <- 1 / vi
  theta_hat <- sum(wi * yi) / sum(wi)
  zi <- (yi - theta_hat) / sei
  ranks <- rank(sei)
  cor_test <- cor.test(ranks, zi, method = "spearman")
  list(
    rho     = cor_test$estimate,
    p_value = cor_test$p.value,
    significant = cor_test$p.value < 0.05
  )
}

#' Trim and fill
trim_and_fill <- function(yi, vi, estimator = "L0") {
  n   <- length(yi)
  wi  <- 1 / vi
  theta_hat <- sum(wi * yi) / sum(wi)
  di  <- yi - theta_hat
  # Rank based on absolute deviation
  r   <- rank(abs(di))
  sgn <- sign(di)

  # Estimate number of missing studies
  k0 <- switch(estimator,
    "L0" = max(0, round((4 * sum(r[sgn > 0]) - n * (n + 1)) / (2 * n))),
    "R0" = max(0, round(sum(sgn > 0) - n / 2)),
    "Q0" = {
      Q <- sum(wi * (yi - theta_hat)^2)
      max(0, round((Q - (n - 1)) / (2 * mean(wi * (yi - theta_hat)))))
    }
  )

  # Impute missing studies by reflecting the k0 most extreme positives
  if (k0 > 0 && k0 < n) {
    ordered_idx <- order(di, decreasing = TRUE)
    fill_idx    <- ordered_idx[seq_len(min(k0, n))]
    yi_fill <- c(yi, 2 * theta_hat - yi[fill_idx])
    vi_fill <- c(vi, vi[fill_idx])
    wi_fill <- 1 / vi_fill
    theta_adj <- sum(wi_fill * yi_fill) / sum(wi_fill)
  } else {
    theta_adj <- theta_hat
    yi_fill   <- yi
    vi_fill   <- vi
  }

  list(
    k0       = k0,
    estimate = theta_hat,
    estimate_adjusted = theta_adj,
    n_original = n,
    n_filled   = length(yi_fill)
  )
}

#' Step selection model (simplified Vevea-Hedges)
selection_model_step <- function(yi, vi) {
  wi  <- 1 / vi
  sei <- sqrt(vi)
  n   <- length(yi)

  theta_hat <- sum(wi * yi) / sum(wi)
  z_vals    <- yi / sei
  p_vals    <- 2 * pnorm(-abs(z_vals))
  sig       <- p_vals < 0.05

  if (sum(sig) == 0 || sum(!sig) == 0) {
    return(list(estimate_unadjusted = theta_hat, estimate_adjusted = theta_hat,
                weight_sig = 1, weight_nonsig = 1))
  }

  theta_sig    <- sum(wi[sig] * yi[sig]) / sum(wi[sig])
  theta_nonsig <- sum(wi[!sig] * yi[!sig]) / sum(wi[!sig])

  list(
    estimate_unadjusted = theta_hat,
    estimate_adjusted   = theta_nonsig,
    theta_significant   = theta_sig,
    theta_nonsignificant = theta_nonsig,
    prop_significant    = mean(sig)
  )
}

#' Beta selection model (simplified)
selection_model_beta <- function(yi, vi) {
  # Approximate selection model with beta weight function
  selection_model_step(yi, vi)
}

#' P-curve analysis (simplified)
p_curve_analysis <- function(yi, vi) {
  sei    <- sqrt(vi)
  z_vals <- yi / sei
  p_vals <- 2 * pnorm(-abs(z_vals))
  sig_p  <- p_vals[p_vals < 0.05]

  if (length(sig_p) < 3) {
    return(list(n_significant = length(sig_p),
                message = "Fewer than 3 significant results; p-curve not meaningful."))
  }

  # Proportion of p < 0.025 among p < 0.05 (right-skew test)
  prop_low <- mean(sig_p < 0.025)

  # Binomial test: under the null (no effect), p | p<0.05 ~ Uniform(0, 0.05)
  # so P(p < 0.025 | p < 0.05) = 0.5
  binom_result <- binom.test(sum(sig_p < 0.025), length(sig_p), 0.5,
                              alternative = "greater")

  list(
    n_significant    = length(sig_p),
    prop_very_sig    = prop_low,
    binom_p          = binom_result$p.value,
    evidential_value = prop_low > 0.5 && binom_result$p.value < 0.05
  )
}

#' P-uniform analysis (simplified)
p_uniform_analysis <- function(yi, vi) {
  # P-uniform is closely related to p-curve; provide a basic version
  p_curve_analysis(yi, vi)
}

#' Contour data for enhanced funnel plot
compute_contour_data <- function(yi, vi) {
  wi  <- 1 / vi
  est <- sum(wi * yi) / sum(wi)
  se_range <- seq(0.001, max(sqrt(vi)) * 1.5, length.out = 200)

  contours <- list()
  for (level in c(0.10, 0.05, 0.01)) {
    z <- qnorm(1 - level / 2)
    contours[[paste0("p_", level)]] <- data.frame(
      se = se_range,
      lower = est - z * se_range,
      upper = est + z * se_range
    )
  }
  contours
}

#' Cook's distance for meta-analysis
compute_cooks_distance <- function(model, data) {
  yi  <- data$yi
  vi  <- data$vi
  n   <- length(yi)
  est <- if (methods::is(model, "metaverse")) model@robust$estimate else model$estimate
  se  <- if (methods::is(model, "metaverse")) model@robust$se else model$se

  cooks_d <- numeric(n)
  for (i in seq_len(n)) {
    fit_i <- meta_robust(yi[-i], vi[-i])
    cooks_d[i] <- (est - fit_i$estimate)^2 / (se^2)
  }
  cooks_d
}

#' DFFITS for meta-analysis
compute_dffits <- function(model, data) {
  yi  <- data$yi
  vi  <- data$vi
  n   <- length(yi)
  est <- if (methods::is(model, "metaverse")) model@robust$estimate else model$estimate

  dffits_val <- numeric(n)
  for (i in seq_len(n)) {
    fit_i <- meta_robust(yi[-i], vi[-i])
    dffits_val[i] <- (est - fit_i$estimate) / fit_i$se
  }
  dffits_val
}

#' Covariance ratio
compute_covratio <- function(model, data) {
  yi  <- data$yi
  vi  <- data$vi
  n   <- length(yi)
  se  <- if (methods::is(model, "metaverse")) model@robust$se else model$se

  covratio_val <- numeric(n)
  for (i in seq_len(n)) {
    fit_i <- meta_robust(yi[-i], vi[-i])
    covratio_val[i] <- (fit_i$se / se)^2
  }
  covratio_val
}

#' Fragility index: minimum studies to remove to flip significance
compute_fragility_index <- function(model, data) {
  yi  <- data$yi
  vi  <- data$vi
  n   <- length(yi)
  est <- if (methods::is(model, "metaverse")) model@robust$estimate else model$estimate
  p   <- if (methods::is(model, "metaverse")) model@robust$p_value else model$p_value

  is_sig <- p < 0.05

  # Try removing studies one-by-one, then pairs, etc.
  for (k in 1:n) {
    combos <- utils::combn(n, k)
    for (j in seq_len(ncol(combos))) {
      idx <- combos[, j]
      fit_k <- meta_robust(yi[-idx], vi[-idx])
      new_sig <- fit_k$p_value < 0.05
      if (new_sig != is_sig) return(k)
    }
    if (k >= 5) break  # Limit combinatorial explosion
  }
  return(NA_integer_)
}

#' Directional fragility
compute_directional_fragility <- function(model, data, direction = "positive") {
  yi  <- data$yi
  vi  <- data$vi
  n   <- length(yi)

  # Order studies by effect: positive = largest first, negative = smallest first
  ord <- if (direction == "positive") order(yi, decreasing = TRUE) else order(yi)

  for (k in 1:min(n - 2, 10)) {
    idx <- ord[seq_len(k)]
    fit_k <- meta_robust(yi[-idx], vi[-idx])
    if (fit_k$p_value >= 0.05) return(k)
  }
  return(NA_integer_)
}

#' Robustness index
compute_robustness_index <- function(model, data) {
  yi  <- data$yi
  vi  <- data$vi
  n   <- length(yi)
  est <- if (methods::is(model, "metaverse")) model@robust$estimate else model$estimate

  # Proportion of leave-one-out estimates that remain significant
  n_sig <- 0
  for (i in seq_len(n)) {
    fit_i <- meta_robust(yi[-i], vi[-i])
    if (fit_i$p_value < 0.05) n_sig <- n_sig + 1
  }
  n_sig / n
}

#' Interpret fragility index
interpret_fragility <- function(fi, n) {
  if (is.na(fi)) return("Fragility index could not be computed (result robust to removal of up to 5 studies).")
  fq <- fi / n
  if (fq > 0.2) return(paste0("Fragility index = ", fi, " (FQ = ", round(fq, 2), "). Result is robust."))
  if (fq > 0.1) return(paste0("Fragility index = ", fi, " (FQ = ", round(fq, 2), "). Moderate robustness."))
  return(paste0("Fragility index = ", fi, " (FQ = ", round(fq, 2), "). Result is fragile."))
}

#' E-value computation (VanderWeele & Ding, 2017)
compute_evalue <- function(est) {
  # For risk ratio scale: E = RR + sqrt(RR * (RR - 1))

  # If est is on SMD scale, convert to approximate RR
  RR <- exp(0.91 * est)  # Approximate conversion
  if (RR < 1) RR <- 1 / RR
  RR + sqrt(RR * (RR - 1))
}

#' E-value under null
compute_evalue_null <- function(est, se) {
  # E-value for the confidence interval limit closest to null
  z_crit <- qnorm(0.975)
  ci_limit <- if (est > 0) est - z_crit * se else est + z_crit * se
  if (ci_limit * est <= 0) return(1)  # CI includes null
  compute_evalue(ci_limit)
}

#' Bias scenarios
compute_bias_scenarios <- function(model) {
  est <- if (methods::is(model, "metaverse")) model@robust$estimate else model$estimate
  se  <- if (methods::is(model, "metaverse")) model@robust$se else model$se

  # Compute how large confounding would need to be to shift estimate to null
  bias_factors <- c(1.5, 2, 3, 5)
  scenarios <- data.frame(
    confound_RR = bias_factors,
    adjusted_est = est / bias_factors,
    still_significant = (est / bias_factors) / se > qnorm(0.975)
  )
  scenarios
}

#' Interpret E-values
interpret_evalues <- function(e_val_est, e_val_ci) {
  msg <- paste0("E-value for estimate = ", round(e_val_est, 2),
                "; E-value for CI = ", round(e_val_ci, 2), ". ")
  if (e_val_ci > 3) {
    msg <- paste0(msg, "Strong evidence: unmeasured confounding would need to be very large to explain away the result.")
  } else if (e_val_ci > 1.5) {
    msg <- paste0(msg, "Moderate evidence against unmeasured confounding.")
  } else {
    msg <- paste0(msg, "Weak evidence: modest unmeasured confounding could explain away the result.")
  }
  msg
}

#' Outlier detection
outliers_complete <- function(model, control = list(), verbose = FALSE) {
  data <- model@data
  yi   <- data$yi
  vi   <- data$vi
  sei  <- sqrt(vi)
  n    <- length(yi)

  est <- model@robust$estimate
  se  <- model@robust$se

  # Standardized residuals
  zi <- (yi - est) / sei

  # Outlier thresholds
  outliers_2sd <- which(abs(zi) > 2)
  outliers_3sd <- which(abs(zi) > 3)

  # Grubbs-like test based on max |z|
  max_z <- max(abs(zi))
  grubbs_p <- 2 * n * pnorm(-max_z)  # Bonferroni-corrected

  list(
    standardized_residuals = zi,
    outliers_2sd = outliers_2sd,
    outliers_3sd = outliers_3sd,
    max_z        = max_z,
    grubbs_p     = grubbs_p,
    n_outliers   = length(outliers_2sd)
  )
}

#' Cumulative meta-analysis
cumulative_analysis <- function(model, control = list(), verbose = FALSE) {
  data <- model@data
  yi   <- data$yi
  vi   <- data$vi
  n    <- length(yi)

  cum_est <- cum_se <- numeric(n)
  for (i in seq_len(n)) {
    wi_i <- 1 / vi[1:i]
    mu_i <- sum(wi_i * yi[1:i]) / sum(wi_i)
    se_i <- sqrt(1 / sum(wi_i))
    cum_est[i] <- mu_i
    cum_se[i]  <- se_i
  }

  z_crit <- qnorm(0.975)
  data.frame(
    k        = seq_len(n),
    estimate = cum_est,
    se       = cum_se,
    ci_lower = cum_est - z_crit * cum_se,
    ci_upper = cum_est + z_crit * cum_se,
    p_value  = 2 * pnorm(-abs(cum_est / cum_se))
  )
}

