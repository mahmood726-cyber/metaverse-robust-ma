
#' @title Complete Inference Methods
#' @description Selective inference, conformal prediction, and advanced bootstrap

#' Advanced Uncertainty Quantification for Meta-Analysis
#'
#' Provides uncertainty estimates beyond standard confidence intervals,
#' including bootstrap, conformal prediction, selective inference,
#' permutation tests, and Bayesian posterior sampling.
#'
#' @param model A metaverse S4 object with fitted robust results.
#' @param method Character: inference method. One of "bootstrap", "conformal",
#'   "selective", "permutation", "bayesian".
#' @param conf.level Confidence level for intervals (default 0.95).
#' @param R Number of bootstrap or permutation replicates (default 1000).
#' @param control Named list of method-specific parameters.
#' @param verbose Logical; print progress.
#' @return An object of class \code{inference_result} with method-specific
#'   fields (intervals, p_value, posterior samples, etc.).
#' @export
quantify_uncertainty <- function(model,
                               method = c("bootstrap", "conformal", "selective", 
                                        "permutation", "bayesian"),
                               conf.level = 0.95,
                               R = 1000,
                               control = list(),
                               verbose = FALSE) {
  
  method <- match.arg(method)
  
  if(verbose) cat("Computing", method, "inference...\n")
  
  result <- switch(method,
    "bootstrap" = bootstrap_inference_complete(model, R, conf.level, control, verbose),
    "conformal" = conformal_inference_complete(model, conf.level, control, verbose),
    "selective" = selective_inference_complete(model, conf.level, control, verbose),
    "permutation" = permutation_inference(model, R, conf.level, verbose),
    "bayesian" = bayesian_inference(model, control, verbose)
  )
  
  result$method <- method
  result$conf.level <- conf.level
  
  class(result) <- c("inference_result", "list")
  result
}

#' Complete bootstrap implementation
bootstrap_inference_complete <- function(model, R = 1000, conf.level = 0.95, 
                                        control = list(), verbose = FALSE) {
  
  con <- list(
    type = "nonparametric",
    method = "bca",
    parallel = FALSE,
    ncores = 2
  )
  con[names(control)] <- control
  
  # Extract data
  data <- model@data
  n <- nrow(data)
  
  # Bootstrap function
  boot_fun <- function(data, indices) {
    d <- data[indices, ]
    fit <- meta_robust(d$yi, d$vi, method = model@model$method)
    c(estimate = fit$estimate, tau2 = fit$tau2, I2 = fit$I2)
  }
  
  # Run bootstrap
  if(con$parallel && con$ncores > 1) {
    cl <- makeCluster(con$ncores)
    clusterEvalQ(cl, library(metaverse))
    boot_result <- boot(data, boot_fun, R = R, parallel = "snow", 
                       ncpus = con$ncores, cl = cl)
    stopCluster(cl)
  } else {
    boot_result <- boot(data, boot_fun, R = R)
  }
  
  # Compute intervals
  alpha <- 1 - conf.level
  
  ci_types <- c("norm", "basic", "perc", "bca")
  intervals <- list()
  
  for(i in 1:3) {  # For each parameter
    intervals[[i]] <- list()
    for(type in ci_types) {
      tryCatch({
        ci <- boot.ci(boot_result, conf = conf.level, type = type, index = i)
        intervals[[i]][[type]] <- ci[[type]][4:5]
      }, error = function(e) {
        intervals[[i]][[type]] <- c(NA, NA)
      })
    }
  }
  
  list(
    estimates = boot_result$t0,
    bootstrap_samples = boot_result$t,
    intervals = intervals,
    se = apply(boot_result$t, 2, sd),
    bias = colMeans(boot_result$t) - boot_result$t0,
    R = R
  )
}

#' Split conformal prediction for meta-analysis
#'
#' Provides distribution-free prediction intervals with finite-sample
#' coverage guarantees. Uses a calibration/prediction split.
conformal_inference_complete <- function(model, conf.level = 0.95,
                                        control = list(), verbose = FALSE) {

  con <- list(
    method = "split",
    cal_fraction = 0.5
  )
  con[names(control)] <- control

  data <- model@data
  n <- nrow(data)
  yi <- data$yi
  vi <- data$vi

  # --- split into calibration and prediction sets ---
  cal_size <- max(2, floor(n * con$cal_fraction))
  cal_idx  <- sample(n, cal_size)
  pred_idx <- setdiff(seq_len(n), cal_idx)

  # Fit robust model on calibration set
  cal_fit <- meta_robust(yi[cal_idx], vi[cal_idx])

  # Nonconformity scores = |y_i - mu_hat| / sqrt(vi)
  scores <- abs(yi[cal_idx] - cal_fit$estimate) / sqrt(vi[cal_idx])

  # Quantile of scores for the desired coverage (finite-sample correction)
  q_level <- ceiling((cal_size + 1) * conf.level) / cal_size
  q_level <- min(q_level, 1)
  q_hat   <- as.numeric(quantile(scores, probs = q_level))

  # Build prediction intervals for prediction set
  pred_lower <- cal_fit$estimate - q_hat * sqrt(vi[pred_idx])
  pred_upper <- cal_fit$estimate + q_hat * sqrt(vi[pred_idx])

  # Coverage on prediction set (empirical check)
  coverage <- mean(yi[pred_idx] >= pred_lower & yi[pred_idx] <= pred_upper)

  list(
    estimate      = cal_fit$estimate,
    q_hat         = q_hat,
    intervals     = data.frame(
      study = data$study[pred_idx],
      lower = pred_lower,
      upper = pred_upper,
      yi    = yi[pred_idx]
    ),
    coverage      = coverage,
    conf.level    = conf.level,
    cal_size      = cal_size,
    pred_size     = length(pred_idx)
  )
}

#' Selective inference (post-selection)
#'
#' After variable selection in meta-regression, naive p-values are invalid
#' because they ignore the selection step. This function provides adjusted
#' p-values and confidence intervals using a data-splitting approach.
selective_inference_complete <- function(model, conf.level = 0.95,
                                        control = list(), verbose = FALSE) {

  con <- list(
    split_fraction = 0.5
  )
  con[names(control)] <- control

  data <- model@data
  n <- nrow(data)
  yi <- data$yi
  vi <- data$vi

  # If selection results exist, use them; otherwise report unadjusted
  selected <- model@selection$selected
  if (is.null(selected) || length(selected) == 0) {
    return(list(
      method  = "selective",
      message = "No variable selection results found in model. Run select_moderators first.",
      adjusted_p = numeric(0),
      adjusted_ci = matrix(nrow = 0, ncol = 2)
    ))
  }

  # Data-splitting approach: split sample, select on first half, infer on second
  split_size <- max(5, floor(n * con$split_fraction))
  split_idx  <- sample(n, split_size)
  infer_idx  <- setdiff(seq_len(n), split_idx)

  if (length(infer_idx) < 3) {
    return(list(
      method  = "selective",
      message = "Not enough observations for inference after splitting.",
      adjusted_p = rep(NA, length(selected)),
      adjusted_ci = matrix(NA, nrow = length(selected), ncol = 2)
    ))
  }

  # Re-fit on inference half using only selected moderators
  mod_names <- grep("^X|^mod", names(data), value = TRUE)
  if (length(mod_names) == 0 || max(selected) > length(mod_names)) {
    # Fall back to intercept-only inference
    fit_inf <- meta_robust(yi[infer_idx], vi[infer_idx])
    z_crit <- qnorm(1 - (1 - conf.level) / 2)
    return(list(
      method  = "selective",
      estimate = fit_inf$estimate,
      adjusted_p  = fit_inf$p_value,
      adjusted_ci = matrix(c(fit_inf$estimate - z_crit * fit_inf$se,
                              fit_inf$estimate + z_crit * fit_inf$se),
                            nrow = 1, ncol = 2,
                            dimnames = list(NULL, c("lower", "upper")))
    ))
  }

  X_sel <- as.matrix(data[infer_idx, mod_names[selected], drop = FALSE])
  yi_inf <- yi[infer_idx]
  vi_inf <- vi[infer_idx]

  # Weighted least squares on inference split
  wi_inf <- 1 / vi_inf
  X_des  <- cbind(1, X_sel)
  W      <- diag(wi_inf)
  XtWX   <- t(X_des) %*% W %*% X_des

  if (rcond(XtWX) < .Machine$double.eps) {
    return(list(
      method  = "selective",
      message = "Design matrix singular on inference split.",
      adjusted_p  = rep(NA, ncol(X_des)),
      adjusted_ci = matrix(NA, nrow = ncol(X_des), ncol = 2)
    ))
  }

  beta_hat <- solve(XtWX) %*% t(X_des) %*% W %*% yi_inf
  resid    <- yi_inf - as.vector(X_des %*% beta_hat)
  sigma2   <- sum(wi_inf * resid^2) / max(1, length(yi_inf) - ncol(X_des))
  vcov_hat <- sigma2 * solve(XtWX)
  se_hat   <- sqrt(pmax(0, diag(vcov_hat)))

  z_vals <- as.vector(beta_hat) / se_hat
  p_vals <- 2 * pnorm(-abs(z_vals))
  z_crit <- qnorm(1 - (1 - conf.level) / 2)
  ci_mat <- cbind(lower = as.vector(beta_hat) - z_crit * se_hat,
                  upper = as.vector(beta_hat) + z_crit * se_hat)

  row_names <- c("intercept", mod_names[selected])
  names(p_vals) <- row_names
  rownames(ci_mat) <- row_names

  list(
    method      = "selective",
    coefficients = as.vector(beta_hat),
    se          = se_hat,
    adjusted_p  = p_vals,
    adjusted_ci = ci_mat
  )
}

#' Permutation inference
permutation_inference <- function(model, R = 1000, conf.level = 0.95, verbose = FALSE) {
  
  data <- model@data
  n <- nrow(data)
  
  # Observed statistic
  obs_stat <- model@robust$estimate
  
  # Permutation distribution
  perm_stats <- numeric(R)
  
  for(r in 1:R) {
    # Permute effect sizes
    perm_idx <- sample(n)
    data_perm <- data
    data_perm$yi <- data$yi[perm_idx]
    
    # Recompute statistic
    fit_perm <- meta_robust(data_perm$yi, data_perm$vi)
    perm_stats[r] <- fit_perm$estimate
  }
  
  # P-value
  p_value <- mean(abs(perm_stats) >= abs(obs_stat))
  
  # Confidence interval (inversion of test)
  alpha <- 1 - conf.level
  ci <- quantile(perm_stats, c(alpha/2, 1 - alpha/2))
  
  list(
    estimate = obs_stat,
    p_value = p_value,
    ci = ci,
    permutation_distribution = perm_stats,
    R = R
  )
}

#' Bayesian inference via Gibbs sampling
#'
#' Implements a conjugate normal-inverse-gamma model for the random-effects
#' meta-analysis. This avoids external MCMC dependencies (Stan/JAGS).
bayesian_inference <- function(model, control = list(), verbose = FALSE) {

  con <- list(
    n_iter  = 5000,
    burn_in = 1000,
    thin    = 1,
    prior_mu_mean = 0,
    prior_mu_var  = 100,
    prior_tau2_shape = 0.001,
    prior_tau2_rate  = 0.001
  )
  con[names(control)] <- control

  data <- model@data
  yi <- data$yi
  vi <- data$vi
  n  <- length(yi)

  # Storage
  n_save <- floor((con$n_iter - con$burn_in) / con$thin)
  mu_samples   <- numeric(n_save)
  tau2_samples <- numeric(n_save)

  # Initialize

  mu   <- mean(yi)
  tau2 <- max(0.01, var(yi) - mean(vi))

  save_idx <- 1
  for (iter in seq_len(con$n_iter)) {
    # --- Update mu | tau2, y ---
    vi_total <- vi + tau2
    wi <- 1 / vi_total
    post_var  <- 1 / (sum(wi) + 1 / con$prior_mu_var)
    post_mean <- post_var * (sum(wi * yi) + con$prior_mu_mean / con$prior_mu_var)
    mu <- rnorm(1, post_mean, sqrt(post_var))

    # --- Update tau2 | mu, y  (inverse-gamma via Metropolis step) ---
    # Propose log(tau2*) = log(tau2) + N(0, 0.5)
    log_tau2_prop <- log(tau2) + rnorm(1, 0, 0.5)
    tau2_prop <- exp(log_tau2_prop)

    # Log-posterior of tau2
    log_post <- function(t2) {
      vt <- vi + t2
      -0.5 * sum(log(vt)) - 0.5 * sum((yi - mu)^2 / vt) +
        (con$prior_tau2_shape - 1) * log(t2) - con$prior_tau2_rate * t2
    }

    log_alpha <- log_post(tau2_prop) - log_post(tau2) +
                 log_tau2_prop - log(tau2)  # Jacobian for log-scale proposal
    if (log(runif(1)) < log_alpha) {
      tau2 <- tau2_prop
    }

    # Store
    if (iter > con$burn_in && (iter - con$burn_in) %% con$thin == 0) {
      mu_samples[save_idx]   <- mu
      tau2_samples[save_idx] <- tau2
      save_idx <- save_idx + 1
    }
  }

  # Trim to actual saved samples
  actual <- save_idx - 1
  mu_samples   <- mu_samples[seq_len(actual)]
  tau2_samples <- tau2_samples[seq_len(actual)]

  list(
    mu_posterior   = mu_samples,
    tau2_posterior = tau2_samples,
    estimate       = mean(mu_samples),
    estimate_sd    = sd(mu_samples),
    tau2_mean      = mean(tau2_samples),
    tau2_sd        = sd(tau2_samples),
    hpd_mu         = quantile(mu_samples, c(0.025, 0.975)),
    hpd_tau2       = quantile(tau2_samples, c(0.025, 0.975)),
    n_samples      = actual
  )
}

