
#' @title Comprehensive Robust Estimation Methods
#' @description Advanced robust estimators with contamination models
#' @import robustbase
#' @import MASS
#' @importFrom stats pnorm qnorm dnorm optim pchisq median mad quantile

#' @export
meta_robust <- function(yi, vi, data = NULL, 
                       mods = NULL,
                       method = c("MM", "M", "tau", "S", "ML", "REML"),
                       contamination = c("none", "t-mixture", "slash", "box-cox"),
                       weights = NULL,
                       control = list(),
                       verbose = FALSE) {
  
  # Argument matching
  method <- match.arg(method)
  contamination <- match.arg(contamination)
  
  # Extract data if provided
  if (!is.null(data)) {
    yi <- eval(substitute(yi), data, parent.frame())
    vi <- eval(substitute(vi), data, parent.frame())
    if(!is.null(substitute(mods))) {
      mods <- eval(substitute(mods), data, parent.frame())
    }
  }
  
  # Input validation
  check_inputs(yi, vi, mods, weights)
  
  # Setup control parameters
  con <- list(
    k = 4.685,           # Tukey biweight tuning
    tol = 1e-6,          # Convergence tolerance  
    maxiter = 1000,      # Maximum iterations
    efficiency = 0.95,   # MM-estimator efficiency
    df = 4,              # Degrees of freedom for t-mixture
    lambda = NULL,       # Box-Cox parameter
    init = "median"      # Initial estimate
  )
  con[names(control)] <- control
  
  n <- length(yi)
  wi <- if(!is.null(weights)) weights else 1/vi
  
  # Initial estimate
  init_val <- switch(con$init,
    "median" = median(yi),
    "weighted" = sum(wi * yi) / sum(wi),
    "trimmed" = mean(yi, trim = 0.1)
  )
  
  if(verbose) cat("Initial estimate:", init_val, "\n")
  
  # Choose main estimation method
  result <- switch(method,
    "MM" = mm_estimator(yi, vi, wi, mods, con, verbose),
    "M" = m_estimator(yi, vi, wi, mods, con, verbose),
    "tau" = tau_estimator(yi, vi, wi, mods, con, verbose),
    "S" = s_estimator(yi, vi, wi, mods, con, verbose),
    "ML" = ml_estimator(yi, vi, wi, mods, con, verbose),
    "REML" = reml_estimator(yi, vi, wi, mods, con, verbose)
  )
  
  # Apply contamination model if specified
  if (contamination != "none") {
    result <- apply_contamination_model(result, yi, vi, mods, 
                                       type = contamination, 
                                       control = con, 
                                       verbose = verbose)
  }
  
  # Add additional statistics
  result <- compute_statistics(result, yi, vi, wi, mods)
  
  # Add method information
  result$method <- method
  result$contamination <- contamination
  result$call <- match.call()
  
  class(result) <- c("meta_robust", "list")
  result
}

#' M-estimator with various psi functions
#' @export

#' M-estimator with various psi functions
#' @export
m_estimator <- function(yi, vi, wi, mods = NULL, control = list(), verbose = FALSE) {
  
  n <- length(yi)
  
  # Design matrix
  if(!is.null(mods)) {
    if(!is.matrix(mods)) mods <- as.matrix(mods)
    X <- cbind(1, mods)
    p <- ncol(X)
  } else {
    X <- matrix(1, n, 1)
    p <- 1
  }
  
  # Initial estimate via weighted least squares
  W <- diag(wi)
  beta_old <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% yi
  
  # Iteratively reweighted least squares
  for(iter in 1:control$maxiter) {
    
    # Fitted values and residuals
    fitted <- as.vector(X %*% beta_old)  # Ensure vector
    resid <- yi - fitted
    
    # Standardized residuals
    ri <- resid / sqrt(vi)
    
    # Robust weights (Tukey biweight)
    w_rob <- psi_tukey(ri, control$k) / ri
    w_rob[is.nan(w_rob) | is.infinite(w_rob)] <- 1
    
    # Combined weights - ensure proper dimensions
    w_combined <- wi * w_rob
    
    # Update estimate - use weighted least squares formula
    XtW <- t(X * w_combined)  # Efficient computation of t(X) %*% W
    XtWX <- XtW %*% X
    XtWy <- as.vector(XtW %*% yi)
    
    # Check if XtWX is singular
    if(rcond(XtWX) < .Machine$double.eps) {
      warning("Singular matrix in iteration ", iter)
      break
    }
    
    beta_new <- solve(XtWX, XtWy)
    
    # Check convergence
    if(max(abs(beta_new - beta_old)) < control$tol) {
      if(verbose) cat("M-estimator converged in", iter, "iterations\n")
      break
    }
    
    beta_old <- beta_new
  }
  
  # Calculate robust covariance matrix
  V <- robust_vcov(X, resid, vi, w_rob, type = "HC3")
  se <- sqrt(diag(V))
  
  # Prepare output
  list(
    coefficients = as.vector(beta_new),
    estimate = beta_new[1],
    se = se[1],
    vcov = V,
    fitted = as.vector(X %*% beta_new),
    residuals = resid,
    weights = w_rob,
    iterations = iter,
    converged = (iter < control$maxiter)
  )
}
mm_estimator <- function(yi, vi, wi, mods = NULL, control = list(), verbose = FALSE) {
  
  if(verbose) cat("Step 1: Computing S-estimate...\n")
  
  # Step 1: S-estimation for high breakdown point
  s_est <- s_estimator(yi, vi, wi, mods, control, verbose = FALSE)
  
  if(verbose) cat("Step 2: Computing M-estimate...\n")
  
  # Step 2: M-estimation for efficiency
  control$init <- s_est$estimate  # Use S-estimate as starting value
  control$k <- tuning_constant(control$efficiency)
  m_est <- m_estimator(yi, vi, wi, mods, control, verbose = FALSE)
  
  # Combine results
  list(
    coefficients = m_est$coefficients,
    estimate = m_est$estimate,
    se = m_est$se,
    vcov = m_est$vcov,
    fitted = m_est$fitted,
    residuals = m_est$residuals,
    weights = m_est$weights,
    s_estimate = s_est$estimate,
    s_scale = s_est$scale,
    iterations = s_est$iterations + m_est$iterations,
    converged = s_est$converged & m_est$converged
  )
}

#' S-estimator for high breakdown point
#' @export
s_estimator <- function(yi, vi, wi, mods = NULL, control = list(), verbose = FALSE) {
  
  n <- length(yi)
  
  # Design matrix
  if(!is.null(mods)) {
    if(!is.matrix(mods)) mods <- as.matrix(mods)
    X <- cbind(1, mods)
    p <- ncol(X)
  } else {
    X <- matrix(1, n, 1)
    p <- 1
  }
  
  # Use median as initial estimate for intercept model
  if(p == 1) {
    theta_init <- median(yi)
    ri <- (yi - theta_init) / sqrt(vi)
    scale <- median(abs(ri)) * 1.4826
    
    return(list(
      coefficients = theta_init,
      estimate = theta_init,
      scale = scale,
      iterations = 1,
      converged = TRUE
    ))
  }
  
  # For models with covariates, use robust regression
  rob_fit <- robustbase::lmrob(yi ~ X - 1, weights = wi)
  
  list(
    coefficients = coef(rob_fit),
    estimate = coef(rob_fit)[1],
    scale = rob_fit$scale,
    iterations = rob_fit$control$iterations,
    converged = rob_fit$converged
  )
}

#' Tau-squared estimator with robust methods
#' @export
tau_estimator <- function(yi, vi, wi, mods = NULL, control = list(), verbose = FALSE) {
  
  n <- length(yi)
  
  # Get robust point estimate first
  m_est <- m_estimator(yi, vi, wi, mods, control, verbose = FALSE)
  theta_robust <- m_est$estimate
  
  # Q statistic with robust estimate
  Q_robust <- sum(wi * (yi - theta_robust)^2)
  
  # Multiple tau2 estimators
  tau2_methods <- list()
  
  # DerSimonian-Laird (robust version)
  c_dl <- sum(wi) - sum(wi^2) / sum(wi)
  tau2_methods$DL <- max(0, (Q_robust - (n - 1)) / c_dl)
  
  # Paule-Mandel estimator
  tau2_methods$PM <- tau2_pm(yi, vi)
  
  # REML estimator
  tau2_methods$REML <- tau2_reml(yi, vi)
  
  # Empirical Bayes
  tau2_methods$EB <- tau2_eb(yi, vi)
  
  # Choose primary estimate
  tau2 <- tau2_methods[[control$tau_method %||% "DL"]]
  
  # Standard error of pooled estimate with tau2
  vi_total <- vi + tau2
  wi_new <- 1 / vi_total
  theta_final <- sum(wi_new * yi) / sum(wi_new)
  se <- sqrt(1 / sum(wi_new))
  
  # Heterogeneity statistics
  I2 <- max(0, (Q_robust - (n - 1)) / Q_robust * 100)
  H2 <- Q_robust / (n - 1)
  
  # Prediction interval
  t_crit <- qt(0.975, n - 1)
  pi_lower <- theta_final - t_crit * sqrt(tau2 + se^2)
  pi_upper <- theta_final + t_crit * sqrt(tau2 + se^2)
  
  list(
    coefficients = theta_final,
    estimate = theta_final,
    se = se,
    tau2 = tau2,
    tau = sqrt(tau2),
    tau2_methods = tau2_methods,
    Q = Q_robust,
    Q_p = pchisq(Q_robust, df = n - 1, lower.tail = FALSE),
    I2 = I2,
    H2 = H2,
    prediction_interval = c(pi_lower, pi_upper),
    iterations = m_est$iterations,
    converged = m_est$converged
  )
}

#' Contamination models
#' @export
apply_contamination_model <- function(result, yi, vi, mods = NULL,
                                     type = "t-mixture", 
                                     control = list(), 
                                     verbose = FALSE) {
  
  if(verbose) cat("Applying", type, "contamination model...\n")
  
  contam_result <- switch(type,
    "t-mixture" = fit_t_mixture(yi, vi, df = control$df %||% 4, verbose),
    "slash" = fit_slash_mixture(yi, vi, verbose),
    "box-cox" = fit_boxcox_transform(yi, vi, lambda = control$lambda, verbose),
    result  # Return unchanged if "none"
  )
  
  # Merge contamination results
  result$contamination_fit <- contam_result
  result$outlier_prob <- contam_result$outlier_prob
  result$transformed_data <- contam_result$transformed
  
  result
}

#' t-mixture contamination model
fit_t_mixture <- function(yi, vi, df = 4, verbose = FALSE) {
  
  n <- length(yi)
  
  # EM algorithm for t-mixture
  # Initial values
  mu <- mean(yi)
  sigma2 <- var(yi)
  weights <- rep(1, n)
  
  for(iter in 1:100) {
    mu_old <- mu
    
    # E-step: compute weights
    for(i in 1:n) {
      weights[i] <- (df + 1) / (df + ((yi[i] - mu)^2) / (sigma2 + vi[i]))
    }
    
    # M-step: update parameters
    mu <- sum(weights * yi / (sigma2 + vi)) / sum(weights / (sigma2 + vi))
    sigma2 <- sum(weights * (yi - mu)^2) / n
    
    # Check convergence
    if(abs(mu - mu_old) < 1e-6) break
  }
  
  # Compute outlier probabilities
  outlier_prob <- 1 - weights
  
  list(
    estimate = mu,
    sigma2 = sigma2,
    df = df,
    weights = weights,
    outlier_prob = outlier_prob,
    iterations = iter
  )
}

#' Slash distribution mixture
fit_slash_mixture <- function(yi, vi, verbose = FALSE) {
  
  # Simplified slash distribution fitting
  # The slash distribution is t/U where U ~ Uniform(0,1)
  
  n <- length(yi)
  
  # Use robust location and scale
  location <- median(yi)
  scale <- mad(yi)
  
  # Compute weights based on slash density approximation
  weights <- numeric(n)
  for(i in 1:n) {
    z <- (yi[i] - location) / scale
    if(abs(z) < 1) {
      weights[i] <- 1
    } else {
      weights[i] <- 1 / abs(z)
    }
  }
  
  # Recompute estimate with weights
  estimate <- sum(weights * yi) / sum(weights)
  
  list(
    estimate = estimate,
    location = location,
    scale = scale,
    weights = weights,
    outlier_prob = 1 - weights
  )
}

#' Box-Cox transformation
fit_boxcox_transform <- function(yi, vi, lambda = NULL, verbose = FALSE) {
  
  # Shift data to be positive if necessary
  shift <- 0
  if(any(yi <= 0)) {
    shift <- abs(min(yi)) + 1
    yi_shifted <- yi + shift
  } else {
    yi_shifted <- yi
  }
  
  # Estimate lambda if not provided
  if(is.null(lambda)) {
    lambda <- optimize_boxcox_lambda(yi_shifted, vi)
  }
  
  # Apply transformation
  if(abs(lambda) < 0.001) {
    yi_transformed <- log(yi_shifted)
  } else {
    yi_transformed <- (yi_shifted^lambda - 1) / lambda
  }
  
  # Fit model on transformed scale
  wi <- 1/vi
  estimate_transformed <- sum(wi * yi_transformed) / sum(wi)
  
  # Back-transform
  if(abs(lambda) < 0.001) {
    estimate <- exp(estimate_transformed) - shift
  } else {
    estimate <- (lambda * estimate_transformed + 1)^(1/lambda) - shift
  }
  
  list(
    estimate = estimate,
    lambda = lambda,
    shift = shift,
    transformed = yi_transformed,
    outlier_prob = rep(0, length(yi))  # Placeholder
  )
}

#' Optimize Box-Cox lambda parameter
optimize_boxcox_lambda <- function(yi, vi) {
  
  loglik <- function(lambda) {
    if(abs(lambda) < 0.001) {
      yt <- log(yi)
    } else {
      yt <- (yi^lambda - 1) / lambda
    }
    
    wi <- 1/vi
    mu <- sum(wi * yt) / sum(wi)
    -sum(dnorm(yt, mu, sqrt(vi), log = TRUE)) - (lambda - 1) * sum(log(yi))
  }
  
  result <- optimize(loglik, interval = c(-2, 2))
  result$minimum
}

#' Additional tau-squared estimators
tau2_pm <- function(yi, vi) {
  # Paule-Mandel estimator
  Q <- function(tau2) {
    wi <- 1/(vi + tau2)
    theta <- sum(wi * yi) / sum(wi)
    sum(wi * (yi - theta)^2) - (length(yi) - 1)
  }
  
  if(Q(0) <= 0) return(0)
  
  result <- uniroot(Q, interval = c(0, var(yi)), extendInt = "yes")
  max(0, result$root)
}

tau2_reml <- function(yi, vi) {
  # Simplified REML estimator
  n <- length(yi)
  
  reml_loglik <- function(tau2) {
    vi_total <- vi + tau2
    wi <- 1/vi_total
    theta <- sum(wi * yi) / sum(wi)
    
    -0.5 * (sum(log(vi_total)) + sum(wi * (yi - theta)^2) + log(sum(wi)))
  }
  
  upper <- max(var(yi), 0.01)  # Ensure valid interval even for homogeneous data
  result <- optimize(reml_loglik, interval = c(0, upper), maximum = TRUE)
  max(0, result$maximum)
}

tau2_eb <- function(yi, vi) {
  # Empirical Bayes estimator
  n <- length(yi)
  
  # Method of moments
  mu <- mean(yi)
  s2 <- var(yi)
  v_avg <- mean(vi)
  
  max(0, s2 - v_avg)
}

#' ML and REML estimators
ml_estimator <- function(yi, vi, wi, mods = NULL, control = list(), verbose = FALSE) {
  # Maximum likelihood estimation
  # (Simplified - would use optimization in practice)
  tau_est <- tau_estimator(yi, vi, wi, mods, control, verbose)
  tau_est$method <- "ML"
  tau_est
}

reml_estimator <- function(yi, vi, wi, mods = NULL, control = list(), verbose = FALSE) {
  # Restricted maximum likelihood
  tau_est <- tau_estimator(yi, vi, wi, mods, control, verbose)
  tau_est$tau2 <- tau2_reml(yi, vi)
  tau_est$method <- "REML"
  tau_est
}

#' Compute additional statistics
compute_statistics <- function(result, yi, vi, wi, mods = NULL) {
  
  n <- length(yi)
  
  # Confidence intervals
  z_crit <- qnorm(0.975)
  result$ci.lower <- result$estimate - z_crit * result$se
  result$ci.upper <- result$estimate + z_crit * result$se
  
  # P-value
  z_stat <- result$estimate / result$se
  result$p_value <- 2 * pnorm(-abs(z_stat))
  
  # Effect size measures
  result$z_value <- z_stat
  
  # Influence measures
  if(!is.null(result$weights)) {
    result$influential <- which(result$weights < 0.5)
  }
  
  result
}

#' Robust variance-covariance matrix
robust_vcov <- function(X, resid, vi, weights, type = "HC3") {
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Hat matrix
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  h <- diag(H)
  
  # Compute robust variance based on type
  if(type == "HC0") {
    # White standard errors
    meat <- t(X) %*% diag(resid^2 * weights) %*% X
  } else if(type == "HC1") {
    # With finite sample correction
    meat <- (n/(n-p)) * t(X) %*% diag(resid^2 * weights) %*% X
  } else if(type == "HC2") {
    # Weighted by leverage
    meat <- t(X) %*% diag(resid^2 / (1 - h) * weights) %*% X
  } else if(type == "HC3") {
    # Default - most conservative
    meat <- t(X) %*% diag(resid^2 / (1 - h)^2 * weights) %*% X
  }
  
  bread <- solve(t(X) %*% X)
  vcov <- bread %*% meat %*% bread
  
  vcov
}

#' Helper functions
psi_tukey <- function(x, k = 4.685) {
  x * (1 - (x/k)^2)^2 * (abs(x) <= k)
}

psi_huber <- function(x, k = 1.345) {
  pmin(pmax(-k, x), k)
}

psi_hampel <- function(x, a = 2, b = 4, c = 8) {
  ifelse(abs(x) <= a, x,
         ifelse(abs(x) <= b, a * sign(x),
                ifelse(abs(x) <= c, a * sign(x) * (c - abs(x))/(c - b), 0)))
}

tuning_constant <- function(efficiency) {
  # Tukey biweight tuning constants for given efficiency
  if(efficiency >= 0.95) return(4.685)
  if(efficiency >= 0.90) return(4.0)
  if(efficiency >= 0.85) return(3.5)
  if(efficiency >= 0.80) return(3.0)
  return(2.5)
}

check_inputs <- function(yi, vi, mods = NULL, weights = NULL) {
  
  if(length(yi) != length(vi)) {
    stop("yi and vi must have the same length")
  }
  
  if(any(vi <= 0, na.rm = TRUE)) {
    stop("All variances must be positive")
  }
  
  if(!is.null(mods)) {
    if(is.vector(mods)) mods <- as.matrix(mods)
    if(nrow(mods) != length(yi)) {
      stop("Number of moderator rows must equal length of yi")
    }
  }
  
  if(!is.null(weights)) {
    if(length(weights) != length(yi)) {
      stop("weights must have the same length as yi")
    }
    if(any(weights < 0, na.rm = TRUE)) {
      stop("All weights must be non-negative")
    }
  }
  
  invisible(TRUE)
}

#' Null coalescing operator
`%||%` <- function(x, y) if(is.null(x)) y else x

#' @export
print.meta_robust <- function(x, digits = 3, ...) {
  cat("\nRobust Meta-Analysis Results\n")
  cat("Method:", x$method, "\n")
  if(x$contamination != "none") {
    cat("Contamination model:", x$contamination, "\n")
  }
  cat("\nEstimate:", round(x$estimate, digits), "\n")
  cat("Std. Error:", round(x$se, digits), "\n")
  cat("95% CI: [", round(x$ci.lower, digits), ", ", 
      round(x$ci.upper, digits), "]\n", sep = "")
  cat("p-value:", format.pval(x$p_value, digits = digits), "\n")
  
  if(!is.null(x$tau2)) {
    cat("\nHeterogeneity:\n")
    cat("  tau^2 =", round(x$tau2, digits + 1), "\n")
    cat("  I^2 =", round(x$I2, 1), "%\n")
  }
  
  invisible(x)
}

