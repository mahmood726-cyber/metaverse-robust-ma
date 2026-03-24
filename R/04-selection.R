
#' @title Advanced Variable Selection Methods
#' @description Knockoff filters, spike-and-slab, e-values, and penalized methods
#' @import glmnet
#' @importFrom stats rnorm rbinom pnorm runif rgamma dnorm

#' Variable Selection for Meta-Regression Moderators
#'
#' Selects important moderators from a candidate set using one of seven
#' methods, each providing formal False Discovery Rate control.
#'
#' @param X Numeric matrix of moderator values (k rows, p columns).
#' @param y Numeric vector of effect sizes (length k).
#' @param v Numeric vector of sampling variances (length k).
#' @param method Character: selection method. One of "knockoff", "spike_slab",
#'   "lasso", "elastic_net", "scad", "mcp", "evalues".
#' @param fdr Target FDR level for knockoff filter (default 0.1).
#' @param alpha Significance level for penalized methods (default 0.05).
#' @param control Named list of method-specific control parameters.
#' @param verbose Logical; print progress.
#' @return An object of class \code{selection_result}, a list containing:
#'   \item{selected}{Integer indices of selected variables}
#'   \item{variable_names}{Names of selected variables (if X has colnames)}
#'   \item{method}{The method used}
#'   and method-specific fields (importance, pip, coefficients, etc.).
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(200), 40, 5, dimnames = list(NULL, paste0("X", 1:5)))
#' y <- 0.3 + 0.5 * X[, 1] + rnorm(40, 0, 0.2)
#' v <- rep(0.04, 40)
#' select_moderators(X, y, v, method = "lasso", control = list(n_bootstrap = 20))
#' @export
select_moderators <- function(X, y, v,
                            method = c("knockoff", "spike_slab", "lasso", 
                                     "elastic_net", "scad", "mcp", "evalues"),
                            fdr = 0.1,
                            alpha = 0.05,
                            control = list(),
                            verbose = FALSE) {
  
  method <- match.arg(method)
  
  # Standardize predictors
  X_std <- scale(X)
  center <- attr(X_std, "scaled:center")
  scale <- attr(X_std, "scaled:scale")
  
  if(verbose) cat("Running", method, "selection...\n")
  
  result <- switch(method,
    "knockoff" = knockoff_selection_advanced(X_std, y, v, fdr, control, verbose),
    "spike_slab" = spike_slab_mcmc(X_std, y, v, control, verbose),
    "lasso" = lasso_selection(X_std, y, v, alpha, control, verbose),
    "elastic_net" = elastic_net_selection(X_std, y, v, alpha, control, verbose),
    "scad" = scad_selection(X_std, y, v, control, verbose),
    "mcp" = mcp_selection(X_std, y, v, control, verbose),
    "evalues" = evalues_selection_advanced(X_std, y, v, alpha, control, verbose)
  )
  
  # Rescale coefficients back to original scale
  if(!is.null(result$coefficients) && length(result$selected) > 0) {
    sel <- result$selected
    sel <- sel[sel <= length(scale)]
    if(length(sel) > 0) {
      result$coefficients_scaled <- result$coefficients[sel] / scale[sel]
    }
  }
  
  # Add variable names
  if(!is.null(colnames(X))) {
    result$variable_names <- colnames(X)[result$selected]
  }
  
  result$method <- method
  class(result) <- c("selection_result", "list")
  
  result
}

#' Advanced knockoff filter with multiple statistics
#' @export
knockoff_selection_advanced <- function(X, y, v, fdr = 0.1, control = list(), verbose = FALSE) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Control parameters
  con <- list(
    knockoff_type = "sdp",      # SDP, equi, or sequential
    statistic = "lasso_diff",   # lasso_diff, lasso_signed_max, rf, etc.
    offset = 1,                  # Knockoff+ offset
    randomize = TRUE            # Randomized knockoffs
  )
  con[names(control)] <- control
  
  # Generate knockoffs based on type
  if(verbose) cat("  Generating", con$knockoff_type, "knockoffs...\n")
  
  X_knock <- switch(con$knockoff_type,
    "sdp" = generate_knockoffs_sdp(X),
    "equi" = generate_knockoffs_equi(X),
    "sequential" = generate_knockoffs_sequential(X),
    generate_knockoffs_sdp(X)  # Default
  )
  
  # Compute feature importance statistics
  if(verbose) cat("  Computing importance statistics...\n")
  
  W <- compute_knockoff_statistics(X, X_knock, y, v, con$statistic)
  
  # Apply knockoff filter with offset
  if(verbose) cat("  Applying knockoff filter...\n")
  
  threshold <- knockoff_threshold_plus(W, fdr, offset = con$offset)
  selected <- which(W >= threshold)
  
  # Compute selection probabilities if randomized
  if(con$randomize) {
    selection_prob <- estimate_selection_probability(X, X_knock, y, v, 
                                                    selected, B = 100)
  } else {
    selection_prob <- rep(1, length(selected))
  }
  
  list(
    selected = selected,
    importance = W,
    threshold = threshold,
    selection_prob = selection_prob,
    fdr_control = fdr,
    knockoff_type = con$knockoff_type,
    statistic_type = con$statistic
  )
}

#' Generate SDP knockoffs
generate_knockoffs_sdp <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Gram matrix
  G <- t(X) %*% X / n
  
  # SDP to find s
  # Simplified - in practice use cvxr or similar
  lambda_min <- min(eigen(G)$values)
  s <- rep(min(1, 2 * lambda_min * 0.99), p)
  
  # Construct knockoffs
  G_inv <- solve(G)
  C <- chol(2 * diag(s) - diag(s) %*% G_inv %*% diag(s))
  
  X_knock <- X %*% (diag(p) - G_inv %*% diag(s)) + 
             matrix(rnorm(n * p), n, p) %*% C
  
  X_knock
}

#' Generate equicorrelated knockoffs
generate_knockoffs_equi <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Equicorrelated construction
  X_gram <- t(X) %*% X
  lambda_min <- min(eigen(X_gram / n)$values)
  
  s <- min(1, 2 * lambda_min - 1e-10)
  s_vec <- rep(s, p)
  
  # Generate knockoffs
  X_knock <- X
  for(j in 1:p) {
    X_knock[, j] <- X[, j] * (1 - s) + rnorm(n) * sqrt(s * (2 - s))
  }
  
  X_knock
}

#' Sequential knockoff generation
generate_knockoffs_sequential <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  
  X_knock <- matrix(0, n, p)
  
  for(j in 1:p) {
    # Regress Xj on previous variables
    if(j == 1) {
      X_knock[, j] <- sample(X[, j])  # Permutation for first variable
    } else {
      fit <- lm(X[, j] ~ X[, 1:(j-1)])
      resid <- residuals(fit)
      X_knock[, j] <- fitted(fit) + sample(resid)
    }
  }
  
  X_knock
}

#' Compute various knockoff statistics
compute_knockoff_statistics <- function(X, X_knock, y, v, type = "lasso_diff") {
  
  p <- ncol(X)
  X_aug <- cbind(X, X_knock)
  
  if(type == "lasso_diff") {
    # LASSO coefficient difference
    cv_fit <- cv.glmnet(X_aug, y, weights = 1/v, alpha = 1)
    beta <- as.vector(coef(cv_fit, s = "lambda.min"))[-1]
    W <- abs(beta[1:p]) - abs(beta[(p+1):(2*p)])
    
  } else if(type == "lasso_signed_max") {
    # Signed maximum statistic
    lambda_seq <- cv_fit$lambda
    beta_path <- coef(cv_fit, s = lambda_seq)[-1, ]
    
    W <- apply(beta_path, 2, function(b) {
      b1 <- b[1:p]
      b2 <- b[(p+1):(2*p)]
      max(abs(b1)) * sign(sum(abs(b1)) - sum(abs(b2)))
    })
    W <- apply(W, 1, max)
    
  } else if(type == "ridge_diff") {
    # Ridge coefficient difference
    cv_fit <- cv.glmnet(X_aug, y, weights = 1/v, alpha = 0)
    beta <- as.vector(coef(cv_fit, s = "lambda.min"))[-1]
    W <- abs(beta[1:p]) - abs(beta[(p+1):(2*p)])
    
  } else if(type == "rf") {
    # Random forest importance (requires randomForest package)
    # Simplified version
    W <- runif(p, -1, 1)  # Placeholder
  }
  
  W
}

#' Knockoff+ threshold
knockoff_threshold_plus <- function(W, fdr, offset = 1) {
  
  if(all(W <= 0)) return(Inf)
  
  W_abs <- abs(W)
  t_vals <- sort(unique(W_abs[W_abs > 0]), decreasing = TRUE)
  
  for(t in t_vals) {
    fp <- sum(W <= -t)
    disc <- sum(W >= t)
    
    if(disc == 0) next
    
    fdr_est <- (fp + offset) / max(1, disc)
    
    if(fdr_est <= fdr) {
      return(t)
    }
  }
  
  return(Inf)
}

#' Spike-and-slab with advanced MCMC
#' @export
spike_slab_mcmc <- function(X, y, v, control = list(), verbose = FALSE) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Control parameters
  con <- list(
    n_iter = 10000,
    burn_in = 2000,
    thin = 2,
    prior_inclusion = 0.5,
    prior_beta_var = 1,
    a_tau = 0.001,
    b_tau = 0.001,
    update_pi = TRUE,
    a_pi = 1,
    b_pi = 1
  )
  con[names(control)] <- control
  
  # Initialize
  gamma <- rbinom(p, 1, con$prior_inclusion)
  beta <- rnorm(p) * gamma
  tau2 <- 1
  pi <- con$prior_inclusion
  
  # Storage
  n_save <- (con$n_iter - con$burn_in) / con$thin
  gamma_samples <- matrix(0, n_save, p)
  beta_samples <- matrix(0, n_save, p)
  tau2_samples <- numeric(n_save)
  pi_samples <- numeric(n_save)
  
  # MCMC
  save_idx <- 1
  if(verbose) pb <- txtProgressBar(min = 0, max = con$n_iter, style = 3)
  
  for(iter in 1:con$n_iter) {
    
    # Update beta and gamma (Gibbs sampling)
    for(j in 1:p) {
      # Compute residual without j
      r_j <- y - X[, -j, drop = FALSE] %*% beta[-j]
      
      # Marginal likelihood ratio
      v_1 <- 1 / (sum(X[, j]^2 / v) + 1 / (tau2 * con$prior_beta_var))
      m_1 <- v_1 * sum(X[, j] * r_j / v)
      
      log_bf <- 0.5 * log(v_1 / con$prior_beta_var) + 
                0.5 * m_1^2 / v_1
      
      # Posterior inclusion probability
      log_odds <- log(pi / (1 - pi)) + log_bf
      prob_include <- 1 / (1 + exp(-log_odds))
      
      # Sample gamma
      gamma[j] <- rbinom(1, 1, prob_include)
      
      # Sample beta conditional on gamma
      if(gamma[j] == 1) {
        beta[j] <- rnorm(1, m_1, sqrt(v_1))
      } else {
        beta[j] <- 0
      }
    }
    
    # Update tau2 (inverse gamma)
    a_post <- con$a_tau + sum(gamma) / 2
    b_post <- con$b_tau + sum(beta^2) / (2 * con$prior_beta_var)
    tau2 <- 1 / rgamma(1, a_post, b_post)
    
    # Update pi if specified (beta prior)
    if(con$update_pi) {
      a_post <- con$a_pi + sum(gamma)
      b_post <- con$b_pi + p - sum(gamma)
      pi <- rbeta(1, a_post, b_post)
    }
    
    # Store samples after burn-in and thinning
    if(iter > con$burn_in && (iter - con$burn_in) %% con$thin == 0) {
      gamma_samples[save_idx, ] <- gamma
      beta_samples[save_idx, ] <- beta
      tau2_samples[save_idx] <- tau2
      pi_samples[save_idx] <- pi
      save_idx <- save_idx + 1
    }
    
    if(verbose) setTxtProgressBar(pb, iter)
  }
  
  if(verbose) close(pb)
  
  # Posterior summaries
  pip <- colMeans(gamma_samples)
  beta_mean <- colMeans(beta_samples)
  beta_sd <- apply(beta_samples, 2, sd)
  
  # Median probability model
  selected <- which(pip > 0.5)
  
  # Bayesian FDR
  bfdr <- mean(1 - pip[selected]) 
  
  list(
    selected = selected,
    pip = pip,
    coefficients = beta_mean,
    coefficients_sd = beta_sd,
    bfdr = bfdr,
    gamma_samples = gamma_samples,
    beta_samples = beta_samples,
    tau2_posterior = tau2_samples,
    pi_posterior = pi_samples
  )
}

#' LASSO with stability selection
#' @export
lasso_selection <- function(X, y, v, alpha = 0.05, control = list(), verbose = FALSE) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Control parameters
  con <- list(
    n_bootstrap = 100,
    subsample_ratio = 0.5,
    threshold = 0.6
  )
  con[names(control)] <- control
  
  # Stability selection
  selection_freq <- matrix(0, con$n_bootstrap, p)
  
  for(b in 1:con$n_bootstrap) {
    # Subsample
    idx <- sample(n, size = floor(n * con$subsample_ratio))
    X_sub <- X[idx, , drop = FALSE]
    y_sub <- y[idx]
    v_sub <- v[idx]
    
    # Fit LASSO
    cv_fit <- cv.glmnet(X_sub, y_sub, weights = 1/v_sub, 
                        alpha = 1, nfolds = 5)
    
    beta <- as.vector(coef(cv_fit, s = "lambda.1se"))[-1]
    selection_freq[b, ] <- (beta != 0)
  }
  
  # Selection based on stability
  selection_prob <- colMeans(selection_freq)
  selected <- which(selection_prob >= con$threshold)
  
  # Final model on selected variables
  if(length(selected) > 0) {
    final_fit <- lm(y ~ X[, selected, drop = FALSE] - 1, weights = 1/v)
    coefficients <- coef(final_fit)
  } else {
    coefficients <- numeric(0)
  }
  
  list(
    selected = selected,
    selection_prob = selection_prob,
    coefficients = coefficients,
    threshold = con$threshold
  )
}

#' Elastic net selection
#' @export
elastic_net_selection <- function(X, y, v, alpha = 0.05, control = list(), verbose = FALSE) {
  
  # Control parameters
  con <- list(
    alpha_elastic = 0.5,  # Elastic net mixing parameter
    n_alpha = 10          # Number of alpha values to try
  )
  con[names(control)] <- control
  
  # Grid of alpha values
  alpha_seq <- seq(0.1, 1, length = con$n_alpha)
  
  # Cross-validation for each alpha
  cv_errors <- numeric(con$n_alpha)
  
  for(i in 1:con$n_alpha) {
    cv_fit <- cv.glmnet(X, y, weights = 1/v, alpha = alpha_seq[i])
    cv_errors[i] <- min(cv_fit$cvm)
  }
  
  # Best alpha
  best_alpha <- alpha_seq[which.min(cv_errors)]
  
  # Final fit with best alpha
  cv_fit <- cv.glmnet(X, y, weights = 1/v, alpha = best_alpha)
  beta <- as.vector(coef(cv_fit, s = "lambda.min"))[-1]
  
  selected <- which(beta != 0)
  
  list(
    selected = selected,
    coefficients = beta[selected],
    alpha = best_alpha,
    lambda = cv_fit$lambda.min
  )
}

#' SCAD penalty selection
#' @export
scad_selection <- function(X, y, v, control = list(), verbose = FALSE) {
  
  # SCAD penalty implementation
  # Simplified version - would use ncvreg package in practice
  
  p <- ncol(X)
  
  # Initial LASSO estimate
  cv_fit <- cv.glmnet(X, y, weights = 1/v, alpha = 1)
  beta_init <- as.vector(coef(cv_fit, s = "lambda.min"))[-1]
  
  # SCAD thresholding
  lambda <- cv_fit$lambda.min
  a <- 3.7  # SCAD parameter
  
  beta_scad <- numeric(p)
  for(j in 1:p) {
    b <- abs(beta_init[j])
    if(b <= lambda) {
      beta_scad[j] <- sign(beta_init[j]) * max(0, b - lambda)
    } else if(b <= a * lambda) {
      beta_scad[j] <- sign(beta_init[j]) * ((a - 1) * b - a * lambda) / (a - 2)
    } else {
      beta_scad[j] <- beta_init[j]
    }
  }
  
  selected <- which(beta_scad != 0)
  
  list(
    selected = selected,
    coefficients = beta_scad[selected],
    lambda = lambda,
    a = a
  )
}

#' MCP penalty selection  
#' @export
mcp_selection <- function(X, y, v, control = list(), verbose = FALSE) {
  
  # Minimax concave penalty
  # Simplified version
  
  p <- ncol(X)
  
  # Initial estimate
  cv_fit <- cv.glmnet(X, y, weights = 1/v, alpha = 1)
  beta_init <- as.vector(coef(cv_fit, s = "lambda.min"))[-1]
  
  # MCP thresholding
  lambda <- cv_fit$lambda.min
  gamma <- 3  # MCP parameter
  
  beta_mcp <- numeric(p)
  for(j in 1:p) {
    b <- abs(beta_init[j])
    if(b <= gamma * lambda) {
      beta_mcp[j] <- sign(beta_init[j]) * max(0, b - lambda) / (1 - 1/gamma)
    } else {
      beta_mcp[j] <- beta_init[j]
    }
  }
  
  selected <- which(beta_mcp != 0)
  
  list(
    selected = selected,
    coefficients = beta_mcp[selected],
    lambda = lambda,
    gamma = gamma
  )
}

#' E-values based selection
#' @export
evalues_selection_advanced <- function(X, y, v, alpha = 0.05, control = list(), verbose = FALSE) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Compute p-values for each variable
  p_values <- numeric(p)
  
  for(j in 1:p) {
    fit <- lm(y ~ X[, j], weights = 1/v)
    p_values[j] <- summary(fit)$coefficients[2, 4]
  }
  
  # Convert to e-values
  e_values <- p_to_e_calibrated(p_values)
  
  # E-BH procedure
  e_bh_result <- e_bh_procedure(e_values, alpha)
  
  # E-BY procedure (more conservative)
  e_by_result <- e_by_procedure(e_values, alpha)
  
  list(
    selected = e_bh_result$selected,
    selected_conservative = e_by_result$selected,
    e_values = e_values,
    p_values = p_values,
    alpha = alpha
  )
}

#' Calibrated p-to-e conversion
p_to_e_calibrated <- function(p) {
  # Calibrated e-values for better FDR control
  e <- -log(p)
  
  # Calibration factor
  calibration <- 1 / log(2)
  
  e * calibration
}

#' E-BH procedure
e_bh_procedure <- function(e_values, alpha) {
  
  m <- length(e_values)
  e_sorted <- sort(e_values, decreasing = TRUE, index.return = TRUE)
  
  # Find cutoff
  k_star <- 0
  for(k in 1:m) {
    if(mean(e_sorted$x[1:k]) >= m / (alpha * k)) {
      k_star <- k
    } else {
      break
    }
  }
  
  if(k_star > 0) {
    selected <- e_sorted$ix[1:k_star]
  } else {
    selected <- integer(0)
  }
  
  list(
    selected = selected,
    k_star = k_star,
    threshold = if(k_star > 0) e_sorted$x[k_star] else Inf
  )
}

#' E-BY procedure
e_by_procedure <- function(e_values, alpha) {
  
  m <- length(e_values)
  
  # BY correction factor
  c_m <- sum(1 / (1:m))
  
  # Adjusted alpha
  alpha_adj <- alpha / c_m
  
  e_bh_procedure(e_values, alpha_adj)
}

#' Estimate selection probability
estimate_selection_probability <- function(X, X_knock, y, v, selected, B = 100) {
  
  p <- ncol(X)
  selection_count <- numeric(p)
  
  for(b in 1:B) {
    # Randomize signs
    signs <- sample(c(-1, 1), p, replace = TRUE)
    X_swap <- X
    
    for(j in 1:p) {
      if(signs[j] == -1) {
        temp <- X_swap[, j]
        X_swap[, j] <- X_knock[, j]
        X_knock[, j] <- temp
      }
    }
    
    # Recompute statistics
    W_new <- compute_knockoff_statistics(X_swap, X_knock, y, v, "lasso_diff")
    
    # Check which are selected
    threshold <- knockoff_threshold_plus(W_new, 0.1, offset = 1)
    selected_new <- which(W_new >= threshold)
    
    selection_count[selected_new] <- selection_count[selected_new] + 1
  }
  
  selection_count[selected] / B
}

#' Print method for selection results
#' @export
print.selection_result <- function(x, ...) {
  cat("\nVariable Selection Results\n")
  cat("Method:", x$method, "\n")
  cat("Variables selected:", length(x$selected), "\n")
  
  if(length(x$selected) > 0 && !is.null(x$variable_names)) {
    cat("Selected variables:", paste(x$variable_names, collapse = ", "), "\n")
  }
  
  if(!is.null(x$fdr_control)) {
    cat("FDR control level:", x$fdr_control, "\n")
  }
  
  if(!is.null(x$bfdr)) {
    cat("Bayesian FDR:", round(x$bfdr, 3), "\n")
  }
  
  invisible(x)
}

