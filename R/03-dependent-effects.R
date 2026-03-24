
#' @title Handling Dependent Effect Sizes
#' @description Functions for multilevel and multivariate meta-analysis
#' @import lme4
#' @import Matrix
#' @importFrom stats vcov model.matrix residuals mad cor

#' @export
meta_multilevel <- function(formula, data, 
                           correlation = c("estimate", "CS", "AR1", "unstructured"),
                           robust = TRUE,
                           cluster = NULL) {
  
  correlation <- match.arg(correlation)
  
  # Parse formula to extract random effects
  # This is a simplified version
  if(robust) {
    model <- fit_robust_multilevel(formula, data, correlation)
  } else {
    model <- lme4::lmer(formula, data = data, REML = TRUE)
  }
  
  # Add cluster-robust SEs if requested
  if(!is.null(cluster)) {
    model$vcov_robust <- cluster_robust_vcov(model, cluster)
  }
  
  model
}

#' @export
fit_robust_multilevel <- function(formula, data, correlation) {
  # Simplified robust multilevel fitting
  # In practice, this would implement robust mixed models
  
  # For now, use lmer with robust weights
  initial_fit <- lme4::lmer(formula, data = data, REML = TRUE)
  
  # Calculate robust weights based on residuals
  resid <- residuals(initial_fit)
  weights <- 1 / (1 + abs(resid/mad(resid)))
  
  # Refit with weights
  robust_fit <- lme4::lmer(formula, data = data, weights = weights, REML = TRUE)
  
  robust_fit
}

#' @export
cluster_robust_vcov <- function(model, cluster, type = "CR2") {
  # Cluster-robust variance-covariance matrix
  # CR2 is recommended for small samples
  
  X <- model.matrix(model)
  n <- nrow(X)
  k <- ncol(X)
  
  clusters <- as.factor(cluster)
  n_clusters <- length(unique(clusters))
  
  # Get residuals
  e <- residuals(model)
  
  # Calculate meat of sandwich
  meat <- matrix(0, k, k)
  
  for(j in unique(clusters)) {
    idx <- which(clusters == j)
    Xj <- X[idx, , drop = FALSE]
    ej <- e[idx]
    
    if(type == "CR2") {
      # CR2 adjustment
      Hj <- Xj %*% solve(t(X) %*% X) %*% t(Xj)
      Aj <- diag(length(idx)) - Hj
      Aj_inv <- solve(Aj)
      meat <- meat + t(Xj) %*% Aj_inv %*% (ej %*% t(ej)) %*% t(Aj_inv) %*% Xj
    } else {
      # CR1
      meat <- meat + t(Xj) %*% (ej %*% t(ej)) %*% Xj
    }
  }
  
  # Bread
  bread <- solve(t(X) %*% X)
  
  # Sandwich
  vcov_robust <- bread %*% meat %*% bread * n_clusters / (n_clusters - 1)
  
  vcov_robust
}

#' @export
working_correlation_matrix <- function(data, type = "CS", rho = NULL) {
  # Generate working correlation matrix
  
  n <- nrow(data)
  
  R <- switch(type,
    "CS" = {
      # Compound symmetry
      if(is.null(rho)) rho <- 0.5
      matrix(rho, n, n) + diag(1 - rho, n)
    },
    "AR1" = {
      # Autoregressive order 1
      if(is.null(rho)) rho <- 0.5
      R <- matrix(0, n, n)
      for(i in 1:n) {
        for(j in 1:n) {
          R[i,j] <- rho^abs(i-j)
        }
      }
      R
    },
    "unstructured" = {
      # Estimate from data
      cor(data)
    }
  )
  
  R
}

