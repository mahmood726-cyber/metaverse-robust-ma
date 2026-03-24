
#' Simulate Meta-Analysis Data
#'
#' Generates synthetic meta-analysis data with configurable true effect,
#' heterogeneity, contamination, study size distribution, and moderators.
#'
#' @param k Number of studies (default 50).
#' @param theta True overall effect size (default 0.5).
#' @param tau2 True between-study variance (default 0.1).
#' @param contamination Proportion of outlier studies (default 0, range 0-1).
#' @param study_size_dist Distribution of study sizes: "uniform", "exponential", or "bimodal".
#' @param effect_size_type Type of effect size: "SMD", "OR", "RR", "MD".
#' @param moderators Either a positive integer (number of moderators to generate)
#'   or a matrix of moderator values with k rows.
#' @param seed Optional random seed for reproducibility.
#' @return A data frame with columns study, yi, vi, ni, sei, and any moderator columns.
#' @examples
#' data <- simulate_meta_data(k = 30, theta = 0.5, tau2 = 0.1, seed = 42)
#' head(data)
#' @export
simulate_meta_data <- function(k = 50,
                              theta = 0.5, 
                              tau2 = 0.1,
                              contamination = 0,
                              study_size_dist = "uniform",
                              effect_size_type = "SMD",
                              moderators = NULL,
                              seed = NULL) {
  
  if(!is.null(seed)) set.seed(seed)
  
  # Generate study sizes
  n <- switch(study_size_dist,
    "uniform" = round(runif(k, 20, 200)),
    "exponential" = round(rexp(k, 1/100) + 20),
    "bimodal" = c(round(runif(k/2, 20, 50)), round(runif(k/2, 150, 200))),
    round(runif(k, 20, 200))
  )
  
  # Generate moderators if specified
  if(!is.null(moderators)) {
    if(is.numeric(moderators)) {
      # Number of moderators to generate
      X <- matrix(rnorm(k * moderators), k, moderators)
      colnames(X) <- paste0("mod", 1:moderators)
    } else {
      X <- moderators
    }
    
    # Effect depends on moderators
    beta <- rnorm(ncol(X), 0, 0.2)
    theta_i <- rnorm(k, theta + X %*% beta, sqrt(tau2))
  } else {
    # Random effects
    theta_i <- rnorm(k, theta, sqrt(tau2))
    X <- NULL
  }
  
  # Add contamination
  if(contamination > 0) {
    n_contam <- round(k * contamination)
    contam_idx <- sample(k, n_contam)
    theta_i[contam_idx] <- theta_i[contam_idx] + 
                           rt(n_contam, df = 2) * sqrt(tau2) * 3
  }
  
  # Generate observed effects based on type
  vi <- switch(effect_size_type,
    "SMD" = 1/n + theta_i^2/(2*n),
    "OR" = 1/n,
    "RR" = 1/n,
    "MD" = runif(k, 0.5, 2) / n,
    1/n
  )
  
  yi <- rnorm(k, theta_i, sqrt(vi))
  
  # Create data frame
  data <- data.frame(
    study = paste0("Study_", 1:k),
    yi = yi,
    vi = vi,
    ni = n,
    sei = sqrt(vi)
  )
  
  if(!is.null(X)) {
    data <- cbind(data, X)
  }
  
  # Add attributes
  attr(data, "true_theta") <- theta
  attr(data, "true_tau2") <- tau2
  attr(data, "contamination") <- contamination
  attr(data, "effect_size_type") <- effect_size_type
  
  data
}
