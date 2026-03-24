# metaverse Evaluation Dataset: Contaminated Meta-Analysis with Moderators
# Simulates a meta-analysis of 60 studies with severe outliers and 10 covariates
set.seed(123)

k <- 60
# True covariates
X1 <- rnorm(k, 0, 1)
X2 <- rbinom(k, 1, 0.5)
# Noise covariates
X_noise <- matrix(rnorm(k * 8), nrow = k, ncol = 8)

# True effect sizes (Baseline effect 0.3, X1 adds 0.2, X2 subtracts 0.15)
true_theta <- 0.3 + 0.2 * X1 - 0.15 * X2
tau2_true <- 0.05
vi <- runif(k, 0.01, 0.08)

# Add random effects and sampling error
yi <- rnorm(k, mean = true_theta, sd = sqrt(vi + tau2_true))

# Contaminate 10% of the studies (Outliers)
outlier_idx <- sample(1:k, size = 6)
yi[outlier_idx] <- yi[outlier_idx] + runif(6, 1.5, 3.0) * sample(c(-1, 1), 6, replace=TRUE)

# Combine into a dataframe
metaverse_eval_data <- data.frame(
  study_id = paste0("Study_", 1:k),
  yi = yi,
  vi = vi,
  se = sqrt(vi),
  X1 = X1,
  X2 = as.numeric(X2)
)

# Add the noise covariates
colnames(X_noise) <- paste0("X", 3:10)
metaverse_eval_data <- cbind(metaverse_eval_data, X_noise)
metaverse_eval_data$is_outlier <- FALSE
metaverse_eval_data$is_outlier[outlier_idx] <- TRUE

# Save to data directory
save(metaverse_eval_data, file = "data/metaverse_eval_data.rda")
cat("✓ metaverse_eval_data.rda created successfully.
")
