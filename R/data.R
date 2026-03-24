#' Metaverse Evaluation Dataset: Contaminated Meta-Analysis
#'
#' A simulated dataset containing effect sizes from 60 studies, explicitly designed
#' to evaluate robust meta-analysis methods. It includes severe data contamination
#' (10% extreme outliers) and 10 potential covariates to test variable selection.
#'
#' @format A data frame with 60 rows and 15 variables:
#' \describe{
#'   \item{study_id}{Character: Unique study identifier}
#'   \item{yi}{Numeric: Effect size estimate (contaminated with outliers)}
#'   \item{vi}{Numeric: Variance of the effect size}
#'   \item{se}{Numeric: Standard error}
#'   \item{X1}{Numeric: True continuous moderator (adds 0.2 to effect)}
#'   \item{X2}{Numeric: True binary moderator (subtracts 0.15 from effect)}
#'   \item{X3-X10}{Numeric: Noise covariates with no true effect}
#'   \item{is_outlier}{Logical: Flag indicating if the study is a true outlier}
#' }
#'
#' @source Simulated evaluation data for the metaverse package demonstration.
#'
#' @examples
#' data(metaverse_eval_data)
#' # Identify outliers
#' sum(metaverse_eval_data$is_outlier)
"metaverse_eval_data"
