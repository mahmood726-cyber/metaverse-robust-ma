
#' @title S4 Class Definitions for metaverse
#' @description Core class structures with validation
#' @importFrom methods setClass setGeneric setMethod setValidity new validObject
#' @import Matrix

# Main metaverse class with validation
setClassUnion("callOrNULL", c("call", "NULL"))

setClass("metaverse",
  slots = list(
    data = "data.frame",
    formula = "formula",
    model = "list",
    robust = "list",
    classical = "list",
    sensitivity = "list",
    inference = "list",
    diagnostics = "list",
    selection = "list",
    plots = "list",
    options = "list",
    call = "callOrNULL",
    version = "character"
  ),
  prototype = list(
    data = data.frame(),
    formula = ~ 1,
    model = list(),
    robust = list(),
    classical = list(),
    sensitivity = list(),
    inference = list(),
    diagnostics = list(),
    selection = list(),
    plots = list(),
    options = list(),
    call = NULL,
    version = "0.1.0"
  )
)

# Validity check for metaverse objects
setValidity("metaverse", function(object) {
  errors <- character()
  
  # Check data requirements
  if(nrow(object@data) > 0) {
    required_cols <- c("yi", "vi")
    if(!all(required_cols %in% names(object@data))) {
      errors <- c(errors, "Data must contain yi and vi columns")
    }
    
    if("vi" %in% names(object@data) && any(object@data$vi <= 0, na.rm = TRUE)) {
      errors <- c(errors, "All variances (vi) must be positive")
    }
  }
  
  if(length(errors) == 0) TRUE else errors
})

# Result class for robust estimation
setClass("robust_result",
  slots = list(
    estimate = "numeric",
    se = "numeric",
    ci = "matrix",
    p_value = "numeric",
    weights = "numeric",
    tau2 = "numeric",
    I2 = "numeric",
    H2 = "numeric",
    Q = "numeric",
    convergence = "list",
    method = "character",
    contamination = "character"
  )
)

# Selection result class
setClass("selection_result",
  slots = list(
    selected = "integer",
    importance = "numeric",
    threshold = "numeric",
    method = "character",
    fdr = "numeric",
    pip = "numeric",
    coefficients = "matrix"
  )
)

# Sensitivity result class
setClass("sensitivity_result",
  slots = list(
    publication_bias = "list",
    influence = "list",
    fragility = "list",
    evalues = "list",
    leave_k_out = "list"
  )
)

# Generic methods
setGeneric("fit_robust", function(object, ...) standardGeneric("fit_robust"))
setGeneric("select_moderators", function(object, ...) standardGeneric("select_moderators"))
setGeneric("quantify_uncertainty", function(object, ...) standardGeneric("quantify_uncertainty"))
setGeneric("assess_sensitivity", function(object, ...) standardGeneric("assess_sensitivity"))
# Note: diagnose(), visualize(), and transport() are regular functions
# that handle dispatch internally via methods::is() checks.
# report() is reserved for future implementation.
setGeneric("report", function(object, ...) standardGeneric("report"))

# Enhanced show method (ASCII-safe for Windows cp1252)
setMethod("show", "metaverse",
  function(object) {
    cat("\n")
    cat("===========================================================\n")
    cat("  METAVERSE Meta-Analysis Results (v", as.character(object@version), ")\n", sep = "")
    cat("===========================================================\n\n")

    # Data summary
    cat("Data Summary:\n")
    cat("   - Studies: ", nrow(object@data), "\n")
    if(nrow(object@data) > 0) {
      cat("   - Total N: ", sum(object@data$ni, na.rm = TRUE), "\n")
      cat("   - Effect size range: [",
          round(min(object@data$yi, na.rm = TRUE), 3), ", ",
          round(max(object@data$yi, na.rm = TRUE), 3), "]\n")
    }

    # Model info
    if(length(object@model) > 0) {
      cat("\nModel Specification:\n")
      cat("   - Effect type: ", object@model$effect_type, "\n")
      cat("   - Method: ", object@model$method, "\n")
      if(!is.null(object@model$moderators)) {
        cat("   - Moderators: ", paste(object@model$moderators, collapse = ", "), "\n")
      }
    }

    # Robust results
    if(length(object@robust) > 0) {
      cat("\nRobust Estimates:\n")
      cat("   - Effect: ", sprintf("%.3f", object@robust$estimate),
          " (SE = ", sprintf("%.3f", object@robust$se), ")\n", sep = "")
      cat("   - 95%% CI: [", sprintf("%.3f", object@robust$ci.lower),
          ", ", sprintf("%.3f", object@robust$ci.upper), "]\n", sep = "")
      cat("   - p-value: ", format.pval(object@robust$p_value, digits = 3), "\n")

      if(!is.null(object@robust$tau2)) {
        cat("\nHeterogeneity:\n")
        cat("   - tau^2 = ", sprintf("%.4f", object@robust$tau2), "\n", sep = "")
        cat("   - I^2 = ", sprintf("%.1f%%", object@robust$I2), "\n", sep = "")
        cat("   - H^2 = ", sprintf("%.2f", object@robust$H2), "\n", sep = "")
      }
    }

    # Selection results
    if(length(object@selection) > 0) {
      cat("\nVariable Selection:\n")
      cat("   - Method: ", object@selection$method, "\n")
      cat("   - Selected: ", length(object@selection$selected), " variables\n")
      if(length(object@selection$selected) > 0) {
        cat("   - Variables: ", paste(object@selection$selected, collapse = ", "), "\n")
      }
    }

    # Sensitivity summary
    if(length(object@sensitivity) > 0) {
      cat("\nSensitivity Analyses Performed:\n")
      analyses <- names(object@sensitivity)
      for(analysis in analyses) {
        cat("   [x] ", analysis, "\n")
      }
    }

    cat("\n-----------------------------------------------------------\n")
    cat("Use summary() for detailed results, plot() for visualizations\n")
    cat("===========================================================\n")
  }
)

# Summary method (use existing generic from methods package)
setMethod("summary", "metaverse",
  function(object, digits = 3) {
    cat("\nDetailed Meta-Analysis Summary\n")
    cat("=====================================\n\n")
    
    # Print all components with formatting
    if(length(object@robust) > 0) {
      print_robust_summary(object@robust, digits)
    }
    
    if(length(object@sensitivity) > 0) {
      print_sensitivity_summary(object@sensitivity, digits)
    }
    
    invisible(object)
  }
)

# Helper functions for pretty printing
print_robust_summary <- function(robust, digits = 3) {
  cat("Robust Meta-Analysis Results\n")
  cat("----------------------------\n")
  
  res_df <- data.frame(
    Estimate = robust$estimate,
    SE = robust$se,
    `CI.lower` = robust$ci.lower,
    `CI.upper` = robust$ci.upper,
    `p-value` = robust$p_value,
    check.names = FALSE
  )
  
  print(round(res_df, digits), row.names = FALSE)
  cat("\n")
}

print_sensitivity_summary <- function(sensitivity, digits = 3) {
  cat("\nSensitivity Analysis Results\n")
  cat("----------------------------\n")
  
  for(name in names(sensitivity)) {
    cat("\n", toupper(name), ":\n", sep = "")
    print(sensitivity[[name]])
  }
}

