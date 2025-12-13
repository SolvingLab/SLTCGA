# ==============================================================================
# Utilities - Validation Layer
# ==============================================================================
# Single validation function for cancer ID and modal type availability
# NOT exported - internal use only
# ==============================================================================


#' Validate Cancer and Modal Availability
#'
#' @description
#' Internal function to validate cancer IDs and check if specified modal types
#' are available for those cancers. Provides helpful messages if data is missing.
#'
#' @param cancer_ids Character vector of cancer IDs (e.g., c("BRCA", "LUAD"))
#' @param modal_type Character. Modal type: "RNAseq", "Protein", "Phospho",
#'   "Methylation", "logCNA", "Mutation", "Clinical", "Survival"
#' @param surv_type Character or NULL. If modal_type="Survival", specify "OS" or "PFS"
#' @param stop_on_error Logical. If TRUE, stops execution. If FALSE, returns warning.
#'
#' @return Invisible TRUE if validation passes
#'
#' @keywords internal
.validate_cancer_modal <- function(cancer_ids,
                                   modal_type,
                                   surv_type = NULL,
                                   stop_on_error = TRUE) {
  # ============================================================================
  # 1. Validate cancer IDs
  # ============================================================================

  valid_cancers <- c(
    "BRCA", "CCRCC", "COAD", "GBM", "HNSCC",
    "LUAD", "LUSC", "OV", "PDAC", "UCEC"
  )

  invalid_cancers <- setdiff(cancer_ids, valid_cancers)

  if (length(invalid_cancers) > 0) {
    msg <- sprintf(
      "Invalid cancer ID(s): %s\nValid options: %s",
      paste(invalid_cancers, collapse = ", "),
      paste(valid_cancers, collapse = ", ")
    )

    if (stop_on_error) {
      stop(msg, call. = FALSE)
    } else {
      warning(msg, call. = FALSE)
      return(invisible(FALSE))
    }
  }

  # ============================================================================
  # 2. Validate modal type availability
  # ============================================================================

  # Modal availability matrix (based on data exploration)
  modal_availability <- list(
    RNAseq      = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LUAD", "LUSC", "OV", "PDAC", "UCEC"),
    Protein     = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LUAD", "LUSC", "OV", "PDAC", "UCEC"),
    Phospho     = c("BRCA", "CCRCC", "GBM", "HNSCC", "LUAD", "LUSC", "PDAC", "UCEC"),
    Methylation = c("CCRCC", "GBM", "HNSCC", "LUAD", "LUSC", "PDAC", "UCEC"),
    logCNA      = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LUAD", "LUSC", "OV", "PDAC", "UCEC"),
    Mutation    = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LUAD", "LUSC", "OV", "PDAC", "UCEC"),
    Clinical    = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LUAD", "LUSC", "OV", "PDAC", "UCEC"),
    Survival    = c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LUAD", "LUSC", "OV", "PDAC", "UCEC")
  )

  # Check if modal type exists (handle scalar only)
  if (length(modal_type) != 1) {
    stop("modal_type must be a single value, not a vector", call. = FALSE)
  }

  if (!(modal_type %in% names(modal_availability))) {
    msg <- sprintf(
      "Invalid modal type: %s\nValid options: %s",
      modal_type,
      paste(names(modal_availability), collapse = ", ")
    )

    if (stop_on_error) {
      stop(msg, call. = FALSE)
    } else {
      warning(msg, call. = FALSE)
      return(invisible(FALSE))
    }
  }

  # Check availability for each cancer
  available_for_modal <- modal_availability[[modal_type]]
  unavailable_cancers <- setdiff(cancer_ids, available_for_modal)

  if (length(unavailable_cancers) > 0) {
    # Provide helpful message
    available_str <- paste(intersect(cancer_ids, available_for_modal), collapse = ", ")
    unavailable_str <- paste(unavailable_cancers, collapse = ", ")

    msg <- sprintf(
      "%s data is NOT available for: %s\n%s data IS available for: %s\n\nSuggestion: Remove unavailable cancer types or choose a different modal type.",
      modal_type,
      unavailable_str,
      modal_type,
      ifelse(available_str == "", "None of your selected cancers", available_str)
    )

    if (stop_on_error) {
      message("\n", msg)
      stop("Modal availability check failed", call. = FALSE)
    } else {
      warning(msg, call. = FALSE)
      return(invisible(FALSE))
    }
  }

  # ============================================================================
  # 3. Special validation for Survival
  # ============================================================================

  if (modal_type == "Survival") {
    if (is.null(surv_type)) {
      msg <- "When modal_type='Survival', you must specify surv_type='OS' or 'PFS'"

      if (stop_on_error) {
        stop(msg, call. = FALSE)
      } else {
        warning(msg, call. = FALSE)
        return(invisible(FALSE))
      }
    }

    if (!surv_type %in% c("OS", "PFS")) {
      msg <- sprintf(
        "Invalid surv_type: %s\nValid options: OS, PFS",
        surv_type
      )

      if (stop_on_error) {
        stop(msg, call. = FALSE)
      } else {
        warning(msg, call. = FALSE)
        return(invisible(FALSE))
      }
    }

    # All cancers have both OS and PFS
    message(sprintf(
      "✓ %s survival data is available for all %d cancer type(s)",
      surv_type,
      length(cancer_ids)
    ))
  } else {
    message(sprintf(
      "✓ %s data validated for %d cancer type(s): %s",
      modal_type,
      length(cancer_ids),
      paste(cancer_ids, collapse = ", ")
    ))
  }

  return(invisible(TRUE))
}
