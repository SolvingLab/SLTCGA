# ==============================================================================
# Data Loading Layer
# ==============================================================================
# Unified data loading function that handles all modal types
# Returns wide dataframe with standardized column naming: CancerType_Feature_Modal
# Feature labels in format: THBS2 (Protein, COAD)
# ==============================================================================


#' Load Multi-Omics Data for CPTAC Analysis
#'
#' @description
#' Core data loading function that retrieves and merges multi-omics data from CPTAC database.
#' Supports all omics layers and automatically handles phosphorylation sites, missing data,
#' and multi-cancer merging. Returns standardized data structure ready for downstream analysis.
#'
#' @param var1 Character vector. Gene names or clinical variable names for variable 1.
#'   \itemize{
#'     \item For genes: Use standard gene symbols (e.g., "TP53", "EGFR", "AKT1")
#'     \item For clinical: Use variable names (e.g., "Age", "Tumor_Stage", "Gender")
#'     \item Can be single gene or multiple genes: c("TP53", "EGFR", "KRAS")
#'     \item Phospho sites are automatically detected (no need to specify sites)
#'   }
#'
#' @param var1_modal Character. Omics layer for var1. Must be one of:
#'   \itemize{
#'     \item "RNAseq" - Gene expression (mRNA level)
#'     \item "Protein" - Protein abundance (most stable for enrichment)
#'     \item "Phospho" - Phosphorylation sites (auto-detects all sites for gene)
#'     \item "Mutation" - Somatic mutations (binary: WildType/Mutation)
#'     \item "Clinical" - Clinical variables (Age, Stage, Gender, BMI, etc.)
#'     \item "logCNA" - Copy number alterations (log-transformed)
#'     \item "Methylation" - DNA methylation levels
#'     \item "Survival" - Survival data (requires surv_type parameter)
#'   }
#'
#' @param var1_cancers Character vector. Cancer type(s) for var1. Options:
#'   "BRCA" (breast), "LUAD" (lung adeno), "COAD" (colon), "CCRCC" (kidney),
#'   "GBM" (glioblastoma), "HNSCC" (head-neck), "LUSC" (lung squamous),
#'   "OV" (ovarian), "PDAC" (pancreatic), "UCEC" (endometrial).
#'   Can be single or multiple: c("BRCA", "LUAD", "COAD")
#'
#' @param var2 Character vector or NULL. Gene names or clinical variables for variable 2.
#'   Same format as var1. Set to NULL for single-variable analysis (survival/enrichment).
#'
#' @param var2_modal Character or NULL. Omics layer for var2. Same options as var1_modal.
#'   Set to NULL for single-variable analysis.
#'
#' @param var2_cancers Character vector or NULL. Cancer types for var2.
#'   Can be same or different from var1_cancers. Set to NULL for single-variable analysis.
#'
#' @param surv_type Character or NULL. Survival type when loading survival data:
#'   \itemize{
#'     \item "OS" - Overall Survival
#'     \item "PFS" - Progression-Free Survival
#'     \item NULL - Not loading survival data (default)
#'   }
#'
#' @return A list with 5 components:
#'   \describe{
#'     \item{data}{Wide data.frame with standardized column names (CancerType_Gene_Modal).
#'                 Rows are samples, columns are features + cancer_type column.
#'                 Example: BRCA_TP53_RNAseq, LUAD_AKT1_Protein, BRCA_S124_AKT1_Phospho}
#'     \item{var1_features}{Character vector of feature labels for var1.
#'                          Format: "GENE (Modal, Cancer)" or "SITE_GENE (Modal, Cancer)".
#'                          Example: "TP53 (RNAseq, BRCA)", "S124_AKT1 (Phospho, LUAD)"}
#'     \item{var2_features}{Character vector of feature labels for var2. Empty if var2 is NULL}
#'     \item{var1_types}{Named character vector of variable types for var1.
#'                       Values: "continuous" (numeric), "categorical" (factor), "survival".
#'                       Names: feature labels matching var1_features}
#'     \item{var2_types}{Named character vector of variable types for var2}
#'   }
#'
#' @details
#' **Key Features**:
#' \itemize{
#'   \item Automatic phosphorylation site detection: Input "AKT1" loads all phospho sites
#'   \item Multi-cancer merging: Samples from different cancers are merged with cancer_type column
#'   \item Missing data handling: Features with >90\% missing values are automatically filtered
#'   \item Clinical variable mapping: "Age" maps to "age" column in clinical data
#'   \item Survival data merging: Time and event columns are named as CancerType_SurvType_time/event
#' }
#'
#' **Data Availability**:
#' \itemize{
#'   \item RNAseq, Protein, Mutation, Clinical: All 10 cancer types
#'   \item Phospho: BRCA, CCRCC, GBM, HNSCC, LUAD, LUSC, PDAC, UCEC (8 types)
#'   \item Methylation: CCRCC, GBM, HNSCC, LUAD, LUSC, PDAC, UCEC (7 types)
#'   \item Survival (OS/PFS): All 10 cancer types
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Load single gene, single cancer
#' data <- cptac_load_modality(
#'   var1 = "TP53",
#'   var1_modal = "RNAseq",
#'   var1_cancers = "BRCA"
#' )
#' # Returns: 1 feature, ~120 samples
#'
#' # Example 2: Load phosphorylation data (auto-detects sites)
#' data <- cptac_load_modality(
#'   var1 = "AKT1",
#'   var1_modal = "Phospho",
#'   var1_cancers = "BRCA"
#' )
#' # Returns: ~9 phospho sites for AKT1
#'
#' # Example 3: Load multiple genes, multiple cancers
#' data <- cptac_load_modality(
#'   var1 = c("TP53", "EGFR", "KRAS"),
#'   var1_modal = "RNAseq",
#'   var1_cancers = c("BRCA", "LUAD")
#' )
#' # Returns: 6 features (3 genes × 2 cancers)
#'
#' # Example 4: Load two variables for correlation
#' data <- cptac_load_modality(
#'   var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
#'   var2 = "TP53", var2_modal = "Protein", var2_cancers = "BRCA"
#' )
#' # Returns: 2 features, ready for correlation analysis
#'
#' # Example 5: Load mutation data
#' data <- cptac_load_modality(
#'   var1 = c("KRAS", "EGFR", "TP53"),
#'   var1_modal = "Mutation",
#'   var1_cancers = "LUAD"
#' )
#' # Returns: 3 mutation features (WildType/Mutation)
#'
#' # Example 6: Load clinical variables
#' data <- cptac_load_modality(
#'   var1 = c("Age", "Tumor_Stage", "Gender"),
#'   var1_modal = "Clinical",
#'   var1_cancers = "BRCA"
#' )
#' # Returns: Clinical variables (categorical or continuous)
#' }
#'
#' @export
cptac_load_modality <- function(var1,
                                var1_modal,
                                var1_cancers,
                                var2 = NULL,
                                var2_modal = NULL,
                                var2_cancers = NULL,
                                surv_type = NULL) {
  message("\n========================================")
  message("CPTAC Data Loading")
  message("========================================")

  # Validate inputs (handle vectors of modal types)
  for (modal in unique(var1_modal)) {
    modal_cancers <- unique(var1_cancers[var1_modal == modal])
    .validate_cancer_modal(modal_cancers, modal, surv_type)
  }

  if (!is.null(var2)) {
    for (modal in unique(var2_modal)) {
      modal_cancers <- unique(var2_cancers[var2_modal == modal])
      .validate_cancer_modal(modal_cancers, modal, surv_type)
    }
  }

  # Load var1 (handle multiple vars with different modals)
  message(
    "\n[Var1] Loading ", paste(unique(var1_modal), collapse = "/"), " data for ",
    paste(var1, collapse = ", "), "..."
  )

  # If all modals are the same, load together; otherwise load separately
  if (length(unique(var1_modal)) == 1) {
    var1_result <- .load_modal_data(
      vars = var1,
      modal = var1_modal[1],
      cancers = var1_cancers,
      surv_type = surv_type,
      var_name = "var1"
    )
  } else {
    # Load each var separately and merge
    all_data <- NULL
    all_features <- character(0)
    all_types <- character(0)

    for (i in seq_along(var1)) {
      result_i <- .load_modal_data(
        vars = var1[i],
        modal = var1_modal[i],
        cancers = var1_cancers[i],
        surv_type = surv_type,
        var_name = "var1"
      )

      if (is.null(all_data)) {
        all_data <- result_i$data
      } else {
        all_data <- merge(all_data, result_i$data, by = "row.names", all = FALSE)
        rownames(all_data) <- all_data$Row.names
        all_data$Row.names <- NULL
      }

      all_features <- c(all_features, result_i$features)
      all_types <- c(all_types, result_i$types)
    }

    var1_result <- list(
      data = all_data,
      features = all_features,
      types = all_types
    )
  }

  # Load var2 if provided
  if (!is.null(var2)) {
    message(
      "\n[Var2] Loading ", paste(unique(var2_modal), collapse = "/"), " data for ",
      paste(var2, collapse = ", "), "..."
    )
    var2_result <- .load_modal_data(
      vars = var2,
      modal = var2_modal,
      cancers = var2_cancers,
      surv_type = surv_type,
      var_name = "var2"
    )
  } else {
    # Placeholder for genome-wide scan or single variable analysis
    var2_result <- list(
      data = data.frame(row.names = rownames(var1_result$data)),
      features = character(0),
      types = character(0)
    )
  }

  # Merge data
  message("\n[Merge] Combining datasets...")
  merged_data <- .merge_modal_data(var1_result$data, var2_result$data)

  # Add cancer_type column (extract from sample IDs)
  merged_data$cancer_type <- sapply(strsplit(rownames(merged_data), "_"), `[`, 1)

  # Also add var1_cancers and var2_cancers for reference
  all_cancers <- unique(c(var1_cancers, var2_cancers))

  message(sprintf(
    "✓ Final dataset: %d samples × %d features",
    nrow(merged_data),
    ncol(merged_data)
  ))

  return(list(
    data = merged_data,
    var1_features = var1_result$features,
    var2_features = var2_result$features,
    var1_types = var1_result$types,
    var2_types = var2_result$types,
    var1_cancers = var1_cancers,
    var2_cancers = var2_cancers
  ))
}


# ==============================================================================
# Internal: Load modal data dispatcher
# ==============================================================================

#' Load data for a single modal type
#' @keywords internal
.load_modal_data <- function(vars, modal, cancers, surv_type, var_name) {
  # Route to appropriate loader
  if (modal == "Clinical") {
    return(.load_clinical_modal(vars, cancers))
  } else if (modal == "Survival") {
    return(.load_survival_modal(cancers, surv_type))
  } else if (modal == "Mutation") {
    return(.load_mutation_modal(vars, cancers))
  } else if (modal == "Phospho") {
    return(.load_phospho_modal(vars, cancers))
  } else {
    # Standard omics: RNAseq, Protein, Methylation, logCNA
    return(.load_standard_modal(vars, modal, cancers))
  }
}


# ==============================================================================
# Internal: Standard omics loader (RNAseq, Protein, Methylation, logCNA)
# ==============================================================================

#' Load standard omics data
#' @keywords internal
.load_standard_modal <- function(vars, modal, cancers) {
  gene_path <- Sys.getenv("SL_BULK_DATA")
  if (gene_path == "") {
    stop("SL_BULK_DATA environment variable not set. Please run:\n",
      "Sys.setenv(SL_BULK_DATA = '/path/to/bulk_data')",
      call. = FALSE
    )
  }

  gene_path <- file.path(gene_path, "CPTAC_Omics_Split")

  all_data <- list()
  all_features <- character(0)

  for (cancer in cancers) {
    for (gene in vars) {
      gene_file <- file.path(gene_path, paste0(gene, "_cptac.qs"))

      if (!file.exists(gene_file)) {
        warning(sprintf("Gene file not found: %s (skipping)", gene))
        next
      }

      gene_data <- qs::qread(gene_file)
      gene_data$cancer_type <- sapply(strsplit(gene_data$sampleid, "_"), `[`, 1)

      # Filter cancer and tumor type
      gene_data <- gene_data[gene_data$cancer_type == cancer &
        gene_data$type == "Tumor", ]

      if (nrow(gene_data) == 0) next

      # Extract feature column
      feature_col <- paste0(gene, "_", modal)

      if (!feature_col %in% colnames(gene_data)) {
        warning(sprintf("%s not found for %s in %s (skipping)", feature_col, gene, cancer))
        next
      }

      # Create standardized column name: CancerType_Gene_Modal
      new_col_name <- paste0(cancer, "_", gene, "_", modal)
      feature_label <- paste0(gene, " (", modal, ", ", cancer, ")")

      # Extract data
      extracted <- data.frame(
        value = gene_data[[feature_col]],
        row.names = gene_data$sampleid
      )
      colnames(extracted) <- new_col_name

      # Convert to numeric (ensure no factors)
      extracted[[new_col_name]] <- as.numeric(as.character(extracted[[new_col_name]]))

      # Check if all values are NA
      if (all(is.na(extracted[[new_col_name]]))) {
        warning(sprintf("%s has all NA values in %s (skipping)", gene, cancer))
        next
      }

      # Check if variance is zero (all same value)
      valid_vals <- extracted[[new_col_name]][!is.na(extracted[[new_col_name]])]
      if (length(valid_vals) > 0 && length(unique(valid_vals)) == 1) {
        warning(sprintf(
          "%s has no variance in %s (all values are %.3f, skipping)",
          gene, cancer, valid_vals[1]
        ))
        next
      }

      all_data[[feature_label]] <- extracted
      all_features <- c(all_features, feature_label)
    }
  }

  if (length(all_data) == 0) {
    stop("No valid data loaded for ", modal, call. = FALSE)
  }

  # Merge all features
  merged <- all_data[[1]]

  # Check if merged has data
  if (nrow(merged) == 0) {
    stop(sprintf(
      "No data available for %s in %s",
      vars[1], cancers[1]
    ), call. = FALSE)
  }

  if (length(all_data) > 1) {
    for (i in 2:length(all_data)) {
      # Check before merging
      if (nrow(all_data[[i]]) == 0) {
        warning(sprintf("Skipping empty dataset: %s", names(all_data)[i]))
        next
      }

      merged <- merge(merged, all_data[[i]], by = "row.names", all = TRUE)
      rownames(merged) <- merged$Row.names
      merged$Row.names <- NULL
    }
  }

  # Final check
  if (nrow(merged) == 0) {
    stop(sprintf(
      "No overlapping samples found for %s",
      paste(vars, collapse = ", ")
    ), call. = FALSE)
  }

  # All standard omics are continuous
  types <- rep("continuous", length(all_features))
  names(types) <- all_features

  message(sprintf(
    "  Loaded %d feature(s) from %d cancer(s)",
    length(all_features), length(cancers)
  ))

  return(list(
    data = merged,
    features = all_features,
    types = types
  ))
}


# ==============================================================================
# Internal: Phospho loader (auto-expands to all sites)
# ==============================================================================

#' Load phosphorylation data
#' @keywords internal
.load_phospho_modal <- function(vars, cancers) {
  merge_path <- Sys.getenv("SL_BULK_DATA")
  if (merge_path == "") {
    stop("SL_BULK_DATA environment variable not set", call. = FALSE)
  }

  phospho_file <- file.path(merge_path, "LinkedOmicsKB_PanCancer_Phospho_Quantification.qs")

  if (!file.exists(phospho_file)) {
    stop("Phospho data file not found: ", phospho_file, call. = FALSE)
  }

  phospho <- qs::qread(phospho_file, nthreads = min(parallel::detectCores(), 6))
  phospho$cancer_type <- sapply(strsplit(rownames(phospho), "_"), `[`, 1)

  # Load clinical to filter tumor samples
  clin_file <- file.path(merge_path, "LinkedOmicsKB_PanCancer_Clin.qs")
  clin <- qs::qread(clin_file)
  clin$type <- ifelse(is.na(clin$type), "Tumor", clin$type)
  clin <- clin[clin$type == "Tumor", c("sampleid", "type")]

  phospho <- phospho[rownames(phospho) %in% clin$sampleid, ]

  # Extract phospho sites for target genes
  all_data <- list()
  all_features <- character(0)

  for (cancer in cancers) {
    cancer_phospho <- phospho[phospho$cancer_type == cancer, ]

    for (gene in vars) {
      # Find all phospho sites for this gene
      gene_pattern <- paste0("_", gene, "$")
      site_cols <- grep(gene_pattern, colnames(cancer_phospho), value = TRUE)

      if (length(site_cols) == 0) {
        warning(sprintf("No phospho sites found for %s in %s (skipping)", gene, cancer))
        next
      }

      # Filter by NA rate (< 90%)
      valid_sites <- character(0)
      for (site in site_cols) {
        na_rate <- sum(is.na(cancer_phospho[[site]])) / nrow(cancer_phospho)
        if (na_rate < 0.9) {
          valid_sites <- c(valid_sites, site)
        }
      }

      if (length(valid_sites) == 0) {
        warning(sprintf("No valid phospho sites for %s in %s (all > 90%% NA, skipping)", gene, cancer))
        next
      }

      # Add each site as separate feature
      for (site in valid_sites) {
        new_col_name <- paste0(cancer, "_", site, "_Phospho")
        feature_label <- paste0(site, " (Phospho, ", cancer, ")")

        extracted <- data.frame(
          value = as.numeric(as.character(cancer_phospho[[site]])),
          row.names = rownames(cancer_phospho)
        )
        colnames(extracted) <- new_col_name

        all_data[[feature_label]] <- extracted
        all_features <- c(all_features, feature_label)
      }
    }
  }

  if (length(all_data) == 0) {
    stop("No valid phospho sites found", call. = FALSE)
  }

  # Merge
  merged <- all_data[[1]]
  if (length(all_data) > 1) {
    for (i in 2:length(all_data)) {
      merged <- merge(merged, all_data[[i]], by = "row.names", all = TRUE)
      rownames(merged) <- merged$Row.names
      merged$Row.names <- NULL
    }
  }

  types <- rep("continuous", length(all_features))
  names(types) <- all_features

  message(sprintf("  Loaded %d phospho site(s)", length(all_features)))

  return(list(
    data = merged,
    features = all_features,
    types = types
  ))
}


# ==============================================================================
# Internal: Mutation loader (binary categorical)
# ==============================================================================

#' Load mutation data
#' @keywords internal
.load_mutation_modal <- function(vars, cancers) {
  gene_path <- Sys.getenv("SL_BULK_DATA")
  if (gene_path == "") {
    stop("SL_BULK_DATA environment variable not set", call. = FALSE)
  }

  gene_path <- file.path(gene_path, "CPTAC_Omics_Split")

  all_data <- list()
  all_features <- character(0)

  for (cancer in cancers) {
    for (gene in vars) {
      gene_file <- file.path(gene_path, paste0(gene, "_cptac.qs"))

      if (!file.exists(gene_file)) {
        warning(sprintf("Gene file not found: %s (skipping)", gene))
        next
      }

      gene_data <- qs::qread(gene_file)
      gene_data$cancer_type <- sapply(strsplit(gene_data$sampleid, "_"), `[`, 1)

      gene_data <- gene_data[gene_data$cancer_type == cancer &
        gene_data$type == "Tumor", ]

      if (nrow(gene_data) == 0) next

      # Extract mutation column
      mut_col <- paste0(gene, "_Mutation")

      if (!mut_col %in% colnames(gene_data)) {
        warning(sprintf("%s not found for %s in %s (skipping)", mut_col, gene, cancer))
        next
      }

      new_col_name <- paste0(cancer, "_", gene, "_Mutation")
      feature_label <- paste0(gene, " (Mutation, ", cancer, ")")

      extracted <- data.frame(
        value = factor(gene_data[[mut_col]], levels = c("WildType", "Mutation")),
        row.names = gene_data$sampleid
      )
      colnames(extracted) <- new_col_name

      all_data[[feature_label]] <- extracted
      all_features <- c(all_features, feature_label)
    }
  }

  if (length(all_data) == 0) {
    stop("No valid mutation data loaded", call. = FALSE)
  }

  # Merge
  merged <- all_data[[1]]

  # Check if first dataset has data
  if (nrow(merged) == 0) {
    stop("No mutation data available for specified genes/cancers", call. = FALSE)
  }

  if (length(all_data) > 1) {
    for (i in 2:length(all_data)) {
      # Skip empty datasets
      if (nrow(all_data[[i]]) == 0) {
        warning(sprintf("Skipping empty mutation data: %s", names(all_data)[i]))
        next
      }

      merged <- merge(merged, all_data[[i]], by = "row.names", all = TRUE)
      rownames(merged) <- merged$Row.names
      merged$Row.names <- NULL
    }
  }

  # Final check
  if (nrow(merged) == 0) {
    stop("No overlapping samples found for mutation data", call. = FALSE)
  }

  types <- rep("categorical", length(all_features))
  names(types) <- all_features

  message(sprintf("  Loaded %d mutation feature(s)", length(all_features)))

  return(list(
    data = merged,
    features = all_features,
    types = types
  ))
}


# ==============================================================================
# Internal: Clinical loader (auto-categorizes Age/BMI)
# ==============================================================================

#' Load clinical data
#' @keywords internal
.load_clinical_modal <- function(vars, cancers) {
  gene_path <- Sys.getenv("SL_BULK_DATA")
  if (gene_path == "") {
    stop("SL_BULK_DATA environment variable not set", call. = FALSE)
  }

  gene_path <- file.path(gene_path, "CPTAC_Omics_Split")

  # Use TP53 as reference (contains all clinical data)
  gene_file <- file.path(gene_path, "TP53_cptac.qs")

  if (!file.exists(gene_file)) {
    stop("Reference gene file not found for clinical data", call. = FALSE)
  }

  gene_data <- qs::qread(gene_file)
  gene_data$cancer_type <- sapply(strsplit(gene_data$sampleid, "_"), `[`, 1)
  gene_data <- gene_data[gene_data$type == "Tumor", ]

  # Get all available clinical columns (exclude mutation columns and omics columns)
  exclude_patterns <- c(
    "_mutation$", # All mutation columns
    "_RNAseq$", "_Protein$", "_Methylation$", "_logCNA$", "_Mutation$", # Omics columns
    "^os_", "^pfs_", # Survival columns
    "^sampleid$", "^type$", "^cancer_type$", # Meta columns
    "^mutation_table$" # Mutation table
  )

  all_cols <- colnames(gene_data)
  clinical_cols <- all_cols

  for (pattern in exclude_patterns) {
    clinical_cols <- grep(pattern, clinical_cols, value = TRUE, invert = TRUE)
  }

  message(sprintf("  Available clinical variables: %s", paste(head(clinical_cols, 10), collapse = ", ")))
  if (length(clinical_cols) > 10) {
    message(sprintf("    ... and %d more", length(clinical_cols) - 10))
  }

  all_data <- list()
  all_features <- character(0)
  all_types <- character(0)

  for (cancer in cancers) {
    cancer_data <- gene_data[gene_data$cancer_type == cancer, ]

    if (nrow(cancer_data) == 0) {
      warning(sprintf("No clinical data for %s (skipping)", cancer))
      next
    }

    for (var in vars) {
      # Try to match clinical column
      matched_col <- .match_clinical_col(var, clinical_cols)

      if (is.null(matched_col)) {
        warning(sprintf("Clinical variable '%s' not found (skipping)", var))
        message(sprintf("  Available columns: %s", paste(head(clinical_cols, 20), collapse = ", ")))
        next
      }

      # Get values
      values <- cancer_data[[matched_col]]

      # Remove NA values
      if (all(is.na(values))) {
        warning(sprintf("All NA values for %s in %s (skipping)", var, cancer))
        next
      }

      # Categorize
      categorized <- .categorize_clinical(values, tolower(matched_col))

      new_col_name <- paste0(cancer, "_", var, "_Clinical")
      feature_label <- paste0(var, " (Clinical, ", cancer, ")")

      extracted <- data.frame(
        value = categorized$values,
        row.names = cancer_data$sampleid
      )
      colnames(extracted) <- new_col_name

      all_data[[feature_label]] <- extracted
      all_features <- c(all_features, feature_label)
      all_types <- c(all_types, categorized$type)
    }
  }

  if (length(all_data) == 0) {
    stop("No valid clinical data loaded", call. = FALSE)
  }

  # Merge
  merged <- all_data[[1]]
  if (length(all_data) > 1) {
    for (i in 2:length(all_data)) {
      merged <- merge(merged, all_data[[i]], by = "row.names", all = TRUE)
      rownames(merged) <- merged$Row.names
      merged$Row.names <- NULL
    }
  }

  names(all_types) <- all_features

  message(sprintf("  Loaded %d clinical feature(s)", length(all_features)))

  return(list(
    data = merged,
    features = all_features,
    types = all_types
  ))
}


#' Match clinical column name
#' @keywords internal
.match_clinical_col <- function(search, available) {
  search_lower <- tolower(search)
  available_lower <- tolower(available)

  # Direct match
  idx <- which(available_lower == search_lower)
  if (length(idx) > 0) {
    return(available[idx[1]])
  }

  # Partial match (for convenience)
  idx <- grep(paste0("^", search_lower, "$"), available_lower)
  if (length(idx) > 0) {
    return(available[idx[1]])
  }

  # Aliases for common variables
  aliases <- list(
    age = c("age"),
    bmi = c("bmi", "body_mass_index"),
    sex = c("sex", "gender"),
    stage = c("stage", "tumor_stage", "ajcc_pathologic_stage"),
    grade = c("grade", "histologic_grade", "histologic_grade"),
    path_stage_pt = c("path_stage_pt", "pt_stage", "t_stage"),
    path_stage_pn = c("path_stage_pn", "pn_stage", "n_stage"),
    tobacco = c("tobacco_smoking_history", "tobacco", "smoking"),
    tumor_necrosis = c("tumor_necrosis", "necrosis"),
    tumor_size = c("tumor_size_cm", "tumor_size", "size")
  )

  for (std in names(aliases)) {
    if (search_lower %in% aliases[[std]]) {
      for (alias in aliases[[std]]) {
        idx <- which(available_lower == alias)
        if (length(idx) > 0) {
          return(available[idx[1]])
        }
      }
    }
  }

  return(NULL)
}


#' Categorize clinical values
#' @keywords internal
.categorize_clinical <- function(values, var_name) {
  # Age: < 60 vs >= 60
  if (var_name %in% c("age")) {
    cutoff <- 60
    categorized <- ifelse(values < cutoff,
      paste0("<", cutoff),
      paste0(">=", cutoff)
    )
    return(list(values = factor(categorized), type = "categorical"))
  }

  # BMI: < 25 vs >= 25
  if (var_name %in% c("bmi", "body_mass_index")) {
    cutoff <- 25
    categorized <- ifelse(values < cutoff,
      paste0("<", cutoff),
      paste0(">=", cutoff)
    )
    return(list(values = factor(categorized), type = "categorical"))
  }

  # All others: already categorical
  return(list(values = factor(values), type = "categorical"))
}


# ==============================================================================
# Internal: Survival loader
# ==============================================================================

#' Load survival data
#' @keywords internal
.load_survival_modal <- function(cancers, surv_type) {
  gene_path <- Sys.getenv("SL_BULK_DATA")
  if (gene_path == "") {
    stop("SL_BULK_DATA environment variable not set", call. = FALSE)
  }

  gene_path <- file.path(gene_path, "CPTAC_Omics_Split")

  # Use TP53 as reference (all genes have same survival data)
  gene_file <- file.path(gene_path, "TP53_cptac.qs")

  if (!file.exists(gene_file)) {
    stop("Reference gene file not found for survival data", call. = FALSE)
  }

  gene_data <- qs::qread(gene_file)
  gene_data$cancer_type <- sapply(strsplit(gene_data$sampleid, "_"), `[`, 1)
  gene_data <- gene_data[gene_data$type == "Tumor", ]

  # Select survival columns based on type
  if (surv_type == "OS") {
    time_col <- "os_time"
    event_col <- "os_event"
  } else {
    time_col <- "pfs_time"
    event_col <- "pfs_event"
  }

  all_data <- list()
  all_features <- character(0)

  for (cancer in cancers) {
    cancer_data <- gene_data[gene_data$cancer_type == cancer, ]

    if (nrow(cancer_data) == 0) next

    # Extract survival columns
    time_col_name <- paste0(cancer, "_", surv_type, "_time")
    event_col_name <- paste0(cancer, "_", surv_type, "_event")

    surv_df <- data.frame(
      time = cancer_data[[time_col]],
      event = cancer_data[[event_col]],
      row.names = cancer_data$sampleid
    )
    colnames(surv_df) <- c(time_col_name, event_col_name)

    all_data[[cancer]] <- surv_df
    all_features <- c(
      all_features,
      paste0(surv_type, " (Survival, ", cancer, ")")
    )
  }

  if (length(all_data) == 0) {
    stop("No survival data loaded", call. = FALSE)
  }

  # Merge
  merged <- all_data[[1]]
  if (length(all_data) > 1) {
    for (i in 2:length(all_data)) {
      merged <- merge(merged, all_data[[i]], by = "row.names", all = TRUE)
      rownames(merged) <- merged$Row.names
      merged$Row.names <- NULL
    }
  }

  # Survival is special type
  types <- rep("survival", length(all_features))
  names(types) <- all_features

  message(sprintf(
    "  Loaded %s survival data for %d cancer(s)",
    surv_type, length(cancers)
  ))

  return(list(
    data = merged,
    features = all_features,
    types = types
  ))
}


# ==============================================================================
# Internal: Merge modal data
# ==============================================================================

#' Merge two modal datasets
#' @keywords internal
.merge_modal_data <- function(data1, data2) {
  # Check for empty datasets
  if (nrow(data1) == 0) {
    stop("No data available for var1", call. = FALSE)
  }

  if (ncol(data2) == 0) {
    return(data1)
  }

  if (nrow(data2) == 0) {
    stop("No data available for var2", call. = FALSE)
  }

  # Check for overlapping column names
  overlap_cols <- intersect(colnames(data1), colnames(data2))

  if (length(overlap_cols) > 0) {
    # Get all samples (union of rownames)
    all_samples <- union(rownames(data1), rownames(data2))

    # Create empty merged data frame
    merged <- data.frame(row.names = all_samples)

    # Add data1 columns
    for (col in colnames(data1)) {
      merged[[col]] <- data1[rownames(merged), col]
    }

    # Add data2 columns (non-overlapping only)
    for (col in setdiff(colnames(data2), overlap_cols)) {
      merged[[col]] <- data2[rownames(merged), col]
    }

    return(merged)
  } else {
    # No overlap, use standard merge
    merged <- merge(data1, data2, by = "row.names", all = TRUE)
    rownames(merged) <- merged$Row.names
    merged$Row.names <- NULL
    return(merged)
  }
}
