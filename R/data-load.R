# ==============================================================================
# Data Loading Layer
# ==============================================================================
# Unified data loading function for 8 TCGA modalities
# Returns standardized wide dataframe ready for analysis
# ==============================================================================

#' Load Multi-Omics Data for TCGA Analysis
#'
#' @description
#' Core data loading function that retrieves and merges multi-omics data from TCGA database.
#' Supports 8 omics layers, 33 main cancer types plus 32 molecular subtypes (65 total).
#' Returns standardized data structure ready for downstream analysis.
#'
#' @param var1 Character vector. Gene names, clinical variables, miRNA names, or immune cells.
#' @param var1_modal Character. Omics layer for var1. Must be one of:
#'   "RNAseq", "Mutation", "CNV", "Methylation", "miRNA",
#'   "Clinical", "Signature", "ImmuneCell", "Survival"
#' @param var1_cancers Character vector. Cancer type(s). Case insensitive.
#'   33 main types + 32 subtypes. Examples: "BRCA", "brca", "BRCA_IDC", "LUAD"
#' @param var2 Character vector or NULL. Second variable (for correlation).
#' @param var2_modal Character or NULL. Modal for var2.
#' @param var2_cancers Character vector or NULL. Cancers for var2.
#' @param surv_type Character or NULL. Survival type: "OS", "DSS", "PFI", "DFI"
#' @param rnaseq_type Character. RNAseq normalization method (default: "log2TPM")
#'   Options: "log2TPM" (recommended), "log2RSEM", "log2FPKM", "log2Counts"
#' @param cnv_type Character. CNV algorithm (default: "SNP6_Array")
#'   Options: "SNP6_Array", "Gistic2_Score", "ABSOLUTE"
#' @param methylation_region Character. Methylation region (default: "Promoter_mean")
#'   Options: "Promoter_mean", "Global_mean", "Enhancer_mean", "TSS200_mean"
#' @param immune_algorithm Character or NULL. Immune deconvolution algorithm.
#'   Options: "cibersort", "xcell", "quantiseq", "mcpcounter", "timer", "epic", "ips", "estimate"
#'   If NULL, auto-matches or loads all algorithms for a cell type.
#'
#' @return A list with 5 components:
#'   \item{data}{Wide data.frame with standardized column names}
#'   \item{var1_features}{Feature labels for var1 (format: "GENE (Modal, Cancer)")}
#'   \item{var2_features}{Feature labels for var2}
#'   \item{var1_types}{Variable types for var1 ("continuous" or "categorical")}
#'   \item{var2_types}{Variable types for var2}
#'
#' @export
tcga_load_modality <- function(var1,
                               var1_modal,
                               var1_cancers,
                               var2 = NULL,
                               var2_modal = NULL,
                               var2_cancers = NULL,
                               surv_type = NULL,
                               rnaseq_type = "log2TPM",
                               cnv_type = "SNP6_Array",
                               methylation_region = "Promoter_mean",
                               immune_algorithm = NULL) {
  message("\n========================================")
  message("TCGA Data Loading")
  message("========================================")

  # Validate inputs
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

  # Load var1
  message(
    "\n[Var1] Loading ", paste(unique(var1_modal), collapse = "/"), " data for ",
    paste(var1, collapse = ", "), "..."
  )

  var1_result <- .load_modal_data_tcga(
    vars = var1,
    modal = var1_modal[1],
    cancers = var1_cancers,
    surv_type = surv_type,
    var_name = "var1",
    rnaseq_type = rnaseq_type,
    cnv_type = cnv_type,
    methylation_region = methylation_region,
    immune_algorithm = immune_algorithm
  )

  # Load var2 if provided
  if (!is.null(var2)) {
    message(
      "\n[Var2] Loading ", paste(unique(var2_modal), collapse = "/"), " data for ",
      paste(var2, collapse = ", "), "..."
    )
    var2_result <- .load_modal_data_tcga(
      vars = var2,
      modal = var2_modal[1],
      cancers = var2_cancers,
      surv_type = surv_type,
      var_name = "var2",
      rnaseq_type = rnaseq_type,
      cnv_type = cnv_type,
      methylation_region = methylation_region,
      immune_algorithm = immune_algorithm
    )
  } else {
    var2_result <- list(
      data = data.frame(row.names = rownames(var1_result$data)),
      features = character(0),
      types = character(0)
    )
  }

  # Merge data
  message("\n[Merge] Combining datasets...")
  merged_data <- .merge_modal_data_tcga(var1_result$data, var2_result$data)

  # Add cancer_type column (with subtype support)
  all_cancers <- unique(c(var1_cancers, var2_cancers))
  merged_data$cancer_type <- .assign_cancer_type(rownames(merged_data), all_cancers)

  message(sprintf(
    "✓ Final dataset: %d samples × %d features",
    nrow(merged_data),
    ncol(merged_data) - 1 # Exclude cancer_type
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
# Modal Data Dispatcher
# ==============================================================================

#' Load data for a single modal type
#' @keywords internal
.load_modal_data_tcga <- function(vars, modal, cancers, surv_type, var_name,
                                  rnaseq_type, cnv_type, methylation_region, immune_algorithm) {
  # Route to appropriate loader
  switch(modal,
    "RNAseq" = .load_rnaseq_modal_tcga(vars, cancers, rnaseq_type),
    "Mutation" = .load_mutation_modal_tcga(vars, cancers),
    "CNV" = .load_cnv_modal_tcga(vars, cancers, cnv_type),
    "Methylation" = .load_methylation_modal_tcga(vars, cancers, methylation_region),
    "miRNA" = .load_mirna_modal_tcga(vars, cancers),
    "Clinical" = .load_clinical_modal_tcga(vars, cancers),
    "Signature" = .load_signature_modal_tcga(vars, cancers),
    "ImmuneCell" = .load_immunecell_modal_tcga(vars, cancers, immune_algorithm),
    "Survival" = .load_survival_modal_tcga(cancers, surv_type),
    stop("Unknown modal type: ", modal, call. = FALSE)
  )
}

# ==============================================================================
# RNAseq Loader (Multiple normalization methods)
# ==============================================================================

#' Load RNAseq data
#' @keywords internal
.load_rnaseq_modal_tcga <- function(vars, cancers, rnaseq_type = "log2TPM") {
  base_path <- Sys.getenv("SL_BULK_DATA")
  if (base_path == "") {
    stop("SL_BULK_DATA environment variable not set", call. = FALSE)
  }

  # Map rnaseq_type to column suffix
  suffix_map <- list(
    "log2TPM" = "RNAseq_TPM_unstranded_log2",
    "log2RSEM" = "RNAseq_RSEM_Batch_Normalized_log2",
    "log2FPKM" = "RNAseq_FPKM_unstranded_log2",
    "log2Counts" = "RNAseq_Counts_unstranded_log2"
  )

  col_suffix <- suffix_map[[rnaseq_type]]
  if (is.null(col_suffix)) {
    stop("Unknown rnaseq_type: ", rnaseq_type, call. = FALSE)
  }

  # Map rnaseq_type to readable label
  label_map <- list(
    "log2TPM" = "RNAseq",
    "log2RSEM" = "RNAseq_RSEM",
    "log2FPKM" = "RNAseq_FPKM",
    "log2Counts" = "RNAseq_Counts"
  )
  modal_label <- label_map[[rnaseq_type]]

  gene_path <- file.path(base_path, "TCGA_Omics_Split")

  all_data <- list()
  all_features <- character(0)

  for (cancer in cancers) {
    cancer_info <- .parse_cancer_type(cancer)

    for (gene in vars) {
      gene_file <- file.path(gene_path, paste0(gene, ".qs"))

      if (!file.exists(gene_file)) {
        warning(sprintf("Gene file not found: %s (skipping)", gene))
        next
      }

      gene_data <- qs::qread(gene_file, nthreads = min(parallel::detectCores(), 4))

      # Filter by cancer type
      parent_mask <- .match_cancer_parent(gene_data$cancertype_tcga_clinical, cancer_info)
      gene_data <- gene_data[parent_mask & gene_data$tissue_type_clinical == "Tumor", ]

      if (nrow(gene_data) == 0) next

      # Subtype filtering
      if (cancer_info$type == "subtype") {
        subtype_samples <- .filter_subtype_samples(cancer_info)
        gene_data <- gene_data[gene_data$patient_barcode_16_clinical %in% subtype_samples, ]
      }

      if (nrow(gene_data) == 0) next

      # Extract feature column
      feature_col <- paste0(gene, "_", col_suffix)

      if (!feature_col %in% colnames(gene_data)) {
        warning(sprintf("%s not found for %s in %s (skipping)", feature_col, gene, cancer))
        next
      }

      # Standardized column name and label
      new_col_name <- paste0(cancer, "_", gene, "_", modal_label)
      feature_label <- paste0(gene, " (", modal_label, ", ", cancer, ")")

      extracted <- data.frame(
        value = as.numeric(gene_data[[feature_col]]),
        row.names = gene_data$patient_barcode_16_clinical
      )
      colnames(extracted) <- new_col_name

      # Filter all NA or zero variance
      if (all(is.na(extracted[[new_col_name]]))) {
        warning(sprintf("%s has all NA values in %s (skipping)", gene, cancer))
        next
      }

      valid_vals <- extracted[[new_col_name]][!is.na(extracted[[new_col_name]])]
      if (length(valid_vals) > 0 && length(unique(valid_vals)) == 1) {
        warning(sprintf("%s has no variance in %s (skipping)", gene, cancer))
        next
      }

      all_data[[feature_label]] <- extracted
      all_features <- c(all_features, feature_label)
    }
  }

  if (length(all_data) == 0) {
    stop("No valid RNAseq data loaded", call. = FALSE)
  }

  merged <- .merge_list_data(all_data)

  types <- rep("continuous", length(all_features))
  names(types) <- all_features

  message(sprintf("  Loaded %d RNAseq feature(s) [%s]", length(all_features), rnaseq_type))

  return(list(data = merged, features = all_features, types = types))
}

# ==============================================================================
# Mutation Loader
# ==============================================================================

#' Load mutation data
#' @keywords internal
.load_mutation_modal_tcga <- function(vars, cancers) {
  base_path <- Sys.getenv("SL_BULK_DATA")
  gene_path <- file.path(base_path, "TCGA_Omics_Split")

  all_data <- list()
  all_features <- character(0)

  for (cancer in cancers) {
    cancer_info <- .parse_cancer_type(cancer)

    for (gene in vars) {
      gene_file <- file.path(gene_path, paste0(gene, ".qs"))

      if (!file.exists(gene_file)) {
        warning(sprintf("Gene file not found: %s (skipping)", gene))
        next
      }

      gene_data <- qs::qread(gene_file, nthreads = 4)

      parent_mask <- .match_cancer_parent(gene_data$cancertype_tcga_clinical, cancer_info)
      gene_data <- gene_data[parent_mask & gene_data$tissue_type_clinical == "Tumor", ]

      if (nrow(gene_data) == 0) next

      # Subtype filtering
      if (cancer_info$type == "subtype") {
        subtype_samples <- .filter_subtype_samples(cancer_info)
        gene_data <- gene_data[gene_data$patient_barcode_16_clinical %in% subtype_samples, ]
      }

      if (nrow(gene_data) == 0) next

      # Extract mutation column (prefer excludeSyn)
      mut_col <- paste0(gene, "_SomaticMutation_WXS_excludeSyn_wm")

      if (!mut_col %in% colnames(gene_data)) {
        warning(sprintf("%s not found for %s in %s (skipping)", mut_col, gene, cancer))
        next
      }

      # Convert to factor
      # TCGA values: "Wild", "Mutant", "SynMutant" (synonymous mutation)
      # Convert to: WildType (Wild + SynMutant) vs Mutation (Mutant)
      mut_values <- gene_data[[mut_col]]
      mut_values <- ifelse(mut_values %in% c("Wild", "SynMutant"), "WildType",
        ifelse(mut_values == "Mutant", "Mutation", NA)
      )
      mut_values <- factor(mut_values, levels = c("WildType", "Mutation"))

      new_col_name <- paste0(cancer, "_", gene, "_Mutation")
      feature_label <- paste0(gene, " (Mutation, ", cancer, ")")

      extracted <- data.frame(
        value = mut_values,
        row.names = gene_data$patient_barcode_16_clinical
      )
      colnames(extracted) <- new_col_name

      all_data[[feature_label]] <- extracted
      all_features <- c(all_features, feature_label)
    }
  }

  if (length(all_data) == 0) {
    stop("No valid Mutation data loaded", call. = FALSE)
  }

  merged <- .merge_list_data(all_data)

  types <- rep("categorical", length(all_features))
  names(types) <- all_features

  message(sprintf("  Loaded %d Mutation feature(s)", length(all_features)))

  return(list(data = merged, features = all_features, types = types))
}

# ==============================================================================
# CNV Loader
# ==============================================================================

#' Load CNV data
#' @keywords internal
.load_cnv_modal_tcga <- function(vars, cancers, cnv_type = "SNP6_Array") {
  base_path <- Sys.getenv("SL_BULK_DATA")
  gene_path <- file.path(base_path, "TCGA_Omics_Split")

  # Map cnv_type to column suffix
  suffix_map <- list(
    "SNP6_Array" = "CNV_SNP6_Array",
    "Gistic2_Score" = "CNV_Gistic2_Score",
    "ABSOLUTE" = "CNV_ABSOLUTE_log2",
    "ASCAT2" = "CNV_ASCAT2_log2"
  )

  col_suffix <- suffix_map[[cnv_type]]
  if (is.null(col_suffix)) {
    stop("Unknown cnv_type: ", cnv_type, call. = FALSE)
  }

  # Modal label for plots
  label_map <- list(
    "SNP6_Array" = "CNV",
    "Gistic2_Score" = "CNV_Gistic2",
    "ABSOLUTE" = "CNV_ABSOLUTE",
    "ASCAT2" = "CNV_ASCAT2"
  )
  modal_label <- label_map[[cnv_type]]

  all_data <- list()
  all_features <- character(0)

  for (cancer in cancers) {
    cancer_info <- .parse_cancer_type(cancer)

    for (gene in vars) {
      gene_file <- file.path(gene_path, paste0(gene, ".qs"))

      if (!file.exists(gene_file)) {
        warning(sprintf("Gene file not found: %s (skipping)", gene))
        next
      }

      gene_data <- qs::qread(gene_file, nthreads = 4)
      parent_mask <- .match_cancer_parent(gene_data$cancertype_tcga_clinical, cancer_info)
      gene_data <- gene_data[parent_mask & gene_data$tissue_type_clinical == "Tumor", ]

      if (nrow(gene_data) == 0) next

      if (cancer_info$type == "subtype") {
        subtype_samples <- .filter_subtype_samples(cancer_info)
        gene_data <- gene_data[gene_data$patient_barcode_16_clinical %in% subtype_samples, ]
      }

      if (nrow(gene_data) == 0) next

      feature_col <- paste0(gene, "_", col_suffix)

      if (!feature_col %in% colnames(gene_data)) {
        warning(sprintf("%s not found for %s (skipping)", feature_col, gene))
        next
      }

      new_col_name <- paste0(cancer, "_", gene, "_", modal_label)
      feature_label <- paste0(gene, " (", modal_label, ", ", cancer, ")")

      extracted <- data.frame(
        value = as.numeric(gene_data[[feature_col]]),
        row.names = gene_data$patient_barcode_16_clinical
      )
      colnames(extracted) <- new_col_name

      all_data[[feature_label]] <- extracted
      all_features <- c(all_features, feature_label)
    }
  }

  if (length(all_data) == 0) {
    stop("No valid CNV data loaded", call. = FALSE)
  }

  merged <- .merge_list_data(all_data)

  # Add cancer type info (with subtype support)
  merged$cancer_type <- .assign_cancer_type(rownames(merged), cancers)

  types <- rep("continuous", length(all_features))
  names(types) <- all_features

  message(sprintf("  Loaded %d CNV feature(s) [%s]", length(all_features), cnv_type))

  return(list(data = merged, features = all_features, types = types))
}

# ==============================================================================
# Methylation Loader
# ==============================================================================

#' Load methylation data
#' @keywords internal
.load_methylation_modal_tcga <- function(vars, cancers, methylation_region = "Promoter_mean") {
  base_path <- Sys.getenv("SL_BULK_DATA")
  gene_path <- file.path(base_path, "TCGA_Omics_Split")

  # Column suffix
  col_suffix <- paste0("Methylation450_", methylation_region)

  # Modal label
  region_label <- gsub("_mean|_median", "", methylation_region)
  modal_label <- paste0("Methylation_", region_label)

  all_data <- list()
  all_features <- character(0)

  for (cancer in cancers) {
    cancer_info <- .parse_cancer_type(cancer)

    for (gene in vars) {
      gene_file <- file.path(gene_path, paste0(gene, ".qs"))

      if (!file.exists(gene_file)) {
        warning(sprintf("Gene file not found: %s (skipping)", gene))
        next
      }

      gene_data <- qs::qread(gene_file, nthreads = 4)
      parent_mask <- .match_cancer_parent(gene_data$cancertype_tcga_clinical, cancer_info)
      gene_data <- gene_data[parent_mask & gene_data$tissue_type_clinical == "Tumor", ]

      if (nrow(gene_data) == 0) next

      if (cancer_info$type == "subtype") {
        subtype_samples <- .filter_subtype_samples(cancer_info)
        gene_data <- gene_data[gene_data$patient_barcode_16_clinical %in% subtype_samples, ]
      }

      if (nrow(gene_data) == 0) next

      feature_col <- paste0(gene, "_", col_suffix)

      if (!feature_col %in% colnames(gene_data)) {
        warning(sprintf("%s not found for %s (skipping)", feature_col, gene))
        next
      }

      new_col_name <- paste0(cancer, "_", gene, "_", modal_label)
      feature_label <- paste0(gene, " (", modal_label, ", ", cancer, ")")

      extracted <- data.frame(
        value = as.numeric(gene_data[[feature_col]]),
        row.names = gene_data$patient_barcode_16_clinical
      )
      colnames(extracted) <- new_col_name

      all_data[[feature_label]] <- extracted
      all_features <- c(all_features, feature_label)
    }
  }

  if (length(all_data) == 0) {
    stop("No valid Methylation data loaded", call. = FALSE)
  }

  merged <- .merge_list_data(all_data)

  # Add cancer type info (with subtype support)
  merged$cancer_type <- .assign_cancer_type(rownames(merged), cancers)

  types <- rep("continuous", length(all_features))
  names(types) <- all_features

  message(sprintf("  Loaded %d Methylation feature(s) [%s]", length(all_features), methylation_region))

  return(list(data = merged, features = all_features, types = types))
}

# ==============================================================================
# miRNA Loader
# ==============================================================================

#' Load miRNA data
#' @keywords internal
.load_mirna_modal_tcga <- function(vars, cancers) {
  base_path <- Sys.getenv("SL_BULK_DATA")
  mirna_file <- file.path(base_path, "TCGA_miRNA_pancancer_scores.qs")

  if (!file.exists(mirna_file)) {
    stop("miRNA file not found: ", mirna_file, call. = FALSE)
  }

  mirna_data <- qs::qread(mirna_file, nthreads = min(parallel::detectCores(), 6))

  # Add cancer type info (with subtype support)
  # Note: miRNA data doesn't have subtype filtering, so use parent cancer
  clin <- .load_clinical_basic()
  mirna_data$cancer_type <- clin[rownames(mirna_data), "cancertype_tcga"]

  all_data <- list()
  all_features <- character(0)

  for (cancer in cancers) {
    cancer_info <- .parse_cancer_type(cancer)
    cancer_mask <- .match_cancer_parent(mirna_data$cancer_type, cancer_info)
    cancer_mirna <- mirna_data[cancer_mask, ]

    if (nrow(cancer_mirna) == 0) {
      warning(sprintf("No miRNA data for %s (skipping)", cancer))
      next
    }

    # Subtype filtering
    if (cancer_info$type == "subtype") {
      subtype_samples <- .filter_subtype_samples(cancer_info)
      cancer_mirna <- cancer_mirna[rownames(cancer_mirna) %in% subtype_samples, ]
    }

    if (nrow(cancer_mirna) == 0) next

    for (mirna in vars) {
      if (!mirna %in% colnames(mirna_data)) {
        warning(sprintf("miRNA %s not found (skipping)", mirna))
        next
      }

      new_col_name <- paste0(cancer, "_", mirna, "_miRNA")
      feature_label <- paste0(mirna, " (miRNA, ", cancer, ")")

      extracted <- data.frame(
        value = as.numeric(cancer_mirna[[mirna]]),
        row.names = rownames(cancer_mirna)
      )
      colnames(extracted) <- new_col_name

      all_data[[feature_label]] <- extracted
      all_features <- c(all_features, feature_label)
    }
  }

  if (length(all_data) == 0) {
    stop("No valid miRNA data loaded", call. = FALSE)
  }

  merged <- .merge_list_data(all_data)

  types <- rep("continuous", length(all_features))
  names(types) <- all_features

  message(sprintf("  Loaded %d miRNA feature(s)", length(all_features)))

  return(list(data = merged, features = all_features, types = types))
}

# ==============================================================================
# Clinical Loader (66 traditional clinical variables)
# ==============================================================================

#' Clinical variable aliases (smart matching)
#' @keywords internal
CLINICAL_ALIASES <- list(
  # Demographic
  "Age" = "demographic_age",
  "Gender" = "demographic_gender",
  "Sex" = "demographic_gender",
  "Race" = "demographic_race",
  "BMI" = "demographic_bmi",
  "Height" = "demographic_height",
  "Weight" = "demographic_weight",
  "VitalStatus" = "demographic_vital_status",
  "Smoking" = "demographic_tobacco_smoking_history",
  "Alcohol" = "demographic_alcohol_history_documented",

  # Stage
  "Stage" = "ajcc_pathologic_stage_clean",
  "T" = "ajcc_pathologic_t_clean",
  "N" = "ajcc_pathologic_n_clean",
  "M" = "ajcc_pathologic_m_clean",
  "StageNum" = "ajcc_pathologic_stage_num",
  "TNum" = "ajcc_pathologic_t_num",
  "NNum" = "ajcc_pathologic_n_num",
  "MNum" = "ajcc_pathologic_m_num",

  # Tumor
  "Grade" = "histological_grade",
  "Type" = "histological_type",
  "ER" = "tumor_er_status_nature2012",
  "PR" = "tumor_pr_status_nature2012",
  "HER2" = "tumor_her2_status_nature2012",
  "MSI" = "tumor_msi_status",
  "TIL" = "tumor_infiltrating_lymphocytes",
  "Residual" = "tumor_residual",

  # Treatment
  "Radiation" = "treatment_radiation_therapy",
  "Outcome" = "treatment_outcome_first_course_clean",
  "Neoadjuvant" = "treatment_history_neoadjuvant_clean"
)

#' Clinical variable groups
#' @keywords internal
CLINICAL_GROUPS <- list(
  "basic" = c("Age", "Gender", "Race", "BMI"),
  "stage" = c("Stage", "T", "N", "M"),
  "tumor" = c("Grade", "ER", "PR", "HER2", "MSI", "TIL"),
  "treatment" = c("Radiation", "Outcome", "Neoadjuvant")
)

#' Load clinical data
#' @keywords internal
.load_clinical_modal_tcga <- function(vars, cancers) {
  clin_data <- .load_clinical_basic()

  # Define clinical columns (exclude signature columns)
  CLINICAL_COLUMNS <- c(
    grep("^demographic_", colnames(clin_data), value = TRUE),
    grep("^ajcc_", colnames(clin_data), value = TRUE),
    grep("^tumor_(residual|laterality|tissue|infiltrating|msi|er|pr|her2)",
      colnames(clin_data),
      value = TRUE
    ),
    grep("^treatment_", colnames(clin_data), value = TRUE),
    grep("^invasion_", colnames(clin_data), value = TRUE),
    grep("^score_(karnofsky|gleason)", colnames(clin_data), value = TRUE),
    grep("^lymphnodes_", colnames(clin_data), value = TRUE),
    "histological_type", "histological_grade"
  )

  CLINICAL_COLUMNS <- unique(CLINICAL_COLUMNS)

  all_data <- list()
  all_features <- character(0)
  all_types <- character(0)

  for (cancer in cancers) {
    cancer_info <- .parse_cancer_type(cancer)
    cancer_mask <- .match_cancer_parent(clin_data$cancertype_tcga, cancer_info)
    cancer_clin <- clin_data[cancer_mask, ]

    # Subtype filtering
    if (cancer_info$type == "subtype") {
      # Use %in% to avoid NA issues (== returns NA when comparing with NA)
      filter_match <- cancer_clin[[cancer_info$filter_col]] %in% cancer_info$filter_value
      cancer_clin <- cancer_clin[filter_match & !is.na(filter_match), ]
    }

    if (nrow(cancer_clin) == 0) {
      warning(sprintf("No clinical data for %s (skipping)", cancer))
      next
    }

    for (var in vars) {
      # Match column (support alias)
      matched_col <- .match_clinical_col_tcga(var, CLINICAL_COLUMNS)

      if (is.null(matched_col)) {
        warning(sprintf("Clinical variable '%s' not found (skipping)", var))
        similar <- .suggest_similar(var, names(CLINICAL_ALIASES))
        if (length(similar) > 0) {
          message(sprintf("  Did you mean: %s?", paste(similar, collapse = ", ")))
        }
        next
      }

      values <- cancer_clin[[matched_col]]

      # Categorize
      categorized <- .categorize_clinical_tcga(values, matched_col)

      new_col_name <- paste0(cancer, "_", var, "_Clinical")
      feature_label <- paste0(var, " (Clinical, ", cancer, ")")

      extracted <- data.frame(
        value = categorized$values,
        row.names = cancer_clin$patient_barcode_16
      )
      colnames(extracted) <- new_col_name

      all_data[[feature_label]] <- extracted
      all_features <- c(all_features, feature_label)
      all_types <- c(all_types, categorized$type)
    }
  }

  if (length(all_data) == 0) {
    stop("No valid Clinical data loaded", call. = FALSE)
  }

  merged <- .merge_list_data(all_data)
  names(all_types) <- all_features

  message(sprintf("  Loaded %d Clinical feature(s)", length(all_features)))

  return(list(data = merged, features = all_features, types = all_types))
}

#' Match clinical column name
#' @keywords internal
.match_clinical_col_tcga <- function(search, available_cols) {
  # First check if it's an alias
  if (search %in% names(CLINICAL_ALIASES)) {
    return(CLINICAL_ALIASES[[search]])
  }

  # Then check direct match (case insensitive)
  search_lower <- tolower(search)
  available_lower <- tolower(available_cols)

  idx <- which(available_lower == search_lower)
  if (length(idx) > 0) {
    return(available_cols[idx[1]])
  }

  # Partial match
  idx <- grep(paste0("^", search_lower), available_lower)
  if (length(idx) > 0) {
    return(available_cols[idx[1]])
  }

  return(NULL)
}

#' Categorize clinical values
#' @keywords internal
.categorize_clinical_tcga <- function(values, var_name) {
  var_lower <- tolower(var_name)

  # Age: < 60 vs >= 60 (match "_age" at end or "age" as whole word to avoid "stage")
  if (grepl("_age$|^age$", var_lower)) {
    if (is.numeric(values)) {
      cutoff <- 60
      categorized <- ifelse(values < cutoff, paste0("<", cutoff), paste0(">=", cutoff))
      return(list(values = factor(categorized), type = "categorical"))
    }
  }

  # BMI: < 25 vs >= 25 (match "_bmi" at end or "bmi" as whole word)
  if (grepl("_bmi$|^bmi$", var_lower)) {
    if (is.numeric(values)) {
      cutoff <- 25
      categorized <- ifelse(values < cutoff, paste0("<", cutoff), paste0(">=", cutoff))
      return(list(values = factor(categorized), type = "categorical"))
    }
  }

  # Numeric variables (keep as continuous)
  if (is.numeric(values)) {
    return(list(values = as.numeric(values), type = "continuous"))
  }

  # Categorical variables
  return(list(values = factor(values), type = "categorical"))
}

# ==============================================================================
# Signature Loader (58 molecular signature variables)
# ==============================================================================

#' Signature variable aliases
#' @keywords internal
SIGNATURE_ALIASES <- list(
  # Immunity2018 (most commonly used)
  "ImmuneSubtype" = "immunity2018_immune_subtype",
  "Leukocyte" = "immunity2018_leukocyte_fraction",
  "Stromal" = "immunity2018_stromal_fraction",
  "TIL_Score" = "immunity2018_lymphocyte_infiltration_signature_score",
  "IFN_Gamma" = "immunity2018_ifn_gamma_response",
  "TGF_Beta" = "immunity2018_tgf_beta_response",
  "Proliferation" = "immunity2018_proliferation",
  "Neoantigen_SNV" = "immunity2018_snv_neoantigens",
  "Neoantigen_Indel" = "immunity2018_indel_neoantigens",
  "Aneuploidy" = "immunity2018_aneuploidy_score",
  "TCR_Diversity" = "immunity2018_tcr_shannon",
  "BCR_Diversity" = "immunity2018_bcr_shannon",
  "Heterogeneity" = "immunity2018_intratumor_heterogeneity",

  # Cell2019
  "Purity" = "cell2019_purity",
  "Ploidy" = "cell2019_ploidy",
  "GenomeDoublings" = "cell2019_genome_doublings",

  # PanCanAtlas2018
  "TMB" = "pancanatlas2018_tmb_nonsynonymous",
  "MSI" = "pancanatlas2018_msi_score_mantis", # Default: Mantis
  "MSI_Mantis" = "pancanatlas2018_msi_score_mantis",
  "MSI_Sensor" = "pancanatlas2018_msi_score_sensor",
  "MutationCount" = "pancanatlas2018_mutation_count",
  "Hypoxia" = "pancanatlas2018_buffa_hypoxia_score", # Default: Buffa
  "Hypoxia_Buffa" = "pancanatlas2018_buffa_hypoxia_score",
  "Hypoxia_Winter" = "pancanatlas2018_winter_hypoxia_score",
  "Hypoxia_Ragnum" = "pancanatlas2018_ragnum_hypoxia_score",

  # UCSCXena - HRD
  "HRD" = "ucscxena_hrd", # Combined HRD score
  "TAI" = "ucscxena_tai", # Telomeric allelic imbalance
  "LST" = "ucscxena_lst", # Large-scale state transitions
  "LOH" = "ucscxena_loh", # Loss of heterozygosity

  # UCSCXena - Stemness
  "Stemness" = "ucscxena_stemness_score_rnass", # Default: RNA-based
  "Stemness_RNA" = "ucscxena_stemness_score_rnass",
  "Stemness_DNA" = "ucscxena_stemness_score_dnass",
  "Stemness_DMP" = "ucscxena_stemness_score_dmpss",
  "Stemness_ENH" = "ucscxena_stemness_score_enhss",

  # Firehose
  "FGA" = "firehose_fga",
  "FGG" = "firehose_fgg",
  "FGL" = "firehose_fgl"
)

#' Signature variable groups
#' @keywords internal
SIGNATURE_GROUPS <- list(
  "immune" = c("Leukocyte", "Stromal", "TIL_Score", "IFN_Gamma", "TGF_Beta"),
  "purity" = c("Purity", "Ploidy", "GenomeDoublings"),
  "mutation" = c("TMB", "MSI_Mantis", "Neoantigen_SNV", "MutationCount"),
  "stemness" = c("Stemness_RNA", "Stemness_DNA"),
  "hrd" = c("HRD", "TAI", "LST", "LOH"),
  "hypoxia" = c("Hypoxia_Buffa", "Hypoxia_Winter", "Hypoxia_Ragnum"),
  "genomic" = c("FGA", "Aneuploidy", "Heterogeneity"),
  "diversity" = c("TCR_Diversity", "BCR_Diversity")
)

#' Load signature data
#' @keywords internal
.load_signature_modal_tcga <- function(vars, cancers) {
  clin_data <- .load_clinical_basic()

  # Define signature columns (58 total)
  SIGNATURE_COLUMNS <- c(
    grep("^immunity2018_", colnames(clin_data), value = TRUE),
    grep("^cell2019_", colnames(clin_data), value = TRUE),
    grep("^pancanatlas2018_", colnames(clin_data), value = TRUE),
    grep("^ucscxena_(tai|lst|loh|hrd|stemness)", colnames(clin_data), value = TRUE),
    grep("^firehose_", colnames(clin_data), value = TRUE)
  )

  all_data <- list()
  all_features <- character(0)
  all_types <- character(0)

  for (cancer in cancers) {
    cancer_info <- .parse_cancer_type(cancer)
    cancer_mask <- .match_cancer_parent(clin_data$cancertype_tcga, cancer_info)
    cancer_clin <- clin_data[cancer_mask, ]

    # Subtype filtering
    if (cancer_info$type == "subtype") {
      # Use %in% to avoid NA issues (== returns NA when comparing with NA)
      filter_match <- cancer_clin[[cancer_info$filter_col]] %in% cancer_info$filter_value
      cancer_clin <- cancer_clin[filter_match & !is.na(filter_match), ]
    }

    if (nrow(cancer_clin) == 0) {
      warning(sprintf("No data for %s (skipping)", cancer))
      next
    }

    for (var in vars) {
      # Match column (support alias)
      matched_col <- .match_signature_col_tcga(var, SIGNATURE_COLUMNS)

      if (is.null(matched_col)) {
        warning(sprintf("Signature '%s' not found (skipping)", var))
        similar <- .suggest_similar(var, names(SIGNATURE_ALIASES))
        if (length(similar) > 0) {
          message(sprintf("  Did you mean: %s?", paste(similar, collapse = ", ")))
        }
        message("  Use list_variables(modal='Signature') to see all available signatures.")
        next
      }

      values <- cancer_clin[[matched_col]]

      # Determine type
      if (is.numeric(values)) {
        var_type <- "continuous"
        processed_values <- as.numeric(values)
      } else {
        var_type <- "categorical"
        processed_values <- factor(values)
      }

      new_col_name <- paste0(cancer, "_", var, "_Signature")
      feature_label <- paste0(var, " (Signature, ", cancer, ")")

      extracted <- data.frame(
        value = processed_values,
        row.names = cancer_clin$patient_barcode_16
      )
      colnames(extracted) <- new_col_name

      all_data[[feature_label]] <- extracted
      all_features <- c(all_features, feature_label)
      all_types <- c(all_types, var_type)
    }
  }

  if (length(all_data) == 0) {
    stop("No valid Signature data loaded", call. = FALSE)
  }

  merged <- .merge_list_data(all_data)
  names(all_types) <- all_features

  message(sprintf("  Loaded %d Signature feature(s)", length(all_features)))

  return(list(data = merged, features = all_features, types = all_types))
}

#' Match signature column name
#' @keywords internal
.match_signature_col_tcga <- function(search, available_cols) {
  # First check if it's an alias
  if (search %in% names(SIGNATURE_ALIASES)) {
    return(SIGNATURE_ALIASES[[search]])
  }

  # Then check direct match
  if (search %in% available_cols) {
    return(search)
  }

  # Case insensitive match
  search_lower <- tolower(search)
  available_lower <- tolower(available_cols)

  idx <- which(available_lower == search_lower)
  if (length(idx) > 0) {
    return(available_cols[idx[1]])
  }

  return(NULL)
}

# ==============================================================================
# ImmuneCell Loader (99 cell types from 8 algorithms)
# ==============================================================================

#' Immune cell categories
#' @keywords internal
IMMUNE_CELL_CATEGORIES <- list(
  "B_cells" = c(
    "B_cells_naive_cibersort", "B_cells_memory_cibersort",
    "B_cells_mcpcounter", "B_cells_quantiseq", "B_cells_timer",
    "B_cells_xcell", "B_cells_Class_switched_memory_xcell",
    "B_cells_Memory_xcell", "B_cells_naive_xcell", "B_cells_pro_xcell"
  ),
  "T_cells_CD4" = c(
    "T_cells_CD4_naive_cibersort", "T_cells_CD4_memory_resting_cibersort",
    "T_cells_CD4_memory_activated_cibersort", "T_cells_follicular_helper_cibersort",
    "T_cells_CD4_epic", "T_cells_CD4_quantiseq", "T_cell_CD4_timer",
    "T_cells_CD4_memory_xcell", "T_cells_CD4_naive_xcell", "T_cells_CD4_xcell",
    "T_cells_CD4_Tcm_xcell", "T_cells_CD4_Tem_xcell",
    "T_cells_CD4_Th1_xcell", "T_cells_CD4_Th2_xcell"
  ),
  "T_cells_CD8" = c(
    "T_cells_CD8_cibersort", "T_cells_CD8_epic", "T_cells_CD8_mcpcounter",
    "T_cells_CD8_quantiseq", "T_cells_CD8_timer",
    "T_cells_CD8_naive_xcell", "T_cells_CD8_xcell",
    "T_cells_CD8_Tcm_xcell", "T_cells_CD8_Tem_xcell"
  ),
  "Tregs" = c("Tregs_cibersort", "Tregs_quantiseq", "Tregs_xcell"),
  "NK_cells" = c(
    "NK_cells_resting_cibersort", "NK_cells_activated_cibersort",
    "NK_cells_epic", "NK_cells_mcpcounter", "NK_cells_quantiseq",
    "NK_cells_xcell", "NKT_cells_xcell"
  ),
  "Macrophages" = c(
    "Macrophages_M0_cibersort", "Macrophages_M1_cibersort", "Macrophages_M2_cibersort",
    "Macrophages_epic", "Macrophages_M1_quantiseq", "Macrophages_M2_quantiseq",
    "Macrophages_timer", "Macrophages_xcell",
    "Macrophages_M1_xcell", "Macrophages_M2_xcell"
  ),
  "DC" = c(
    "DC_resting_cibersort", "DC_activated_cibersort",
    "DC_myeloid_mcpcounter", "DC_quantiseq", "DC_timer",
    "aDC_xcell", "cDC_xcell", "DC_xcell", "iDC_xcell", "pDC_xcell"
  ),
  "Neutrophils" = c(
    "Neutrophils_cibersort", "Neutrophils_mcpcounter",
    "Neutrophils_quantiseq", "Neutrophils_timer", "Neutrophils_xcell"
  ),
  "Monocytes" = c(
    "Monocytes_cibersort", "Monocytes_mcpcounter",
    "Monocytes_quantiseq", "Monocytes_xcell"
  ),
  "Microenvironment" = c(
    "Endothelial_cells_epic", "Endothelial_cells_mcpcounter", "Endothelial_cells_xcell",
    "Endothelial_cells_ly_xcell", "Endothelial_cells_mv_xcell",
    "Fibroblasts_CAFs_epic", "Fibroblasts_mcpcounter", "Fibroblasts_xcell"
  ),
  "ESTIMATE" = c("StromalScore_estimate", "ImmuneScore_estimate", "ESTIMATEScore_estimate")
)

#' Load immune cell data
#' @keywords internal
.load_immunecell_modal_tcga <- function(vars, cancers, immune_algorithm = NULL) {
  base_path <- Sys.getenv("SL_BULK_DATA")
  immune_file <- file.path(base_path, "TCGA_DeconvCell_All_scores.qs")

  if (!file.exists(immune_file)) {
    stop("ImmuneCell file not found: ", immune_file, call. = FALSE)
  }

  immune_data <- qs::qread(immune_file, nthreads = min(parallel::detectCores(), 6))

  # Add cancer type info (with subtype support)
  # Note: ImmuneCell data doesn't have subtype info yet, so use parent cancer
  clin <- .load_clinical_basic()
  immune_data$cancer_type <- clin[rownames(immune_data), "cancertype_tcga"]

  all_data <- list()
  all_features <- character(0)

  for (cancer in cancers) {
    cancer_info <- .parse_cancer_type(cancer)
    cancer_mask <- .match_cancer_parent(immune_data$cancer_type, cancer_info)
    cancer_immune <- immune_data[cancer_mask, ]

    if (nrow(cancer_immune) == 0) {
      warning(sprintf("No ImmuneCell data for %s (skipping)", cancer))
      next
    }

    # Subtype filtering
    if (cancer_info$type == "subtype") {
      subtype_samples <- .filter_subtype_samples(cancer_info)
      cancer_immune <- cancer_immune[rownames(cancer_immune) %in% subtype_samples, ]
    }

    if (nrow(cancer_immune) == 0) next

    # Handle special keyword AFTER filtering by cancer
    vars_to_load <- vars
    if (length(vars) == 1 && vars == "ALL_IMMUNE_CELLS") {
      if (!is.null(immune_algorithm)) {
        vars_to_load <- grep(paste0("_", immune_algorithm, "$"),
          colnames(cancer_immune),
          value = TRUE
        )
        vars_to_load <- vars_to_load[vars_to_load != "cancer_type"]
        if (cancer == cancers[1]) { # Only message once
          message(sprintf("  Using ALL %s cells (%d cells)", immune_algorithm, length(vars_to_load)))
        }
      } else {
        vars_to_load <- colnames(cancer_immune)[colnames(cancer_immune) != "cancer_type"]
        if (cancer == cancers[1]) {
          message(sprintf("  Using ALL immune cells (%d cells)", length(vars_to_load)))
        }
      }
    }

    for (var in vars_to_load) {
      # Smart matching (use cancer_immune columns, not full immune_data)
      matched_cols <- .match_immune_cell_col(var, colnames(cancer_immune), immune_algorithm)

      if (length(matched_cols) == 0) {
        warning(sprintf("Immune cell '%s' not found (skipping)", var))
        next
      }

      # Load all matched columns (e.g., "B_cells" matches multiple algorithms)
      for (col in matched_cols) {
        new_col_name <- paste0(cancer, "_", col, "_ImmuneCell")
        feature_label <- paste0(col, " (ImmuneCell, ", cancer, ")")

        extracted <- data.frame(
          value = as.numeric(cancer_immune[[col]]),
          row.names = rownames(cancer_immune)
        )
        colnames(extracted) <- new_col_name

        all_data[[feature_label]] <- extracted
        all_features <- c(all_features, feature_label)
      }
    }
  }

  if (length(all_data) == 0) {
    stop("No valid ImmuneCell data loaded", call. = FALSE)
  }

  merged <- .merge_list_data(all_data)

  types <- rep("continuous", length(all_features))
  names(types) <- all_features

  message(sprintf("  Loaded %d ImmuneCell feature(s)", length(all_features)))

  return(list(data = merged, features = all_features, types = types))
}

#' Match immune cell column name
#' @keywords internal
.match_immune_cell_col <- function(search, available, algorithm = NULL) {
  # Complete match first
  if (search %in% available) {
    return(search)
  }

  # Prefix match + algorithm filter
  pattern <- paste0("^", gsub("_", ".*", search))
  matches <- grep(pattern, available, value = TRUE, ignore.case = TRUE)

  # Filter by algorithm if specified
  if (!is.null(algorithm)) {
    matches <- grep(paste0("_", algorithm, "$"), matches, value = TRUE)
  }

  # If still no matches, try fuzzy search
  if (length(matches) == 0) {
    search_clean <- gsub("[_-]", "", tolower(search))
    available_clean <- gsub("[_-]", "", tolower(available))
    idx <- grep(search_clean, available_clean)
    matches <- available[idx]

    if (!is.null(algorithm)) {
      matches <- grep(paste0("_", algorithm, "$"), matches, value = TRUE)
    }
  }

  return(matches)
}

# ==============================================================================
# Survival Loader (4 survival types)
# ==============================================================================

#' Load survival data
#' @keywords internal
.load_survival_modal_tcga <- function(cancers, surv_type) {
  clin_data <- .load_clinical_basic()

  # Map survival type to column names
  surv_map <- list(
    "OS" = c(time = "os_time_years", event = "os_event"),
    "DSS" = c(time = "dss_time_years", event = "dss_event"),
    "PFI" = c(time = "pfi_time_years", event = "pfi_event"),
    "DFI" = c(time = "dfi_time_years", event = "dfi_event"),
    # Aliases for common usage
    "PFS" = c(time = "pfi_time_years", event = "pfi_event"), # PFS -> PFI
    "RFS" = c(time = "dfi_time_years", event = "dfi_event") # RFS -> DFI
  )

  if (!surv_type %in% names(surv_map)) {
    stop("Unknown surv_type: ", surv_type, ". Use: OS, DSS, PFI (or PFS), DFI (or RFS)", call. = FALSE)
  }

  time_col <- surv_map[[surv_type]]["time"]
  event_col <- surv_map[[surv_type]]["event"]

  all_data <- list()
  all_features <- character(0)

  for (cancer in cancers) {
    cancer_info <- .parse_cancer_type(cancer)
    cancer_mask <- .match_cancer_parent(clin_data$cancertype_tcga, cancer_info)
    cancer_clin <- clin_data[cancer_mask, ]

    # Subtype filtering
    if (cancer_info$type == "subtype") {
      # Use %in% to avoid NA issues (== returns NA when comparing with NA)
      filter_match <- cancer_clin[[cancer_info$filter_col]] %in% cancer_info$filter_value
      cancer_clin <- cancer_clin[filter_match & !is.na(filter_match), ]
    }

    if (nrow(cancer_clin) == 0) next

    # Extract survival columns
    time_col_name <- paste0(cancer, "_", surv_type, "_time")
    event_col_name <- paste0(cancer, "_", surv_type, "_event")

    surv_df <- data.frame(
      time = cancer_clin[[time_col]],
      event = cancer_clin[[event_col]],
      row.names = cancer_clin$patient_barcode_16
    )
    colnames(surv_df) <- c(time_col_name, event_col_name)

    all_data[[cancer]] <- surv_df
    all_features <- c(all_features, paste0(surv_type, " (Survival, ", cancer, ")"))
  }

  if (length(all_data) == 0) {
    stop("No survival data loaded", call. = FALSE)
  }

  merged <- .merge_list_data(all_data)

  types <- rep("survival", length(all_features))
  names(types) <- all_features

  message(sprintf("  Loaded %s survival data for %d cancer(s)", surv_type, length(cancers)))

  return(list(data = merged, features = all_features, types = types))
}

# ==============================================================================
# Data Merging
# ==============================================================================

#' Merge two modal datasets
#' @keywords internal
.merge_modal_data_tcga <- function(data1, data2) {
  if (nrow(data1) == 0) {
    stop("No data available for var1", call. = FALSE)
  }

  if (ncol(data2) == 0) {
    return(data1)
  }

  if (nrow(data2) == 0) {
    stop("No data available for var2", call. = FALSE)
  }

  # Check for overlapping columns
  overlap_cols <- intersect(colnames(data1), colnames(data2))

  if (length(overlap_cols) > 0) {
    # Union merge (keep all samples, exclude NA rownames)
    all_samples <- union(rownames(data1), rownames(data2))
    all_samples <- all_samples[!is.na(all_samples) & all_samples != "NA"] # Remove NA and "NA" string
    merged <- data.frame(row.names = all_samples)

    for (col in colnames(data1)) {
      merged[[col]] <- data1[rownames(merged), col]
    }

    for (col in setdiff(colnames(data2), overlap_cols)) {
      merged[[col]] <- data2[rownames(merged), col]
    }

    return(merged)
  } else {
    # Standard merge (exclude NA rownames)
    merged <- merge(data1, data2, by = "row.names", all = TRUE)
    rownames(merged) <- merged$Row.names
    merged$Row.names <- NULL

    # Remove rows with NA or "NA" string rownames
    valid_rows <- !is.na(rownames(merged)) & rownames(merged) != "NA"
    merged <- merged[valid_rows, ]

    return(merged)
  }
}

# ==============================================================================
# Genome-wide Data Loader (for enrichment analysis)
# ==============================================================================

#' Load genome-wide RNAseq data for enrichment analysis
#' @keywords internal
.load_genome_data_tcga <- function(cancers, rnaseq_type = "log2TPM") {
  base_path <- Sys.getenv("SL_BULK_DATA")

  # Use SplitCancer files (faster for full genome loading)
  split_path <- file.path(base_path, "TCGA_RNAseq_log2TPM_SplitCancer")

  if (!dir.exists(split_path)) {
    stop("RNAseq_log2TPM_SplitCancer directory not found", call. = FALSE)
  }

  all_data <- list()

  for (cancer in cancers) {
    cancer_info <- .parse_cancer_type(cancer)
    parent_targets <- if (cancer_info$type == "combined") cancer_info$parents else cancer_info$cancer

    for (parent_cancer in parent_targets) {
      cancer_file <- file.path(split_path, paste0(parent_cancer, ".qs"))

      if (!file.exists(cancer_file)) {
        warning(sprintf("RNAseq file not found for %s (component: %s)", cancer, parent_cancer))
        next
      }

      cancer_data <- qs::qread(cancer_file, nthreads = min(parallel::detectCores(), 6))

      cancer_data_t <- t(cancer_data)

      if (cancer_info$type == "subtype") {
        subtype_samples <- .filter_subtype_samples(cancer_info)
        cancer_data_t <- cancer_data_t[rownames(cancer_data_t) %in% subtype_samples, , drop = FALSE]
      }

      if (nrow(cancer_data_t) == 0) next

      list_key <- paste(cancer, parent_cancer, sep = "__")
      all_data[[list_key]] <- cancer_data_t
    }
  }

  if (length(all_data) == 0) {
    stop("No genome-wide data loaded", call. = FALSE)
  }

  # Combine all cancers (sample × gene)
  combined <- do.call(rbind, all_data)

  # Transpose back to gene × sample matrix (for DEA/correlation)
  genome_matrix <- t(combined)

  message(sprintf("  Genome-wide: %d genes × %d samples", nrow(genome_matrix), ncol(genome_matrix)))

  return(genome_matrix)
}
