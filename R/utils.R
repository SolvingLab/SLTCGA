# ==============================================================================
# Utility Functions Layer
# ==============================================================================
# Cancer type management, data availability validation, helper functions
# ==============================================================================

# ==============================================================================
# Cancer Type Management
# ==============================================================================

#' Main cancer types (33)
#' @keywords internal
TCGA_MAIN_CANCERS <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
  "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

#' Subtype mapping (32 subtypes)
#' @keywords internal
TCGA_SUBTYPE_MAP <- list(
  # BRCA subtypes (3)
  "BRCA_IDC" = list(parent = "BRCA", clinical_col = "tumor_subtype_brca_idc", value = "IDC"),
  "BRCA_ILC" = list(parent = "BRCA", clinical_col = "tumor_subtype_brca_ilc", value = "ILC"),
  "BRCA_TNBC" = list(parent = "BRCA", clinical_col = "tumor_subtype_brca_tnbc", value = "TNBC"),

  # CESC subtypes (1)
  "CESC_CSCC" = list(parent = "CESC", clinical_col = "tumor_subtype_cesc_cscc", value = "CSCC"),

  # COAD subtypes (3)
  "COAD_LCC" = list(parent = "COAD", clinical_col = "tumor_subtype_coad_lcc", value = "LCC"),
  "COAD_RCC" = list(parent = "COAD", clinical_col = "tumor_subtype_coad_rcc", value = "RCC"),
  "COAD_MAC" = list(parent = "COAD", clinical_col = "tumor_subtype_coad_mac", value = "MAC"),

  # ESCA subtypes (2)
  "ESCA_ESCC" = list(parent = "ESCA", clinical_col = "tumor_subtype_esca_escc", value = "ESCC"),
  "ESCA_EAC" = list(parent = "ESCA", clinical_col = "tumor_subtype_esca_eac", value = "EAC"),

  # HNSC subtypes (2)
  "HNSC_OSCC" = list(parent = "HNSC", clinical_col = "tumor_subtype_hnsc_oscc", value = "OSCC"),
  "HNSC_LSCC" = list(parent = "HNSC", clinical_col = "tumor_subtype_hnsc_lscc", value = "LSCC"),

  # LGG subtypes (3)
  "LGG_ASTROCYTOMA" = list(parent = "LGG", clinical_col = "tumor_subtype_lgg_astrocytoma", value = "Astrocytoma"),
  "LGG_OLIGOASTROCYTOMA" = list(parent = "LGG", clinical_col = "tumor_subtype_lgg_oligoastrocytoma", value = "Oligoastrocytoma"),
  "LGG_OLIGODENDROGLIOMA" = list(parent = "LGG", clinical_col = "tumor_subtype_lgg_oligodendroglioma", value = "Oligodendroglioma"),

  # PAAD subtypes (1)
  "PAAD_PADC" = list(parent = "PAAD", clinical_col = "tumor_subtype_paad_padc", value = "PADC"),

  # PCPG subtypes (2)
  "PCPG_PHEOCHROMOCYTOMA" = list(parent = "PCPG", clinical_col = "tumor_subtype_pcpg_pheochromocytoma", value = "Pheochromocytoma"),
  "PCPG_PARAGANGLIOMA" = list(parent = "PCPG", clinical_col = "tumor_subtype_pcpg_paraganglioma", value = "Paraganglioma"),

  # SARC subtypes (2)
  "SARC_DDLPS" = list(parent = "SARC", clinical_col = "tumor_subtype_sarc_ddlps", value = "DDLPS"),
  "SARC_LMS" = list(parent = "SARC", clinical_col = "tumor_subtype_sarc_lms", value = "LMS"),

  # SKCM subtypes (2)
  "SKCM_MCM" = list(parent = "SKCM", clinical_col = "tumor_subtype_skcm_mcm", value = "MCM"),
  "SKCM_PCM" = list(parent = "SKCM", clinical_col = "tumor_subtype_skcm_pcm", value = "PCM"),

  # STAD subtypes (2)
  "STAD_DGA" = list(parent = "STAD", clinical_col = "tumor_subtype_stad_dga", value = "DGA"),
  "STAD_SRCC" = list(parent = "STAD", clinical_col = "tumor_subtype_stad_srcc", value = "SRCC"),

  # TGCT subtypes (1)
  "TGCT_SEMINOMA" = list(parent = "TGCT", clinical_col = "tumor_subtype_tgct_seminoma", value = "Seminoma"),

  # THCA subtypes (1)
  "THCA_CPTC" = list(parent = "THCA", clinical_col = "tumor_subtype_thca_cptc", value = "cPTC"),

  # UCEC subtypes (2)
  "UCEC_EEA" = list(parent = "UCEC", clinical_col = "tumor_subtype_ucec_eea", value = "EEA"),
  "UCEC_SEA" = list(parent = "UCEC", clinical_col = "tumor_subtype_ucec_sea", value = "SEA")
)

# Combined cancer groups (alias for multiple parent cancers)
TCGA_COMBINED_MAP <- list(
  "NSCLC" = c("LUAD", "LUSC"),
  "CRC" = c("COAD", "READ"),
  "GLIOMA" = c("GBM", "LGG"),
  "KRCC" = c("KIRC", "KIRP", "KICH"),
  "GIC" = c("ESCA", "STAD", "COAD", "READ")
)

#' Parse cancer type input
#' @param cancer_input User input cancer type (case insensitive)
#' @return List with type, cancer, subtype info
#' @keywords internal
.parse_cancer_type <- function(cancer_input) {
  cancer_upper <- toupper(cancer_input)

  # Check if main cancer type
  if (cancer_upper %in% TCGA_MAIN_CANCERS) {
    return(list(
      type = "main",
      cancer = cancer_upper,
      subtype = NULL,
      filter_col = NULL,
      filter_value = NULL
    ))
  }

  # Check if subtype
  if (cancer_upper %in% names(TCGA_SUBTYPE_MAP)) {
    subtype_info <- TCGA_SUBTYPE_MAP[[cancer_upper]]
    return(list(
      type = "subtype",
      cancer = subtype_info$parent,
      subtype = cancer_upper,
      filter_col = subtype_info$clinical_col,
      filter_value = subtype_info$value
    ))
  }

  # Check combined cancer groups
  if (cancer_upper %in% names(TCGA_COMBINED_MAP)) {
    return(list(
      type = "combined",
      cancer = cancer_upper,
      parents = TCGA_COMBINED_MAP[[cancer_upper]],
      subtype = NULL,
      filter_col = NULL,
      filter_value = NULL
    ))
  }

  # Not found
  stop(
    "Unknown cancer type: ", cancer_input, "\n",
    "Use list_cancer_types() to see all available cancer types and subtypes.",
    call. = FALSE
  )
}

#' Filter samples by subtype
#' @param cancer_info Cancer info from .parse_cancer_type()
#' @return Vector of sample IDs
#' @keywords internal
.filter_subtype_samples <- function(cancer_info) {
  if (cancer_info$type != "subtype") {
    stop("This function is only for subtypes", call. = FALSE)
  }

  clin <- .load_clinical_basic()

  # Filter by parent cancer and subtype column
  subtype_samples <- clin[
    clin$cancertype_tcga == cancer_info$cancer &
      clin[[cancer_info$filter_col]] == cancer_info$filter_value,
    "patient_barcode_16"
  ]

  subtype_samples <- subtype_samples[!is.na(subtype_samples)]

  if (length(subtype_samples) == 0) {
    warning(sprintf(
      "No samples found for subtype %s (parent: %s, filter: %s = %s)",
      cancer_info$subtype, cancer_info$cancer,
      cancer_info$filter_col, cancer_info$filter_value
    ))
  }

  return(subtype_samples)
}

# ==============================================================================
# Data Availability Validation
# ==============================================================================

#' Validate cancer type and modal availability
#' @param cancers Character vector of cancer types
#' @param modal Character, modal type
#' @param surv_type Character, survival type (for Survival modal)
#' @param stop_on_error Logical, whether to stop on error
#' @return TRUE if all valid
#' @keywords internal
.validate_cancer_modal <- function(cancers, modal, surv_type = NULL, stop_on_error = TRUE) {
  # Parse all cancer types
  cancer_parsed <- lapply(cancers, .parse_cancer_type)

  # Modal availability (all have same availability)
  # RNAseq, Mutation, CNV, Methylation, miRNA, Clinical, Signature, ImmuneCell: ALL
  # (No special restrictions for TCGA)

  # All modals available for all cancer types
  return(TRUE)
}

# ==============================================================================
# Clinical Data Loading (Basic)
# ==============================================================================

#' Load basic clinical data (used by other loaders)
#' @return Data frame with basic clinical info
#' @keywords internal
.load_clinical_basic <- function() {
  base_path <- Sys.getenv("SL_BULK_DATA")
  if (base_path == "") {
    stop("SL_BULK_DATA environment variable not set", call. = FALSE)
  }

  clin_file <- file.path(base_path, "TCGA_Clinical_Final_SovingLab.qs")

  if (!file.exists(clin_file)) {
    stop("Clinical file not found: ", clin_file, call. = FALSE)
  }

  # Cache in environment to avoid repeated loading
  if (!exists(".clinical_cache", envir = .GlobalEnv)) {
    clin <- qs::qread(clin_file, nthreads = min(parallel::detectCores(), 6))
    # Set rownames to patient_barcode_16 for matching with other datasets
    rownames(clin) <- clin$patient_barcode_16
    assign(".clinical_cache", clin, envir = .GlobalEnv)
  }

  return(get(".clinical_cache", envir = .GlobalEnv))
}

# ==============================================================================
# Column Name Extraction
# ==============================================================================

#' Extract column name from feature label
#' @param labels Feature labels in format "GENE (Modal, Cancer)"
#' @param data Data frame to search columns in
#' @return Vector of column names
#' @keywords internal
.extract_colname_from_label <- function(labels, data) {
  colnames_out <- character(length(labels))

  for (i in seq_along(labels)) {
    label <- labels[i]

    # Extract components using regex
    # Format: "GENE (Modal, Cancer)" or "SITE_GENE (Modal, Cancer)"
    pattern <- "^(.+?)\\s+\\((.+?),\\s*(.+?)\\)$"
    matches <- regmatches(label, regexec(pattern, label))[[1]]

    if (length(matches) == 4) {
      gene <- matches[2]
      modal <- matches[3]
      cancer <- matches[4]

      # Construct column name: Cancer_Gene_Modal
      col_name <- paste0(cancer, "_", gene, "_", modal)

      if (col_name %in% colnames(data)) {
        colnames_out[i] <- col_name
      } else {
        warning(sprintf("Column not found for label: %s (expected: %s)", label, col_name))
        colnames_out[i] <- NA
      }
    } else {
      warning(sprintf("Cannot parse label: %s", label))
      colnames_out[i] <- NA
    }
  }

  return(colnames_out)
}

#' Map column names to feature labels (helper version, not used in stats)
#' @param colnames Column names
#' @param modal Modal type
#' @return Vector of feature labels
#' @keywords internal
.map_colname_to_label_simple <- function(colnames, modal = NULL) {
  labels <- character(length(colnames))

  for (i in seq_along(colnames)) {
    col <- colnames[i]

    # Parse column name: Cancer_Gene_Modal
    parts <- strsplit(col, "_")[[1]]

    if (length(parts) >= 3) {
      cancer <- parts[1]
      gene <- paste(parts[2:(length(parts) - 1)], collapse = "_")
      modal_part <- parts[length(parts)]

      # Format: "GENE (Modal, Cancer)"
      labels[i] <- paste0(gene, " (", modal_part, ", ", cancer, ")")
    } else {
      labels[i] <- col
    }
  }

  return(labels)
}

#' Map column name back to feature label (used in stats.R)
#' @keywords internal
.map_colname_to_label <- function(colnames, colname_vector, label_vector) {
  mapped <- character(length(colnames))

  for (i in seq_along(colnames)) {
    idx <- which(colname_vector == colnames[i])
    if (length(idx) > 0) {
      mapped[i] <- label_vector[idx[1]]
    } else {
      mapped[i] <- colnames[i]
    }
  }

  return(mapped)
}

# ==============================================================================
# Data Merging
# ==============================================================================

#' Merge list of data frames
#' @param data_list List of data frames
#' @return Merged data frame
#' @keywords internal
.merge_list_data <- function(data_list) {
  if (length(data_list) == 0) {
    stop("No data to merge", call. = FALSE)
  }

  if (length(data_list) == 1) {
    return(data_list[[1]])
  }

  # Get all unique sample IDs (exclude NA and "NA" string)
  all_samples <- unique(unlist(lapply(data_list, rownames)))
  all_samples <- all_samples[!is.na(all_samples) & all_samples != "NA"] # Remove NA and "NA" string

  # Create empty result data frame
  merged <- data.frame(row.names = all_samples)

  # Add each data frame's columns
  for (df in data_list) {
    for (col in colnames(df)) {
      merged[[col]] <- df[rownames(merged), col]
    }
  }

  return(merged)
}

# ==============================================================================
# Cancer Type Assignment for Subtypes
# ==============================================================================

#' Assign cancer type labels (including subtypes) to samples
#' @keywords internal
.assign_cancer_type <- function(sample_ids, requested_cancers) {
  clin <- .load_clinical_basic()

  cancer_labels <- rep(NA_character_, length(sample_ids))
  names(cancer_labels) <- sample_ids

  for (cancer in requested_cancers) {
    cancer_info <- .parse_cancer_type(cancer)

    if (cancer_info$type == "subtype") {
      subtype_samples <- .filter_subtype_samples(cancer_info)
      matching_samples <- sample_ids[sample_ids %in% subtype_samples]
      cancer_labels[matching_samples] <- cancer
    } else if (cancer_info$type == "combined") {
      present_ids <- intersect(sample_ids, rownames(clin))
      if (length(present_ids) > 0) {
        matches <- present_ids[clin[present_ids, "cancertype_tcga"] %in% cancer_info$parents]
        cancer_labels[matches] <- cancer
      }
    } else {
      # Main cancer: retain parent cancer name
      present_ids <- intersect(sample_ids, rownames(clin))
      if (length(present_ids) > 0) {
        matches <- present_ids[clin[present_ids, "cancertype_tcga"] == cancer_info$cancer]
        cancer_labels[matches] <- cancer_info$cancer
      }
    }
  }

  # For unmatched samples, use parent cancer from clinical
  unmatched <- is.na(cancer_labels)
  if (any(unmatched)) {
    cancer_labels[unmatched] <- clin[sample_ids[unmatched], "cancertype_tcga"]
  }

  return(cancer_labels)
}

#' Helper to match samples to requested cancers (main/subtype/combined)
#' @keywords internal
.match_cancer_parent <- function(values, cancer_info) {
  if (cancer_info$type == "combined") {
    matches <- values %in% cancer_info$parents
  } else {
    matches <- values == cancer_info$cancer
  }
  matches[is.na(matches)] <- FALSE
  matches
}

# ==============================================================================
# Gene Filtering
# ==============================================================================

#' Get protein-coding gene symbols (embedded list)
#' @keywords internal
.get_protein_coding_genes <- function() {
  tryCatch(
    {
      gene_file <- system.file("extdata", "protein_coding_genes.rds", package = "SLTCGA")
      if (gene_file == "") {
        warning("Protein-coding gene list not found in package data")
        return(NULL)
      }
      protein_genes <- readRDS(gene_file)
      message(sprintf("  Filtering to %d protein-coding genes", length(protein_genes)))
      return(protein_genes)
    },
    error = function(e) {
      warning("Failed to load protein-coding gene list: ", e$message)
      return(NULL)
    }
  )
}

# ==============================================================================
# Filename Generation
# ==============================================================================

#' Generate filename for plots
#' @param analysis Analysis type (correlation, enrichment, survival)
#' @param var1 Variable 1
#' @param var1_modal Modal 1
#' @param var1_cancers Cancers 1
#' @param var2 Variable 2 (optional)
#' @param var2_modal Modal 2 (optional)
#' @param var2_cancers Cancers 2 (optional)
#' @return Filename
#' @keywords internal
.generate_filename <- function(analysis, var1, var1_modal, var1_cancers,
                               var2 = NULL, var2_modal = NULL, var2_cancers = NULL) {
  cancer_str <- paste(unique(c(var1_cancers, var2_cancers)), collapse = "-")
  var1_str <- paste(var1, collapse = "-")

  if (!is.null(var2)) {
    # Limit var2_str length to avoid too long filenames
    if (length(var2) > 3) {
      var2_str <- sprintf("%d_%s_genes", length(var2), var2_modal)
    } else {
      var2_str <- paste(var2, collapse = "-")
    }

    filename <- sprintf(
      "%s_%s_%s_%s_vs_%s_%s.png",
      analysis, cancer_str, var1_str, var1_modal, var2_str, var2_modal
    )
  } else {
    filename <- sprintf(
      "%s_%s_%s_%s.png",
      analysis, cancer_str, var1_str, var1_modal
    )
  }

  # Clean filename (remove special characters)
  filename <- gsub("[^a-zA-Z0-9_.-]", "_", filename)

  # Ensure filename is not too long (max 200 characters)
  if (nchar(filename) > 200) {
    # Truncate and add hash
    hash_suffix <- sprintf("_%s.png", substr(digest::digest(filename), 1, 8))
    filename <- paste0(substr(filename, 1, 200 - nchar(hash_suffix)), hash_suffix)
  }

  return(filename)
}

# ==============================================================================
# String Similarity (for suggestions)
# ==============================================================================

#' Calculate string similarity (Levenshtein distance)
#' @param s1 String 1
#' @param s2 String 2
#' @return Similarity score (0-1)
#' @keywords internal
.string_similarity <- function(s1, s2) {
  s1 <- tolower(s1)
  s2 <- tolower(s2)

  # Simple substring check
  if (grepl(s1, s2, fixed = TRUE) || grepl(s2, s1, fixed = TRUE)) {
    return(0.8)
  }

  # Levenshtein distance
  d <- adist(s1, s2)[1, 1]
  max_len <- max(nchar(s1), nchar(s2))

  return(1 - d / max_len)
}

#' Suggest similar strings
#' @param query Query string
#' @param candidates Candidate strings
#' @param n Number of suggestions
#' @return Top n similar strings
#' @keywords internal
.suggest_similar <- function(query, candidates, n = 3) {
  scores <- sapply(candidates, function(x) .string_similarity(query, x))
  top_idx <- order(scores, decreasing = TRUE)[1:min(n, length(scores))]
  return(candidates[top_idx][scores[top_idx] > 0.3])
}
