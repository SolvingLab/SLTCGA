# ==============================================================================
# Helper Functions for Users
# ==============================================================================
# User-friendly functions to explore available data and variables
# ==============================================================================


#' List All Available Modalities
#'
#' @description
#' Display all omics layers available in SLTCGA with descriptions.
#'
#' @return Data frame with modality information (invisible)
#'
#' @examples
#' \dontrun{
#' list_modalities()
#' }
#'
#' @export
list_modalities <- function() {
  modalities <- data.frame(
    Modal = c(
      "RNAseq", "Mutation", "CNV", "Methylation", "miRNA",
      "Clinical", "Signature", "ImmuneCell"
    ),
    Description = c(
      "Gene expression (multiple normalizations)",
      "Somatic mutations (WildType/Mutation)",
      "Copy number variations (multiple algorithms)",
      "DNA methylation (450K array)",
      "MicroRNA expression (1881 miRNAs)",
      "Traditional clinical information (66 variables)",
      "Molecular signatures and scores (58 variables)",
      "Immune cell infiltration (99 cell types, 8 algorithms)"
    ),
    N_Variables = c(
      "~20,000", "~20,000", "~20,000", "~20,000", "1,881",
      "66", "58", "99"
    ),
    Data_Type = c(
      "Continuous", "Categorical", "Continuous", "Continuous",
      "Continuous", "Mixed", "Mixed", "Continuous"
    ),
    stringsAsFactors = FALSE
  )

  cat("\n╔════════════════════════════════════════════════════════════╗\n")
  cat("║          SLTCGA Available Modalities                      ║\n")
  cat("╚════════════════════════════════════════════════════════════╝\n\n")

  print(modalities, row.names = FALSE, right = FALSE)

  cat("\n")
  cat("Quick tips:\n")
  cat("  • Use list_variables(modal='Clinical') to see Clinical variables\n")
  cat("  • Use list_variables(modal='Signature') to see Signature variables\n")
  cat("  • Use list_immune_cells() to see all 99 immune cell types\n")
  cat("  • Use list_cancer_types() to see all 65 cancer types\n")
  cat("\n")

  invisible(modalities)
}


#' List Variables for a Modal
#'
#' @description
#' Display available variables for Clinical, Signature, or ImmuneCell modalities.
#'
#' @param modal Character. Modal type: "Clinical", "Signature", or "ImmuneCell"
#' @param pattern Character. Optional pattern to filter variables
#' @param show_groups Logical. Show variable groups (default: TRUE)
#'
#' @return Named vector of variables (invisible)
#'
#' @examples
#' \dontrun{
#' # List all clinical variables
#' list_variables(modal = "Clinical")
#'
#' # List signature variables
#' list_variables(modal = "Signature")
#'
#' # Search for TMB-related signatures
#' list_variables(modal = "Signature", pattern = "TMB")
#' }
#'
#' @export
list_variables <- function(modal, pattern = NULL, show_groups = TRUE) {
  modal <- match.arg(modal, c("Clinical", "Signature", "ImmuneCell"))

  if (modal == "Clinical") {
    vars <- CLINICAL_ALIASES
    groups <- CLINICAL_GROUPS
  } else if (modal == "Signature") {
    vars <- SIGNATURE_ALIASES
    groups <- SIGNATURE_GROUPS
  } else {
    # ImmuneCell: load from file
    base_path <- Sys.getenv("SL_BULK_DATA")
    if (base_path == "") {
      stop("SL_BULK_DATA not set. Use: Sys.setenv(SL_BULK_DATA='/path/to/data')", call. = FALSE)
    }

    immune_file <- file.path(base_path, "TCGA_DeconvCell_All_scores.qs")
    if (!file.exists(immune_file)) {
      stop("ImmuneCell file not found: ", immune_file, call. = FALSE)
    }

    immune_data <- qs::qread(immune_file, nthreads = 2)
    all_cols <- colnames(immune_data)
    immune_cols <- all_cols[all_cols != "cancer_type"]

    vars <- setNames(immune_cols, immune_cols)
    groups <- IMMUNE_CELL_CATEGORIES
  }

  # Filter by pattern
  if (!is.null(pattern)) {
    vars <- vars[grepl(pattern, names(vars), ignore.case = TRUE) |
      grepl(pattern, vars, ignore.case = TRUE)]
  }

  cat(sprintf("\n╔════════════════════════════════════════════════════════════╗\n"))
  cat(sprintf(
    "║  %s Modal: %d variables%s\n",
    modal,
    length(vars),
    strrep(" ", max(0, 42 - nchar(modal) - nchar(as.character(length(vars)))))
  ))
  cat(sprintf("╚════════════════════════════════════════════════════════════╝\n\n"))

  if (show_groups && length(groups) > 0) {
    cat("Available groups:\n")
    for (group_name in names(groups)) {
      cat(sprintf("  • %s: %d variables\n", group_name, length(groups[[group_name]])))
    }
    cat("\n")
  }

  cat("Variables (alias → full_name):\n")
  n_display <- min(30, length(vars))
  for (i in 1:n_display) {
    if (names(vars)[i] == vars[i]) {
      cat(sprintf("  %s\n", vars[i]))
    } else {
      cat(sprintf("  %-20s → %s\n", names(vars)[i], vars[i]))
    }
  }

  if (length(vars) > n_display) {
    cat(sprintf("\n  ... and %d more variables\n", length(vars) - n_display))
    cat("  Use pattern='keyword' to filter\n")
  }
  cat("\n")

  invisible(vars)
}


#' Search Variables Across Modalities
#'
#' @description
#' Search for variables matching a keyword across Clinical, Signature, and ImmuneCell.
#'
#' @param keyword Character. Keyword to search
#' @param modal Character or NULL. Optional modal to restrict search
#'
#' @return List of matched variables (invisible)
#'
#' @examples
#' \dontrun{
#' search_variables("TMB")
#' search_variables("T_cells", modal = "ImmuneCell")
#' }
#'
#' @export
search_variables <- function(keyword, modal = NULL) {
  keyword_lower <- tolower(keyword)
  results <- list()

  modals_to_search <- if (!is.null(modal)) modal else c("Clinical", "Signature", "ImmuneCell")

  for (m in modals_to_search) {
    if (m == "Clinical") {
      vars <- CLINICAL_ALIASES
    } else if (m == "Signature") {
      vars <- SIGNATURE_ALIASES
    } else {
      base_path <- Sys.getenv("SL_BULK_DATA")
      if (base_path == "") next

      immune_file <- file.path(base_path, "TCGA_DeconvCell_All_scores.qs")
      if (!file.exists(immune_file)) next

      immune_data <- qs::qread(immune_file, nthreads = 2)
      immune_cols <- colnames(immune_data)[colnames(immune_data) != "cancer_type"]
      vars <- setNames(immune_cols, immune_cols)
    }

    matches <- vars[grepl(keyword_lower, tolower(names(vars))) |
      grepl(keyword_lower, tolower(vars))]

    if (length(matches) > 0) {
      results[[m]] <- matches
    }
  }

  if (length(results) == 0) {
    cat(sprintf("No variables found matching '%s'\n", keyword))
    return(invisible(NULL))
  }

  cat(sprintf("\n╔════════════════════════════════════════════════════════════╗\n"))
  cat(sprintf(
    "║  Search results for '%s'%s\n", keyword,
    strrep(" ", max(0, 47 - nchar(keyword)))
  ))
  cat(sprintf("╚════════════════════════════════════════════════════════════╝\n\n"))

  for (m in names(results)) {
    cat(sprintf("%s (%d matches):\n", m, length(results[[m]])))
    n_show <- min(5, length(results[[m]]))
    for (i in 1:n_show) {
      if (names(results[[m]])[i] == results[[m]][i]) {
        cat(sprintf("  • %s\n", results[[m]][i]))
      } else {
        cat(sprintf("  • %-20s → %s\n", names(results[[m]])[i], results[[m]][i]))
      }
    }
    if (length(results[[m]]) > n_show) {
      cat(sprintf("  ... and %d more\n", length(results[[m]]) - n_show))
    }
    cat("\n")
  }

  invisible(results)
}


#' List Immune Cell Types
#'
#' @description
#' Display all available immune cell types with algorithm and category information.
#'
#' @param algorithm Character or NULL. Filter by algorithm:
#'   "cibersort", "xcell", "quantiseq", "mcpcounter", "timer", "epic", "ips", "estimate"
#' @param category Character or NULL. Filter by cell category:
#'   "B_cells", "T_cells_CD4", "T_cells_CD8", "Tregs", "NK_cells",
#'   "Macrophages", "DC", "Neutrophils", "Monocytes", "Microenvironment", "ESTIMATE"
#'
#' @return Data frame with cell information (invisible)
#'
#' @examples
#' \dontrun{
#' # List all immune cells
#' list_immune_cells()
#'
#' # List only CIBERSORT cells
#' list_immune_cells(algorithm = "cibersort")
#'
#' # List only CD8+ T cells
#' list_immune_cells(category = "T_cells_CD8")
#' }
#'
#' @export
list_immune_cells <- function(algorithm = NULL, category = NULL) {
  base_path <- Sys.getenv("SL_BULK_DATA")
  if (base_path == "") {
    stop("SL_BULK_DATA not set", call. = FALSE)
  }

  immune_file <- file.path(base_path, "TCGA_DeconvCell_All_scores.qs")
  if (!file.exists(immune_file)) {
    stop("ImmuneCell file not found", call. = FALSE)
  }

  immune_data <- qs::qread(immune_file, nthreads = 2)
  cell_names <- colnames(immune_data)[colnames(immune_data) != "cancer_type"]

  # Extract algorithm
  algos <- sapply(strsplit(cell_names, "_"), function(x) tail(x, 1))

  # Match category
  cats <- sapply(cell_names, function(cell) {
    for (cat_name in names(IMMUNE_CELL_CATEGORIES)) {
      if (cell %in% IMMUNE_CELL_CATEGORIES[[cat_name]]) {
        return(cat_name)
      }
    }
    return("Other")
  })

  df <- data.frame(
    Cell_Name = cell_names,
    Algorithm = algos,
    Category = cats,
    stringsAsFactors = FALSE
  )

  # Filter
  if (!is.null(algorithm)) {
    df <- df[df$Algorithm == algorithm, ]
  }
  if (!is.null(category)) {
    df <- df[df$Category == category, ]
  }

  cat(sprintf("\n╔════════════════════════════════════════════════════════════╗\n"))
  cat(sprintf(
    "║  Immune Cells: %d found%s\n", nrow(df),
    strrep(" ", max(0, 44 - nchar(as.character(nrow(df)))))
  ))
  cat(sprintf("╚════════════════════════════════════════════════════════════╝\n\n"))

  if (!is.null(algorithm)) {
    cat(sprintf("Algorithm: %s\n", algorithm))
  }
  if (!is.null(category)) {
    cat(sprintf("Category: %s\n", category))
  }
  if (!is.null(algorithm) || !is.null(category)) {
    cat("\n")
  }

  n_show <- min(30, nrow(df))
  print(head(df, n_show), row.names = FALSE, right = FALSE)

  if (nrow(df) > n_show) {
    cat(sprintf("\n... and %d more\n", nrow(df) - n_show))
  }
  cat("\n")

  invisible(df)
}


#' Search Immune Cell Types
#'
#' @description
#' Search for immune cell types matching a keyword.
#'
#' @param keyword Character. Keyword to search
#'
#' @return Data frame with matching cells (invisible)
#'
#' @examples
#' \dontrun{
#' search_immune_cells("B_cells")
#' search_immune_cells("macrophage")
#' }
#'
#' @export
search_immune_cells <- function(keyword) {
  df <- list_immune_cells()
  keyword_lower <- tolower(keyword)

  matches <- df[grepl(keyword_lower, tolower(df$Cell_Name)), ]

  if (nrow(matches) == 0) {
    cat(sprintf("No immune cells found matching '%s'\n", keyword))
    return(invisible(NULL))
  }

  cat(sprintf("\n╔════════════════════════════════════════════════════════════╗\n"))
  cat(sprintf(
    "║  Search '%s': %d matches%s\n", keyword, nrow(matches),
    strrep(" ", max(0, 43 - nchar(keyword) - nchar(as.character(nrow(matches)))))
  ))
  cat(sprintf("╚════════════════════════════════════════════════════════════╝\n\n"))

  print(matches, row.names = FALSE, right = FALSE)
  cat("\n")

  invisible(matches)
}


#' List Cancer Types
#'
#' @description
#' Display all available cancer types (33 main + 32 subtypes = 65 total).
#'
#' @param show_subtypes Logical. Show subtypes (default: TRUE)
#'
#' @return Data frame with cancer information (invisible)
#'
#' @examples
#' \dontrun{
#' list_cancer_types()
#' list_cancer_types(show_subtypes = FALSE)
#' }
#'
#' @export
list_cancer_types <- function(show_subtypes = TRUE) {
  cat("\n╔════════════════════════════════════════════════════════════╗\n")
  cat(sprintf(
    "║  TCGA Cancer Types: %d Main + %d Subtypes + %d Combined = %d Total ║\n",
    length(TCGA_MAIN_CANCERS),
    length(TCGA_SUBTYPE_MAP),
    length(TCGA_COMBINED_MAP),
    length(TCGA_MAIN_CANCERS) + length(TCGA_SUBTYPE_MAP) + length(TCGA_COMBINED_MAP)
  ))
  cat("╚════════════════════════════════════════════════════════════╝\n\n")

  cat("Main cancer types (33):\n")
  for (i in seq_along(TCGA_MAIN_CANCERS)) {
    cat(sprintf("  %s", TCGA_MAIN_CANCERS[i]))
    if (i %% 5 == 0) cat("\n")
  }
  if (length(TCGA_MAIN_CANCERS) %% 5 != 0) cat("\n")

  if (show_subtypes) {
    cat(sprintf("\nMolecular subtypes (%d):\n", length(TCGA_SUBTYPE_MAP)))

    subtype_df <- data.frame(
      Subtype = names(TCGA_SUBTYPE_MAP),
      Parent = sapply(TCGA_SUBTYPE_MAP, function(x) x$parent),
      stringsAsFactors = FALSE
    )

    # Group by parent
    parents <- unique(subtype_df$Parent)
    for (parent in parents) {
      subtypes <- subtype_df$Subtype[subtype_df$Parent == parent]
      cat(sprintf("  %s: %s\n", parent, paste(subtypes, collapse = ", ")))
    }
  }

  if (length(TCGA_COMBINED_MAP) > 0) {
    cat("\nCombined cancer groups:\n")
    for (name in names(TCGA_COMBINED_MAP)) {
      cat(sprintf("  %s: %s\n", name, paste(TCGA_COMBINED_MAP[[name]], collapse = " + ")))
    }
  }

  cat("\nNote: Cancer type input is case-insensitive\n")
  cat("      Examples: 'BRCA', 'brca', 'Brca' all work\n\n")

  # Create full list
  all_cancers <- data.frame(
    Cancer_Type = c(TCGA_MAIN_CANCERS, names(TCGA_SUBTYPE_MAP), names(TCGA_COMBINED_MAP)),
    Type = c(
      rep("Main", length(TCGA_MAIN_CANCERS)),
      rep("Subtype", length(TCGA_SUBTYPE_MAP)),
      rep("Combined", length(TCGA_COMBINED_MAP))
    ),
    stringsAsFactors = FALSE
  )

  invisible(all_cancers)
}


#' Get Variable Groups
#'
#' @description
#' Retrieve predefined variable groups for quick access.
#'
#' @param modal Character. Modal type: "Clinical" or "Signature"
#' @param group Character or NULL. Specific group name, or NULL to list all
#'
#' @return Named list of variable groups or character vector for specific group
#'
#' @examples
#' \dontrun{
#' # List all clinical groups
#' get_variable_groups(modal = "Clinical")
#'
#' # Get basic clinical variables
#' basic_clin <- get_variable_groups(modal = "Clinical", group = "basic")
#'
#' # Get immune signatures
#' immune_sig <- get_variable_groups(modal = "Signature", group = "immune")
#' }
#'
#' @export
get_variable_groups <- function(modal, group = NULL) {
  modal <- match.arg(modal, c("Clinical", "Signature"))

  groups <- if (modal == "Clinical") {
    CLINICAL_GROUPS
  } else {
    SIGNATURE_GROUPS
  }

  if (is.null(group)) {
    cat(sprintf("\nAvailable %s groups:\n", modal))
    for (group_name in names(groups)) {
      cat(sprintf("  • %s: %d variables\n", group_name, length(groups[[group_name]])))
    }
    cat("\nUse get_variable_groups(modal='", modal, "', group='name') to retrieve specific group\n\n", sep = "")
    return(invisible(groups))
  } else {
    if (!group %in% names(groups)) {
      stop("Group '", group, "' not found. Use get_variable_groups(modal='", modal, "') to see all groups.", call. = FALSE)
    }
    return(groups[[group]])
  }
}


#' Quick Summary of SLTCGA Package
#'
#' @description
#' Display a quick reference guide for SLTCGA usage.
#'
#' @export
sltcga_guide <- function() {
  cat("\n╔════════════════════════════════════════════════════════════╗\n")
  cat("║              SLTCGA Quick Reference Guide                 ║\n")
  cat("╚════════════════════════════════════════════════════════════╝\n\n")

  cat("Main Functions:\n")
  cat("  • tcga_correlation()  - Correlation/association analysis (Scenarios 1-7)\n")
  cat("  • tcga_enrichment()   - Pathway enrichment analysis (Scenarios 8-15)\n")
  cat("  • tcga_survival()     - Survival analysis (Scenarios 16-17)\n\n")

  cat("Helper Functions:\n")
  cat("  • list_modalities()              - Show all 8 omics layers\n")
  cat("  • list_variables(modal='...')    - Show variables for a modal\n")
  cat("  • search_variables('keyword')    - Search across all modalities\n")
  cat("  • list_immune_cells()            - Show all 99 immune cell types\n")
  cat("  • list_cancer_types()            - Show all 65 cancer types\n\n")

  cat("Data Setup:\n")
  cat("  Sys.setenv(SL_BULK_DATA = '/path/to/bulk_data')\n\n")

  cat("Example Usage:\n")
  cat("  # Correlation analysis\n")
  cat("  result <- tcga_correlation(\n")
  cat("    var1 = 'TP53', var1_modal = 'RNAseq', var1_cancers = 'BRCA',\n")
  cat("    var2 = 'TMB', var2_modal = 'Signature', var2_cancers = 'BRCA'\n")
  cat("  )\n\n")

  cat("  # Enrichment analysis\n")
  cat("  result <- tcga_enrichment(\n")
  cat("    var1 = 'TP53', var1_modal = 'Mutation', var1_cancers = 'LUAD',\n")
  cat("    analysis_type = 'enrichment', top_n = 20\n")
  cat("  )\n\n")

  cat("  # Survival analysis\n")
  cat("  result <- tcga_survival(\n")
  cat("    var1 = 'TP53', var1_modal = 'RNAseq', var1_cancers = 'BRCA',\n")
  cat("    surv_type = 'OS', cutoff_type = 'optimal'\n")
  cat("  )\n\n")

  cat("For detailed documentation, see:\n")
  cat("  ?tcga_correlation\n")
  cat("  ?tcga_enrichment\n")
  cat("  ?tcga_survival\n\n")
}
