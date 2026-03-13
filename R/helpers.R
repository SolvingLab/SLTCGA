# ==============================================================================
# Helper Functions for Users
# ==============================================================================
# User-friendly functions to explore available data and variables
# ==============================================================================


#' List All Available Data Modalities in TCGA Database
#'
#' @description
#' Displays comprehensive overview of all 8 data modalities available in SLTCGA:
#' 5 omics layers (RNAseq, Mutation, CNV, Methylation, miRNA), clinical data
#' (Clinical), molecular signatures (Signature), and immune infiltration scores
#' (ImmuneCell). Shows data types, variable counts, and descriptions. Essential
#' for understanding data structure before using \code{tcga_correlation()},
#' \code{tcga_enrichment()}, or \code{tcga_survival()}. Returns invisible data
#' frame for programmatic access.
#'
#' @return Data frame with 4 columns (invisible):
#'   \describe{
#'     \item{Modal}{Modality name (e.g., "RNAseq", "Mutation")}
#'     \item{Description}{Brief description of data type}
#'     \item{N_Variables}{Approximate number of variables}
#'     \item{Data_Type}{Variable type ("Continuous", "Categorical", "Mixed")}
#'   }
#'
#' @details
#' **8 Data Modalities** (3 categories):
#'
#' **Multi-Omics Layers** (5 modalities, ~20,000 genes each):
#' \itemize{
#'   \item **RNAseq**: Gene expression (continuous, multiple normalizations)
#'   \item **Mutation**: Somatic mutations (categorical: WildType/Mutation)
#'   \item **CNV**: Copy number variations (continuous, multiple algorithms)
#'   \item **Methylation**: DNA methylation (continuous, 450K array)
#'   \item **miRNA**: MicroRNA expression (1,881 miRNAs, continuous)
#' }
#'
#' **Clinical Data** (1 modality):
#' \itemize{
#'   \item **Clinical**: Traditional clinical variables (66 variables, mixed types)
#'     - Demographics, treatment, outcome, histology, molecular features
#' }
#'
#' **Derived Data** (2 modalities):
#' \itemize{
#'   \item **Signature**: Molecular signatures and scores (58 variables, mixed types)
#'     - Immune, metabolic, pathway, stemness, clinical scores (TMB, MSI, etc.)
#'   \item **ImmuneCell**: Immune cell infiltration (99 cell types, continuous)
#'     - Deconvolution results from 8 algorithms applied to RNAseq data
#' }
#'
#' **Next Steps**:
#' \itemize{
#'   \item **Omics layers**: Use gene symbols directly (e.g., "TP53", "BRCA1")
#'   \item **Clinical**: Use \code{\link{list_variables}(modal = "Clinical")} to see 66 variables
#'   \item **Signature**: Use \code{\link{list_variables}(modal = "Signature")} to see 58 scores
#'   \item **ImmuneCell**: Use \code{\link{list_immune_cells}()} to see 99 cell types
#'   \item **Cancer types**: Use \code{\link{list_cancer_types}()} to see 65 cancer types
#' }
#'
#' @section User Queries:
#' **Data Exploration**:
#' \itemize{
#'   \item What data types are available in SLTCGA?
#'   \item What omics layers can I analyze?
#'   \item How many genes are covered in TCGA?
#'   \item What is the difference between RNAseq and Mutation modalities?
#'   \item Which modalities are continuous vs categorical?
#'   \item Can I analyze immune cell infiltration?
#'   \item What molecular signatures are available?
#'   \item How many clinical variables are there?
#'   \item What is the difference between omics data and clinical data?
#'   \item Are Signature and ImmuneCell original omics or derived data?
#' }
#'
#' **Method Selection**:
#' \itemize{
#'   \item Which modality should I use for gene expression analysis?
#'   \item How do I analyze DNA methylation?
#'   \item What modality contains tumor mutation burden (TMB)?
#'   \item Which modality has immune cell data?
#'   \item Can I analyze miRNA expression?
#'   \item What is the ImmuneCell modality?
#' }
#'
#' @references
#' **TCGA Database**:
#'
#' The Cancer Genome Atlas Research Network (2013). The Cancer Genome Atlas
#' Pan-Cancer analysis project. Nature Genetics, 45(10):1113-1120.
#' \doi{10.1038/ng.2764}
#'
#' Database portal: \url{https://www.cancer.gov/tcga}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{list_variables}} - View variables for specific modality
#'   \item \code{\link{search_variables}} - Search variables by keyword
#'   \item \code{\link{list_immune_cells}} - View all immune cell types
#'   \item \code{\link{list_cancer_types}} - View all cancer types
#'   \item \code{\link{tcga_correlation}} - Correlation analysis using modalities
#'   \item \code{\link{sltcga_guide}} - Quick reference guide
#' }
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: View all available modalities
#' # ===========================================================================
#' # Research Question: What data types can I analyze in TCGA?
#'
#' list_modalities()
#'
#' # Output shows 8 modalities in 3 categories:
#' # - Multi-Omics: RNAseq, Mutation, CNV, Methylation, miRNA
#' # - Clinical: Clinical variables
#' # - Derived: Signature scores, ImmuneCell infiltration
#'
#' # ===========================================================================
#' # Example 2: Programmatic access to modality information
#' # ===========================================================================
#'
#' modals <- list_modalities()
#'
#' # Check data types
#' continuous_modals <- modals$Modal[modals$Data_Type == "Continuous"]
#' # Returns: "RNAseq", "CNV", "Methylation", "miRNA", "ImmuneCell"
#'
#' categorical_modals <- modals$Modal[modals$Data_Type == "Categorical"]
#' # Returns: "Mutation"
#'
#' # ===========================================================================
#' # Next Steps
#' # ===========================================================================
#' # After viewing modalities:
#' # 1. Use list_variables(modal = "Clinical") to explore clinical variables
#' # 2. Use list_immune_cells() to see immune cell types
#' # 3. Use list_cancer_types() to see available cancer types
#' # 4. Start analysis with tcga_correlation(), tcga_enrichment(), or tcga_survival()
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
  cat("║          SLTCGA Available Data Modalities                 ║\n")
  cat("╚════════════════════════════════════════════════════════════╝\n\n")

  print(modalities, row.names = FALSE, right = FALSE)

  cat("\n")
  cat("Data categories:\n")
  cat("  • Multi-Omics (5): RNAseq, Mutation, CNV, Methylation, miRNA\n")
  cat("  • Clinical (1): Clinical variables\n")
  cat("  • Derived (2): Signature scores, ImmuneCell infiltration\n\n")
  cat("Quick tips:\n")
  cat("  • For omics: Use gene symbols directly (e.g., 'TP53', 'BRCA1')\n")
  cat("  • For Clinical: Use list_variables(modal='Clinical')\n")
  cat("  • For Signature: Use list_variables(modal='Signature')\n")
  cat("  • For ImmuneCell: Use list_immune_cells()\n")
  cat("  • For cancer types: Use list_cancer_types()\n")
  cat("\n")

  invisible(modalities)
}


#' List Variables for Clinical, Signature, or ImmuneCell Modalities
#'
#' @description
#' Displays available variables for Clinical (66 variables), Signature (58 variables),
#' or ImmuneCell (99 cell types) modalities with optional pattern filtering and
#' grouping. For RNAseq/Mutation/CNV/Methylation/miRNA, use gene symbols directly
#' (~20,000 genes each). Essential for discovering variable names before analysis.
#' Returns invisible named vector for programmatic access.
#'
#' @param modal Character. Modality type to display (required).
#'   Options: "Clinical", "Signature", "ImmuneCell".
#'   Note: RNAseq, Mutation, CNV, Methylation, miRNA use standard gene symbols.
#' @param pattern Character or NULL. Optional regex pattern to filter variables (default: NULL).
#'   Case-insensitive matching on both variable names and aliases.
#'   Examples: "TMB", "T_cells", "Age", "Stage".
#' @param show_groups Logical. Display variable groups/categories (default: TRUE).
#'   For Clinical: basic, treatment, outcome, histology, molecular.
#'   For Signature: immune, metabolic, pathway, stemness, clinical_scores.
#'   For ImmuneCell: algorithm-based grouping.
#'
#' @return Named character vector (invisible):
#'   Names = user-friendly aliases, Values = full variable names.
#'   Access programmatically: \code{vars <- list_variables(modal = "Clinical")}.
#'
#' @details
#' **Clinical Variables (66)**:
#' \itemize{
#'   \item **Basic** (11): Age, Gender, Race, Ethnicity, Height, Weight, BMI, etc.
#'   \item **Treatment** (8): Radiation, Chemotherapy, Immunotherapy, etc.
#'   \item **Outcome** (12): Vital_status, Recurrence, Progression, etc.
#'   \item **Histology** (18): Histology_type, Grade, Stage, T/N/M, etc.
#'   \item **Molecular** (17): MSI_status, TMB_category, HRD_score, etc.
#' }
#'
#' **Signature Variables (58)**:
#' \itemize{
#'   \item **Immune** (15): CYT, IFNG, TIS, Tcell_inflamed, etc.
#'   \item **Metabolic** (8): Glycolysis, OXPHOS, Fatty_acid_metabolism, etc.
#'   \item **Pathway** (12): EMT, Angiogenesis, Apoptosis, DNA_repair, etc.
#'   \item **Stemness** (6): OCLR, DMPss, ENHss, EREG_METHss, etc.
#'   \item **Clinical Scores** (17): TMB, MSI_score, HRD_score, Aneuploidy, etc.
#' }
#'
#' **ImmuneCell Variables (99)**:
#' \itemize{
#'   \item 8 algorithms: CIBERSORT, xCell, quanTIseq, MCPcounter, TIMER, EPIC, IPS, ESTIMATE
#'   \item 11 categories: B_cells, T_cells_CD4, T_cells_CD8, Tregs, NK_cells,
#'         Macrophages, DC, Neutrophils, Monocytes, Microenvironment, ESTIMATE
#'   \item Use \code{\link{list_immune_cells}()} for detailed cell information
#' }
#'
#' **Filtering with Pattern**:
#' \itemize{
#'   \item \code{pattern = "TMB"} → Finds TMB, TMB_NonSynonymous, TMB_Nonsilent
#'   \item \code{pattern = "T_cells"} → Finds all T cell types
#'   \item \code{pattern = "Stage"} → Finds Stage, Stage_T, Stage_N, Stage_M
#' }
#'
#' @section User Queries:
#' **Clinical Variables**:
#' \itemize{
#'   \item What clinical variables are available?
#'   \item How do I access patient age, gender, or race?
#'   \item What treatment information is available?
#'   \item Can I analyze tumor stage or grade?
#'   \item Is MSI status available?
#'   \item What survival endpoints are there?
#'   \item How do I find histology information?
#' }
#'
#' **Signature Variables**:
#' \itemize{
#'   \item What molecular signatures can I analyze?
#'   \item How do I access tumor mutation burden (TMB)?
#'   \item Is there an EMT score?
#'   \item What immune signatures are available?
#'   \item Can I analyze hypoxia or angiogenesis?
#'   \item Is there a stemness score?
#'   \item What metabolic signatures exist?
#'   \item How do I find DNA repair signatures?
#' }
#'
#' **ImmuneCell Variables**:
#' \itemize{
#'   \item What immune cell types are available?
#'   \item How many deconvolution algorithms are included?
#'   \item Can I analyze CD8+ T cell infiltration?
#'   \item Is macrophage infiltration data available?
#'   \item What is the difference between CIBERSORT and xCell?
#'   \item How do I find regulatory T cells (Tregs)?
#' }
#'
#' **Pattern Filtering**:
#' \itemize{
#'   \item How do I search for TMB-related variables?
#'   \item Can I filter signatures by keyword?
#'   \item How do I find all T cell types?
#'   \item Can I search for stage-related clinical variables?
#' }
#'
#' @references
#' **TCGA Database**:
#'
#' The Cancer Genome Atlas Research Network (2013). The Cancer Genome Atlas
#' Pan-Cancer analysis project. Nature Genetics, 45(10):1113-1120.
#' \doi{10.1038/ng.2764}
#'
#' Database portal: \url{https://www.cancer.gov/tcga}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{list_modalities}} - View all omics modalities
#'   \item \code{\link{search_variables}} - Search across multiple modalities
#'   \item \code{\link{list_immune_cells}} - Detailed immune cell information
#'   \item \code{\link{get_variable_groups}} - Retrieve predefined variable groups
#'   \item \code{\link{tcga_correlation}} - Use variables in correlation analysis
#' }
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: List all clinical variables
#' # ===========================================================================
#' # Research Question: What clinical data can I analyze?
#'
#' list_variables(modal = "Clinical")
#'
#' # Output shows 66 clinical variables grouped by category
#' # Use these in tcga_correlation() or tcga_survival()
#'
#' # ===========================================================================
#' # Example 2: List molecular signatures
#' # ===========================================================================
#' # Research Question: What molecular signatures are available?
#'
#' list_variables(modal = "Signature")
#'
#' # Output shows 58 signatures: immune, metabolic, pathway, stemness, clinical scores
#'
#' # ===========================================================================
#' # Example 3: Search for TMB-related signatures
#' # ===========================================================================
#' # Research Question: How do I find tumor mutation burden variables?
#'
#' list_variables(modal = "Signature", pattern = "TMB")
#'
#' # Returns: TMB, TMB_NonSynonymous, TMB_Nonsilent
#'
#' # ===========================================================================
#' # Example 4: Programmatic access to variables
#' # ===========================================================================
#'
#' # Get all clinical variables
#' clin_vars <- list_variables(modal = "Clinical")
#'
#' # Get signature names
#' sig_names <- names(list_variables(modal = "Signature"))
#'
#' # Filter immune signatures
#' immune_sigs <- list_variables(modal = "Signature", pattern = "immune|CYT|IFNG")
#'
#' # ===========================================================================
#' # Example 5: Explore immune cells (use dedicated function)
#' # ===========================================================================
#'
#' # For ImmuneCell, use list_immune_cells() for better display
#' list_immune_cells(algorithm = "cibersort")
#' list_immune_cells(category = "T_cells_CD8")
#'
#' # ===========================================================================
#' # Next Steps
#' # ===========================================================================
#' # After finding variable names:
#' # 1. Use tcga_correlation(var1 = "Age", var1_modal = "Clinical", ...)
#' # 2. Use tcga_correlation(var1 = "TMB", var1_modal = "Signature", ...)
#' # 3. Use tcga_survival(var1 = "Stage", var1_modal = "Clinical", ...)
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


#' Search Variables Across Clinical, Signature, and ImmuneCell Modalities
#'
#' @description
#' Performs case-insensitive keyword search across Clinical (66 variables),
#' Signature (58 variables), and ImmuneCell (99 cell types) modalities, with
#' optional restriction to specific modality. Searches both variable names and
#' aliases. Returns grouped results showing matches per modality. Essential for
#' discovering relevant variables when variable name is uncertain.
#'
#' @param keyword Character. Search keyword or pattern (required).
#'   Case-insensitive regex matching on variable names and aliases.
#'   Examples: "TMB", "T_cells", "Stage", "immune", "stemness".
#' @param modal Character or NULL. Optional modality restriction (default: NULL).
#'   Options: "Clinical", "Signature", "ImmuneCell", or NULL for all.
#'   NULL searches across all three modalities.
#'
#' @return Named list of matched variables (invisible):
#'   \describe{
#'     \item{Clinical}{Named vector of Clinical matches (if any)}
#'     \item{Signature}{Named vector of Signature matches (if any)}
#'     \item{ImmuneCell}{Named vector of ImmuneCell matches (if any)}
#'   }
#'   Returns NULL if no matches found.
#'
#' @details
#' **Search Scope**:
#' \itemize{
#'   \item **Clinical** (66 variables): Age, Gender, Stage, Grade, MSI_status, etc.
#'   \item **Signature** (58 variables): TMB, EMT, Hypoxia, CYT, Stemness, etc.
#'   \item **ImmuneCell** (99 cell types): CD8_T_cells, Macrophages, B_cells, etc.
#'   \item **Not searched**: RNAseq, Mutation, CNV, Methylation, miRNA (use gene symbols directly)
#' }
#'
#' **Search Behavior**:
#' \itemize{
#'   \item Case-insensitive: "tmb", "TMB", "Tmb" all match
#'   \item Partial matching: "T_cells" matches "CD8_T_cells", "CD4_T_cells", etc.
#'   \item Regex support: "Stage|Grade" matches both Stage and Grade variables
#'   \item Alias matching: Searches both user-friendly names and full variable names
#' }
#'
#' **Common Search Terms**:
#' \itemize{
#'   \item Immune: "T_cells", "Macrophage", "B_cells", "NK", "Treg", "CD8"
#'   \item Clinical: "Age", "Gender", "Stage", "Grade", "MSI", "Race"
#'   \item Signatures: "TMB", "EMT", "Hypoxia", "Angiogenesis", "Stemness", "CYT"
#'   \item Metabolic: "Glycolysis", "OXPHOS", "Fatty_acid"
#' }
#'
#' @section User Queries:
#' **General Search**:
#' \itemize{
#'   \item How do I find variables related to tumor mutation burden?
#'   \item Can I search for immune-related variables?
#'   \item What variables contain "Stage" in their name?
#'   \item How do I find all T cell types?
#'   \item Can I search for EMT-related signatures?
#'   \item Is there a variable for hypoxia?
#' }
#'
#' **Immune Cell Search**:
#' \itemize{
#'   \item How do I find CD8+ T cells?
#'   \item What macrophage variables are available?
#'   \item Can I search for B cell infiltration data?
#'   \item How do I find regulatory T cells (Tregs)?
#'   \item What NK cell types exist?
#' }
#'
#' **Clinical Search**:
#' \itemize{
#'   \item How do I find patient age variable?
#'   \item What is the variable name for tumor stage?
#'   \item Can I search for MSI status?
#'   \item How do I find treatment-related variables?
#'   \item What survival endpoints are available?
#' }
#'
#' **Signature Search**:
#' \itemize{
#'   \item How do I find immune signature scores?
#'   \item What stemness signatures exist?
#'   \item Can I search for metabolic pathway scores?
#'   \item Is there a DNA repair signature?
#'   \item How do I find angiogenesis scores?
#' }
#'
#' **Workflow Questions**:
#' \itemize{
#'   \item I want to analyze immune infiltration but don't know exact variable names
#'   \item How do I discover what variables are available for my research topic?
#'   \item Can I search across multiple modalities at once?
#'   \item How do I find alternative names for the same variable?
#' }
#'
#' @references
#' **TCGA Database**:
#'
#' The Cancer Genome Atlas Research Network (2013). The Cancer Genome Atlas
#' Pan-Cancer analysis project. Nature Genetics, 45(10):1113-1120.
#' \doi{10.1038/ng.2764}
#'
#' Database portal: \url{https://www.cancer.gov/tcga}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{list_variables}} - List all variables for a modality
#'   \item \code{\link{list_modalities}} - View all omics modalities
#'   \item \code{\link{list_immune_cells}} - Detailed immune cell search
#'   \item \code{\link{tcga_correlation}} - Use found variables in analysis
#' }
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Search for tumor mutation burden (TMB)
#' # ===========================================================================
#' # Research Question: What TMB-related variables exist?
#'
#' search_variables("TMB")
#'
#' # Returns matches in Signature:
#' # - TMB
#' # - TMB_NonSynonymous
#' # - TMB_Nonsilent
#'
#' # ===========================================================================
#' # Example 2: Search for T cell types across modalities
#' # ===========================================================================
#' # Research Question: What T cell infiltration data is available?
#'
#' search_variables("T_cells")
#'
#' # Returns matches in ImmuneCell:
#' # - CD8_T_cells_cibersort
#' # - CD4_T_cells_cibersort
#' # - Tregs_cibersort
#' # - T_cells_CD8_xcell
#' # ... (shows first 5 per modality)
#'
#' # ===========================================================================
#' # Example 3: Search immune cells only
#' # ===========================================================================
#'
#' search_variables("Macrophage", modal = "ImmuneCell")
#'
#' # Returns only ImmuneCell matches:
#' # - Macrophages_M0_cibersort
#' # - Macrophages_M1_cibersort
#' # - Macrophages_M2_cibersort
#' # - Macrophages_xcell
#'
#' # ===========================================================================
#' # Example 4: Search for stage-related clinical variables
#' # ===========================================================================
#'
#' search_variables("Stage", modal = "Clinical")
#'
#' # Returns Clinical matches:
#' # - Stage
#' # - Stage_T
#' # - Stage_N
#' # - Stage_M
#' # - Pathologic_stage
#'
#' # ===========================================================================
#' # Example 5: Search immune signatures
#' # ===========================================================================
#'
#' search_variables("immune|CYT|IFNG", modal = "Signature")
#'
#' # Returns Signature matches:
#' # - CYT (Cytolytic activity)
#' # - IFNG (Interferon gamma)
#' # - TIS (T cell inflamed signature)
#'
#' # ===========================================================================
#' # Example 6: Programmatic use
#' # ===========================================================================
#'
#' # Get search results
#' results <- search_variables("EMT")
#'
#' # Extract Signature matches
#' emt_vars <- results$Signature
#'
#' # Use in analysis
#' # tcga_correlation(
#' #   var1 = names(emt_vars)[1], var1_modal = "Signature", ...
#' # )
#'
#' # ===========================================================================
#' # Next Steps
#' # ===========================================================================
#' # After finding variable names:
#' # 1. Use list_variables() to see full variable lists
#' # 2. Use found variables in tcga_correlation(), tcga_enrichment(), tcga_survival()
#' # 3. For genes, use gene symbols directly with RNAseq/Mutation/CNV modality
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


#' List All Immune Cell Types from 8 Deconvolution Algorithms
#'
#' @description
#' Displays comprehensive catalog of all 99 immune cell types available in SLTCGA,
#' derived from 8 deconvolution algorithms (CIBERSORT, xCell, quanTIseq, MCPcounter,
#' TIMER, EPIC, IPS, ESTIMATE) with category classification (11 categories). Supports
#' filtering by algorithm or cell category. Essential for immune infiltration analysis.
#' Returns invisible data frame for programmatic access.
#'
#' @param algorithm Character or NULL. Filter by deconvolution algorithm (default: NULL).
#'   Options: "cibersort", "xcell", "quantiseq", "mcpcounter", "timer", "epic", "ips", "estimate".
#'   NULL shows all 99 cell types across all 8 algorithms.
#'   Each algorithm provides different cell type resolution and granularity.
#' @param category Character or NULL. Filter by functional cell category (default: NULL).
#'   Options: "B_cells", "T_cells_CD4", "T_cells_CD8", "Tregs", "NK_cells",
#'            "Macrophages", "DC", "Neutrophils", "Monocytes", "Microenvironment", "ESTIMATE".
#'   NULL shows all categories.
#'   Use to focus on specific immune cell lineages.
#'
#' @return Data frame with 3 columns (invisible):
#'   \describe{
#'     \item{Cell_Name}{Full cell type name (e.g., "CD8_T_cells_cibersort")}
#'     \item{Algorithm}{Deconvolution algorithm (e.g., "cibersort")}
#'     \item{Category}{Functional category (e.g., "T_cells_CD8")}
#'   }
#'
#' @details
#' **8 Deconvolution Algorithms**:
#' \itemize{
#'   \item **CIBERSORT**: 22 cell types, leukocyte-focused, high resolution
#'   \item **xCell**: 64 cell types, tissue-based, includes stromal cells
#'   \item **quanTIseq**: 10 cell types, RNA-seq optimized, absolute quantification
#'   \item **MCPcounter**: 10 cell types, microenvironment profiling
#'   \item **TIMER**: 6 cell types, tumor-infiltrating cells
#'   \item **EPIC**: 8 cell types, bulk RNA-seq deconvolution
#'   \item **IPS**: Immunophenoscore, 4 components
#'   \item **ESTIMATE**: Immune/Stromal/Tumor purity scores
#' }
#'
#' **11 Functional Categories**:
#' \itemize{
#'   \item **B_cells**: Naive B, Memory B, Plasma cells
#'   \item **T_cells_CD4**: Naive CD4+, Memory CD4+, Helper T cells (Th1, Th2, Th17)
#'   \item **T_cells_CD8**: CD8+ T cells, Cytotoxic T cells
#'   \item **Tregs**: Regulatory T cells (CD4+ CD25+ FOXP3+)
#'   \item **NK_cells**: Natural Killer cells (activated/resting)
#'   \item **Macrophages**: M0, M1, M2 macrophages
#'   \item **DC**: Dendritic cells (activated/resting, myeloid/plasmacytoid)
#'   \item **Neutrophils**: Neutrophils
#'   \item **Monocytes**: Monocytes (classical/non-classical)
#'   \item **Microenvironment**: Fibroblasts, Endothelial, MSC, Adipocytes
#'   \item **ESTIMATE**: ImmuneScore, StromalScore, TumorPurity
#' }
#'
#' **Algorithm Comparison**:
#' \itemize{
#'   \item **Resolution**: xCell (64 types) > CIBERSORT (22) > quanTIseq (10)
#'   \item **Specificity**: CIBERSORT best for lymphocytes, xCell for stromal cells
#'   \item **Quantification**: quanTIseq provides absolute fractions, others are relative
#'   \item **Speed**: MCPcounter fastest, CIBERSORT slowest
#' }
#'
#' @section User Queries:
#' **General Exploration**:
#' \itemize{
#'   \item What immune cell types can I analyze?
#'   \item How many immune cells are available in TCGA?
#'   \item What deconvolution algorithms are included?
#'   \item Which algorithm should I use for my analysis?
#'   \item What is the difference between CIBERSORT and xCell?
#' }
#'
#' **Specific Cell Types**:
#' \itemize{
#'   \item How do I find CD8+ T cells?
#'   \item What macrophage subtypes are available (M0/M1/M2)?
#'   \item Can I analyze regulatory T cells (Tregs)?
#'   \item Is B cell infiltration data available?
#'   \item What NK cell types exist?
#'   \item Can I analyze dendritic cells?
#'   \item Is neutrophil infiltration available?
#' }
#'
#' **Algorithm Selection**:
#' \itemize{
#'   \item Which CIBERSORT cell types are available?
#'   \item What does xCell provide?
#'   \item How many cell types does quanTIseq have?
#'   \item What is ESTIMATE algorithm?
#'   \item Can I compare results across different algorithms?
#' }
#'
#' **Category Filtering**:
#' \itemize{
#'   \item How do I list all T cell types?
#'   \item What macrophage types exist?
#'   \item Can I see all CD8+ T cell variants?
#'   \item How do I find stromal/microenvironment cells?
#'   \item What B cell subtypes are available?
#' }
#'
#' **Research Applications**:
#' \itemize{
#'   \item Can I correlate immune infiltration with gene expression?
#'   \item How do I analyze immune-clinical associations?
#'   \item Can I compare immune profiles across cancer types?
#'   \item Is tumor purity information available?
#'   \item How do I study tumor microenvironment?
#' }
#'
#' @references
#' **TCGA Database**:
#'
#' The Cancer Genome Atlas Research Network (2013). The Cancer Genome Atlas
#' Pan-Cancer analysis project. Nature Genetics, 45(10):1113-1120.
#' \doi{10.1038/ng.2764}
#'
#' Database portal: \url{https://www.cancer.gov/tcga}
#'
#' **Deconvolution Methods**:
#'
#' Newman AM, et al. (2015). Robust enumeration of cell subsets from tissue
#' expression profiles. Nature Methods, 12(5):453-457. (CIBERSORT)
#'
#' Aran D, et al. (2017). xCell: digitally portraying the tissue cellular
#' heterogeneity landscape. Genome Biology, 18:220. (xCell)
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{search_immune_cells}} - Search immune cells by keyword
#'   \item \code{\link{list_variables}} - View ImmuneCell as modality
#'   \item \code{\link{tcga_correlation}} - Correlate immune cells with genes/clinical
#'   \item \code{\link{tcga_survival}} - Immune cell prognostic analysis
#' }
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: List all 99 immune cell types
#' # ===========================================================================
#' # Research Question: What immune cell data is available?
#'
#' list_immune_cells()
#'
#' # Output shows all 99 cell types with algorithm and category
#'
#' # ===========================================================================
#' # Example 2: List CIBERSORT cell types (22 types)
#' # ===========================================================================
#' # Research Question: What cell types does CIBERSORT provide?
#'
#' list_immune_cells(algorithm = "cibersort")
#'
#' # Returns 22 CIBERSORT cell types:
#' # - B_cells_naive, B_cells_memory, Plasma_cells
#' # - CD8_T_cells, CD4_T_cells_naive, CD4_T_cells_memory_resting, etc.
#' # - Macrophages_M0, Macrophages_M1, Macrophages_M2
#' # - NK_cells_activated, NK_cells_resting
#' # - Dendritic_cells_activated, Dendritic_cells_resting
#'
#' # ===========================================================================
#' # Example 3: List all CD8+ T cell variants
#' # ===========================================================================
#' # Research Question: What CD8+ T cell types are available?
#'
#' list_immune_cells(category = "T_cells_CD8")
#'
#' # Returns CD8+ T cells from multiple algorithms:
#' # - CD8_T_cells_cibersort
#' # - T_cells_CD8_xcell
#' # - T_cells_CD8_quantiseq
#' # - CD8_T_cells_mcpcounter
#'
#' # ===========================================================================
#' # Example 4: Compare macrophage types across algorithms
#' # ===========================================================================
#'
#' list_immune_cells(category = "Macrophages")
#'
#' # Shows macrophage variants:
#' # - CIBERSORT: M0, M1, M2 (3 subtypes)
#' # - xCell: Macrophages, Macrophages_M1, Macrophages_M2
#' # - quanTIseq: Macrophages_M1, Macrophages_M2
#' # - MCPcounter: Monocytic_lineage
#'
#' # ===========================================================================
#' # Example 5: Explore tumor microenvironment
#' # ===========================================================================
#'
#' list_immune_cells(category = "Microenvironment")
#'
#' # Returns stromal cells:
#' # - Fibroblasts_xcell
#' # - Endothelial_cells_xcell
#' # - MSC_xcell (Mesenchymal stem cells)
#' # - Adipocytes_xcell
#'
#' # ===========================================================================
#' # Example 6: Get tumor purity scores
#' # ===========================================================================
#'
#' list_immune_cells(algorithm = "estimate")
#'
#' # Returns ESTIMATE scores:
#' # - ImmuneScore_estimate
#' # - StromalScore_estimate
#' # - TumorPurity_estimate
#'
#' # ===========================================================================
#' # Example 7: Programmatic access
#' # ===========================================================================
#'
#' # Get all immune cells
#' cells <- list_immune_cells()
#'
#' # Filter CD8 cells
#' cd8_cells <- cells$Cell_Name[cells$Category == "T_cells_CD8"]
#'
#' # Use in analysis
#' # tcga_correlation(
#' #   var1 = "TP53", var1_modal = "Mutation",
#' #   var2 = cd8_cells[1], var2_modal = "ImmuneCell", ...
#' # )
#'
#' # ===========================================================================
#' # Next Steps
#' # ===========================================================================
#' # After identifying cell types:
#' # 1. Use tcga_correlation() to correlate with genes/mutations/clinical
#' # 2. Use tcga_survival() to test prognostic significance
#' # 3. Use tcga_enrichment() to find associated pathways
#' # 4. Compare across algorithms to validate findings
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


#' Search Immune Cell Types by Keyword
#'
#' @description
#' Performs case-insensitive keyword search across all 99 immune cell types from
#' 8 deconvolution algorithms. Searches cell names with partial matching support.
#' Returns matched cells with algorithm and category information. Convenient wrapper
#' around \code{list_immune_cells()} for quick cell type discovery.
#'
#' @param keyword Character. Search keyword or pattern (required).
#'   Case-insensitive partial matching on cell names.
#'   Examples: "CD8", "Macrophage", "B_cells", "Treg", "NK".
#'
#' @return Data frame with 3 columns (invisible):
#'   \describe{
#'     \item{Cell_Name}{Full cell type name}
#'     \item{Algorithm}{Deconvolution algorithm}
#'     \item{Category}{Functional category}
#'   }
#'   Returns NULL if no matches found.
#'
#' @details
#' **Common Search Terms**:
#' \itemize{
#'   \item **T cells**: "CD8", "CD4", "T_cells", "Treg", "Cytotoxic"
#'   \item **B cells**: "B_cells", "Plasma", "Naive_B", "Memory_B"
#'   \item **Myeloid**: "Macrophage", "Monocyte", "Dendritic", "Neutrophil", "M1", "M2"
#'   \item **NK cells**: "NK", "Natural_Killer"
#'   \item **Stromal**: "Fibroblast", "Endothelial", "MSC", "Adipocyte"
#'   \item **Scores**: "Immune", "Stromal", "Purity"
#' }
#'
#' **Search Tips**:
#' \itemize{
#'   \item Partial matching: "Macro" matches "Macrophages_M0", "Macrophages_M1", etc.
#'   \item Case-insensitive: "cd8", "CD8", "Cd8" all work
#'   \item Use specific terms for fewer results: "CD8" better than "T"
#'   \item Use algorithm suffix to filter: "_cibersort", "_xcell"
#' }
#'
#' @section User Queries:
#' **Quick Search**:
#' \itemize{
#'   \item How do I quickly find CD8+ T cells?
#'   \item What macrophage types match "M1"?
#'   \item Can I search for B cells?
#'   \item How do I find regulatory T cells?
#'   \item What cells match "NK"?
#' }
#'
#' **Use Cases**:
#' \itemize{
#'   \item I know the cell type but not the exact variable name
#'   \item I want to see all variants of a cell type across algorithms
#'   \item I need to quickly check if a cell type exists
#'   \item I want to compare similar cells from different algorithms
#' }
#'
#' @references
#' **TCGA Database**:
#'
#' The Cancer Genome Atlas Research Network (2013). The Cancer Genome Atlas
#' Pan-Cancer analysis project. Nature Genetics, 45(10):1113-1120.
#' \doi{10.1038/ng.2764}
#'
#' Database portal: \url{https://www.cancer.gov/tcga}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{list_immune_cells}} - List all immune cells with filtering
#'   \item \code{\link{search_variables}} - Search across Clinical/Signature/ImmuneCell
#'   \item \code{\link{tcga_correlation}} - Use found cells in analysis
#' }
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Search for CD8+ T cells
#' # ===========================================================================
#'
#' search_immune_cells("CD8")
#'
#' # Returns CD8 variants:
#' # - CD8_T_cells_cibersort
#' # - T_cells_CD8_xcell
#' # - T_cells_CD8_quantiseq
#' # - CD8_T_cells_mcpcounter
#'
#' # ===========================================================================
#' # Example 2: Search for macrophages
#' # ===========================================================================
#'
#' search_immune_cells("Macrophage")
#'
#' # Returns all macrophage types:
#' # - Macrophages_M0_cibersort
#' # - Macrophages_M1_cibersort
#' # - Macrophages_M2_cibersort
#' # - Macrophages_xcell
#' # - Macrophages_M1_quantiseq
#' # - Macrophages_M2_quantiseq
#'
#' # ===========================================================================
#' # Example 3: Search for B cells
#' # ===========================================================================
#'
#' search_immune_cells("B_cells")
#'
#' # Returns B cell variants:
#' # - B_cells_naive_cibersort
#' # - B_cells_memory_cibersort
#' # - Plasma_cells_cibersort
#' # - B_cells_xcell
#'
#' # ===========================================================================
#' # Example 4: Quick check if cell type exists
#' # ===========================================================================
#'
#' # Check for regulatory T cells
#' treg_results <- search_immune_cells("Treg")
#'
#' if (!is.null(treg_results)) {
#'   cat("Found", nrow(treg_results), "Treg variants\n")
#' }
#'
#' # ===========================================================================
#' # Next Steps
#' # ===========================================================================
#' # After finding cell names:
#' # 1. Use in tcga_correlation(var2 = "CD8_T_cells_cibersort", var2_modal = "ImmuneCell")
#' # 2. Use in tcga_survival(var1 = "Macrophages_M1_cibersort", var1_modal = "ImmuneCell")
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


#' List All TCGA Cancer Types and Molecular Subtypes
#'
#' @description
#' Displays comprehensive catalog of all TCGA cancer types: 33 main cancer types,
#' 32 molecular subtypes across 14 cancers, and combined cancer groups. Cancer
#' type input is case-insensitive (BRCA, brca, Brca all work). Essential for
#' identifying available cancer types before analysis. Returns invisible data frame
#' for programmatic access.
#'
#' @param show_subtypes Logical. Display molecular subtypes (default: TRUE).
#'   TRUE shows main types + subtypes + combined groups.
#'   FALSE shows only 33 main cancer types.
#'   Subtypes provide finer molecular classification (e.g., BRCA-Basal, BRCA-LumA).
#'
#' @return Data frame with 2 columns (invisible):
#'   \describe{
#'     \item{Cancer_Type}{Cancer type code (e.g., "BRCA", "BRCA-Basal")}
#'     \item{Type}{Classification: "Main", "Subtype", or "Combined"}
#'   }
#'
#' @details
#' **33 Main Cancer Types**:
#' \itemize{
#'   \item **Breast**: BRCA (Breast invasive carcinoma)
#'   \item **Lung**: LUAD (Lung adenocarcinoma), LUSC (Lung squamous cell), SCLC (Small cell lung)
#'   \item **Colorectal**: COAD (Colon), READ (Rectum), COADREAD (Combined)
#'   \item **Brain**: GBM (Glioblastoma), LGG (Low-grade glioma), GBMLGG (Combined)
#'   \item **Kidney**: KIRC (Clear cell), KIRP (Papillary), KICH (Chromophobe)
#'   \item **Liver**: LIHC (Hepatocellular carcinoma)
#'   \item **Stomach**: STAD (Stomach adenocarcinoma)
#'   \item **Prostate**: PRAD (Prostate adenocarcinoma)
#'   \item **Skin**: SKCM (Skin cutaneous melanoma)
#'   \item **Thyroid**: THCA (Thyroid carcinoma)
#'   \item **Bladder**: BLCA (Bladder urothelial carcinoma)
#'   \item **Ovarian**: OV (Ovarian serous cystadenocarcinoma)
#'   \item **Cervical**: CESC (Cervical squamous cell)
#'   \item **Pancreas**: PAAD (Pancreatic adenocarcinoma)
#'   \item **Endometrial**: UCEC (Uterine corpus endometrial)
#'   \item **Head/Neck**: HNSC (Head and neck squamous cell)
#'   \item **Others**: ACC, CHOL, DLBC, ESCA, LAML, MESO, PCPG, SARC, TGCT, THYM, UCS, UVM
#' }
#'
#' **32 Molecular Subtypes** (14 cancers):
#' \itemize{
#'   \item **BRCA**: Basal, Her2, LumA, LumB (4 subtypes)
#'   \item **GBM**: Classical, Mesenchymal, Neural, Proneural (4 subtypes)
#'   \item **LGG**: IDHmut-codel, IDHmut-non-codel, IDHwt (3 subtypes)
#'   \item **LUAD**: Proximal-inflammatory, Proximal-proliferative, Terminal respiratory unit (3 subtypes)
#'   \item **LUSC**: Basal, Classical, Primitive, Secretory (4 subtypes)
#'   \item **STAD**: CIN, EBV, GS, MSI (4 subtypes)
#'   \item **COAD**: CIN, GS, HM-SNV, HM-indel, POLE (5 subtypes)
#'   \item **OV**: Differentiated, Immunoreactive, Mesenchymal, Proliferative (4 subtypes)
#'   \item **KIRC**: m1, m2, m3, m4 (4 subtypes)
#'   \item **UCEC**: CN_HIGH, CN_LOW, MSI, POLE (4 subtypes)
#'   \item **PRAD**: ERG, ETS, Triple-negative (3 subtypes)
#'   \item **HNSC**: Atypical, Basal, Classical, Mesenchymal (4 subtypes)
#'   \item **BLCA**: Basal, Luminal (2 subtypes)
#'   \item **LIHC**: iCluster1, iCluster2, iCluster3 (3 subtypes)
#' }
#'
#' **Usage Notes**:
#' \itemize{
#'   \item Case-insensitive: "BRCA", "brca", "Brca" all accepted
#'   \item Subtypes use hyphen: "BRCA-Basal", "BRCA-LumA"
#'   \item Combined groups: "COADREAD" (COAD + READ), "GBMLGG" (GBM + LGG)
#'   \item Use list_cancer_types() output to verify exact spelling
#' }
#'
#' @section User Queries:
#' **General Exploration**:
#' \itemize{
#'   \item What cancer types are available in TCGA?
#'   \item How many cancer types does SLTCGA support?
#'   \item What is the difference between main types and subtypes?
#'   \item Can I analyze breast cancer subtypes?
#'   \item What lung cancer types are available?
#' }
#'
#' **Specific Cancer Types**:
#' \itemize{
#'   \item Is breast cancer (BRCA) data available?
#'   \item Can I analyze lung adenocarcinoma (LUAD)?
#'   \item Is glioblastoma (GBM) included?
#'   \item What kidney cancer types exist?
#'   \item Is pancreatic cancer (PAAD) available?
#'   \item Can I analyze melanoma (SKCM)?
#'   \item Is liver cancer (LIHC) included?
#' }
#'
#' **Molecular Subtypes**:
#' \itemize{
#'   \item What breast cancer subtypes are available?
#'   \item Can I analyze basal vs luminal breast cancer?
#'   \item What are BRCA-Basal, BRCA-LumA, BRCA-LumB, BRCA-Her2?
#'   \item Are glioblastoma molecular subtypes available?
#'   \item Can I analyze IDH-mutant vs IDH-wildtype gliomas?
#'   \item What gastric cancer (STAD) subtypes exist?
#'   \item Are colon cancer molecular subtypes available?
#' }
#'
#' **Combined Groups**:
#' \itemize{
#'   \item Can I analyze colorectal cancer (COAD + READ) together?
#'   \item Is there a combined glioma group (GBM + LGG)?
#'   \item How do I use combined cancer groups?
#' }
#'
#' **Input Format**:
#' \itemize{
#'   \item Is cancer type case-sensitive?
#'   \item Can I use lowercase cancer names?
#'   \item What is the correct format for subtypes?
#'   \item How do I specify molecular subtypes in analysis functions?
#' }
#'
#' @references
#' **TCGA Database**:
#'
#' The Cancer Genome Atlas Research Network (2013). The Cancer Genome Atlas
#' Pan-Cancer analysis project. Nature Genetics, 45(10):1113-1120.
#' \doi{10.1038/ng.2764}
#'
#' Database portal: \url{https://www.cancer.gov/tcga}
#'
#' **Molecular Subtypes**:
#'
#' Hoadley KA, et al. (2018). Cell-of-Origin Patterns Dominate the Molecular
#' Classification of 10,000 Tumors from 33 Types of Cancer. Cell, 173(2):291-304.
#' \doi{10.1016/j.cell.2018.03.022}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{list_modalities}} - View available omics layers
#'   \item \code{\link{tcga_correlation}} - Use cancer types in correlation analysis
#'   \item \code{\link{tcga_enrichment}} - Pan-cancer enrichment analysis
#'   \item \code{\link{tcga_survival}} - Cancer-specific survival analysis
#' }
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: List all cancer types (main + subtypes)
#' # ===========================================================================
#' # Research Question: What cancer types can I analyze?
#'
#' list_cancer_types()
#'
#' # Output shows:
#' # - 33 main cancer types
#' # - 32 molecular subtypes (grouped by parent)
#' # - Combined groups (COADREAD, GBMLGG)
#'
#' # ===========================================================================
#' # Example 2: List main cancer types only
#' # ===========================================================================
#'
#' list_cancer_types(show_subtypes = FALSE)
#'
#' # Shows only 33 main cancer types (no subtypes)
#'
#' # ===========================================================================
#' # Example 3: Identify breast cancer subtypes
#' # ===========================================================================
#' # Research Question: What BRCA subtypes are available?
#'
#' list_cancer_types(show_subtypes = TRUE)
#'
#' # BRCA subtypes shown:
#' # - BRCA-Basal (Triple-negative, aggressive)
#' # - BRCA-Her2 (HER2-enriched)
#' # - BRCA-LumA (Luminal A, ER+, best prognosis)
#' # - BRCA-LumB (Luminal B, ER+, intermediate prognosis)
#'
#' # ===========================================================================
#' # Example 4: Programmatic access
#' # ===========================================================================
#'
#' # Get all cancer types
#' cancers <- list_cancer_types()
#'
#' # Filter main types only
#' main_types <- cancers$Cancer_Type[cancers$Type == "Main"]
#'
#' # Filter subtypes only
#' subtypes <- cancers$Cancer_Type[cancers$Type == "Subtype"]
#'
#' # Filter breast cancer subtypes
#' brca_subtypes <- subtypes[grepl("^BRCA-", subtypes)]
#' # Returns: "BRCA-Basal", "BRCA-Her2", "BRCA-LumA", "BRCA-LumB"
#'
#' # ===========================================================================
#' # Example 5: Use cancer types in analysis
#' # ===========================================================================
#'
#' # Main cancer type
#' # tcga_correlation(
#' #   var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA", ...
#' # )
#'
#' # Molecular subtype
#' # tcga_correlation(
#' #   var1 = "ESR1", var1_modal = "RNAseq", var1_cancers = "BRCA-LumA", ...
#' # )
#'
#' # Multiple cancers
#' # tcga_correlation(
#' #   var1 = "TP53", var1_modal = "Mutation",
#' #   var1_cancers = c("BRCA", "LUAD", "COAD"), ...
#' # )
#'
#' # ===========================================================================
#' # Example 6: Pan-cancer analysis
#' # ===========================================================================
#'
#' # Get all main cancer types for pan-cancer study
#' all_cancers <- list_cancer_types(show_subtypes = FALSE)
#' main_cancers <- all_cancers$Cancer_Type
#'
#' # Use in pan-cancer correlation
#' # tcga_correlation(
#' #   var1 = "TP53", var1_modal = "Mutation",
#' #   var1_cancers = main_cancers[1:10], # First 10 cancers
#' #   var2 = "TMB", var2_modal = "Signature",
#' #   var2_cancers = main_cancers[1:10]
#' # )
#'
#' # ===========================================================================
#' # Next Steps
#' # ===========================================================================
#' # After identifying cancer types:
#' # 1. Use in tcga_correlation(var1_cancers = "BRCA", var2_cancers = "BRCA")
#' # 2. Use in tcga_enrichment(var1_cancers = "LUAD")
#' # 3. Use in tcga_survival(var1_cancers = "BRCA-Basal")
#' # 4. Compare subtypes: var1_cancers = c("BRCA-Basal", "BRCA-LumA")
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
  cat("  • list_modalities()              - Show all 8 data modalities\n")
  cat("  • list_variables(modal='...')    - Show variables for a modality\n")
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
