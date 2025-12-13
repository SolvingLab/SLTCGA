# ==============================================================================
# Main Functions Layer
# ==============================================================================
# Three main exported functions: cptac_correlation, cptac_enrichment, cptac_survival
# All functions automatically save plots following the same pattern
# ==============================================================================


#' CPTAC Correlation and Association Analysis
#'
#' @description
#' Comprehensive correlation and association analysis covering 7 scenarios:
#' \itemize{
#'   \item Scenario 1: 1 continuous vs 1 continuous → Scatter plot (CorPlot) or Lollipop
#'   \item Scenario 2: 1 vs multiple continuous → LollipopPlot or DotPlot
#'   \item Scenario 3: Multiple vs multiple continuous → DotPlot (correlation matrix)
#'   \item Scenario 4: 1 categorical vs 1 continuous → BoxPlot
#'   \item Scenario 5-6: Multiple BoxPlots
#'   \item Scenario 7: Categorical vs categorical → Percentage BarPlot or Heatmap
#' }
#' Automatically detects scenario and applies appropriate statistical test and visualization.
#'
#' @param var1 Character vector. Gene names or clinical variables for variable 1.
#'   Examples: "TP53", c("TP53", "EGFR"), c("KRAS", "EGFR", "ALK")
#'
#' @param var1_modal Character. Omics layer for var1. Options:
#'   "RNAseq", "Protein", "Phospho", "Mutation", "Clinical", "logCNA", "Methylation"
#'
#' @param var1_cancers Character vector. Cancer types for var1.
#'   Options: "BRCA", "LUAD", "COAD", "CCRCC", "GBM", "HNSCC", "LUSC", "OV", "PDAC", "UCEC"
#'   Can be single or multiple: "BRCA" or c("BRCA", "LUAD", "COAD")
#'
#' @param var2 Character vector. Gene names or clinical variables for variable 2.
#'   Same format as var1. Required for correlation analysis.
#'
#' @param var2_modal Character. Omics layer for var2. Same options as var1_modal.
#'
#' @param var2_cancers Character vector. Cancer types for var2.
#'   Can be same as or different from var1_cancers.
#'
#' @param method Character. Correlation method for continuous variables (default: "pearson").
#'   Options: "pearson", "spearman", "kendall"
#'   Note: Only used for continuous vs continuous scenarios. Ignored for categorical variables.
#'
#' @param use Character. Handling of missing values (default: "pairwise.complete.obs").
#'   Options: "everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs"
#'
#' @param p_adjust_method Character. Multiple testing correction method (default: "BH").
#'   Options: "BH" (Benjamini-Hochberg), "bonferroni", "holm", "hochberg", "hommel", "BY", "fdr", "none"
#'
#' @param alpha Numeric. Significance threshold for marking significant results (default: 0.05).
#'   Used for categorical vs continuous tests (Wilcoxon, Kruskal-Wallis).
#'
#' @return A list with 3 components:
#'   \describe{
#'     \item{stats}{Data frame with statistical results. Columns vary by scenario:
#'       \itemize{
#'         \item Continuous: var1_feature, var2_feature, r, p, p_adjusted, method
#'         \item Categorical: var1_feature, var2_feature, p_value, test_method, effect_size, odds_ratio, log2_or
#'         \item Mixed: categorical, continuous, p_value, test_method, effect_size, n_groups
#'       }}
#'     \item{plot}{Plot object (ggplot, patchwork, or ComplexHeatmap).
#'       Direct access: result$plot. Size info: attr(result$plot, "width/height")}
#'     \item{raw_data}{Data frame with merged input data (all samples and features)}
#'   }
#'
#' @details
#' **Statistical Methods**:
#' \itemize{
#'   \item Continuous vs Continuous: Pearson/Spearman correlation, filters diagonal (self-correlation)
#'   \item Categorical vs Categorical: Chi-square or Fisher's exact test, calculates Odds Ratio and log2(OR)
#'   \item Categorical vs Continuous: Wilcoxon (2 groups) or Kruskal-Wallis (3+ groups)
#' }
#'
#' **Visualization Features**:
#' \itemize{
#'   \item Single cancer: Titles show "CPTAC-BRCA", axis labels exclude cancer
#'   \item Multi-cancer: Titles show "CPTAC-Database", axis labels include cancer
#'   \item Heatmap: Shows log2(OR) with red (co-occurrence) and blue (mutual exclusivity)
#'   \item BoxPlot: Categorical variable on x-axis, continuous on y-axis (smart swap)
#' }
#'
#' @examples
#' \dontrun{
#' # Scenario 1: mRNA-Protein correlation
#' result <- cptac_correlation(
#'   var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
#'   var2 = "TP53", var2_modal = "Protein", var2_cancers = "BRCA"
#' )
#'
#' # Scenario 2: Protein vs multiple Phospho sites
#' result <- cptac_correlation(
#'   var1 = "AKT1", var1_modal = "Protein", var1_cancers = "BRCA",
#'   var2 = c("AKT1", "MTOR", "RPS6"), var2_modal = "Phospho", var2_cancers = "BRCA"
#' )
#'
#' # Scenario 3: Phospho correlation matrix (removes diagonal)
#' result <- cptac_correlation(
#'   var1 = "AKT1", var1_modal = "Phospho", var1_cancers = "BRCA",
#'   var2 = "AKT1", var2_modal = "Phospho", var2_cancers = "BRCA"
#' )
#'
#' # Scenario 4: Mutation impact on expression
#' result <- cptac_correlation(
#'   var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
#'   var2 = "EGFR", var2_modal = "RNAseq", var2_cancers = "LUAD"
#' )
#'
#' # Scenario 5: Multiple mutations vs protein
#' result <- cptac_correlation(
#'   var1 = "AKT1", var1_modal = "Protein", var1_cancers = "BRCA",
#'   var2 = c("PIK3CA", "TP53"), var2_modal = "Mutation", var2_cancers = "BRCA"
#' )
#'
#' # Scenario 6: Clinical vs Phospho
#' result <- cptac_correlation(
#'   var1 = "Tumor_Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
#'   var2 = "AKT1", var2_modal = "Phospho", var2_cancers = "LUAD"
#' )
#'
#' # Scenario 7: Co-mutation analysis (log2(OR) heatmap)
#' result <- cptac_correlation(
#'   var1 = c("KRAS", "EGFR", "ALK"), var1_modal = "Mutation", var1_cancers = "LUAD",
#'   var2 = c("TP53", "STK11"), var2_modal = "Mutation", var2_cancers = "LUAD"
#' )
#'
#' # Multi-cancer comparison
#' result <- cptac_correlation(
#'   var1 = "TP53", var1_modal = "RNAseq", var1_cancers = c("BRCA", "LUAD", "COAD"),
#'   var2 = "TP53", var2_modal = "Protein", var2_cancers = c("BRCA", "LUAD", "COAD")
#' )
#' }
#'
#' @export
cptac_correlation <- function(var1,
                              var1_modal,
                              var1_cancers,
                              var2,
                              var2_modal,
                              var2_cancers,
                              method = "pearson",
                              use = "pairwise.complete.obs",
                              p_adjust_method = "BH",
                              alpha = 0.05) {
  message("\n========================================")
  message("CPTAC Correlation Analysis")
  message("========================================")

  # Load data
  loaded <- cptac_load_modality(
    var1 = var1,
    var1_modal = var1_modal,
    var1_cancers = var1_cancers,
    var2 = var2,
    var2_modal = var2_modal,
    var2_cancers = var2_cancers,
    surv_type = NULL
  )

  # Detect scenario
  scenario_info <- .detect_correlation_scenario(
    var1_features = loaded$var1_features,
    var2_features = loaded$var2_features,
    var1_types = loaded$var1_types,
    var2_types = loaded$var2_types,
    n_cancers = length(unique(c(var1_cancers, var2_cancers)))
  )

  # Perform statistics
  message("\n[Statistics] Running analysis...")

  if (scenario_info$var1_class == "continuous" && scenario_info$var2_class == "continuous") {
    # Correlation
    stats <- .stats_correlation(
      data = loaded$data,
      var1_features = loaded$var1_features,
      var2_features = loaded$var2_features,
      method = method,
      use = use,
      p_adjust_method = p_adjust_method
    )
  } else if (scenario_info$var1_class == "categorical" && scenario_info$var2_class == "categorical") {
    # Association
    stats <- .stats_association(
      data = loaded$data,
      var1_features = loaded$var1_features,
      var2_features = loaded$var2_features,
      alpha = alpha,
      p_adjust_method = p_adjust_method
    )
  } else {
    # Group difference
    if (scenario_info$var1_class == "categorical") {
      cat_features <- loaded$var1_features
      con_features <- loaded$var2_features
    } else {
      cat_features <- loaded$var2_features
      con_features <- loaded$var1_features
    }

    stats <- .stats_group_difference(
      data = loaded$data,
      cat_features = cat_features,
      con_features = con_features,
      alpha = alpha,
      p_adjust_method = p_adjust_method
    )
  }

  message(sprintf("  Completed: %d pairwise comparison(s)", nrow(stats)))

  # Create plot
  message("\n[Visualization] Generating plot...")

  plot_result <- .dispatch_correlation_plot(
    scenario_info = scenario_info,
    data = loaded$data,
    stats = stats,
    var1_features = loaded$var1_features,
    var2_features = loaded$var2_features
  )

  message(sprintf(
    "  Plot: %s (%.1f × %.1f inches)",
    scenario_info$plot_type,
    attr(plot_result, "width"),
    attr(plot_result, "height")
  ))

  # Save plot
  output_dir <- file.path(getwd(), "slcptac_output")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  filename <- .generate_filename(
    analysis = "correlation",
    var1 = var1,
    var1_modal = var1_modal,
    var1_cancers = var1_cancers,
    var2 = var2,
    var2_modal = var2_modal,
    var2_cancers = var2_cancers
  )

  filepath <- file.path(output_dir, filename)

  # Check if it's a ComplexHeatmap object (needs special handling)
  if (!is.null(attr(plot_result, "plot_type")) && attr(plot_result, "plot_type") == "heatmap") {
    # Save ComplexHeatmap using png device
    png(
      filename = filepath,
      width = attr(plot_result, "width"),
      height = attr(plot_result, "height"),
      units = "in",
      res = 300
    )
    ComplexHeatmap::draw(plot_result)
    dev.off()
  } else {
    # Save ggplot using ggsave
    ggplot2::ggsave(
      filename = filepath,
      plot = plot_result,
      width = attr(plot_result, "width"),
      height = attr(plot_result, "height"),
      dpi = 300, limitsize = FALSE
    )
  }

  message(sprintf("\n✓ Plot saved: %s", filepath))
  message("========================================\n")

  return(list(
    stats = stats,
    plot = plot_result,
    raw_data = loaded$data
  ))
}


#' CPTAC Enrichment Analysis
#'
#' @description
#' Perform enrichment analysis (Scenarios 8-15):
#' - Categorical variable:
#'   * Genome-wide (Scenario 8): DEA → NetworkPlot
#'   * Enrichment (Scenario 9): DEA → GSEA → Paired DotPlot
#' - Multiple categorical variables:
#'   * Genome-wide (Scenario 10): Multi-DEA → DotPlot
#'   * Enrichment (Scenario 11): Multi-DEA → GSEA → Matrix DotPlot
#' - Continuous variable:
#'   * Genome-wide (Scenario 12): Correlation → NetworkPlot
#'   * Enrichment (Scenario 13): Correlation → GSEA → Paired DotPlot
#' - Multiple continuous variables:
#'   * Genome-wide (Scenario 14): Multi-Correlation → DotPlot
#'   * Enrichment (Scenario 15): Multi-Correlation → GSEA → Matrix DotPlot
#'
#' @param var1 Character vector. Variable names (genes or clinical variables)
#' @param var1_modal Character. Modal type for var1. Options: "RNAseq", "Protein", "Phospho", "Mutation", "Clinical", "logCNA", "Methylation"
#' @param var1_cancers Character vector. Cancer types. Options: "BRCA", "LUAD", "COAD", "CCRCC", "GBM", "HNSCC", "LUSC", "OV", "PDAC", "UCEC"
#'
#' @param analysis_type Character. Type of enrichment analysis (default: "enrichment")
#'   - "genome": Genome-wide scan (DEA or correlation) → NetworkPlot or DotPlot
#'   - "enrichment": Pathway enrichment (GSEA) → Paired DotPlot or Matrix
#'
#' @param enrich_database Character. Database for pathway enrichment (default: "MsigDB")
#'   Options: "MsigDB" (recommended), "GO", "KEGG", "Wiki", "Reactome", "Mesh", "HgDisease", "Enrichrdb"
#'   Note: Different databases have different numbers of gene sets (MsigDB Hallmark: 50, GO BP: ~15000, KEGG: ~300)
#'
#' @param enrich_ont Character. Gene Ontology sub-ontology, only used when enrich_database = "GO" (default: "BP")
#'   Options: "BP" (Biological Process), "CC" (Cellular Component), "MF" (Molecular Function), "all"
#'
#' @param genome_modal Character. Omics layer to scan in genome-wide analysis (default: "Protein")
#'   Options: "Protein", "RNAseq", "Phospho", "Methylation", "logCNA"
#'   **IMPORTANT**: For analysis_type = "enrichment" (GSEA), genome_modal is automatically set to "Protein" regardless of input
#'   Reason: GSEA requires stable gene-level expression, and Protein is the most suitable omics layer
#'
#' @param method Character. Correlation method for continuous variables (default: "pearson")
#'   Options: "pearson", "spearman", "kendall"
#'   Note: Only used for continuous variables (RNAseq, Protein, Phospho). Ignored for categorical variables (Mutation, Clinical)
#'
#' @param top_n Integer. Number of top pathways to display in plot (default: 50)
#'   Note: stats will return ALL pathways, but plot only shows top N most significant pathways for clarity
#'
#' @param n_workers Integer. Number of parallel workers for GSEA computation (default: 6)
#'   Tip: Increase for faster computation on multi-core systems, decrease if memory is limited
#'
#' @param kegg_category Character. KEGG database category (default: "pathway")
#'   Options: "pathway", "module", "enzyme", "disease", "drug", "network"
#'   Only used when enrich_database = "KEGG"
#'
#' @param msigdb_category Character. MsigDB collection (default: "H" for Hallmark)
#'   Options: "H" (Hallmark, 50 gene sets), "C1" (Positional), "C2-CGP" (Chemical/Genetic),
#'            "C2-CP" (Canonical Pathways), "C5-GO-BP" (GO Biological Process), etc.
#'   Only used when enrich_database = "MsigDB"
#'
#' @param hgdisease_source Character. Human disease database source (default: "do")
#'   Options: "do" (Disease Ontology), "ncg_v7", "ncg_v6", "disgenet", "covid19"
#'   Only used when enrich_database = "HgDisease"
#'
#' @param mesh_method Character. MeSH mapping method (default: "gendoo")
#'   Options: "gendoo", "gene2pubmed", "RBBH"
#'   Only used when enrich_database = "Mesh"
#'
#' @param mesh_category Character. MeSH descriptor category (default: "A")
#'   Only used when enrich_database = "Mesh"
#'
#' @param enrichrdb_library Character. Enrichr library name (default: "Cancer_Cell_Line_Encyclopedia")
#'   Only used when enrich_database = "Enrichrdb"
#'
#' @return List with three components:
#'   \item{stats}{Data frame with enrichment results
#'     - For genome scan: Top genes (controlled by top_n for NetworkPlot, or all significant for DotPlot)
#'     - For GSEA enrichment: ALL pathways with NES, p-value, q-value, etc.}
#'   \item{plot}{Plot object (patchwork or ggplot)
#'     - Direct access: result$plot (no need for result$plot$plot)
#'     - Width/height stored as attributes: attr(result$plot, "width"), attr(result$plot, "height")
#'     - Plot types: NetworkPlot, DotPlot Paired, GSEA Paired, or GSEA Matrix}
#'   \item{raw_data}{Complete genome-wide analysis results (all genes, not just top N)
#'     - For categorical variables: Full DEA results (data.frame with logFC, p-value for all genes)
#'     - For continuous variables: Full correlation results (data.frame with r, p-value for all genes)
#'     - For multiple variables: List of results for each variable}
#'
#' @examples
#' \dontrun{
#' # Example 1: Mutation vs Protein genome scan (Scenario 8)
#' result <- cptac_enrichment(
#'   var1 = "KRAS",
#'   var1_modal = "Mutation",
#'   var1_cancers = "LUAD",
#'   analysis_type = "genome",
#'   genome_modal = "Protein",
#'   top_n = 30
#' )
#'
#' # Example 2: Mutation vs GSEA enrichment (Scenario 9, default MsigDB Hallmark)
#' result <- cptac_enrichment(
#'   var1 = "PIK3CA",
#'   var1_modal = "Mutation",
#'   var1_cancers = "BRCA",
#'   analysis_type = "enrichment",
#'   top_n = 20 # genome_modal自动设为Protein
#' )
#'
#' # Example 3: Use different databases
#' # GO Biological Process
#' result <- cptac_enrichment(
#'   var1 = "TP53",
#'   var1_modal = "RNAseq",
#'   var1_cancers = "BRCA",
#'   analysis_type = "enrichment",
#'   enrich_database = "GO",
#'   enrich_ont = "BP"
#' )
#'
#' # KEGG Pathway
#' result <- cptac_enrichment(
#'   var1 = "EGFR",
#'   var1_modal = "Protein",
#'   var1_cancers = "LUAD",
#'   analysis_type = "enrichment",
#'   enrich_database = "KEGG",
#'   kegg_category = "pathway"
#' )
#'
#' # Reactome Pathways
#' result <- cptac_enrichment(
#'   var1 = "AKT1",
#'   var1_modal = "Phospho",
#'   var1_cancers = c("BRCA", "LUAD"),
#'   analysis_type = "enrichment",
#'   enrich_database = "Reactome"
#' )
#'
#' # Access results
#' head(result$stats) # All pathways
#' result$plot # View plot
#' head(result$raw_data) # Full DEA/correlation results
#' }
#'
#' @export
cptac_enrichment <- function(var1,
                             var1_modal,
                             var1_cancers,
                             analysis_type = "enrichment",
                             enrich_database = "MsigDB",
                             enrich_ont = "BP",
                             genome_modal = "Protein",
                             method = "pearson",
                             top_n = 50,
                             n_workers = 6,
                             kegg_category = "pathway",
                             msigdb_category = "H",
                             hgdisease_source = "do",
                             mesh_method = "gendoo",
                             mesh_category = "A",
                             enrichrdb_library = "Cancer_Cell_Line_Encyclopedia") {
  message("\n========================================")
  message("CPTAC Enrichment Analysis")
  message("========================================")

  # Validate: Clinical variables with genome scan
  if (any(var1_modal == "Clinical") && analysis_type == "genome") {
    stop(
      "Clinical variables cannot be used for genome-wide scans.\n",
      "Reason: Genome-wide DEA requires exactly 2 groups, but clinical variables often have >2 categories (e.g., Stage I/II/III/IV).\n",
      "\n",
      "Suggested alternatives:\n",
      "  1. Use cptac_correlation() to compare clinical variables with specific genes/proteins\n",
      "     Example: cptac_correlation(var1='Tumor_Stage', var1_modal='Clinical', var2=c('TP53','EGFR'), var2_modal='Protein')\n",
      "  \n",
      "  2. For binary clinical variables (e.g., Gender: Male/Female), enrichment analysis may work\n",
      "     Note: Ensure your clinical variable has exactly 2 categories before using enrichment",
      call. = FALSE
    )
  }

  # Validate and fix: enrichment分析必须使用Protein
  if (analysis_type == "enrichment" && genome_modal != "Protein") {
    message("\n[Note] For pathway enrichment analysis (GSEA), genome_modal is automatically set to 'Protein'")
    message("       Reason: GSEA requires gene-level expression data, and Protein is the most stable omics layer")
    message("       Your specified genome_modal '", genome_modal, "' has been overridden\n")
    genome_modal <- "Protein"
  }

  # Load variable data
  loaded <- cptac_load_modality(
    var1 = var1,
    var1_modal = var1_modal,
    var1_cancers = var1_cancers,
    var2 = NULL,
    var2_modal = NULL,
    var2_cancers = NULL
  )

  # Detect scenario
  scenario_info <- .detect_enrichment_scenario(
    var_features = loaded$var1_features,
    var_types = loaded$var1_types,
    analysis_type = analysis_type
  )

  # Perform analysis
  message("\n[Analysis] Running enrichment pipeline...")

  if (scenario_info$var_class == "categorical") {
    # DEA-based enrichment (Scenarios 8-11)
    result <- .run_categorical_enrichment(
      data = loaded$data,
      var_features = loaded$var1_features,
      var1_cancers = var1_cancers,
      genome_modal = genome_modal,
      analysis_type = analysis_type,
      enrich_database = enrich_database,
      enrich_ont = enrich_ont,
      top_n = top_n,
      n_workers = n_workers,
      kegg_category = kegg_category,
      msigdb_category = msigdb_category,
      hgdisease_source = hgdisease_source,
      mesh_method = mesh_method,
      mesh_category = mesh_category,
      enrichrdb_library = enrichrdb_library
    )
  } else {
    # Correlation-based enrichment (Scenarios 12-15)
    result <- .run_continuous_enrichment(
      data = loaded$data,
      var_features = loaded$var1_features,
      var1_cancers = var1_cancers,
      genome_modal = genome_modal,
      analysis_type = analysis_type,
      enrich_database = enrich_database,
      enrich_ont = enrich_ont,
      method = method,
      top_n = top_n,
      n_workers = n_workers,
      kegg_category = kegg_category,
      msigdb_category = msigdb_category,
      hgdisease_source = hgdisease_source,
      mesh_method = mesh_method,
      mesh_category = mesh_category,
      enrichrdb_library = enrichrdb_library
    )
  }

  message("\n✓ Enrichment analysis completed")

  # Save plot (following cptac_correlation pattern)
  output_dir <- file.path(getwd(), "slcptac_output")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  var_str <- paste(unique(gsub(" \\(.*\\)", "", loaded$var1_features)), collapse = "-")
  cancer_str <- paste(unique(var1_cancers), collapse = "-")
  modal_str <- paste(unique(var1_modal), collapse = "-")
  analysis_str <- if (analysis_type == "genome") "GenomeScan" else "GSEA"

  filename <- sprintf(
    "enrichment_%s_%s_%s_%s_%s.png",
    analysis_str, cancer_str, var_str, modal_str, genome_modal
  )
  filepath <- file.path(output_dir, filename)

  # Enrichment的plot直接是plot对象（带width/height属性）
  ggplot2::ggsave(
    filename = filepath,
    plot = result$plot,
    width = attr(result$plot, "width"),
    height = attr(result$plot, "height"),
    dpi = 300, limitsize = FALSE
  )

  message(sprintf("✓ Plot saved: %s", filepath))
  message("========================================\n")

  return(result)
}


#' CPTAC Survival Analysis (Kaplan-Meier and Cox Regression)
#'
#' @description
#' Comprehensive survival analysis covering 2 scenarios:
#' \itemize{
#'   \item Scenario 16: Single feature (1 gene in 1 cancer) → Kaplan-Meier curve + Cox regression
#'   \item Scenario 17: Multiple features → Forest plot showing hazard ratios
#'         Multiple features include: multiple genes, multiple cancers, or multiple phospho sites
#' }
#' Automatically detects scenario based on number of features and applies appropriate analysis.
#'
#' @param var1 Character vector. Gene names or clinical variables.
#'   \itemize{
#'     \item Single gene, single cancer → Scenario 16 (KM + Cox)
#'     \item Single gene, multiple cancers → Scenario 17 (Forest plot)
#'     \item Multiple genes → Scenario 17 (Forest plot)
#'     \item Phospho sites: Each site analyzed independently in Forest plot
#'   }
#'   Examples: "TP53", c("TP53", "EGFR"), c("AKT1", "MTOR")
#'
#' @param var1_modal Character. Omics layer. Options:
#'   "RNAseq", "Protein", "Phospho", "Mutation", "Clinical", "logCNA", "Methylation"
#'   Note: Phospho sites will generate multiple features per gene
#'
#' @param var1_cancers Character vector. Cancer types.
#'   Options: "BRCA", "LUAD", "COAD", "CCRCC", "GBM", "HNSCC", "LUSC", "OV", "PDAC", "UCEC"
#'   Multiple cancers will use Forest plot (Scenario 17)
#'
#' @param surv_type Character. Type of survival endpoint (default: "OS").
#'   \itemize{
#'     \item "OS" - Overall Survival (time to death from any cause)
#'     \item "PFS" - Progression-Free Survival (time to disease progression or death)
#'   }
#'
#' @param cutoff_type Character. Method to dichotomize continuous variables (default: "optimal").
#'   Only used for continuous variables (RNAseq, Protein, Phospho, logCNA, Methylation).
#'   Ignored for categorical variables (Mutation, Clinical).
#'   \itemize{
#'     \item "optimal" - Maximizes log-rank test statistic (recommended)
#'     \item "median" - Use median value as cutoff
#'     \item "mean" - Use mean value as cutoff
#'     \item "quantile" - Use specified percentile (see percent parameter)
#'   }
#'
#' @param minprop Numeric. Minimum proportion in each group for optimal cutoff (default: 0.1).
#'   Ensures at least 10\% samples in each group (High/Low). Range: 0.05-0.3
#'
#' @param percent Numeric. Percentile for quantile cutoff (default: 0.25).
#'   Only used when cutoff_type = "quantile". Range: 0-1 (0.25 = first quartile)
#'
#' @param palette Character vector. Colors for survival curves (default: c("#ED6355", "#41A98E", "#EFA63A", "#3a6ea5")).
#'   Typically 2 colors for Scenario 16 (High/Low or Mutation/WildType)
#'
#' @param show_cindex Logical. Show concordance index in plot (default: TRUE).
#'   C-index measures predictive accuracy (0.5 = random, 1.0 = perfect)
#'
#' @return A list with 3 components:
#'   \describe{
#'     \item{stats}{Data frame with survival analysis results:
#'       \itemize{
#'         \item Scenario 16: km_pvalue, cox_hr, cox_hr_lower, cox_hr_upper, cox_pvalue, cox_cindex
#'         \item Scenario 17: variable (feature label), hr, hr_lower, hr_upper, p_value, cindex
#'       }
#'       Note: hr > 1 indicates worse survival (higher risk)}
#'     \item{plot}{Plot object (ggplot or patchwork).
#'       Scenario 16: KM curve + Cox curve side-by-side.
#'       Scenario 17: Forest plot with HR and confidence intervals.
#'       Access: result$plot, attr(result$plot, "width/height")}
#'     \item{raw_data}{Data frame with merged data including survival time and event columns}
#'   }
#'
#' @details
#' **Scenario Detection**:
#' \itemize{
#'   \item n_features == 1 → Scenario 16 (Kaplan-Meier + Cox)
#'   \item n_features > 1 → Scenario 17 (Forest plot)
#'   \item Multiple cancers for 1 gene → Forest plot (each cancer independently)
#'   \item Multiple phospho sites → Forest plot (each site independently)
#' }
#'
#' **Statistical Methods**:
#' \itemize{
#'   \item Kaplan-Meier: Log-rank test for group comparison
#'   \item Cox Regression: Proportional hazards model for HR estimation
#'   \item Optimal cutoff: Maximizes separation using log-rank statistic
#' }
#'
#' **Multi-Cancer Handling**:
#' \itemize{
#'   \item Each cancer uses its own survival time/event columns
#'   \item Forest plot displays: "GENE (Modal, Cancer)" format
#'   \item Allows direct comparison of same gene across cancers
#' }
#'
#' @examples
#' \dontrun{
#' # Scenario 16: Single gene, single cancer (KM + Cox)
#' result <- cptac_survival(
#'   var1 = "TP53",
#'   var1_modal = "RNAseq",
#'   var1_cancers = "BRCA",
#'   surv_type = "OS",
#'   cutoff_type = "optimal"
#' )
#'
#' # Scenario 17: Multiple genes (Forest plot)
#' result <- cptac_survival(
#'   var1 = c("TP53", "EGFR", "KRAS"),
#'   var1_modal = "RNAseq",
#'   var1_cancers = "LUAD",
#'   surv_type = "OS"
#' )
#'
#' # Scenario 17: Single gene, multiple cancers (Forest plot)
#' result <- cptac_survival(
#'   var1 = "TP53",
#'   var1_modal = "RNAseq",
#'   var1_cancers = c("BRCA", "LUAD", "COAD"),
#'   surv_type = "OS"
#' )
#'
#' # Phospho sites survival (each site independent)
#' result <- cptac_survival(
#'   var1 = "AKT1",
#'   var1_modal = "Phospho",
#'   var1_cancers = c("BRCA", "LUAD"),
#'   surv_type = "OS",
#'   cutoff_type = "optimal"
#' )
#'
#' # Mutation survival (categorical, no cutoff needed)
#' result <- cptac_survival(
#'   var1 = "KRAS",
#'   var1_modal = "Mutation",
#'   var1_cancers = "LUAD",
#'   surv_type = "PFS"
#' )
#'
#' # Clinical variables
#' result <- cptac_survival(
#'   var1 = c("Age", "Tumor_Stage"),
#'   var1_modal = "Clinical",
#'   var1_cancers = "BRCA",
#'   surv_type = "OS"
#' )
#' }
#'
#' @export
cptac_survival <- function(var1,
                           var1_modal,
                           var1_cancers,
                           surv_type = "OS",
                           cutoff_type = "optimal",
                           minprop = 0.1,
                           percent = 0.25,
                           palette = c("#ED6355", "#41A98E", "#EFA63A", "#3a6ea5"),
                           show_cindex = TRUE) {
  message("\n========================================")
  message("CPTAC Survival Analysis")
  message("========================================")

  # Load variable data
  loaded_var <- cptac_load_modality(
    var1 = var1,
    var1_modal = var1_modal,
    var1_cancers = var1_cancers,
    var2 = NULL,
    var2_modal = NULL,
    var2_cancers = NULL
  )

  # Load survival data
  loaded_surv <- cptac_load_modality(
    var1 = surv_type,
    var1_modal = "Survival",
    var1_cancers = var1_cancers,
    var2 = NULL,
    var2_modal = NULL,
    var2_cancers = NULL,
    surv_type = surv_type
  )

  # Merge data
  merged_data <- .merge_modal_data(loaded_var$data, loaded_surv$data)

  # Detect scenario
  scenario_info <- .detect_survival_scenario(
    var_features = loaded_var$var1_features,
    var_types = loaded_var$var1_types,
    n_cancers = length(var1_cancers)
  )

  # Perform survival analysis
  message("\n[Analysis] Running survival analysis...")

  time_col <- paste0(var1_cancers[1], "_", surv_type, "_time")
  event_col <- paste0(var1_cancers[1], "_", surv_type, "_event")

  if (scenario_info$scenario_id == 16) {
    # Scenario 16: Single feature → KM + Cox
    result <- .run_survival_single(
      merged_data = merged_data,
      var_feature = loaded_var$var1_features[1],
      var_type = loaded_var$var1_types[1],
      time_col = time_col,
      event_col = event_col,
      surv_type = surv_type,
      cutoff_type = cutoff_type,
      minprop = minprop,
      palette = palette,
      var_cancers = var1_cancers,
      var_col_override = NULL
    )
  } else {
    # Scenario 17: Multiple variables → Forest plot
    result <- .run_survival_forest(
      merged_data = merged_data,
      var_features = loaded_var$var1_features,
      var_types = loaded_var$var1_types,
      time_col = time_col,
      event_col = event_col,
      surv_type = surv_type,
      cutoff_type = cutoff_type,
      minprop = minprop,
      var1_cancers = var1_cancers
    )
  }

  message("\n✓ Survival analysis completed")

  # Save plot (following cptac_correlation pattern)
  output_dir <- file.path(getwd(), "slcptac_output")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  var_str <- paste(unique(gsub(" \\(.*\\)", "", loaded_var$var1_features)), collapse = "-")
  cancer_str <- paste(unique(var1_cancers), collapse = "-")
  modal_str <- paste(unique(var1_modal), collapse = "-")
  analysis_str <- if (scenario_info$n_vars == 1) "KM_Cox" else "Forest"

  filename <- sprintf(
    "survival_%s_%s_%s_%s_%s.png",
    analysis_str, surv_type, cancer_str, var_str, modal_str
  )
  filepath <- file.path(output_dir, filename)

  # 统一使用属性方案
  ggplot2::ggsave(
    filename = filepath,
    plot = result$plot,
    width = attr(result$plot, "width"),
    height = attr(result$plot, "height"),
    dpi = 300,
    limitsize = FALSE
  )

  message(sprintf("✓ Plot saved: %s", filepath))
  message("========================================\n")

  return(result)
}


# ==============================================================================
# Internal Helper Functions
# ==============================================================================

#' Dispatch correlation plot based on scenario
#' @keywords internal
.dispatch_correlation_plot <- function(scenario_info, data, stats, var1_features, var2_features) {
  if (scenario_info$scenario_id == 1) {
    .plot_scenario1(data, stats, var1_features[1], var2_features[1], scenario_info)
  } else if (scenario_info$scenario_id == 2) {
    .plot_scenario2(data, stats, var1_features, var2_features, scenario_info)
  } else if (scenario_info$scenario_id == 3) {
    .plot_scenario3(data, stats, var1_features, var2_features, scenario_info)
  } else if (scenario_info$scenario_id == 4) {
    .plot_scenario4(data, stats, var1_features, var2_features, scenario_info)
  } else if (scenario_info$scenario_id %in% c(5, 6)) {
    .plot_scenario5_6(data, stats, var1_features, var2_features, scenario_info)
  } else {
    # Scenario 7
    .plot_scenario7(data, stats, var1_features, var2_features, scenario_info)
  }
}


#' Run categorical enrichment (DEA-based, Scenarios 8-11)
#' @keywords internal
.run_categorical_enrichment <- function(data, var_features, var1_cancers, genome_modal,
                                        analysis_type, enrich_database, enrich_ont,
                                        top_n, n_workers, kegg_category, msigdb_category,
                                        hgdisease_source, mesh_method, mesh_category, enrichrdb_library) {
  # Load genome-wide data
  message("\n[Step 1] Loading genome-wide data...")
  genome_matrix <- .load_genome_data(var1_cancers, genome_modal)

  if (length(var_features) == 1) {
    # ========== Scenarios 8 & 9: Single variable ==========
    var_label <- var_features[1]
    var_col <- .extract_colname_from_label(c(var_label), data)[1]
    var_data <- factor(data[[var_col]])
    names(var_data) <- rownames(data)

    # Perform DEA
    message("\n[Step 2] Performing DEA...")
    dea_stats <- .stats_dea_genome(var_data, genome_matrix, var1_cancers)

    if (analysis_type == "genome") {
      # Scenario 8: Genome-wide → NetworkPlot
      message("\n[Step 3] Creating NetworkPlot...")

      query_gene <- gsub("\\s*\\(.*\\)", "", var_label)

      # Select top 50 up + top 50 down genes
      dea_up <- dea_stats[dea_stats$logFC > 0, ]
      dea_down <- dea_stats[dea_stats$logFC < 0, ]
      p_col <- if ("P.Value" %in% colnames(dea_up)) "P.Value" else "pvalue"
      dea_up <- dea_up[order(dea_up[[p_col]]), ]
      dea_down <- dea_down[order(dea_down[[p_col]]), ]
      dea_top <- rbind(head(dea_up, 50), head(dea_down, 50))

      plot_result <- .plot_network(
        stats = dea_top,
        var1_name = query_gene,
        edge_metric = "logFC",
        query_omics = "Mutation",
        genome_omics = genome_modal,
        cancer_type = var1_cancers[1],
        analysis_type = "DEA",
        method = NULL
      )

      return(list(
        stats = dea_top,
        plot = plot_result,
        raw_data = dea_stats
      ))
    } else {
      # Scenario 9: Enrichment → GSEA
      message("\n[Step 3] Performing GSEA...")

      ranked_genes <- setNames(dea_stats$logFC, dea_stats$gene)
      ranked_genes <- sort(ranked_genes, decreasing = TRUE)

      gsea_result <- .perform_gsea(
        ranked_genes = ranked_genes,
        enrich_type = enrich_database,
        GO_ont = enrich_ont,
        n_workers = n_workers,
        kegg_category = kegg_category,
        msigdb_category = msigdb_category,
        hgdisease_source = hgdisease_source,
        mesh_method = mesh_method,
        mesh_category = mesh_category,
        enrichrdb_library = enrichrdb_library
      )

      message("\n[Step 4] Creating GSEA plots...")

      gene_name <- gsub("\\s*\\(.*", "", var_features[1])
      modal_match <- regmatches(var_features[1], regexpr("\\(([^,]+),", var_features[1]))
      modal_type <- if (length(modal_match) > 0) gsub("[\\(,]", "", modal_match) else "Mutation"

      plot_result <- .plot_gsea_paired(
        gsea_stats = gsea_result,
        var_name = gene_name,
        omics_type = modal_type,
        cancer_types = var1_cancers,
        enrich_type = enrich_database,
        GO_ont = enrich_ont,
        method = NULL, # DEA doesn't use correlation method
        top_n = top_n
      )

      return(list(
        stats = gsea_result,
        plot = plot_result,
        raw_data = dea_stats
      ))
    }
  } else {
    # ========== Scenarios 10 & 11: Multiple variables ==========
    message("\n[Step 2] Processing multiple variables...")

    all_stats <- list()

    for (var_label in var_features) {
      var_col <- .extract_colname_from_label(c(var_label), data)[1]
      var_data <- factor(data[[var_col]])
      names(var_data) <- rownames(data)

      message(sprintf("  Processing %s...", var_label))
      dea_stats <- .stats_dea_genome(var_data, genome_matrix, var1_cancers)
      dea_stats$var_name <- var_label
      all_stats[[var_label]] <- dea_stats
    }

    combined_stats <- do.call(rbind, all_stats)

    if (analysis_type == "genome") {
      # Scenario 10: Genome-wide → DotPlot
      message("\n[Step 3] Creating DotPlot...")

      plot_result <- .plot_dotplot_paired(
        all_stats = combined_stats,
        analysis_type = "DEA",
        cancer_types = var1_cancers,
        genome_omics = genome_modal,
        is_mutation = TRUE,
        use_mean = FALSE,
        feature_list = NULL,
        method = NULL,
        top_n = top_n
      )

      # 使用筛选后的数据作为stats（每个变量top N × 2）
      filtered_stats <- attr(plot_result, "filtered_stats")

      return(list(
        stats = if (!is.null(filtered_stats)) filtered_stats else combined_stats,
        plot = plot_result,
        raw_data = all_stats
      ))
    } else {
      # Scenario 11: Enrichment → GSEA Matrix
      message("\n[Step 3] Performing GSEA for multiple variables...")

      gsea_combined <- data.frame()

      for (var_label in var_features) {
        var_stats <- all_stats[[var_label]]
        ranked_genes <- setNames(var_stats$logFC, var_stats$gene)
        ranked_genes <- sort(ranked_genes, decreasing = TRUE)

        gsea_result <- .perform_gsea(
          ranked_genes = ranked_genes,
          enrich_type = enrich_database,
          GO_ont = enrich_ont,
          n_workers = n_workers,
          kegg_category = kegg_category,
          msigdb_category = msigdb_category,
          hgdisease_source = hgdisease_source,
          mesh_method = mesh_method,
          mesh_category = mesh_category,
          enrichrdb_library = enrichrdb_library
        )
        gsea_result$var_name <- var_label
        gsea_combined <- rbind(gsea_combined, gsea_result)
      }

      message("\n[Step 4] Creating GSEA matrix plot...")
      plot_result <- .plot_gsea_matrix(
        all_gsea_stats = gsea_combined,
        enrich_type = enrich_database,
        GO_ont = enrich_ont,
        cancer_types = var1_cancers,
        method = NULL, # DEA doesn't use correlation
        use_mean = FALSE,
        feature_list = NULL,
        top_n = top_n,
        omics_type = "Mutation" # Categorical变量都是Mutation或Clinical
      )

      # 使用筛选后的pathways作为stats
      filtered_stats <- attr(plot_result, "filtered_stats")

      return(list(
        stats = if (!is.null(filtered_stats)) filtered_stats else gsea_combined,
        plot = plot_result,
        raw_data = all_stats
      ))
    }
  }
}


#' Run continuous enrichment (Correlation-based, Scenarios 12-15)
#' @keywords internal
.run_continuous_enrichment <- function(data, var_features, var1_cancers, genome_modal,
                                       analysis_type, enrich_database, enrich_ont,
                                       method, top_n, n_workers, kegg_category, msigdb_category,
                                       hgdisease_source, mesh_method, mesh_category, enrichrdb_library) {
  # Load genome-wide data
  message("\n[Step 1] Loading genome-wide data...")
  genome_matrix <- .load_genome_data(var1_cancers, genome_modal)

  if (length(var_features) == 1) {
    # ========== Scenarios 12 & 13: Single variable ==========
    var_label <- var_features[1]
    var_col <- .extract_colname_from_label(c(var_label), data)[1]
    var_data <- as.numeric(data[[var_col]])
    names(var_data) <- rownames(data)

    # Perform correlation
    message("\n[Step 2] Calculating correlations...")
    cor_stats <- .stats_cor_genome(var_data, genome_matrix, var1_cancers, method)

    if (analysis_type == "genome") {
      # Scenario 12: Genome-wide → NetworkPlot
      message("\n[Step 3] Creating NetworkPlot...")

      query_gene <- gsub("\\s*\\(.*\\)", "", var_label)

      # Extract modal type from label
      modal_match <- regmatches(var_label, regexpr("\\(([^,]+),", var_label))
      modal_type <- if (length(modal_match) > 0) gsub("[\\(,]", "", modal_match) else "Protein"

      # Select top 50 positive + top 50 negative correlations
      cor_pos <- cor_stats[cor_stats$r > 0, ]
      cor_neg <- cor_stats[cor_stats$r < 0, ]
      cor_pos <- cor_pos[order(cor_pos$pvalue), ]
      cor_neg <- cor_neg[order(cor_neg$pvalue), ]
      cor_top <- rbind(head(cor_pos, 50), head(cor_neg, 50))

      plot_result <- .plot_network(
        stats = cor_top,
        var1_name = query_gene,
        edge_metric = "r",
        query_omics = modal_type,
        genome_omics = genome_modal,
        cancer_type = var1_cancers[1],
        analysis_type = "correlation",
        method = method
      )

      return(list(
        stats = cor_top,
        plot = plot_result,
        raw_data = cor_stats
      ))
    } else {
      # Scenario 13: Enrichment → GSEA
      message("\n[Step 3] Performing GSEA...")

      ranked_genes <- setNames(cor_stats$r, cor_stats$gene)
      ranked_genes <- sort(ranked_genes, decreasing = TRUE)

      gsea_result <- .perform_gsea(
        ranked_genes = ranked_genes,
        enrich_type = enrich_database,
        GO_ont = enrich_ont,
        n_workers = n_workers,
        kegg_category = kegg_category,
        msigdb_category = msigdb_category,
        hgdisease_source = hgdisease_source,
        mesh_method = mesh_method,
        mesh_category = mesh_category,
        enrichrdb_library = enrichrdb_library
      )

      message("\n[Step 4] Creating GSEA plots...")

      gene_name <- gsub("\\s*\\(.*", "", var_features[1])
      modal_match <- regmatches(var_features[1], regexpr("\\(([^,]+),", var_features[1]))
      modal_type <- if (length(modal_match) > 0) gsub("[\\(,]", "", modal_match) else "Protein"

      plot_result <- .plot_gsea_paired(
        gsea_stats = gsea_result,
        var_name = gene_name,
        omics_type = modal_type,
        cancer_types = var1_cancers,
        enrich_type = enrich_database,
        GO_ont = enrich_ont,
        method = method, # Pass correlation method
        top_n = top_n
      )

      return(list(
        stats = gsea_result,
        plot = plot_result,
        raw_data = cor_stats
      ))
    }
  } else {
    # ========== Scenarios 14 & 15: Multiple variables ==========
    message("\n[Step 2] Processing multiple variables...")

    all_stats <- list()

    for (var_label in var_features) {
      var_col <- .extract_colname_from_label(c(var_label), data)[1]
      var_data <- as.numeric(data[[var_col]])
      names(var_data) <- rownames(data)

      message(sprintf("  Processing %s...", var_label))
      cor_stats <- .stats_cor_genome(var_data, genome_matrix, var1_cancers, method)
      cor_stats$var_name <- var_label
      all_stats[[var_label]] <- cor_stats
    }

    combined_stats <- do.call(rbind, all_stats)

    if (analysis_type == "genome") {
      # Scenario 14: Genome-wide → DotPlot
      message("\n[Step 3] Creating DotPlot...")

      # Determine if we need mean expression
      n_unique_cancers <- length(unique(var1_cancers))
      use_mean_expr <- (n_unique_cancers == 1 && length(var_features) > 1)
      feature_names <- if (use_mean_expr) gsub("\\s*\\(.*", "", var_features) else NULL

      plot_result <- .plot_dotplot_paired(
        all_stats = combined_stats,
        analysis_type = "Correlation",
        cancer_types = var1_cancers,
        genome_omics = genome_modal,
        is_mutation = FALSE,
        use_mean = use_mean_expr,
        feature_list = feature_names,
        method = method,
        top_n = top_n
      )

      # 使用筛选后的数据作为stats（每个变量top N × 2）
      filtered_stats <- attr(plot_result, "filtered_stats")

      return(list(
        stats = if (!is.null(filtered_stats)) filtered_stats else combined_stats,
        plot = plot_result,
        raw_data = all_stats
      ))
    } else {
      # Scenario 15: Enrichment → GSEA Matrix
      message("\n[Step 3] Performing GSEA for multiple variables...")

      gsea_combined <- data.frame()

      for (var_label in var_features) {
        var_stats <- all_stats[[var_label]]
        ranked_genes <- setNames(var_stats$r, var_stats$gene)
        ranked_genes <- sort(ranked_genes, decreasing = TRUE)

        gsea_result <- .perform_gsea(
          ranked_genes = ranked_genes,
          enrich_type = enrich_database,
          GO_ont = enrich_ont,
          n_workers = n_workers,
          kegg_category = kegg_category,
          msigdb_category = msigdb_category,
          hgdisease_source = hgdisease_source,
          mesh_method = mesh_method,
          mesh_category = mesh_category,
          enrichrdb_library = enrichrdb_library
        )
        gsea_result$var_name <- var_label
        gsea_combined <- rbind(gsea_combined, gsea_result)
      }

      message("\n[Step 4] Creating GSEA matrix plot...")

      # Determine if we need mean expression
      n_unique_cancers <- length(unique(var1_cancers))
      use_mean_expr <- (n_unique_cancers == 1 && length(var_features) > 1)
      feature_names <- if (use_mean_expr) gsub("\\s*\\(.*", "", var_features) else NULL

      # Extract omics_type from first feature (format: "GENE (Modal, Cancer)")
      modal_match <- regmatches(var_features[1], regexpr("\\(([^,]+),", var_features[1]))
      modal_type <- if (length(modal_match) > 0) gsub("[\\(,]", "", modal_match) else NULL

      plot_result <- .plot_gsea_matrix(
        all_gsea_stats = gsea_combined,
        enrich_type = enrich_database,
        GO_ont = enrich_ont,
        cancer_types = var1_cancers,
        method = method,
        use_mean = use_mean_expr,
        feature_list = feature_names,
        top_n = top_n,
        omics_type = modal_type
      )

      # 使用筛选后的pathways作为stats（用于图表展示）
      filtered_stats <- attr(plot_result, "filtered_stats")

      return(list(
        stats = if (!is.null(filtered_stats)) filtered_stats else gsea_combined,
        plot = plot_result,
        raw_data = all_stats
      ))
    }
  }
}


#' Run single variable survival analysis (Scenario 16)
#' @keywords internal
.run_survival_single <- function(merged_data, var_feature, var_type, time_col, event_col,
                                 surv_type, cutoff_type, minprop, palette, var_cancers,
                                 var_col_override = NULL) {
  message("\n[Step 1] Preparing survival data...")

  # 使用override列名（用于phospho mean），否则从feature提取
  var_col <- if (!is.null(var_col_override)) {
    var_col_override
  } else {
    .extract_colname_from_label(c(var_feature), merged_data)[1]
  }

  # Validate columns exist
  if (!var_col %in% colnames(merged_data)) {
    stop(sprintf(
      "Variable column '%s' not found in data. Available: %s",
      var_col, paste(colnames(merged_data), collapse = ", ")
    ), call. = FALSE)
  }
  if (!time_col %in% colnames(merged_data)) {
    stop(sprintf("Time column '%s' not found in data", time_col), call. = FALSE)
  }
  if (!event_col %in% colnames(merged_data)) {
    stop(sprintf("Event column '%s' not found in data", event_col), call. = FALSE)
  }

  # Handle continuous vs categorical
  if (var_type == "continuous") {
    message(sprintf("  Calculating %s cutoff...", cutoff_type))

    # Filter to complete cases first
    valid_idx <- complete.cases(merged_data[, c(var_col, time_col, event_col)])
    if (sum(valid_idx) < 10) {
      stop("Too few valid samples for survival analysis (need at least 10)", call. = FALSE)
    }

    cutoff <- .calc_optimal_cutoff(merged_data[valid_idx, ], var_col, time_col, event_col, minprop)

    # Create group on ALL data (not just valid)
    merged_data$group <- ifelse(merged_data[[var_col]] > cutoff, "High", "Low")
    merged_data$group <- factor(merged_data$group, levels = c("Low", "High"))
    group_col <- "group"

    message(sprintf(
      "  Cutoff: %.3f (Low: n=%d, High: n=%d)",
      cutoff,
      sum(merged_data$group == "Low", na.rm = TRUE),
      sum(merged_data$group == "High", na.rm = TRUE)
    ))
  } else {
    group_col <- var_col
  }

  # Perform KM analysis
  message("\n[Step 2] Performing Kaplan-Meier analysis...")
  km_result <- .perform_km_analysis(merged_data, group_col, time_col, event_col)

  # Perform Cox analysis
  message("\n[Step 3] Performing Cox regression...")
  cox_result <- .perform_cox_analysis(merged_data, var_col, time_col, event_col)

  # Create combined plot
  message("\n[Step 4] Creating combined KM + Cox plot...")

  gene_name <- gsub("\\s*\\(.*", "", var_feature)
  modal_match <- regmatches(var_feature, regexpr("\\(([^,]+),", var_feature))
  modal_type <- if (length(modal_match) > 0) gsub("[\\(,]", "", modal_match) else "Unknown"

  cox_stats_df <- data.frame(
    variable = var_feature,
    hr = cox_result$hr,
    hr_lower = cox_result$hr_lower,
    hr_upper = cox_result$hr_upper,
    p_value = cox_result$p_value,
    stringsAsFactors = FALSE
  )

  plot_result <- .plot_km_cox_combined(
    km_fit = km_result$survfit,
    cox_model_stats = cox_stats_df,
    data = merged_data,
    time_col = time_col,
    event_col = event_col,
    group_col = group_col,
    var_name = gene_name,
    omics_type = modal_type,
    cancer_type = var_cancers[1],
    surv_type = surv_type,
    var_col = var_col # Pass the actual feature column
  )

  # Prepare stats
  stats <- data.frame(
    variable = var_feature,
    km_pvalue = km_result$p_value,
    cox_hr = cox_result$hr,
    cox_hr_lower = cox_result$hr_lower,
    cox_hr_upper = cox_result$hr_upper,
    cox_pvalue = cox_result$p_value,
    cox_cindex = cox_result$cindex,
    stringsAsFactors = FALSE
  )

  return(list(
    stats = stats,
    plot = plot_result,
    raw_data = merged_data
  ))
}


#' Run multiple variables survival analysis (Scenario 17)
#' @keywords internal
.run_survival_forest <- function(merged_data, var_features, var_types, time_col, event_col,
                                 surv_type, cutoff_type, minprop, var1_cancers) {
  message("\n[Step 1] Performing Cox regression for multiple variables...")

  forest_data <- data.frame()

  for (i in seq_along(var_features)) {
    var_label <- var_features[i]
    var_type <- var_types[i]
    var_col <- .extract_colname_from_label(c(var_label), merged_data)[1]

    message(sprintf("  Processing %s...", var_label))

    # 从feature label提取癌种
    cancer_match <- regmatches(var_label, regexpr(",\\s*([^)]+)\\)", var_label))
    feature_cancer <- if (length(cancer_match) > 0) {
      gsub(",\\s*|\\)", "", cancer_match)
    } else {
      var1_cancers[1]
    }

    # 使用该feature对应癌种的time/event列
    feature_time_col <- paste0(feature_cancer, "_", surv_type, "_time")
    feature_event_col <- paste0(feature_cancer, "_", surv_type, "_event")

    # 只用该feature对应癌种的数据
    feature_data <- merged_data[merged_data$cancer_type == feature_cancer, ]

    # For continuous variables, calculate cutoff (for reference only)
    if (var_type == "continuous") {
      tryCatch(
        {
          cutoff <- .calc_optimal_cutoff(feature_data, var_col, feature_time_col, feature_event_col, minprop)
        },
        error = function(e) {
          # 如果cutoff计算失败，使用median
          cutoff <- median(feature_data[[var_col]], na.rm = TRUE)
        }
      )
    }

    # Perform Cox regression
    cox_result <- .perform_cox_analysis(feature_data, var_col, feature_time_col, feature_event_col)

    # 使用完整的feature label（包含phospho site、modal、cancer）
    forest_data <- rbind(forest_data, data.frame(
      variable = var_label, # 使用完整label而不是只有基因名
      hr = cox_result$hr,
      hr_lower = cox_result$hr_lower,
      hr_upper = cox_result$hr_upper,
      p_value = cox_result$p_value,
      cindex = cox_result$cindex,
      stringsAsFactors = FALSE
    ))
  }

  # Create forest plot
  message("\n[Step 2] Creating forest plot...")

  # 提取所有涉及的癌种
  all_cancers <- unique(sapply(var_features, function(f) {
    match <- regmatches(f, regexpr(",\\s*([^)]+)\\)", f))
    if (length(match) > 0) gsub(",\\s*|\\)", "", match) else NA
  }))
  all_cancers <- all_cancers[!is.na(all_cancers)]
  cancer_label <- if (length(all_cancers) > 1) "Database" else all_cancers[1]

  plot_result <- .plot_forest(
    cox_stats = forest_data,
    surv_type = surv_type,
    cancer_type = cancer_label
  )

  return(list(
    stats = forest_data,
    plot = plot_result,
    raw_data = merged_data
  ))
}


#' Generate filename for saved plots
#' @keywords internal
.generate_filename <- function(analysis, var1, var1_modal, var1_cancers,
                               var2 = NULL, var2_modal = NULL, var2_cancers = NULL) {
  cancer_str <- paste(unique(c(var1_cancers, var2_cancers)), collapse = "-")
  var1_str <- paste(var1, collapse = "-")

  if (!is.null(var2)) {
    var2_str <- paste(var2, collapse = "-")
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

  return(filename)
}
