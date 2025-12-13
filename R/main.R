# ==============================================================================
# Main Functions Layer
# ==============================================================================
# Three main exported functions: tcga_correlation, tcga_enrichment, tcga_survival
# ==============================================================================


#' TCGA Correlation and Association Analysis
#'
#' @description
#' Comprehensive correlation and association analysis covering 7 scenarios.
#' Automatically detects scenario and applies appropriate statistical test and visualization.
#'
#' @param var1 Character vector. Gene names, clinical variables, miRNA, or immune cells.
#' @param var1_modal Character. Omics layer for var1. Options:
#'   "RNAseq", "Mutation", "CNV", "Methylation", "miRNA",
#'   "Clinical", "Signature", "ImmuneCell"
#' @param var1_cancers Character vector. Cancer types (case insensitive).
#'   33 main types + 32 subtypes. Examples: "BRCA", "LUAD", "BRCA_IDC"
#' @param var2 Character vector. Second variable for correlation.
#' @param var2_modal Character. Omics layer for var2.
#' @param var2_cancers Character vector. Cancer types for var2.
#' @param method Character. Correlation method (default: "pearson").
#'   Options: "pearson", "spearman", "kendall"
#' @param use Character. Missing value handling (default: "pairwise.complete.obs").
#' @param p_adjust_method Character. Multiple testing correction (default: "BH").
#' @param alpha Numeric. Significance threshold (default: 0.05).
#' @param rnaseq_type Character. RNAseq normalization (default: "log2TPM").
#'   Options: "log2TPM", "log2RSEM", "log2FPKM", "log2Counts"
#' @param cnv_type Character. CNV algorithm (default: "SNP6_Array").
#' @param methylation_region Character. Methylation region (default: "Promoter_mean").
#' @param immune_algorithm Character. Immune algorithm (default: NULL for all).
#' @param plot_type Character. Plot type for Scenario 6 (Mutation vs ImmuneCell).
#'   Options: "auto" (default), "boxplot", "heatmap"
#'   - "auto": Use heatmap if >=8 cells, otherwise boxplot
#'   - "boxplot": Force multiple boxplots
#'   - "heatmap": Force heatmap with difference barplot
#'
#' @return A list with 3 components:
#'   \item{stats}{Data frame with statistical results}
#'   \item{plot}{Plot object (ggplot, patchwork, or ComplexHeatmap)}
#'   \item{raw_data}{Data frame with merged input data}
#'
#' @export
tcga_correlation <- function(var1,
                             var1_modal,
                             var1_cancers,
                             var2,
                             var2_modal,
                             var2_cancers,
                             method = "pearson",
                             use = "pairwise.complete.obs",
                             p_adjust_method = "BH",
                             alpha = 0.05,
                             rnaseq_type = "log2TPM",
                             cnv_type = "SNP6_Array",
                             methylation_region = "Promoter_mean",
                             immune_algorithm = NULL,
                             plot_type = "auto") {
  message("\n========================================")
  message("TCGA Correlation Analysis")
  message("========================================")

  # Load data
  loaded <- tcga_load_modality(
    var1 = var1,
    var1_modal = var1_modal,
    var1_cancers = var1_cancers,
    var2 = var2,
    var2_modal = var2_modal,
    var2_cancers = var2_cancers,
    surv_type = NULL,
    rnaseq_type = rnaseq_type,
    cnv_type = cnv_type,
    methylation_region = methylation_region,
    immune_algorithm = immune_algorithm
  )

  # Filter out features with all NA values
  var1_features_filtered <- loaded$var1_features
  var2_features_filtered <- loaded$var2_features

  var1_cols <- .extract_colname_from_label(loaded$var1_features, loaded$data)
  var2_cols <- .extract_colname_from_label(loaded$var2_features, loaded$data)

  # Check var1 for all-NA columns
  var1_valid <- sapply(var1_cols, function(col) {
    !all(is.na(loaded$data[[col]]))
  })

  if (any(!var1_valid)) {
    removed_features <- loaded$var1_features[!var1_valid]
    message("\n[Filter] Removing var1 features with no valid data:")
    for (feat in removed_features) {
      message(sprintf("  - %s (all NA)", feat))
    }
    var1_features_filtered <- loaded$var1_features[var1_valid]
  }

  # Check var2 for all-NA columns
  var2_valid <- sapply(var2_cols, function(col) {
    !all(is.na(loaded$data[[col]]))
  })

  if (any(!var2_valid)) {
    removed_features <- loaded$var2_features[!var2_valid]
    message("\n[Filter] Removing var2 features with no valid data:")
    for (feat in removed_features) {
      message(sprintf("  - %s (all NA)", feat))
    }
    var2_features_filtered <- loaded$var2_features[var2_valid]
  }

  # Check if we have any features left
  if (length(var1_features_filtered) == 0) {
    stop("All var1 features have no valid data (all NA)", call. = FALSE)
  }
  if (length(var2_features_filtered) == 0) {
    stop("All var2 features have no valid data (all NA)", call. = FALSE)
  }

  # Detect scenario (use filtered features)
  var1_types_filtered <- loaded$var1_types[names(loaded$var1_types) %in% var1_features_filtered]
  var2_types_filtered <- loaded$var2_types[names(loaded$var2_types) %in% var2_features_filtered]

  scenario_info <- .detect_correlation_scenario(
    var1_features = var1_features_filtered,
    var2_features = var2_features_filtered,
    var1_types = var1_types_filtered,
    var2_types = var2_types_filtered,
    n_cancers = length(unique(c(var1_cancers, var2_cancers)))
  )

  # Smart pairing: if var1_cancers == var2_cancers, only compare matching cancers
  if (length(var1_cancers) > 1 && identical(var1_cancers, var2_cancers)) {
    message("\n[Smart Pairing] Detected identical cancer lists - filtering to matching cancer pairs only")

    # Extract cancer from feature labels: "VAR (Modal, CANCER)" -> "CANCER"
    extract_cancer <- function(feature) {
      match <- regmatches(feature, regexpr(", ([^)]+)\\)", feature))
      if (length(match) > 0) {
        return(gsub(", |\\)", "", match))
      }
      return(NA)
    }

    var1_cancers_in_features <- sapply(loaded$var1_features, extract_cancer)
    var2_cancers_in_features <- sapply(loaded$var2_features, extract_cancer)

    # Only keep matching cancer pairs
    keep_pairs <- outer(var1_cancers_in_features, var2_cancers_in_features, "==")

    message(sprintf(
      "  Original pairs: %d x %d = %d",
      length(loaded$var1_features), length(loaded$var2_features),
      length(loaded$var1_features) * length(loaded$var2_features)
    ))
    message(sprintf("  Filtered to same-cancer pairs: %d", sum(keep_pairs, na.rm = TRUE)))
  }

  # Perform statistics
  message("\n[Statistics] Running analysis...")

  if (scenario_info$var1_class == "continuous" && scenario_info$var2_class == "continuous") {
    stats <- .stats_correlation(
      data = loaded$data,
      var1_features = var1_features_filtered,
      var2_features = var2_features_filtered,
      method = method,
      use = use,
      p_adjust_method = p_adjust_method
    )

    # Filter stats to matching cancer pairs if needed
    if (length(var1_cancers) > 1 && identical(var1_cancers, var2_cancers)) {
      var1_cancer_in_stats <- sapply(stats$var1_feature, extract_cancer)
      var2_cancer_in_stats <- sapply(stats$var2_feature, extract_cancer)
      same_cancer <- var1_cancer_in_stats == var2_cancer_in_stats & !is.na(var1_cancer_in_stats)
      stats <- stats[same_cancer, ]

      if (nrow(stats) == 0) {
        stop("No valid same-cancer correlations found", call. = FALSE)
      }
    }
  } else if (scenario_info$var1_class == "categorical" && scenario_info$var2_class == "categorical") {
    stats <- .stats_association(
      data = loaded$data,
      var1_features = var1_features_filtered,
      var2_features = var2_features_filtered,
      alpha = alpha,
      p_adjust_method = p_adjust_method
    )

    # Filter stats to matching cancer pairs if needed
    if (length(var1_cancers) > 1 && identical(var1_cancers, var2_cancers)) {
      extract_cancer <- function(feature) {
        match <- regmatches(feature, regexpr(", ([^)]+)\\)", feature))
        if (length(match) > 0) {
          return(gsub(", |\\)", "", match))
        }
        return(NA)
      }
      var1_cancer_in_stats <- sapply(stats$var1_feature, extract_cancer)
      var2_cancer_in_stats <- sapply(stats$var2_feature, extract_cancer)
      same_cancer <- var1_cancer_in_stats == var2_cancer_in_stats & !is.na(var1_cancer_in_stats)
      stats <- stats[same_cancer, ]

      if (nrow(stats) == 0) {
        stop("No valid same-cancer associations found", call. = FALSE)
      }
    }
  } else {
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

    # Filter stats to matching cancer pairs if needed
    if (length(var1_cancers) > 1 && identical(var1_cancers, var2_cancers)) {
      extract_cancer <- function(feature) {
        match <- regmatches(feature, regexpr(", ([^)]+)\\)", feature))
        if (length(match) > 0) {
          return(gsub(", |\\)", "", match))
        }
        return(NA)
      }

      # Check both var1 and var2 features (depends on which is categorical)
      if (scenario_info$var1_class == "categorical") {
        var1_cancer_in_stats <- sapply(stats$var1_feature, extract_cancer)
        var2_cancer_in_stats <- sapply(stats$var2_feature, extract_cancer)
      } else {
        var1_cancer_in_stats <- sapply(stats$var2_feature, extract_cancer)
        var2_cancer_in_stats <- sapply(stats$var1_feature, extract_cancer)
      }

      same_cancer <- var1_cancer_in_stats == var2_cancer_in_stats & !is.na(var1_cancer_in_stats)
      stats <- stats[same_cancer, ]

      if (nrow(stats) == 0) {
        stop("No valid same-cancer comparisons found", call. = FALSE)
      }
    }
  }

  message(sprintf("  Completed: %d pairwise comparison(s)", nrow(stats)))

  # Create plot
  message("\n[Visualization] Generating plot...")

  plot_result <- .dispatch_correlation_plot(
    scenario_info = scenario_info,
    data = loaded$data,
    stats = stats,
    var1_features = var1_features_filtered,
    var2_features = var2_features_filtered,
    plot_type = plot_type
  )

  plot_type_attr <- attr(plot_result, "plot_type")
  plot_label <- scenario_info$plot_type
  if (!is.null(plot_type_attr) && grepl("heatmap", plot_type_attr)) {
    plot_label <- "Heatmap"
  }

  message(sprintf(
    "  Plot: %s (%.1f × %.1f inches)",
    plot_label,
    attr(plot_result, "width"),
    attr(plot_result, "height")
  ))

  # Save plot
  output_dir <- file.path(getwd(), "sltcga_output")
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

  # Check plot type for saving
  plot_type_attr <- attr(plot_result, "plot_type")

  if (!is.null(plot_type_attr) && plot_type_attr %in% c("heatmap", "heatmap_immune")) {
    # ComplexHeatmap object
    png(
      filename = filepath,
      width = attr(plot_result, "width"),
      height = attr(plot_result, "height"),
      units = "in",
      res = 300
    )

    # Draw heatmap (caption is now in column_title)
    ComplexHeatmap::draw(plot_result,
      heatmap_legend_side = "right",
      annotation_legend_side = "right"
    )

    dev.off()
  } else {
    # ggplot object
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


#' TCGA Enrichment Analysis
#'
#' @description
#' Perform enrichment analysis (Scenarios 8-15).
#' TCGA uses RNAseq for enrichment (no Protein data available).
#'
#' @param var1 Character vector. Variable names (genes, clinical, signatures, etc.)
#' @param var1_modal Character. Modal type for var1.
#' @param var1_cancers Character vector. Cancer types.
#' @param analysis_type Character. Type of analysis (default: "enrichment").
#'   - "genome": Genome-wide scan → NetworkPlot or DotPlot
#'   - "enrichment": Pathway enrichment (GSEA) → Paired DotPlot or Matrix
#' @param enrich_database Character. Database for GSEA (default: "MsigDB").
#'   Options: "MsigDB", "GO", "KEGG", "Wiki", "Reactome", "Mesh", "HgDisease", "Enrichrdb"
#' @param enrich_ont Character. GO sub-ontology (default: "BP"). Only for GO.
#' @param genome_modal Character. Omics layer for genome scan (default: "RNAseq").
#'   **Note**: For TCGA, enrichment always uses RNAseq (no Protein data).
#' @param method Character. Correlation method (default: "pearson").
#' @param top_n Integer. Number of top results to display (default: 50).
#' @param n_workers Integer. Parallel workers for GSEA (default: 6).
#' @param rnaseq_type Character. RNAseq normalization (default: "log2TPM").
#' @param kegg_category Character. KEGG category (default: "pathway").
#' @param msigdb_category Character. MsigDB collection (default: "H").
#' @param hgdisease_source Character. Disease database (default: "do").
#' @param mesh_method Character. MeSH method (default: "gendoo").
#' @param mesh_category Character. MeSH category (default: "A").
#' @param enrichrdb_library Character. Enrichr library (default: "Cancer_Cell_Line_Encyclopedia").
#' @param immune_algorithm Character. Immune deconvolution algorithm (default: NULL for all).
#'   Options: "cibersort", "epic", "mcpcounter", "quantiseq", "timer", "xcell"
#'   Only used when var1_modal = "ImmuneCell"
#'
#' @return List with three components:
#'   \item{stats}{Enrichment results}
#'   \item{plot}{Plot object}
#'   \item{raw_data}{Complete genome-wide analysis results}
#'
#' @export
tcga_enrichment <- function(var1,
                            var1_modal,
                            var1_cancers,
                            analysis_type = "enrichment",
                            enrich_database = "MsigDB",
                            enrich_ont = "BP",
                            genome_modal = "RNAseq",
                            method = "pearson",
                            top_n = 50,
                            n_workers = 6,
                            rnaseq_type = "log2TPM",
                            kegg_category = "pathway",
                            msigdb_category = "H",
                            hgdisease_source = "do",
                            mesh_method = "gendoo",
                            mesh_category = "A",
                            enrichrdb_library = "Cancer_Cell_Line_Encyclopedia",
                            immune_algorithm = NULL) {
  message("\n========================================")
  message("TCGA Enrichment Analysis")
  message("========================================")

  # Validate: Clinical variables with genome scan
  if (any(var1_modal == "Clinical") && analysis_type == "genome") {
    stop(
      "Clinical variables cannot be used for genome-wide scans.\n",
      "Use tcga_correlation() instead.",
      call. = FALSE
    )
  }

  # Force genome_modal to RNAseq for enrichment (TCGA has no Protein data)
  if (analysis_type == "enrichment" && genome_modal != "RNAseq") {
    message("\n[Note] TCGA enrichment analysis uses RNAseq (no Protein data available)")
    message("       Your specified genome_modal '", genome_modal, "' has been overridden to 'RNAseq'\n")
    genome_modal <- "RNAseq"
  }

  # Load variable data
  loaded <- tcga_load_modality(
    var1 = var1,
    var1_modal = var1_modal,
    var1_cancers = var1_cancers,
    var2 = NULL,
    var2_modal = NULL,
    var2_cancers = NULL,
    rnaseq_type = rnaseq_type,
    immune_algorithm = immune_algorithm
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
      rnaseq_type = rnaseq_type,
      kegg_category = kegg_category,
      msigdb_category = msigdb_category,
      hgdisease_source = hgdisease_source,
      mesh_method = mesh_method,
      mesh_category = mesh_category,
      enrichrdb_library = enrichrdb_library
    )
  } else {
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
      rnaseq_type = rnaseq_type,
      kegg_category = kegg_category,
      msigdb_category = msigdb_category,
      hgdisease_source = hgdisease_source,
      mesh_method = mesh_method,
      mesh_category = mesh_category,
      enrichrdb_library = enrichrdb_library
    )
  }

  message("\n✓ Enrichment analysis completed")

  # Save plot
  output_dir <- file.path(getwd(), "sltcga_output")
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
  filename <- gsub("[^a-zA-Z0-9_.-]", "_", filename)
  filepath <- file.path(output_dir, filename)

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


#' TCGA Survival Analysis (Kaplan-Meier and Cox Regression)
#'
#' @description
#' Comprehensive survival analysis covering 2 scenarios:
#' - Scenario 16: Single feature → KM + Cox curves
#' - Scenario 17: Multiple features → Forest plot
#'
#' @param var1 Character vector. Gene names, clinical variables, etc.
#' @param var1_modal Character. Omics layer.
#' @param var1_cancers Character vector. Cancer types.
#' @param surv_type Character. Survival endpoint (default: "OS").
#'   Options: "OS", "DSS", "PFI", "DFI"
#' @param cutoff_type Character. Cutoff method for continuous variables (default: "optimal").
#'   Options: "optimal", "median", "mean", "quantile"
#' @param minprop Numeric. Minimum proportion per group (default: 0.1).
#' @param percent Numeric. Percentile for quantile cutoff (default: 0.25).
#' @param palette Character vector. Colors for survival curves.
#' @param show_cindex Logical. Show concordance index (default: TRUE).
#' @param rnaseq_type Character. RNAseq normalization (default: "log2TPM").
#' @param cnv_type Character. CNV algorithm (default: "SNP6_Array").
#' @param methylation_region Character. Methylation region (default: "Promoter_mean").
#' @param immune_algorithm Character. Immune algorithm (default: NULL).
#'
#' @return A list with 3 components:
#'   \item{stats}{Survival analysis results}
#'   \item{plot}{Plot object}
#'   \item{raw_data}{Merged data with survival columns}
#'
#' @export
tcga_survival <- function(var1,
                          var1_modal,
                          var1_cancers,
                          surv_type = "OS",
                          cutoff_type = "optimal",
                          minprop = 0.1,
                          percent = 0.25,
                          palette = c("#ED6355", "#41A98E", "#EFA63A", "#3a6ea5"),
                          show_cindex = TRUE,
                          rnaseq_type = "log2TPM",
                          cnv_type = "SNP6_Array",
                          methylation_region = "Promoter_mean",
                          immune_algorithm = NULL) {
  message("\n========================================")
  message("TCGA Survival Analysis")
  message("========================================")

  # Load variable data
  loaded_var <- tcga_load_modality(
    var1 = var1,
    var1_modal = var1_modal,
    var1_cancers = var1_cancers,
    var2 = NULL,
    var2_modal = NULL,
    var2_cancers = NULL,
    rnaseq_type = rnaseq_type,
    cnv_type = cnv_type,
    methylation_region = methylation_region,
    immune_algorithm = immune_algorithm
  )

  # Load survival data
  loaded_surv <- tcga_load_modality(
    var1 = surv_type,
    var1_modal = "Survival",
    var1_cancers = var1_cancers,
    var2 = NULL,
    var2_modal = NULL,
    var2_cancers = NULL,
    surv_type = surv_type
  )

  # Merge data
  merged_data <- .merge_modal_data_tcga(loaded_var$data, loaded_surv$data)

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

  # Save plot
  output_dir <- file.path(getwd(), "sltcga_output")
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
  filename <- gsub("[^a-zA-Z0-9_.-]", "_", filename)
  filepath <- file.path(output_dir, filename)

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
.dispatch_correlation_plot <- function(scenario_info, data, stats, var1_features, var2_features, plot_type = "auto") {
  if (scenario_info$scenario_id == 1) {
    .plot_scenario1(data, stats, var1_features[1], var2_features[1], scenario_info)
  } else if (scenario_info$scenario_id == 2) {
    .plot_scenario2(data, stats, var1_features, var2_features, scenario_info)
  } else if (scenario_info$scenario_id == 3) {
    .plot_scenario3(data, stats, var1_features, var2_features, scenario_info)
  } else if (scenario_info$scenario_id == 4) {
    .plot_scenario4(data, stats, var1_features, var2_features, scenario_info)
  } else if (scenario_info$scenario_id %in% c(5, 6)) {
    # Check if this is Categorical vs ImmuneCell scenario
    is_categorical_immune <- any(grepl("ImmuneCell", c(var1_features, var2_features)))

    # Determine which are categorical and continuous
    cat_features <- if (scenario_info$var1_class == "categorical") var1_features else var2_features
    con_features <- if (scenario_info$var1_class == "categorical") var2_features else var1_features

    n_continuous <- length(con_features)

    # Decide plot type
    use_heatmap <- FALSE
    if (is_categorical_immune) {
      if (plot_type == "heatmap") {
        use_heatmap <- TRUE
      } else if (plot_type == "auto" && n_continuous >= 8) {
        use_heatmap <- TRUE
        message("  Using heatmap (>=8 immune cells). Use plot_type='boxplot' to force boxplots.")
      }
    }

    if (use_heatmap) {
      scenario_info$plot_type <- "Heatmap"
      .plot_categorical_immune_heatmap(data, stats, cat_features, con_features)
    } else {
      .plot_scenario5_6(data, stats, var1_features, var2_features, scenario_info)
    }
  } else {
    .plot_scenario7(data, stats, var1_features, var2_features, scenario_info)
  }
}


#' Run categorical enrichment (DEA-based, Scenarios 8-11)
#' @keywords internal
.run_categorical_enrichment <- function(data, var_features, var1_cancers, genome_modal,
                                        analysis_type, enrich_database, enrich_ont,
                                        top_n, n_workers, rnaseq_type, kegg_category, msigdb_category,
                                        hgdisease_source, mesh_method, mesh_category, enrichrdb_library) {
  # Load genome-wide data
  message("\n[Step 1] Loading genome-wide RNAseq data...")
  genome_matrix <- .load_genome_data_tcga(var1_cancers, rnaseq_type)

  if (length(var_features) == 1) {
    # Scenarios 8 & 9: Single variable
    var_label <- var_features[1]
    var_col <- .extract_colname_from_label(c(var_label), data)[1]
    var_data <- factor(data[[var_col]])
    names(var_data) <- rownames(data)

    message("\n[Step 2] Performing DEA...")
    is_for_enrichment <- (analysis_type == "enrichment")
    dea_stats <- .stats_dea_genome(var_data, genome_matrix, var1_cancers, for_enrichment = is_for_enrichment)

    if (analysis_type == "genome") {
      # Scenario 8: NetworkPlot
      message("\n[Step 3] Creating NetworkPlot...")

      query_gene <- gsub("\\s*\\(.*\\)", "", var_label)

      dea_up <- dea_stats[dea_stats$logFC > 0, ]
      dea_down <- dea_stats[dea_stats$logFC < 0, ]
      dea_up <- dea_up[order(dea_up$pvalue), ]
      dea_down <- dea_down[order(dea_down$pvalue), ]
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
      # Scenario 9: GSEA
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
        method = NULL,
        top_n = top_n
      )

      return(list(
        stats = gsea_result,
        plot = plot_result,
        raw_data = dea_stats
      ))
    }
  } else {
    # Scenarios 10 & 11: Multiple variables
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
      # Scenario 10: DotPlot
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

      filtered_stats <- attr(plot_result, "filtered_stats")

      return(list(
        stats = if (!is.null(filtered_stats)) filtered_stats else combined_stats,
        plot = plot_result,
        raw_data = all_stats
      ))
    } else {
      # Scenario 11: GSEA Matrix
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
        method = NULL,
        use_mean = FALSE,
        feature_list = NULL,
        top_n = top_n,
        omics_type = "Mutation"
      )

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
                                       method, top_n, n_workers, rnaseq_type, kegg_category, msigdb_category,
                                       hgdisease_source, mesh_method, mesh_category, enrichrdb_library) {
  # Load genome-wide data
  message("\n[Step 1] Loading genome-wide RNAseq data...")
  genome_matrix <- .load_genome_data_tcga(var1_cancers, rnaseq_type)

  if (length(var_features) == 1) {
    # Scenarios 12 & 13: Single variable
    var_label <- var_features[1]
    var_col <- .extract_colname_from_label(c(var_label), data)[1]
    var_data <- as.numeric(data[[var_col]])
    names(var_data) <- rownames(data)

    message("\n[Step 2] Calculating correlations...")
    is_for_enrichment <- (analysis_type == "enrichment")
    cor_stats <- .stats_cor_genome(var_data, genome_matrix, var1_cancers, method, for_enrichment = is_for_enrichment)

    if (analysis_type == "genome") {
      # Scenario 12: NetworkPlot
      message("\n[Step 3] Creating NetworkPlot...")

      query_gene <- gsub("\\s*\\(.*\\)", "", var_label)
      modal_match <- regmatches(var_label, regexpr("\\(([^,]+),", var_label))
      modal_type <- if (length(modal_match) > 0) gsub("[\\(,]", "", modal_match) else "RNAseq"

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
      # Scenario 13: GSEA
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
      modal_type <- if (length(modal_match) > 0) gsub("[\\(,]", "", modal_match) else "RNAseq"

      plot_result <- .plot_gsea_paired(
        gsea_stats = gsea_result,
        var_name = gene_name,
        omics_type = modal_type,
        cancer_types = var1_cancers,
        enrich_type = enrich_database,
        GO_ont = enrich_ont,
        method = method,
        top_n = top_n
      )

      return(list(
        stats = gsea_result,
        plot = plot_result,
        raw_data = cor_stats
      ))
    }
  } else {
    # Scenarios 14 & 15: Multiple variables
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
      # Scenario 14: DotPlot
      message("\n[Step 3] Creating DotPlot...")

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

      filtered_stats <- attr(plot_result, "filtered_stats")

      return(list(
        stats = if (!is.null(filtered_stats)) filtered_stats else combined_stats,
        plot = plot_result,
        raw_data = all_stats
      ))
    } else {
      # Scenario 15: GSEA Matrix
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

      n_unique_cancers <- length(unique(var1_cancers))
      use_mean_expr <- (n_unique_cancers == 1 && length(var_features) > 1)
      feature_names <- if (use_mean_expr) gsub("\\s*\\(.*", "", var_features) else NULL

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

  var_col <- if (!is.null(var_col_override)) {
    var_col_override
  } else {
    .extract_colname_from_label(c(var_feature), merged_data)[1]
  }

  if (!var_col %in% colnames(merged_data)) {
    stop(sprintf("Variable column '%s' not found", var_col), call. = FALSE)
  }
  if (!time_col %in% colnames(merged_data)) {
    stop(sprintf("Time column '%s' not found", time_col), call. = FALSE)
  }
  if (!event_col %in% colnames(merged_data)) {
    stop(sprintf("Event column '%s' not found", event_col), call. = FALSE)
  }

  # Check if multi-level categorical variable (>2 levels)
  # If so, use forest plot instead of KM+Cox
  if (var_type == "categorical") {
    n_levels <- nlevels(merged_data[[var_col]])
    if (n_levels > 2) {
      message(sprintf("  Detected %d-level categorical variable", n_levels))
      message("  Switching to forest plot (KM plot not suitable for >2 groups)")

      # Use forest plot workflow for single multi-level variable
      var_types_named <- setNames(var_type, var_feature)
      return(.run_survival_forest(
        merged_data = merged_data,
        var_features = c(var_feature),
        var_types = var_types_named,
        time_col = time_col,
        event_col = event_col,
        surv_type = surv_type,
        cutoff_type = cutoff_type,
        minprop = minprop,
        var1_cancers = var_cancers
      ))
    }
  }

  if (var_type == "continuous") {
    message(sprintf("  Calculating %s cutoff...", cutoff_type))

    valid_idx <- complete.cases(merged_data[, c(var_col, time_col, event_col)])
    if (sum(valid_idx) < 10) {
      stop("Too few valid samples for survival analysis", call. = FALSE)
    }

    cutoff <- .calc_optimal_cutoff(merged_data[valid_idx, ], var_col, time_col, event_col, minprop)

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

  message("\n[Step 2] Performing Kaplan-Meier analysis...")
  km_result <- .perform_km_analysis(merged_data, group_col, time_col, event_col)

  message("\n[Step 3] Performing Cox regression...")
  cox_result <- .perform_cox_analysis(merged_data, var_col, time_col, event_col)

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
    var_col = var_col
  )

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

    # Extract cancer from feature label
    cancer_match <- regmatches(var_label, regexpr(",\\s*([^)]+)\\)", var_label))
    feature_cancer <- if (length(cancer_match) > 0) {
      gsub(",\\s*|\\)", "", cancer_match)
    } else {
      var1_cancers[1]
    }

    feature_time_col <- paste0(feature_cancer, "_", surv_type, "_time")
    feature_event_col <- paste0(feature_cancer, "_", surv_type, "_event")

    feature_data <- merged_data[merged_data$cancer_type == feature_cancer, ]

    if (var_type == "continuous") {
      tryCatch(
        {
          cutoff <- .calc_optimal_cutoff(feature_data, var_col, feature_time_col, feature_event_col, minprop)
        },
        error = function(e) {
          cutoff <- median(feature_data[[var_col]], na.rm = TRUE)
        }
      )
    }

    cox_result <- .perform_cox_analysis(feature_data, var_col, feature_time_col, feature_event_col)

    # Handle multi-level categorical variables
    if (cox_result$n_coefs == 1) {
      # Single coefficient (binary or continuous)
      forest_data <- rbind(forest_data, data.frame(
        variable = var_label,
        hr = cox_result$hr,
        hr_lower = cox_result$hr_lower,
        hr_upper = cox_result$hr_upper,
        p_value = cox_result$p_value,
        cindex = cox_result$cindex,
        stringsAsFactors = FALSE
      ))
    } else {
      # Multiple coefficients (multi-level categorical)
      var_base_name <- gsub(" \\(.*", "", var_label)

      for (j in seq_along(cox_result$coef_names)) {
        # Extract level name from coefficient name
        # E.g., "BRCA_Race_ClinicalASIAN" -> "ASIAN"
        coef_name <- cox_result$coef_names[j]
        level_name <- gsub(var_col, "", coef_name)

        # Create readable label
        var_label_with_level <- paste0(var_base_name, ": ", level_name)

        forest_data <- rbind(forest_data, data.frame(
          variable = var_label_with_level,
          hr = cox_result$hrs[j],
          hr_lower = cox_result$hr_lowers[j],
          hr_upper = cox_result$hr_uppers[j],
          p_value = cox_result$p_values[j],
          cindex = cox_result$cindex,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  message("\n[Step 2] Creating forest plot...")

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
