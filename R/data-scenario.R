# ==============================================================================
# Scenario Detection Layer
# ==============================================================================
# Detects analysis scenario based on variable types and counts
# Scenarios 1-17 covering correlation, enrichment, and survival analysis
# ==============================================================================


#' Detect Correlation Scenario (Scenarios 1-7)
#'
#' @description
#' Determines which of 7 correlation scenarios applies based on:
#' - Number of continuous variables
#' - Number of categorical variables
#' - Number of cancers
#'
#' Scenarios:
#' 1. 1 continuous vs 1 continuous
#' 2. 1 continuous vs multiple continuous (or vice versa)
#' 3. multiple continuous vs multiple continuous
#' 4. 1 categorical vs 1 categorical
#' 5. 1 categorical vs multiple categorical (or vice versa)
#' 6. multiple categorical vs multiple categorical
#' 7. categorical vs continuous (any combination)
#'
#' @param var1_features Character vector of var1 feature labels
#' @param var2_features Character vector of var2 feature labels
#' @param var1_types Character vector of var1 types
#' @param var2_types Character vector of var2 types
#' @param n_cancers Integer. Number of cancers
#'
#' @return List with scenario information
#'
#' @keywords internal
.detect_correlation_scenario <- function(var1_features,
                                         var2_features,
                                         var1_types,
                                         var2_types,
                                         n_cancers) {
  n_var1 <- length(var1_features)
  n_var2 <- length(var2_features)

  # Helper function to extract gene name from feature label
  # Format: "GENE (Modal, Cancer)" or "SITE_GENE (Modal, Cancer)"
  extract_gene_name <- function(features, remove_phospho_prefix = FALSE) {
    sapply(features, function(feat) {
      # Remove everything from " (" onwards
      gene_with_site <- gsub(" \\(.*", "", feat)
      # 如果需要去掉phospho site前缀（只用于计数基因数）
      if (remove_phospho_prefix) {
        gene_name <- gsub("^[STY][0-9]+_", "", gene_with_site)
        return(gene_name)
      } else {
        return(gene_with_site)
      }
    })
  }

  # 检查是否是phospho sites
  var1_is_phospho <- any(grepl("^[STY][0-9]+_", var1_features))
  var2_is_phospho <- any(grepl("^[STY][0-9]+_", var2_features))

  # Extract unique gene names
  # 对phospho，提取实际基因名（去掉site前缀）用于计数
  var1_genes <- unique(extract_gene_name(var1_features, remove_phospho_prefix = var1_is_phospho))
  var2_genes <- unique(extract_gene_name(var2_features, remove_phospho_prefix = var2_is_phospho))

  # Count unique genes
  n_var1_genes <- length(var1_genes)
  n_var2_genes <- length(var2_genes)

  # Determine if var1 and var2 are continuous or categorical
  var1_class <- if (all(var1_types == "continuous")) "continuous" else "categorical"
  var2_class <- if (all(var2_types == "continuous")) "continuous" else "categorical"

  # ============================================================================
  # Scenario detection logic
  # ============================================================================

  if (var1_class == "continuous" && var2_class == "continuous") {
    # Continuous vs Continuous (Scenarios 1-3)
    # Use unique gene counts (not feature counts across cancers)

    if (n_var1_genes == 1 && n_var2_genes == 1) {
      # Scenario 1: 1con vs 1con
      scenario_id <- 1
      scenario_name <- "1 continuous vs 1 continuous"
      plot_type <- if (n_cancers == 1) "CorPlot" else "LollipopPlot"
    } else if (n_var1_genes == 1 || n_var2_genes == 1) {
      # Scenario 2: 1con vs mcon
      scenario_id <- 2
      scenario_name <- "1 continuous vs multiple continuous"
      plot_type <- if (n_cancers == 1) "LollipopPlot" else "DotPlot"
    } else {
      # Scenario 3: mcon vs mcon
      scenario_id <- 3
      scenario_name <- "multiple continuous vs multiple continuous"
      plot_type <- "DotPlot"
    }
  } else if (var1_class == "categorical" && var2_class == "categorical") {
    # Categorical vs Categorical → Scenario 7 (Heatmap)
    scenario_id <- 7
    scenario_name <- "multiple categorical vs multiple categorical"
    plot_type <- "Heatmap"
  } else {
    # Mixed: Categorical vs Continuous (Scenarios 4-6)
    # S4: 1con vs 1cat → BoxPlot
    # S5: 1con vs mcat (or mcat vs 1con) → Multiple BoxPlots
    # S6: mcon vs 1cat (or 1cat vs mcon) → Multiple BoxPlots
    # S7: mcat vs mcon (or mcon vs mcat) → DotPlot/Heatmap
    # Use unique gene counts (not feature counts across cancers)

    n_con_genes <- if (var1_class == "continuous") n_var1_genes else n_var2_genes
    n_cat_genes <- if (var1_class == "categorical") n_var1_genes else n_var2_genes

    if (n_con_genes == 1 && n_cat_genes == 1) {
      # Scenario 4: 1con vs 1cat
      scenario_id <- 4
      scenario_name <- "1 continuous vs 1 categorical"
      plot_type <- "BoxPlot"
    } else if (n_con_genes == 1 && n_cat_genes > 1) {
      # Scenario 5: 1con vs mcat
      scenario_id <- 5
      scenario_name <- "1 continuous vs multiple categorical"
      plot_type <- "Multiple_BoxPlots"
    } else if (n_con_genes > 1 && n_cat_genes == 1) {
      # Scenario 6: mcon vs 1cat
      scenario_id <- 6
      scenario_name <- "multiple continuous vs 1 categorical"
      plot_type <- "Multiple_BoxPlots"
    } else {
      # Edge case: mcon vs mcat (not in original 7 scenarios, treat as S6)
      scenario_id <- 6
      scenario_name <- "multiple continuous vs multiple categorical (mixed)"
      plot_type <- "Multiple_BoxPlots"
    }
  }

  # ============================================================================
  # Return scenario info
  # ============================================================================

  # 构建更准确的描述
  # 检查是否是phospho sites（从原始features检查）
  var1_is_phospho <- any(grepl("^[STY][0-9]+_", var1_features))
  var2_is_phospho <- any(grepl("^[STY][0-9]+_", var2_features))

  var1_desc <- if (var1_is_phospho) {
    # Phospho sites: 总是显示phospho site数量
    if (n_var1_genes == 1) {
      sprintf("%d phospho site(s) from 1 gene", n_var1)
    } else {
      sprintf("%d phospho site(s) from %d genes", n_var1, n_var1_genes)
    }
  } else if (n_var1 != n_var1_genes) {
    # 多癌种：同一基因跨癌种
    sprintf("%d feature(s) from %d gene(s) across multiple cancers", n_var1, n_var1_genes)
  } else {
    # 普通情况
    sprintf("%d %s feature(s)", n_var1, var1_class)
  }

  var2_desc <- if (var2_is_phospho) {
    # Phospho sites: 总是显示phospho site数量
    if (n_var2_genes == 1) {
      sprintf("%d phospho site(s) from 1 gene", n_var2)
    } else {
      sprintf("%d phospho site(s) from %d genes", n_var2, n_var2_genes)
    }
  } else if (n_var2 != n_var2_genes) {
    # 多癌种：同一基因跨癌种
    sprintf("%d feature(s) from %d gene(s) across multiple cancers", n_var2, n_var2_genes)
  } else {
    # 普通情况
    sprintf("%d %s feature(s)", n_var2, var2_class)
  }

  message(sprintf(
    "\n[Scenario] Detected: Scenario %d - %s",
    scenario_id,
    scenario_name
  ))
  message(sprintf("  Var1: %s", var1_desc))
  message(sprintf("  Var2: %s", var2_desc))
  message(sprintf(
    "  Plot type: %s",
    plot_type
  ))

  return(list(
    scenario_id = scenario_id,
    scenario_name = scenario_name,
    var1_class = var1_class,
    var2_class = var2_class,
    n_var1 = n_var1,
    n_var2 = n_var2,
    n_var1_genes = n_var1_genes,
    n_var2_genes = n_var2_genes,
    n_cancers = n_cancers,
    plot_type = plot_type
  ))
}


#' Detect Enrichment Scenario (Scenarios 8-15)
#'
#' @description
#' Determines enrichment scenario based on:
#' - Variable type (continuous or categorical)
#' - Number of variables
#' - Analysis type (vs genome or vs enrichment)
#'
#' Scenarios:
#' 8. 1 categorical vs genome-wide
#' 9. 1 categorical vs enrichment
#' 10. multiple categorical vs genome-wide
#' 11. multiple categorical vs enrichment
#' 12. 1 continuous vs genome-wide
#' 13. 1 continuous vs enrichment
#' 14. multiple continuous vs genome-wide
#' 15. multiple continuous vs enrichment
#'
#' @param var_features Character vector of variable feature labels
#' @param var_types Character vector of variable types
#' @param analysis_type "genome" or "enrichment"
#'
#' @return List with scenario information
#'
#' @keywords internal
.detect_enrichment_scenario <- function(var_features,
                                        var_types,
                                        analysis_type) {
  n_vars <- length(var_features)

  # Extract unique gene names (without cancer/modal info)
  var_genes <- unique(sapply(var_features, function(feat) {
    gsub(" \\(.*", "", feat)
  }))
  n_unique_genes <- length(var_genes)

  var_class <- if (all(var_types == "continuous")) "continuous" else "categorical"

  # Scenario mapping (use unique gene count)
  if (var_class == "categorical") {
    if (n_unique_genes == 1) {
      if (analysis_type == "genome") {
        scenario_id <- 8
        scenario_name <- "1 categorical vs genome-wide"
        plot_type <- "NetworkPlot_Paired"
      } else {
        scenario_id <- 9
        scenario_name <- "1 categorical vs enrichment"
        plot_type <- "GSEA_DotPlot_Paired"
      }
    } else {
      if (analysis_type == "genome") {
        scenario_id <- 10
        scenario_name <- "multiple categorical vs genome-wide"
        plot_type <- "DotPlot_Paired"
      } else {
        scenario_id <- 11
        scenario_name <- "multiple categorical vs enrichment"
        plot_type <- "GSEA_DotPlot_Paired"
      }
    }
  } else {
    # continuous
    if (n_unique_genes == 1) {
      if (analysis_type == "genome") {
        scenario_id <- 12
        scenario_name <- "1 continuous vs genome-wide"
        plot_type <- "NetworkPlot_Paired"
      } else {
        scenario_id <- 13
        scenario_name <- "1 continuous vs enrichment"
        plot_type <- "GSEA_DotPlot_Paired"
      }
    } else {
      if (analysis_type == "genome") {
        scenario_id <- 14
        scenario_name <- "multiple continuous vs genome-wide"
        plot_type <- "DotPlot_Paired"
      } else {
        scenario_id <- 15
        scenario_name <- "multiple continuous vs enrichment"
        plot_type <- "GSEA_DotPlot_Paired"
      }
    }
  }

  message(sprintf(
    "\n[Scenario] Detected: Scenario %d - %s",
    scenario_id,
    scenario_name
  ))

  return(list(
    scenario_id = scenario_id,
    scenario_name = scenario_name,
    var_class = var_class,
    n_vars = n_vars,
    analysis_type = analysis_type,
    plot_type = plot_type
  ))
}


#' Detect Survival Scenario (Scenarios 16-17)
#'
#' @description
#' Determines survival scenario based on number of variables
#'
#' Scenarios:
#' 16. 1 variable vs survival (KM + Cox)
#' 17. multiple variables vs survival (Forest plot)
#'
#' @param var_features Character vector of variable feature labels
#' @param var_types Character vector of variable types
#' @param n_cancers Integer. Number of cancers
#'
#' @return List with scenario information
#'
#' @keywords internal
.detect_survival_scenario <- function(var_features,
                                      var_types,
                                      n_cancers) {
  n_vars <- length(var_features)

  var_class <- if (all(var_types == "continuous")) "continuous" else if (all(var_types == "categorical")) "categorical" else "mixed"

  # 简化逻辑：只看features数量
  # 1个feature → 场景16（KM+Cox）
  # 多个features → 场景17（Forest plot），无论来源（多癌种/多基因/多phospho sites）
  if (n_vars == 1) {
    scenario_id <- 16
    scenario_name <- "1 variable vs survival"
    plot_type <- "KM_Cox_Combined"
  } else {
    scenario_id <- 17
    scenario_name <- "multiple variables vs survival"
    plot_type <- "ForestPlot"
  }

  message(sprintf(
    "\n[Scenario] Detected: Scenario %d - %s",
    scenario_id,
    scenario_name
  ))
  message(sprintf(
    "  Variables: %d %s feature(s) across %d cancer(s)",
    n_vars,
    var_class,
    n_cancers
  ))

  return(list(
    scenario_id = scenario_id,
    scenario_name = scenario_name,
    var_class = var_class,
    n_vars = n_vars,
    n_cancers = n_cancers,
    plot_type = plot_type
  ))
}
