# ==============================================================================
# Statistical Analysis Layer
# ==============================================================================
# Statistical functions for correlation, DEA, enrichment, survival
# Adapted for TCGA data (uses RNAseq instead of Protein for enrichment)
# ==============================================================================


# ==============================================================================
# Correlation Statistics
# ==============================================================================

#' Calculate Correlation Statistics
#' @keywords internal
.stats_correlation <- function(data,
                               var1_features,
                               var2_features,
                               method = "pearson",
                               use = "pairwise.complete.obs",
                               p_adjust_method = "BH") {
  var1_cols <- .extract_colname_from_label(var1_features, data)
  var2_cols <- .extract_colname_from_label(var2_features, data)

  cor_data <- data[, c(var1_cols, var2_cols), drop = FALSE]

  cor_result <- astat::stat_cor(
    x = cor_data,
    cor.method = method,
    use = use,
    p.adjust = FALSE
  )

  stats <- cor_result[cor_result$x %in% var1_cols & cor_result$y %in% var2_cols, ]

  # Remove diagonal (self-correlation)
  stats <- stats[stats$x != stats$y, ]

  if (nrow(stats) == 0) {
    stop("No valid correlations calculated", call. = FALSE)
  }

  stats$p_adjusted <- p.adjust(stats$p, method = p_adjust_method)
  stats$method <- method
  stats$p_value <- stats$p
  stats$var1_feature <- .map_colname_to_label(stats$x, var1_cols, var1_features)
  stats$var2_feature <- .map_colname_to_label(stats$y, var2_cols, var2_features)
  stats$var1 <- stats$var1_feature
  stats$var2 <- stats$var2_feature

  # Extract cancer type
  stats$cancer_type <- sapply(strsplit(stats$var1_feature, ", "), function(x) {
    if (length(x) >= 2) gsub("\\)", "", x[2]) else NA_character_
  })

  stats$x <- NULL
  stats$y <- NULL
  stats <- stats[, c(
    "var1_feature", "var2_feature", "var1", "var2", "r", "p",
    "p_adjusted", "method", "p_value", "cancer_type"
  )]

  return(stats)
}


#' Calculate Association Statistics (Categorical vs Categorical)
#' @keywords internal
.stats_association <- function(data,
                               var1_features,
                               var2_features,
                               alpha = 0.05,
                               p_adjust_method = "BH") {
  var1_cols <- .extract_colname_from_label(var1_features, data)
  var2_cols <- .extract_colname_from_label(var2_features, data)

  stats_list <- list()

  for (v1 in var1_cols) {
    for (v2 in var2_cols) {
      valid_idx <- complete.cases(data[, c(v1, v2)])
      valid_data <- data[valid_idx, ]

      if (nrow(valid_data) < 3) next

      tbl <- table(valid_data[[v1]], valid_data[[v2]])

      odds_ratio <- NA
      log2_or <- NA
      cramers_v <- NA

      if (nrow(tbl) == 2 && ncol(tbl) == 2) {
        a <- tbl[2, 2]
        b <- tbl[2, 1]
        c <- tbl[1, 2]
        d <- tbl[1, 1]

        if (b == 0 || c == 0) {
          a <- a + 0.5
          b <- b + 0.5
          c <- c + 0.5
          d <- d + 0.5
        }

        odds_ratio <- (a * d) / (b * c)
        log2_or <- log2(odds_ratio)
      }

      if (any(tbl < 5)) {
        test_result <- fisher.test(tbl, simulate.p.value = TRUE)
        test_method <- "Fisher's exact"
        effect_size <- ifelse(!is.na(log2_or), log2_or, NA)
      } else {
        test_result <- chisq.test(tbl)
        test_method <- "Chi-square"
        cramers_v <- sqrt(test_result$statistic / (sum(tbl) * (min(dim(tbl)) - 1)))
        cramers_v <- as.numeric(cramers_v)
        effect_size <- ifelse(!is.na(log2_or), log2_or, cramers_v)
      }

      stats_list[[paste(v1, v2, sep = "_vs_")]] <- data.frame(
        var1 = v1,
        var2 = v2,
        p_value = test_result$p.value,
        test_method = test_method,
        effect_size = effect_size,
        odds_ratio = odds_ratio,
        log2_or = log2_or,
        cramers_v = cramers_v,
        n = nrow(valid_data),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(stats_list) == 0) {
    stop("No valid associations calculated", call. = FALSE)
  }

  stats <- do.call(rbind, stats_list)
  rownames(stats) <- NULL

  stats$p_adjusted <- p.adjust(stats$p_value, method = p_adjust_method)
  stats$significant <- stats$p_value < alpha
  stats$var1_feature <- .map_colname_to_label(stats$var1, var1_cols, var1_features)
  stats$var2_feature <- .map_colname_to_label(stats$var2, var2_cols, var2_features)
  stats$method <- stats$test_method

  return(stats)
}


#' Calculate Group Difference Statistics (Categorical vs Continuous)
#' @keywords internal
.stats_group_difference <- function(data,
                                    cat_features,
                                    con_features,
                                    alpha = 0.05,
                                    p_adjust_method = "BH") {
  cat_cols <- .extract_colname_from_label(cat_features, data)
  con_cols <- .extract_colname_from_label(con_features, data)

  stats_list <- list()

  for (cat_col in cat_cols) {
    for (con_col in con_cols) {
      valid_idx <- complete.cases(data[, c(cat_col, con_col)])
      valid_data <- data[valid_idx, ]

      if (nrow(valid_data) < 3) next

      # Check number of groups after removing NA
      cat_factor <- droplevels(valid_data[[cat_col]])
      n_groups <- nlevels(cat_factor)

      # Skip if only 1 group (need at least 2 for comparison)
      if (n_groups < 2) {
        warning(sprintf(
          "Skipping %s vs %s: only %d group found (need at least 2 groups). Check data overlap.",
          cat_col, con_col, n_groups
        ))
        next
      }

      n_groups <- length(unique(valid_data[[cat_col]]))

      if (n_groups == 2) {
        # Use backticks to handle special characters in column names (e.g., miRNA: hsa-let-7a-1)
        test_result <- wilcox.test(
          as.formula(paste0("`", con_col, "` ~ `", cat_col, "`")),
          data = valid_data
        )
        test_method <- "Wilcoxon"
        effect_size <- abs(qnorm(test_result$p.value / 2)) / sqrt(nrow(valid_data))
      } else {
        # Use backticks to handle special characters in column names
        test_result <- kruskal.test(
          as.formula(paste0("`", con_col, "` ~ `", cat_col, "`")),
          data = valid_data
        )
        test_method <- "Kruskal-Wallis"
        effect_size <- (test_result$statistic - n_groups + 1) / (nrow(valid_data) - n_groups)
        names(effect_size) <- NULL
      }

      stats_list[[paste(cat_col, con_col, sep = "_vs_")]] <- data.frame(
        categorical = cat_col,
        continuous = con_col,
        p_value = test_result$p.value,
        test_method = test_method,
        effect_size = effect_size,
        n_groups = n_groups,
        n = nrow(valid_data),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(stats_list) == 0) {
    stop("No valid group differences calculated", call. = FALSE)
  }

  stats <- do.call(rbind, stats_list)
  rownames(stats) <- NULL

  stats$p_adjusted <- p.adjust(stats$p_value, method = p_adjust_method)
  stats$significant <- stats$p_value < alpha
  stats$var1_feature <- .map_colname_to_label(stats$categorical, cat_cols, cat_features)
  stats$var2_feature <- .map_colname_to_label(stats$continuous, con_cols, con_features)
  stats$method <- stats$test_method

  return(stats)
}


# ==============================================================================
# Genome-wide DEA
# ==============================================================================

#' Perform genome-wide DEA using BioEnricher
#' @keywords internal
.stats_dea_genome <- function(var_data, genome_matrix, var_cancers, for_enrichment = FALSE) {
  # genome_matrix is gene × sample, need to check colnames for samples
  common_samples <- intersect(names(var_data), colnames(genome_matrix))

  if (length(common_samples) < 10) {
    stop("Too few overlapping samples (", length(common_samples), ")", call. = FALSE)
  }

  var_values <- var_data[common_samples]
  genome_values <- genome_matrix[, common_samples, drop = FALSE] # Select columns (samples)

  # Filter to protein-coding genes for enrichment analysis
  if (for_enrichment) {
    protein_genes <- .get_protein_coding_genes()
    if (!is.null(protein_genes)) {
      genes_to_keep <- intersect(rownames(genome_values), protein_genes)
      if (length(genes_to_keep) > 0) {
        genome_values <- genome_values[genes_to_keep, , drop = FALSE]
      } else {
        warning("No protein-coding genes found in genome data, using all genes")
      }
    }
  }

  var_values <- droplevels(var_values)

  group_levels <- levels(var_values)
  if (length(group_levels) == 0) {
    stop("No data available: all values are missing", call. = FALSE)
  }
  if (length(group_levels) == 1) {
    stop("DEA requires 2 groups, but all samples belong to the same group", call. = FALSE)
  }
  if (length(group_levels) != 2) {
    stop("DEA requires exactly 2 groups, found: ", length(group_levels), call. = FALSE)
  }

  group_counts <- table(var_values)
  if (any(group_counts < 3)) {
    stop(sprintf(
      "Too few samples in some groups: %s. Need at least 3 per group.",
      paste(paste0(names(group_counts), "=", group_counts), collapse = ", ")
    ), call. = FALSE)
  }

  message(sprintf(
    "  Performing DEA: %s vs %s (%d samples, %d genes)",
    group_levels[2], group_levels[1], length(common_samples), nrow(genome_values)
  ))

  expr_matrix <- genome_values # Already gene × sample
  groups <- var_values
  names(groups) <- common_samples
  contrasts <- paste0(group_levels[2], "-", group_levels[1])

  dea_result_list <- BioEnricher::DEA(
    expr = expr_matrix,
    groups = groups,
    contrasts = contrasts,
    Select.method = "limma",
    cutoff.P = 1,
    cutoff.logFC = 0,
    verbose = FALSE
  )

  dea_results <- dea_result_list$res

  colnames(dea_results)[colnames(dea_results) == "Gene"] <- "gene"
  colnames(dea_results)[colnames(dea_results) == "P.Value"] <- "pvalue"
  colnames(dea_results)[colnames(dea_results) == "adj.P.Val"] <- "p_adjusted"

  rownames(dea_results) <- NULL

  message(sprintf(
    "  ✓ DEA completed: %d significant genes (p<0.05)",
    sum(dea_results$pvalue < 0.05, na.rm = TRUE)
  ))

  return(dea_results)
}


#' Perform genome-wide correlation
#' @keywords internal
.stats_cor_genome <- function(var_data, genome_matrix, var_cancers, method = "pearson", for_enrichment = FALSE) {
  # genome_matrix is gene × sample, need to check colnames for samples
  common_samples <- intersect(names(var_data), colnames(genome_matrix))

  if (length(common_samples) < 10) {
    stop("Too few overlapping samples (", length(common_samples), ")", call. = FALSE)
  }

  var_values <- var_data[common_samples]
  genome_values <- genome_matrix[, common_samples, drop = FALSE] # Select columns (samples)

  # Filter to protein-coding genes for enrichment analysis
  if (for_enrichment) {
    protein_genes <- .get_protein_coding_genes()
    if (!is.null(protein_genes)) {
      genes_to_keep <- intersect(rownames(genome_values), protein_genes)
      if (length(genes_to_keep) > 0) {
        genome_values <- genome_values[genes_to_keep, , drop = FALSE]
      } else {
        warning("No protein-coding genes found in genome data, using all genes")
      }
    }
  }

  message(sprintf(
    "  Calculating correlations (%d samples, %d genes)",
    length(common_samples), nrow(genome_values)
  ))

  var_df <- data.frame(var = var_values)
  rownames(var_df) <- common_samples

  # Transpose genome_values to sample × gene for correlation
  genome_values_t <- t(genome_values)
  rownames(genome_values_t) <- common_samples

  cor_result <- astat::stat_cor(
    x = var_df,
    y = as.data.frame(genome_values_t),
    cor.method = method,
    use = "pairwise.complete.obs",
    p.adjust = FALSE
  )

  colnames(cor_result)[colnames(cor_result) == "y"] <- "gene"
  colnames(cor_result)[colnames(cor_result) == "p"] <- "pvalue"

  cor_result$p_adjusted <- p.adjust(cor_result$pvalue, method = "BH")

  message(sprintf(
    "  ✓ Correlation completed: %d significant genes (p<0.05)",
    sum(cor_result$pvalue < 0.05)
  ))

  return(cor_result)
}


#' Perform GSEA
#' @keywords internal
.perform_gsea <- function(ranked_genes,
                          enrich_type = "GO",
                          GO_ont = "BP",
                          n_workers = 6,
                          minSize = 10,
                          maxSize = 500,
                          kegg_category = "pathway",
                          msigdb_category = "H",
                          hgdisease_source = "do",
                          mesh_method = "gendoo",
                          mesh_category = "A",
                          enrichrdb_library = "Cancer_Cell_Line_Encyclopedia") {
  ont_str <- if (!is.null(GO_ont) && enrich_type == "GO") paste0(" ", GO_ont) else ""
  message(sprintf("  Running GSEA (%s%s)...", enrich_type, ont_str))

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' required. Install: BiocManager::install('fgsea')", call. = FALSE)
  }
  if (!requireNamespace("BioEnricher", quietly = TRUE)) {
    stop("Package 'BioEnricher' required", call. = FALSE)
  }

  geneset_df <- BioEnricher::get.geneset2(
    type = enrich_type,
    org = "hs",
    GO.ont = tolower(GO_ont),
    KEGG.category = kegg_category,
    Msigdb.category = msigdb_category,
    HgDisease.source = hgdisease_source,
    Mesh.method = mesh_method,
    Mesh.category = mesh_category,
    Enrichrdb.library = enrichrdb_library,
    Return.symbol = TRUE
  )

  pathways <- split(geneset_df$gene, geneset_df$id)
  pathway_descriptions <- setNames(geneset_df$term, geneset_df$id)

  message(sprintf("  Loaded %d pathways", length(pathways)))

  fgsea_result <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranked_genes,
    minSize = minSize,
    maxSize = maxSize,
    nproc = n_workers
  )

  stats <- as.data.frame(fgsea_result)

  if (nrow(stats) > 0) {
    stats$Description <- pathway_descriptions[stats$pathway]
    colnames(stats)[colnames(stats) == "pathway"] <- "ID"
    colnames(stats)[colnames(stats) == "pval"] <- "pvalue"
    colnames(stats)[colnames(stats) == "padj"] <- "qvalue"

    stats$leadingEdge <- sapply(stats$leadingEdge, function(x) {
      paste(unlist(x), collapse = ",")
    })

    cols <- c("ID", "Description", "NES", "pvalue", "qvalue", "size", "leadingEdge")
    stats <- stats[, intersect(cols, colnames(stats)), drop = FALSE]
    stats <- stats[order(stats$pvalue), ]
  }

  message(sprintf("  ✓ GSEA completed: %d pathways", nrow(stats)))

  return(stats)
}


# ==============================================================================
# Survival Analysis Functions
# ==============================================================================

#' Calculate optimal cutoff for continuous variable
#' @keywords internal
.calc_optimal_cutoff <- function(data, var_col, time_col, event_col, minprop = 0.1) {
  if (requireNamespace("survminer", quietly = TRUE)) {
    cutoff_df <- data.frame(
      time = data[[time_col]],
      event = data[[event_col]],
      var = data[[var_col]]
    )

    cutoff_df <- cutoff_df[complete.cases(cutoff_df), ]

    if (nrow(cutoff_df) < 20) {
      warning("Too few samples for optimal cutoff, using median")
      return(median(data[[var_col]], na.rm = TRUE))
    }

    tryCatch(
      {
        cutoff_result <- survminer::surv_cutpoint(
          data = cutoff_df,
          time = "time",
          event = "event",
          variables = "var",
          minprop = minprop
        )
        return(cutoff_result$cutpoint$cutpoint[1])
      },
      error = function(e) {
        warning("Optimal cutoff failed, using median: ", e$message)
        return(median(data[[var_col]], na.rm = TRUE))
      }
    )
  } else {
    return(median(data[[var_col]], na.rm = TRUE))
  }
}


#' Perform Kaplan-Meier analysis
#' @keywords internal
.perform_km_analysis <- function(data, group_col, time_col, event_col) {
  complete_idx <- complete.cases(data[, c(group_col, time_col, event_col)])
  n_complete <- sum(complete_idx)

  if (n_complete == 0) {
    stop("No complete observations for survival analysis", call. = FALSE)
  }

  if (n_complete < 10) {
    stop(sprintf("Too few complete observations (%d). Need at least 10.", n_complete),
      call. = FALSE
    )
  }

  data_complete <- data[complete_idx, ]

  # Use backticks to handle special characters in column names
  formula_str <- paste0("survival::Surv(`", time_col, "`, `", event_col, "`) ~ `", group_col, "`")
  survfit_obj <- survival::survfit(as.formula(formula_str), data = data_complete)
  survdiff_obj <- survival::survdiff(as.formula(formula_str), data = data_complete)

  survfit_obj$call$formula <- as.formula(formula_str)

  p_value <- 1 - pchisq(survdiff_obj$chisq, df = length(survdiff_obj$n) - 1)

  return(list(
    survfit = survfit_obj,
    survdiff = survdiff_obj,
    p_value = p_value
  ))
}


#' Perform Cox regression
#' @keywords internal
.perform_cox_analysis <- function(data, var_col, time_col, event_col) {
  complete_idx <- complete.cases(data[, c(var_col, time_col, event_col)])
  n_complete <- sum(complete_idx)

  if (n_complete == 0) {
    stop("No complete observations for Cox regression", call. = FALSE)
  }

  if (n_complete < 10) {
    stop(sprintf("Too few complete observations (%d). Need at least 10.", n_complete),
      call. = FALSE
    )
  }

  data_complete <- data[complete_idx, ]

  # Use backticks to handle special characters in column names
  cox_formula <- as.formula(paste0("survival::Surv(`", time_col, "`, `", event_col, "`) ~ `", var_col, "`"))
  cox_model <- survival::coxph(cox_formula, data = data_complete)
  cox_summary <- summary(cox_model)

  # Extract all coefficients (for multi-level categorical variables)
  coefs <- coef(cox_model)
  ci <- exp(confint(cox_model))
  pvals <- cox_summary$coefficients[, 5]
  cindex <- cox_summary$concordance[1]

  # For single coefficient (binary or continuous)
  if (length(coefs) == 1) {
    hr <- exp(coefs)
    hr_lower <- ci[1]
    hr_upper <- ci[2]
    p_value <- pvals

    return(list(
      model = cox_model,
      summary = cox_summary,
      hr = hr,
      hr_lower = hr_lower,
      hr_upper = hr_upper,
      p_value = p_value,
      cindex = cindex,
      n_coefs = 1,
      coef_names = NULL
    ))
  } else {
    # For multiple coefficients (multi-level categorical)
    hrs <- exp(coefs)

    if (is.matrix(ci)) {
      hr_lowers <- ci[, 1]
      hr_uppers <- ci[, 2]
    } else {
      hr_lowers <- rep(ci[1], length(coefs))
      hr_uppers <- rep(ci[2], length(coefs))
    }

    coef_names <- names(coefs)

    return(list(
      model = cox_model,
      summary = cox_summary,
      hrs = hrs,
      hr_lowers = hr_lowers,
      hr_uppers = hr_uppers,
      p_values = pvals,
      cindex = cindex,
      n_coefs = length(coefs),
      coef_names = coef_names
    ))
  }
}
