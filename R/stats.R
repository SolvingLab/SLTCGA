# ==============================================================================
# Statistical Analysis Layer
# ==============================================================================
# General-purpose statistical functions for correlation, DEA, enrichment, survival
# ==============================================================================


# ==============================================================================
# Correlation Statistics
# ==============================================================================

#' Calculate Correlation Statistics
#'
#' @description
#' General correlation function for continuous vs continuous
#'
#' @param data Data frame with features
#' @param var1_features Character vector of var1 feature column names
#' @param var2_features Character vector of var2 feature column names
#' @param method Correlation method ("pearson", "spearman", "kendall")
#' @param use Missing value handling
#' @param p_adjust_method P-value adjustment method
#'
#' @return Data frame with correlation statistics
#'
#' @keywords internal
.stats_correlation <- function(data,
                               var1_features,
                               var2_features,
                               method = "pearson",
                               use = "pairwise.complete.obs",
                               p_adjust_method = "BH") {
  # Extract feature columns (remove annotation part)
  var1_cols <- .extract_colname_from_label(var1_features, data)
  var2_cols <- .extract_colname_from_label(var2_features, data)

  # Prepare data
  cor_data <- data[, c(var1_cols, var2_cols), drop = FALSE]

  # Calculate correlations using astat
  cor_result <- astat::stat_cor(
    x = cor_data,
    cor.method = method,
    use = use,
    p.adjust = FALSE
  )

  # Filter to only var1 vs var2 pairs
  stats <- cor_result[cor_result$x %in% var1_cols & cor_result$y %in% var2_cols, ]

  # 过滤对角线（自己和自己的相关性永远是1，无意义）
  stats <- stats[stats$x != stats$y, ]

  if (nrow(stats) == 0) {
    stop("No valid correlations calculated", call. = FALSE)
  }

  # P-value adjustment
  stats$p_adjusted <- p.adjust(stats$p, method = p_adjust_method)

  # Add metadata
  stats$method <- method
  stats$p_adjust_method <- p_adjust_method

  # Add both old and new column names for compatibility
  stats$p_value <- stats$p # Keep both p and p_value
  stats$var1_feature <- .map_colname_to_label(stats$x, var1_cols, var1_features)
  stats$var2_feature <- .map_colname_to_label(stats$y, var2_cols, var2_features)

  # Keep simple column names for plotting (for backwards compatibility)
  stats$var1 <- stats$var1_feature
  stats$var2 <- stats$var2_feature

  # Extract cancer_type from feature labels (format: "GENE (Modal, Cancer)")
  # 优先从var1_feature提取，如果失败则从var2_feature提取
  stats$cancer_type <- sapply(strsplit(stats$var1_feature, ", "), function(x) {
    if (length(x) >= 2) {
      gsub("\\)", "", x[2])
    } else {
      NA_character_
    }
  })

  # 如果var1_feature没有癌种信息，尝试从var2_feature提取
  na_idx <- is.na(stats$cancer_type)
  if (any(na_idx)) {
    stats$cancer_type[na_idx] <- sapply(strsplit(stats$var2_feature[na_idx], ", "), function(x) {
      if (length(x) >= 2) {
        gsub("\\)", "", x[2])
      } else {
        NA_character_
      }
    })
  }

  # Remove only x and y (keep p for backwards compatibility)
  stats$x <- NULL
  stats$y <- NULL
  stats <- stats[, c(
    "var1_feature", "var2_feature", "var1", "var2", "r", "p",
    "p_adjusted", "method", "p_adjust_method", "p_value", "cancer_type"
  )]
  return(stats)
}


#' Calculate Association Statistics (Categorical vs Categorical)
#'
#' @description
#' Chi-square or Fisher's exact test for categorical associations
#'
#' @param data Data frame
#' @param var1_features Character vector of var1 feature labels
#' @param var2_features Character vector of var2 feature labels
#' @param alpha Significance level
#' @param p_adjust_method P-value adjustment method
#'
#' @return Data frame with association statistics
#'
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
      # Remove NAs
      valid_idx <- complete.cases(data[, c(v1, v2)])
      valid_data <- data[valid_idx, ]

      if (nrow(valid_data) < 3) next

      # Contingency table
      tbl <- table(valid_data[[v1]], valid_data[[v2]])

      # Calculate Odds Ratio for 2x2 tables
      odds_ratio <- NA
      log2_or <- NA
      cramers_v <- NA

      if (nrow(tbl) == 2 && ncol(tbl) == 2) {
        # 2x2 table: calculate OR
        # tbl format: rows=var1, cols=var2
        a <- tbl[2, 2] # Both positive
        b <- tbl[2, 1] # var1 positive, var2 negative
        c <- tbl[1, 2] # var1 negative, var2 positive
        d <- tbl[1, 1] # Both negative

        # Add continuity correction to avoid division by zero
        if (b == 0 || c == 0) {
          a <- a + 0.5
          b <- b + 0.5
          c <- c + 0.5
          d <- d + 0.5
        }

        odds_ratio <- (a * d) / (b * c)
        log2_or <- log2(odds_ratio)
      }

      # Choose test
      if (any(tbl < 5)) {
        test_result <- fisher.test(tbl, simulate.p.value = TRUE)
        test_method <- "Fisher's exact"
        # Fisher's test provides OR for 2x2 tables
        if (!is.na(odds_ratio)) {
          effect_size <- log2_or # Use log2(OR) as effect size
        } else {
          effect_size <- NA
        }
      } else {
        test_result <- chisq.test(tbl)
        test_method <- "Chi-square"
        # Cramer's V as alternative measure
        cramers_v <- sqrt(test_result$statistic / (sum(tbl) * (min(dim(tbl)) - 1)))
        cramers_v <- as.numeric(cramers_v)
        # Use log2(OR) if available, otherwise use Cramer's V
        if (!is.na(log2_or)) {
          effect_size <- log2_or
        } else {
          effect_size <- cramers_v
        }
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

  # P-value adjustment
  stats$p_adjusted <- p.adjust(stats$p_value, method = p_adjust_method)

  # Significant判断
  stats$significant <- stats$p_value < alpha

  # Map back to labels (use var1_feature, var2_feature for consistency)
  stats$var1_feature <- .map_colname_to_label(stats$var1, var1_cols, var1_features)
  stats$var2_feature <- .map_colname_to_label(stats$var2, var2_cols, var2_features)

  # Also keep method column for consistency
  stats$method <- stats$test_method

  return(stats)
}


#' Calculate Group Difference Statistics (Categorical vs Continuous)
#'
#' @description
#' Wilcoxon or Kruskal-Wallis test for group differences
#'
#' @param data Data frame
#' @param cat_features Character vector of categorical feature labels
#' @param con_features Character vector of continuous feature labels
#' @param alpha Significance level
#' @param p_adjust_method P-value adjustment method
#'
#' @return Data frame with test statistics
#'
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

      # Determine test
      n_groups <- length(unique(valid_data[[cat_col]]))

      if (n_groups == 2) {
        test_result <- wilcox.test(
          as.formula(paste(con_col, "~", cat_col)),
          data = valid_data
        )
        test_method <- "Wilcoxon"
        # Effect size: r = Z / sqrt(N)
        effect_size <- abs(qnorm(test_result$p.value / 2)) / sqrt(nrow(valid_data))
      } else {
        test_result <- kruskal.test(
          as.formula(paste(con_col, "~", cat_col)),
          data = valid_data
        )
        test_method <- "Kruskal-Wallis"
        # Effect size: epsilon squared
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

  # P-value adjustment
  stats$p_adjusted <- p.adjust(stats$p_value, method = p_adjust_method)

  # Significant判断
  stats$significant <- stats$p_value < alpha

  # Map back to labels (use var1_feature, var2_feature for consistency)
  stats$var1_feature <- .map_colname_to_label(stats$categorical, cat_cols, cat_features)
  stats$var2_feature <- .map_colname_to_label(stats$continuous, con_cols, con_features)

  # Also keep method column for consistency
  stats$method <- stats$test_method

  return(stats)
}


# ==============================================================================
# Genome-wide Analysis Functions
# ==============================================================================

#' Load genome-wide data from pre-processed files
#' @keywords internal
.load_genome_data <- function(cancers, genome_modal) {
  bulk_path <- Sys.getenv("SL_BULK_DATA")
  if (bulk_path == "") {
    stop("SL_BULK_DATA environment variable not set", call. = FALSE)
  }

  # Map modal type to file name
  file_mapping <- list(
    RNAseq = "LinkedOmicsKB_PanCancer_RNAseq_RSEM.qs",
    Protein = "LinkedOmicsKB_PanCancer_Protein_Quantification.qs",
    Phospho = "LinkedOmicsKB_PanCancer_Phospho_Quantification.qs",
    Methylation = "LinkedOmicsKB_PanCancer_Methylation.qs",
    logCNA = "LinkedOmicsKB_PanCancer_CNV_logCNA.qs",
    Mutation = "LinkedOmicsKB_PanCancer_Mutation_Binary.qs"
  )

  if (!genome_modal %in% names(file_mapping)) {
    stop("Unsupported genome_modal: ", genome_modal, call. = FALSE)
  }

  genome_file <- file.path(bulk_path, file_mapping[[genome_modal]])

  if (!file.exists(genome_file)) {
    stop("Genome file not found: ", genome_file, call. = FALSE)
  }

  message(sprintf("  Loading genome-wide %s data from pre-processed file...", genome_modal))

  # Read the entire genome matrix
  genome_data <- qs::qread(genome_file)

  # Add cancer type info
  genome_data$cancer_type <- sapply(strsplit(rownames(genome_data), "_"), `[`, 1)

  # Load clinical to filter tumor samples
  clin_file <- file.path(bulk_path, "LinkedOmicsKB_PanCancer_Clin.qs")
  clin <- qs::qread(clin_file)
  clin$type <- ifelse(is.na(clin$type), "Tumor", clin$type)
  tumor_samples <- clin$sampleid[clin$type == "Tumor"]

  # Filter by cancer type and tumor samples
  genome_data <- genome_data[genome_data$cancer_type %in% cancers &
    rownames(genome_data) %in% tumor_samples, ]

  # Remove cancer_type column
  genome_data$cancer_type <- NULL

  if (nrow(genome_data) == 0) {
    stop("No samples found for specified cancers", call. = FALSE)
  }

  # Remove genes with >90% NA
  na_rates <- colSums(is.na(genome_data)) / nrow(genome_data)
  valid_genes <- na_rates < 0.9
  genome_data <- genome_data[, valid_genes]

  message(sprintf(
    "  ✓ Loaded %d samples × %d genes (filtered <90%% NA)",
    nrow(genome_data), ncol(genome_data)
  ))

  # Convert to matrix
  genome_matrix <- as.matrix(genome_data)

  return(genome_matrix)
}


#' Perform genome-wide DEA using BioEnricher
#' @keywords internal
.stats_dea_genome <- function(var_data, genome_matrix, var_cancers) {
  # Match samples
  common_samples <- intersect(names(var_data), rownames(genome_matrix))

  if (length(common_samples) < 10) {
    stop("Too few overlapping samples (", length(common_samples), ")", call. = FALSE)
  }

  # Use only common samples
  var_values <- var_data[common_samples]
  genome_values <- genome_matrix[common_samples, ]

  # Drop unused factor levels
  var_values <- droplevels(var_values)

  # Get group levels
  group_levels <- levels(var_values)
  if (length(group_levels) == 0) {
    stop("No data available: all values are missing or invalid", call. = FALSE)
  }
  if (length(group_levels) == 1) {
    stop("DEA requires 2 groups, but all samples belong to the same group (",
      group_levels[1], "). Insufficient variation for analysis.",
      call. = FALSE
    )
  }
  if (length(group_levels) != 2) {
    stop("DEA requires exactly 2 groups, found: ", length(group_levels), call. = FALSE)
  }

  # Check sample size per group
  group_counts <- table(var_values)
  if (any(group_counts < 3)) {
    stop(
      sprintf(
        "Too few samples in some groups: %s. Need at least 3 per group.",
        paste(paste0(names(group_counts), "=", group_counts), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  message(sprintf(
    "  Performing DEA: %s vs %s (%d samples, %d genes)",
    group_levels[2], group_levels[1], length(common_samples), ncol(genome_values)
  ))

  # Transpose for BioEnricher::DEA (genes as rows, samples as columns)
  expr_matrix <- t(genome_values)

  # Prepare groups vector (must be factor with proper names)
  groups <- var_values
  names(groups) <- common_samples

  # Prepare contrasts: "Mutation-WildType"
  contrasts <- paste0(group_levels[2], "-", group_levels[1])

  # Call BioEnricher::DEA
  dea_result_list <- BioEnricher::DEA(
    expr = expr_matrix,
    groups = groups,
    contrasts = contrasts,
    Select.method = "limma",
    cutoff.P = 1, # Don't filter, we'll do it later
    cutoff.logFC = 0, # Don't filter by logFC
    verbose = FALSE
  )

  # Extract results data.frame from list
  dea_results <- dea_result_list$res

  # Standardize column names
  colnames(dea_results)[colnames(dea_results) == "Gene"] <- "gene"
  colnames(dea_results)[colnames(dea_results) == "P.Value"] <- "pvalue"
  colnames(dea_results)[colnames(dea_results) == "adj.P.Val"] <- "p_adjusted"

  rownames(dea_results) <- NULL

  message(sprintf("  ✓ DEA completed: %d significant genes (p<0.05)", sum(dea_results$pvalue < 0.05, na.rm = TRUE)))

  return(dea_results)
}


#' Perform genome-wide correlation
#' @keywords internal
.stats_cor_genome <- function(var_data, genome_matrix, var_cancers, method = "pearson") {
  # Match samples
  common_samples <- intersect(names(var_data), rownames(genome_matrix))

  if (length(common_samples) < 10) {
    stop("Too few overlapping samples (", length(common_samples), ")", call. = FALSE)
  }

  # Use only common samples
  var_values <- var_data[common_samples]
  genome_values <- genome_matrix[common_samples, ]

  message(sprintf(
    "  Calculating correlations (%d samples, %d genes)",
    length(common_samples), ncol(genome_values)
  ))

  # Prepare data for astat
  var_df <- data.frame(var = var_values)
  rownames(var_df) <- common_samples

  # Ensure genome_values has proper rownames
  rownames(genome_values) <- common_samples

  # Use astat for fast calculation
  cor_result <- astat::stat_cor(
    x = var_df,
    y = as.data.frame(genome_values),
    cor.method = method,
    use = "pairwise.complete.obs",
    p.adjust = FALSE
  )

  # Rename columns
  colnames(cor_result)[colnames(cor_result) == "y"] <- "gene"
  colnames(cor_result)[colnames(cor_result) == "p"] <- "pvalue"

  # P-value adjustment
  cor_result$p_adjusted <- p.adjust(cor_result$pvalue, method = "BH")

  message(sprintf(
    "  ✓ Correlation completed: %d significant genes (p<0.05)",
    sum(cor_result$pvalue < 0.05)
  ))

  return(cor_result)
}


#' Perform GSEA
#' @keywords internal
#' Perform GSEA using fgsea (following original implementation)
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

  # Check dependencies
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' required. Install: BiocManager::install('fgsea')", call. = FALSE)
  }
  if (!requireNamespace("BioEnricher", quietly = TRUE)) {
    stop("Package 'BioEnricher' required for gene sets", call. = FALSE)
  }

  # Load gene sets using BioEnricher::get.geneset2
  # 所有数据库都使用symbol（GO、KEGG、WIKI、Reactome均支持）
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

  # Run fgseaMultilevel
  fgsea_result <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranked_genes,
    minSize = minSize,
    maxSize = maxSize,
    nproc = n_workers
  )

  # Format results
  stats <- as.data.frame(fgsea_result)

  if (nrow(stats) > 0) {
    stats$Description <- pathway_descriptions[stats$pathway]
    colnames(stats)[colnames(stats) == "pathway"] <- "ID"
    colnames(stats)[colnames(stats) == "pval"] <- "pvalue"
    colnames(stats)[colnames(stats) == "padj"] <- "qvalue"

    # Convert leadingEdge list to character
    stats$leadingEdge <- sapply(stats$leadingEdge, function(x) {
      paste(unlist(x), collapse = ",")
    })

    # Reorder columns
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
        # Correct way to access cutpoint value
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
  # Check data availability
  complete_idx <- complete.cases(data[, c(group_col, time_col, event_col)])
  n_complete <- sum(complete_idx)

  if (n_complete == 0) {
    stop("No complete observations available for survival analysis. All data is missing.",
      call. = FALSE
    )
  }

  if (n_complete < 10) {
    stop(sprintf("Too few complete observations (%d) for survival analysis. Need at least 10.", n_complete),
      call. = FALSE
    )
  }

  # Use only complete cases
  data_complete <- data[complete_idx, ]

  # Create formula directly in survfit/survdiff calls to avoid symbol issue
  formula_str <- paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", group_col)

  survfit_obj <- survival::survfit(as.formula(formula_str), data = data_complete)
  survdiff_obj <- survival::survdiff(as.formula(formula_str), data = data_complete)

  # Fix survfit$call$formula to be actual formula, not symbol
  survfit_obj$call$formula <- as.formula(formula_str)

  # Extract p-value
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
  # Check data availability
  complete_idx <- complete.cases(data[, c(var_col, time_col, event_col)])
  n_complete <- sum(complete_idx)

  if (n_complete == 0) {
    stop("No complete observations available for Cox regression. All data is missing.",
      call. = FALSE
    )
  }

  if (n_complete < 10) {
    stop(sprintf("Too few complete observations (%d) for Cox regression. Need at least 10.", n_complete),
      call. = FALSE
    )
  }

  # Use only complete cases
  data_complete <- data[complete_idx, ]

  cox_formula <- as.formula(paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", var_col))

  cox_model <- survival::coxph(cox_formula, data = data_complete)
  cox_summary <- summary(cox_model)

  # Extract key statistics
  hr <- exp(coef(cox_model))
  ci <- exp(confint(cox_model))
  hr_lower <- ci[1]
  hr_upper <- ci[2]
  p_value <- cox_summary$coefficients[1, 5]
  cindex <- cox_summary$concordance[1]

  return(list(
    model = cox_model,
    summary = cox_summary,
    hr = hr,
    hr_lower = hr_lower,
    hr_upper = hr_upper,
    p_value = p_value,
    cindex = cindex
  ))
}


# ==============================================================================
# Helper Functions
# ==============================================================================

#' Extract column name from feature label
#'
#' @description
#' Converts "THBS2 (Protein, COAD)" → "COAD_THBS2_Protein"
#'
#' @keywords internal
.extract_colname_from_label <- function(labels, data) {
  # Parse labels to extract cancer_feature_modal format
  colnames_out <- character(length(labels))

  for (i in seq_along(labels)) {
    label <- labels[i]

    # Extract gene, modal, and cancer from label
    # Format: "GENE (Modal, Cancer)" or "SITE_GENE (Modal, Cancer)"
    parts <- strsplit(label, " \\(")[[1]]
    gene_or_site <- parts[1]

    modal_cancer <- gsub("\\)", "", parts[2])
    modal_cancer_parts <- strsplit(modal_cancer, ", ")[[1]]
    modal <- modal_cancer_parts[1]
    cancer <- modal_cancer_parts[2]

    # Construct column name with modal: Cancer_Gene_Modal
    colname <- paste0(cancer, "_", gene_or_site, "_", modal)

    # Check if exists in data
    if (colname %in% colnames(data)) {
      colnames_out[i] <- colname
    } else {
      # Try alternative formats
      possible_names <- grep(paste0("^", cancer, "_.*", gsub("_.*", "", gene_or_site)),
        colnames(data),
        value = TRUE
      )
      if (length(possible_names) > 0) {
        colnames_out[i] <- possible_names[1]
      } else {
        stop("Cannot find column for label: ", label, call. = FALSE)
      }
    }
  }

  return(colnames_out)
}


#' Map column name back to feature label
#'
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
