# SLTCGA Test Suite - ImmuneCell Modality
# Setup
library(devtools)
load_all()
Sys.setenv(SL_BULK_DATA = "/Users/liuzaoqu/Desktop/develop/DataMiner_Dev/bulk_data")

# ==============================================================================
# ImmuneCell vs RNAseq
# ==============================================================================

# Case 1: T_cells_CD8 vs CD274 (LUAD)
res <- tcga_correlation(
  var1 = "T_cells_CD8", var1_modal = "ImmuneCell", var1_cancers = "LUAD",
  var2 = "CD274", var2_modal = "RNAseq", var2_cancers = "LUAD",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 2: Macrophages_M1 vs IL6 (BRCA)
res <- tcga_correlation(
  var1 = "Macrophages_M1", var1_modal = "ImmuneCell", var1_cancers = "BRCA",
  var2 = "IL6", var2_modal = "RNAseq", var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 3: NK_cells vs IFNG (SKCM)
res <- tcga_correlation(
  var1 = "NK_cells_activated", var1_modal = "ImmuneCell", var1_cancers = "SKCM",
  var2 = "IFNG", var2_modal = "RNAseq", var2_cancers = "SKCM",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 4: Multiple immune cells vs gene (KIRC)
res <- tcga_correlation(
  var1 = c("T_cells_CD8", "T_cells_CD4_memory_resting", "Macrophages_M1"),
  var1_modal = "ImmuneCell", var1_cancers = "KIRC",
  var2 = "PDCD1", var2_modal = "RNAseq", var2_cancers = "KIRC",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 5: ALL_IMMUNE_CELLS vs CD274 (LUAD)
res <- tcga_correlation(
  var1 = "ALL_IMMUNE_CELLS", var1_modal = "ImmuneCell", var1_cancers = "LUAD",
  var2 = "CD274", var2_modal = "RNAseq", var2_cancers = "LUAD",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# ImmuneCell vs Mutation
# ==============================================================================

# Case 6: T_cells_CD8 vs TP53 mutation (LUAD, boxplot)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "LUAD",
  var2 = "T_cells_CD8", var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 7: ALL_IMMUNE_CELLS vs TP53 mutation (LUAD, heatmap)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "LUAD",
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "cibersort",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 8: ALL_IMMUNE_CELLS vs KRAS mutation (LUAD, heatmap)
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "cibersort",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 9: Immune cells vs PIK3CA mutation (BRCA)
res <- tcga_correlation(
  var1 = "PIK3CA", var1_modal = "Mutation", var1_cancers = "BRCA",
  var2 = c("T_cells_CD8", "Macrophages_M1", "NK_cells_activated"),
  var2_modal = "ImmuneCell", var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 10: Immune cells vs BRAF mutation (SKCM)
res <- tcga_correlation(
  var1 = "BRAF", var1_modal = "Mutation", var1_cancers = "SKCM",
  var2 = c("T_cells_CD8", "T_cells_CD4_memory_resting"),
  var2_modal = "ImmuneCell", var2_cancers = "SKCM",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# ImmuneCell vs Clinical
# ==============================================================================

# Case 11: T_cells_CD8 vs Age (BRCA)
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "T_cells_CD8", var2_modal = "ImmuneCell", var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 12: Immune cells vs Stage (LUAD, heatmap)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = c("T_cells_CD8", "T_cells_CD4_memory_resting", "Macrophages_M1",
           "Macrophages_M2", "B_cells_naive", "NK_cells_activated",
           "DC_resting", "Monocytes", "Plasma_cells", "Mast_cells_resting"),
  var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "cibersort",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 13: Macrophages vs Gender (KIRC)
res <- tcga_correlation(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "KIRC",
  var2 = c("Macrophages_M1", "Macrophages_M2"), var2_modal = "ImmuneCell",
  var2_cancers = "KIRC",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# ImmuneCell vs Signature
# ==============================================================================

# Case 14: T_cells_CD8 vs TMB (LUAD)
res <- tcga_correlation(
  var1 = "T_cells_CD8", var1_modal = "ImmuneCell", var1_cancers = "LUAD",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "LUAD",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 15: Multiple immune cells vs TIL_Score (BRCA)
res <- tcga_correlation(
  var1 = c("T_cells_CD8", "Macrophages_M1", "NK_cells_activated"),
  var1_modal = "ImmuneCell", var1_cancers = "BRCA",
  var2 = "TIL_Score", var2_modal = "Signature", var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 16: Macrophages vs Purity (LIHC)
res <- tcga_correlation(
  var1 = c("Macrophages_M1", "Macrophages_M2"), var1_modal = "ImmuneCell",
  var1_cancers = "LIHC",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "LIHC",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# ImmuneCell vs ImmuneCell
# ==============================================================================

# Case 17: T_cells_CD8 vs Macrophages_M1 (LUAD)
res <- tcga_correlation(
  var1 = "T_cells_CD8", var1_modal = "ImmuneCell", var1_cancers = "LUAD",
  var2 = "Macrophages_M1", var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 18: Multiple T cells vs multiple macrophages (BRCA)
res <- tcga_correlation(
  var1 = c("T_cells_CD8", "T_cells_CD4_memory_resting"),
  var1_modal = "ImmuneCell", var1_cancers = "BRCA",
  var2 = c("Macrophages_M1", "Macrophages_M2"),
  var2_modal = "ImmuneCell", var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 19: B cells vs T cells (SKCM)
res <- tcga_correlation(
  var1 = "B_cells_naive", var1_modal = "ImmuneCell", var1_cancers = "SKCM",
  var2 = "T_cells_CD8", var2_modal = "ImmuneCell", var2_cancers = "SKCM",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# ImmuneCell Survival Analysis
# ==============================================================================

# Case 20: T_cells_CD8 survival (LUAD)
res <- tcga_survival(
  var1 = "T_cells_CD8", var1_modal = "ImmuneCell", var1_cancers = "LUAD",
  surv_type = "OS",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 21: Multiple immune cells survival (BRCA)
res <- tcga_survival(
  var1 = c("T_cells_CD8", "Macrophages_M1", "NK_cells_activated"),
  var1_modal = "ImmuneCell", var1_cancers = "BRCA",
  surv_type = "OS",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 22: Macrophages_M1 survival (LIHC)
res <- tcga_survival(
  var1 = "Macrophages_M1", var1_modal = "ImmuneCell", var1_cancers = "LIHC",
  surv_type = "PFS",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)


