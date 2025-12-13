# SLTCGA Test Suite - Clinical Modality
# Setup
library(devtools)
load_all()
Sys.setenv(SL_BULK_DATA = "/Users/liuzaoqu/Desktop/develop/DataMiner_Dev/bulk_data")

# ==============================================================================
# Clinical vs RNAseq
# ==============================================================================

# Case 1: Age vs TP53 (BRCA)
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 2: Gender vs ESR1 (BRCA)
res <- tcga_correlation(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "ESR1", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 3: Stage vs EGFR (LUAD)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "EGFR", var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 4: Race vs TP53 (BRCA)
res <- tcga_correlation(
  var1 = "Race", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 5: Multiple clinical vs gene (BRCA)
res <- tcga_correlation(
  var1 = c("Age", "Gender", "Race"), var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 6: Stage vs multiple genes (LUAD)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = c("EGFR", "KRAS", "ALK"), var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Clinical vs Signature
# ==============================================================================

# Case 7: Age vs TMB (BRCA)
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 8: Stage vs Purity (LUAD)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 9: Gender vs TIL_Score (BRCA)
res <- tcga_correlation(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TIL_Score", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 10: Stage vs Stemness (COAD)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "COAD",
  var2 = "Stemness", var2_modal = "Signature", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Clinical vs ImmuneCell
# ==============================================================================

# Case 11: Stage vs immune cells (LUAD, heatmap)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "cibersort",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 12: Age vs T cells (BRCA)
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = c("T_cells_CD8", "T_cells_CD4_memory_resting"), var2_modal = "ImmuneCell",
  var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 13: Gender vs macrophages (KIRC)
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
# Clinical vs Mutation
# ==============================================================================

# Case 14: Age vs TP53 mutation (BRCA)
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "Mutation", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 15: Stage vs KRAS mutation (LUAD)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "KRAS", var2_modal = "Mutation", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 16: Gender vs PIK3CA mutation (UCEC)
res <- tcga_correlation(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "UCEC",
  var2 = "PIK3CA", var2_modal = "Mutation", var2_cancers = "UCEC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 17: Stage vs multiple mutations (LUAD)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = c("KRAS", "EGFR", "TP53"), var2_modal = "Mutation", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Clinical vs Clinical
# ==============================================================================

# Case 18: Age vs Gender (BRCA)
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "Gender", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 19: Stage vs Age (LUAD)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "Age", var2_modal = "Clinical", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 20: Gender vs Race (HNSC)
res <- tcga_correlation(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "HNSC",
  var2 = "Race", var2_modal = "Clinical", var2_cancers = "HNSC"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Clinical Survival Analysis
# ==============================================================================

# Case 21: Age survival (BRCA)
res <- tcga_survival(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 22: Stage survival (LUAD)
res <- tcga_survival(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 23: Gender survival (KIRC)
res <- tcga_survival(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "KIRC",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 24: Multiple clinical survival (BRCA)
res <- tcga_survival(
  var1 = c("Age", "Gender", "Race"), var1_modal = "Clinical", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)


