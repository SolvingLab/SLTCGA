# SLTCGA Test Suite - CNV Modality
# Setup
library(devtools)
load_all()
Sys.setenv(SL_BULK_DATA = "/Users/liuzaoqu/Desktop/develop/DataMiner_Dev/bulk_data")

# ==============================================================================
# CNV vs RNAseq
# ==============================================================================

# Case 1: MYC CNV vs MYC expression (BRCA)
res <- tcga_correlation(
  var1 = "MYC", var1_modal = "CNV", var1_cancers = "BRCA",
  var2 = "MYC", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 2: ERBB2 CNV vs ERBB2 expression (BRCA)
res <- tcga_correlation(
  var1 = "ERBB2", var1_modal = "CNV", var1_cancers = "BRCA",
  var2 = "ERBB2", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 3: EGFR CNV vs EGFR expression (LUAD)
res <- tcga_correlation(
  var1 = "EGFR", var1_modal = "CNV", var1_cancers = "LUAD",
  var2 = "EGFR", var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 4: CCND1 CNV vs CCND1 expression (HNSC)
res <- tcga_correlation(
  var1 = "CCND1", var1_modal = "CNV", var1_cancers = "HNSC",
  var2 = "CCND1", var2_modal = "RNAseq", var2_cancers = "HNSC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 5: MYC CNV vs downstream targets (OV)
res <- tcga_correlation(
  var1 = "MYC", var1_modal = "CNV", var1_cancers = "OV",
  var2 = c("CCND1", "CDK4", "E2F1"), var2_modal = "RNAseq", var2_cancers = "OV"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 6: Multiple CNVs vs genes (BRCA)
res <- tcga_correlation(
  var1 = c("MYC", "ERBB2", "CCND1"), var1_modal = "CNV", var1_cancers = "BRCA",
  var2 = c("MYC", "ERBB2", "CCND1"), var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# CNV vs Signature
# ==============================================================================

# Case 7: MYC CNV vs Purity (BRCA)
res <- tcga_correlation(
  var1 = "MYC", var1_modal = "CNV", var1_cancers = "BRCA",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 8: ERBB2 CNV vs TMB (BRCA)
res <- tcga_correlation(
  var1 = "ERBB2", var1_modal = "CNV", var1_cancers = "BRCA",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 9: EGFR CNV vs TIL_Score (LUAD)
res <- tcga_correlation(
  var1 = "EGFR", var1_modal = "CNV", var1_cancers = "LUAD",
  var2 = "TIL_Score", var2_modal = "Signature", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 10: CCND1 CNV vs Stemness (HNSC)
res <- tcga_correlation(
  var1 = "CCND1", var1_modal = "CNV", var1_cancers = "HNSC",
  var2 = "Stemness", var2_modal = "Signature", var2_cancers = "HNSC"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# CNV vs Mutation
# ==============================================================================

# Case 11: MYC CNV vs MYC mutation (BRCA)
res <- tcga_correlation(
  var1 = "MYC", var1_modal = "CNV", var1_cancers = "BRCA",
  var2 = "MYC", var2_modal = "Mutation", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 12: EGFR CNV vs KRAS mutation (LUAD)
res <- tcga_correlation(
  var1 = "EGFR", var1_modal = "CNV", var1_cancers = "LUAD",
  var2 = "KRAS", var2_modal = "Mutation", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 13: ERBB2 CNV vs PIK3CA mutation (BRCA)
res <- tcga_correlation(
  var1 = "ERBB2", var1_modal = "CNV", var1_cancers = "BRCA",
  var2 = "PIK3CA", var2_modal = "Mutation", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# CNV vs ImmuneCell
# ==============================================================================

# Case 14: MYC CNV vs T cells (BRCA)
res <- tcga_correlation(
  var1 = "MYC", var1_modal = "CNV", var1_cancers = "BRCA",
  var2 = c("T_cells_CD8", "T_cells_CD4_memory_resting"), var2_modal = "ImmuneCell",
  var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 15: ERBB2 CNV vs macrophages (BRCA)
res <- tcga_correlation(
  var1 = "ERBB2", var1_modal = "CNV", var1_cancers = "BRCA",
  var2 = c("Macrophages_M1", "Macrophages_M2"), var2_modal = "ImmuneCell",
  var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# CNV vs Clinical
# ==============================================================================

# Case 16: MYC CNV vs Age (BRCA)
res <- tcga_correlation(
  var1 = "MYC", var1_modal = "CNV", var1_cancers = "BRCA",
  var2 = "Age", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 17: EGFR CNV vs Stage (LUAD)
res <- tcga_correlation(
  var1 = "EGFR", var1_modal = "CNV", var1_cancers = "LUAD",
  var2 = "Stage", var2_modal = "Clinical", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 18: ERBB2 CNV vs Gender (BRCA)
res <- tcga_correlation(
  var1 = "ERBB2", var1_modal = "CNV", var1_cancers = "BRCA",
  var2 = "Gender", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# CNV vs CNV
# ==============================================================================

# Case 19: MYC vs CCND1 CNV (BRCA)
res <- tcga_correlation(
  var1 = "MYC", var1_modal = "CNV", var1_cancers = "BRCA",
  var2 = "CCND1", var2_modal = "CNV", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 20: ERBB2 vs PIK3CA CNV (BRCA)
res <- tcga_correlation(
  var1 = "ERBB2", var1_modal = "CNV", var1_cancers = "BRCA",
  var2 = "PIK3CA", var2_modal = "CNV", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 21: EGFR vs KRAS CNV (LUAD)
res <- tcga_correlation(
  var1 = "EGFR", var1_modal = "CNV", var1_cancers = "LUAD",
  var2 = "KRAS", var2_modal = "CNV", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# CNV Survival Analysis
# ==============================================================================

# Case 22: MYC CNV survival (BRCA)
res <- tcga_survival(
  var1 = "MYC", var1_modal = "CNV", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 23: ERBB2 CNV survival (BRCA)
res <- tcga_survival(
  var1 = "ERBB2", var1_modal = "CNV", var1_cancers = "BRCA",
  surv_type = "PFS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 24: Multiple CNVs survival (BRCA)
res <- tcga_survival(
  var1 = c("MYC", "ERBB2", "CCND1"), var1_modal = "CNV", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)


