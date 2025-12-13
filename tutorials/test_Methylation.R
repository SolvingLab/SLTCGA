# SLTCGA Test Suite - Methylation Modality
# Setup
library(devtools)
load_all()
Sys.setenv(SL_BULK_DATA = "/Users/liuzaoqu/Desktop/develop/DataMiner_Dev/bulk_data")

# ==============================================================================
# Methylation vs RNAseq
# ==============================================================================

# Case 1: TP53 methylation vs TP53 expression (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Methylation", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "BRCA",
  methylation_region = "Promoter_mean"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 2: BRCA1 methylation vs BRCA1 expression (OV)
res <- tcga_correlation(
  var1 = "BRCA1", var1_modal = "Methylation", var1_cancers = "OV",
  var2 = "BRCA1", var2_modal = "RNAseq", var2_cancers = "OV"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 3: MLH1 methylation vs MLH1 expression (COAD)
res <- tcga_correlation(
  var1 = "MLH1", var1_modal = "Methylation", var1_cancers = "COAD",
  var2 = "MLH1", var2_modal = "RNAseq", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 4: CDKN2A methylation vs CDKN2A expression (LUAD)
res <- tcga_correlation(
  var1 = "CDKN2A", var1_modal = "Methylation", var1_cancers = "LUAD",
  var2 = "CDKN2A", var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 5: ESR1 methylation vs ESR1 expression (BRCA)
res <- tcga_correlation(
  var1 = "ESR1", var1_modal = "Methylation", var1_cancers = "BRCA",
  var2 = "ESR1", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 6: Multiple methylation vs genes (BRCA)
res <- tcga_correlation(
  var1 = c("TP53", "ESR1", "PGR"), var1_modal = "Methylation", var1_cancers = "BRCA",
  var2 = c("TP53", "ESR1", "PGR"), var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Methylation vs Signature
# ==============================================================================

# Case 7: TP53 methylation vs TMB (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Methylation", var1_cancers = "BRCA",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 8: MLH1 methylation vs MSI (COAD)
res <- tcga_correlation(
  var1 = "MLH1", var1_modal = "Methylation", var1_cancers = "COAD",
  var2 = "MSI", var2_modal = "Signature", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 9: BRCA1 methylation vs HRD (OV)
res <- tcga_correlation(
  var1 = "BRCA1", var1_modal = "Methylation", var1_cancers = "OV",
  var2 = "HRD", var2_modal = "Signature", var2_cancers = "OV"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 10: ESR1 methylation vs Purity (BRCA)
res <- tcga_correlation(
  var1 = "ESR1", var1_modal = "Methylation", var1_cancers = "BRCA",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Methylation vs Mutation
# ==============================================================================

# Case 11: TP53 methylation vs TP53 mutation (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Methylation", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "Mutation", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 12: MLH1 methylation vs MLH1 mutation (COAD)
res <- tcga_correlation(
  var1 = "MLH1", var1_modal = "Methylation", var1_cancers = "COAD",
  var2 = "MLH1", var2_modal = "Mutation", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 13: BRCA1 methylation vs BRCA1 mutation (OV)
res <- tcga_correlation(
  var1 = "THBS2", var1_modal = "Methylation", var1_cancers = "STAD",
  var2 = "THBS2", var2_modal = "Mutation", var2_cancers = "STAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Methylation vs ImmuneCell
# ==============================================================================

# Case 14: TP53 methylation vs T cells (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Methylation", var1_cancers = "BRCA",
  var2 = c("T_cells_CD8", "T_cells_CD4_memory_resting"), var2_modal = "ImmuneCell",
  var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 15: CDKN2A methylation vs macrophages (LUAD)
res <- tcga_correlation(
  var1 = "CDKN2A", var1_modal = "Methylation", var1_cancers = "LUAD",
  var2 = c("Macrophages_M1", "Macrophages_M2"), var2_modal = "ImmuneCell",
  var2_cancers = "LUAD",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Methylation vs Clinical
# ==============================================================================

# Case 16: TP53 methylation vs Age (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Methylation", var1_cancers = "BRCA",
  var2 = "Age", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 17: MLH1 methylation vs Stage (COAD)
res <- tcga_correlation(
  var1 = "MLH1", var1_modal = "Methylation", var1_cancers = "COAD",
  var2 = "Stage", var2_modal = "Clinical", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 18: ESR1 methylation vs Gender (BRCA)
res <- tcga_correlation(
  var1 = "ESR1", var1_modal = "Methylation", var1_cancers = "BRCA",
  var2 = "Gender", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Methylation vs Methylation
# ==============================================================================

# Case 19: TP53 vs ESR1 methylation (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Methylation", var1_cancers = "BRCA",
  var2 = "ESR1", var2_modal = "Methylation", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 20: MLH1 vs MSH2 methylation (COAD)
res <- tcga_correlation(
  var1 = "MLH1", var1_modal = "Methylation", var1_cancers = "COAD",
  var2 = "MSH2", var2_modal = "Methylation", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Methylation Survival Analysis
# ==============================================================================

# Case 21: TP53 methylation survival (BRCA)
res <- tcga_survival(
  var1 = "TP53", var1_modal = "Methylation", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 22: MLH1 methylation survival (COAD)
res <- tcga_survival(
  var1 = "MLH1", var1_modal = "Methylation", var1_cancers = "COAD",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 23: Multiple methylation survival (BRCA)
res <- tcga_survival(
  var1 = c("TP53", "ESR1", "BRCA1"), var1_modal = "Methylation", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)


