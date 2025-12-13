# SLTCGA Test Suite - Signature Modality
# Setup
library(devtools)
load_all()
Sys.setenv(SL_BULK_DATA = "/Users/liuzaoqu/Desktop/develop/DataMiner_Dev/bulk_data")

# ==============================================================================
# Signature vs RNAseq
# ==============================================================================

# Case 1: TMB vs TP53 (BRCA)
res <- tcga_correlation(
  var1 = "TMB", var1_modal = "Signature", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 2: Purity vs ESR1 (BRCA)
res <- tcga_correlation(
  var1 = "Purity", var1_modal = "Signature", var1_cancers = "BRCA",
  var2 = "ESR1", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 3: TIL_Score vs CD274 (LUAD)
res <- tcga_correlation(
  var1 = "TIL_Score", var1_modal = "Signature", var1_cancers = "LUAD",
  var2 = "CD274", var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 4: Stemness vs MYC (COAD)
res <- tcga_correlation(
  var1 = "Stemness", var1_modal = "Signature", var1_cancers = "COAD",
  var2 = "MYC", var2_modal = "RNAseq", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 5: HRD vs BRCA1 (OV)
res <- tcga_correlation(
  var1 = "HRD", var1_modal = "Signature", var1_cancers = "OV",
  var2 = "BRCA1", var2_modal = "RNAseq", var2_cancers = "OV"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 6: Multiple signatures vs gene (LUAD)
res <- tcga_correlation(
  var1 = c("TMB", "TIL_Score", "Purity"), var1_modal = "Signature", var1_cancers = "LUAD",
  var2 = "KRAS", var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 7: Immune signatures vs CD274 (BRCA)
res <- tcga_correlation(
  var1 = c("Leukocyte", "Stromal", "TIL_Score"), var1_modal = "Signature", var1_cancers = "BRCA",
  var2 = "CD274", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Signature vs Signature
# ==============================================================================

# Case 8: TMB vs Purity (BRCA)
res <- tcga_correlation(
  var1 = "TMB", var1_modal = "Signature", var1_cancers = "BRCA",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 9: TIL_Score vs Leukocyte (LUAD)
res <- tcga_correlation(
  var1 = "TIL_Score", var1_modal = "Signature", var1_cancers = "LUAD",
  var2 = "Leukocyte", var2_modal = "Signature", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 10: Stemness vs Ploidy (COAD)
res <- tcga_correlation(
  var1 = "Stemness", var1_modal = "Signature", var1_cancers = "COAD",
  var2 = "Ploidy", var2_modal = "Signature", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 11: HRD vs Aneuploidy (OV)
res <- tcga_correlation(
  var1 = "HRD", var1_modal = "Signature", var1_cancers = "OV",
  var2 = "Aneuploidy", var2_modal = "Signature", var2_cancers = "OV"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 12: Multiple vs multiple signatures (BRCA)
res <- tcga_correlation(
  var1 = c("TMB", "Purity"), var1_modal = "Signature", var1_cancers = "BRCA",
  var2 = c("TIL_Score", "Leukocyte"), var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Signature vs Mutation
# ==============================================================================

# Case 13: TMB vs TP53 mutation (BRCA)
res <- tcga_correlation(
  var1 = "TMB", var1_modal = "Signature", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "Mutation", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 14: TIL_Score vs KRAS mutation (LUAD)
res <- tcga_correlation(
  var1 = "TIL_Score", var1_modal = "Signature", var1_cancers = "LUAD",
  var2 = "KRAS", var2_modal = "Mutation", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 15: HRD vs BRCA1 mutation (OV)
res <- tcga_correlation(
  var1 = "HRD", var1_modal = "Signature", var1_cancers = "OV",
  var2 = "BRCA1", var2_modal = "Mutation", var2_cancers = "OV"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 16: MSI vs MLH1 mutation (COAD)
res <- tcga_correlation(
  var1 = "MSI", var1_modal = "Signature", var1_cancers = "COAD",
  var2 = "MLH1", var2_modal = "Mutation", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Signature vs ImmuneCell
# ==============================================================================

# Case 17: TMB vs immune cells (LUAD)
res <- tcga_correlation(
  var1 = "TMB", var1_modal = "Signature", var1_cancers = "LUAD",
  var2 = c("T_cells_CD8", "T_cells_CD4_memory_resting", "Macrophages_M1"),
  var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 18: Purity vs immune cells (BRCA)
res <- tcga_correlation(
  var1 = "Purity", var1_modal = "Signature", var1_cancers = "BRCA",
  var2 = c("T_cells_CD8", "B_cells_naive", "NK_cells_activated"),
  var2_modal = "ImmuneCell", var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 19: TIL_Score vs macrophages (SKCM)
res <- tcga_correlation(
  var1 = "TIL_Score", var1_modal = "Signature", var1_cancers = "SKCM",
  var2 = c("Macrophages_M1", "Macrophages_M2"), var2_modal = "ImmuneCell",
  var2_cancers = "SKCM",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Signature vs Clinical
# ==============================================================================

# Case 20: TMB vs Age (BRCA)
res <- tcga_correlation(
  var1 = "TMB", var1_modal = "Signature", var1_cancers = "BRCA",
  var2 = "Age", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 21: Purity vs Stage (LUAD)
res <- tcga_correlation(
  var1 = "Purity", var1_modal = "Signature", var1_cancers = "LUAD",
  var2 = "Stage", var2_modal = "Clinical", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 22: TIL_Score vs Gender (KIRC)
res <- tcga_correlation(
  var1 = "TIL_Score", var1_modal = "Signature", var1_cancers = "KIRC",
  var2 = "Gender", var2_modal = "Clinical", var2_cancers = "KIRC"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Signature Survival Analysis
# ==============================================================================

# Case 23: TMB survival (BRCA)
res <- tcga_survival(
  var1 = "TMB", var1_modal = "Signature", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 24: TIL_Score survival (LUAD)
res <- tcga_survival(
  var1 = "TIL_Score", var1_modal = "Signature", var1_cancers = "LUAD",
  surv_type = "PFS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 25: Multiple signatures survival (BRCA)
res <- tcga_survival(
  var1 = c("TMB", "TIL_Score", "IFN_Gamma"), var1_modal = "Signature",
  var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 26: Purity survival (OV)
res <- tcga_survival(
  var1 = "Purity", var1_modal = "Signature", var1_cancers = "OV",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 27: Stemness survival (COAD)
res <- tcga_survival(
  var1 = "Stemness", var1_modal = "Signature", var1_cancers = "COAD",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Signature Subtype Analysis
# ==============================================================================

# Case 28: TMB across BRCA subtypes
res <- tcga_correlation(
  var1 = "TMB", var1_modal = "Signature",
  var1_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC"),
  var2 = "Purity", var2_modal = "Signature",
  var2_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 29: TIL_Score across COAD subtypes
res <- tcga_correlation(
  var1 = "TIL_Score", var1_modal = "Signature",
  var1_cancers = c("COAD_LCC", "COAD_RCC"),
  var2 = "Stemness", var2_modal = "Signature",
  var2_cancers = c("COAD_LCC", "COAD_RCC")
)
res$plot
head(res$stats)
head(res$raw_data)


