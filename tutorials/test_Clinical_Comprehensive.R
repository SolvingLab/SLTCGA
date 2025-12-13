# SLTCGA Test Suite - Clinical Variables Comprehensive Testing
# Testing all 66+ clinical variables
# Setup
library(devtools)
load_all()
Sys.setenv(SL_BULK_DATA = "/Users/liuzaoqu/Desktop/develop/DataMiner_Dev/bulk_data")

# ==============================================================================
# SECTION 1: Demographic Variables
# ==============================================================================

# Case 1: Age vs TP53 expression (BRCA)
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 2: Gender vs ESR1 expression (BRCA)
res <- tcga_correlation(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "ESR1", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 3: Race vs TP53 expression (BRCA)
res <- tcga_correlation(
  var1 = "Race", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 4: VitalStatus vs gene expression (LUAD)
res <- tcga_correlation(
  var1 = "VitalStatus", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "KRAS", var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 5: Smoking vs KRAS mutation (LUAD)
res <- tcga_correlation(
  var1 = "Smoking", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "KRAS", var2_modal = "Mutation", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 6: Alcohol vs gene expression (LIHC)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LIHC",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "LIHC"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 2: Staging Variables
# ==============================================================================

# Case 7: Stage vs gene expression (LUAD)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "EGFR", var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 8: T stage vs gene expression (BRCA)
res <- tcga_correlation(
  var1 = "T", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 9: N stage vs gene expression (COAD)
res <- tcga_correlation(
  var1 = "N", var1_modal = "Clinical", var1_cancers = "COAD",
  var2 = "KRAS", var2_modal = "RNAseq", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 10: M stage vs gene expression (KIRC)
res <- tcga_correlation(
  var1 = "M", var1_modal = "Clinical", var1_cancers = "KIRC",
  var2 = "VHL", var2_modal = "RNAseq", var2_cancers = "KIRC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 11: StageNum (numeric) vs signature (LUAD)
res <- tcga_correlation(
  var1 = "StageNum", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 3: Tumor Characteristics
# ==============================================================================

# Case 12: Grade vs gene expression (BRCA)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "MYC", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 13: Type vs gene expression (BRCA)
res <- tcga_correlation(
  var1 = "Type", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "ESR1", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 14: ER status vs ESR1 expression (BRCA)
res <- tcga_correlation(
  var1 = "ER", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "ESR1", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 15: PR status vs PGR expression (BRCA)
res <- tcga_correlation(
  var1 = "PR", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "PGR", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 16: HER2 status vs ERBB2 expression (BRCA)
res <- tcga_correlation(
  var1 = "HER2", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "ERBB2", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 4: Clinical vs Signature
# ==============================================================================

# Case 17: Age vs TMB (BRCA)
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 18: Stage vs Purity (LUAD)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 19: Gender vs TIL_Score (KIRC)
res <- tcga_correlation(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "KIRC",
  var2 = "TIL_Score", var2_modal = "Signature", var2_cancers = "KIRC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 20: Stage vs Stemness (COAD)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "COAD",
  var2 = "Stemness", var2_modal = "Signature", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 21: Grade vs Proliferation (BRCA)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "Proliferation", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 22: ER status vs IFN_Gamma (BRCA)
res <- tcga_correlation(
  var1 = "ER", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "IFN_Gamma", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 5: Clinical vs Mutation
# ==============================================================================

# Case 23: Age vs TP53 mutation (BRCA)
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "Mutation", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 24: Stage vs KRAS mutation (LUAD)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "KRAS", var2_modal = "Mutation", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 25: Gender vs PIK3CA mutation (UCEC)
res <- tcga_correlation(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "UCEC",
  var2 = "PIK3CA", var2_modal = "Mutation", var2_cancers = "UCEC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 26: Grade vs TP53 mutation (GBM)
res <- tcga_correlation(
  var1 = "Grade", var1_modal = "Clinical", var1_cancers = "GBM",
  var2 = "TP53", var2_modal = "Mutation", var2_cancers = "GBM"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 27: Multiple clinical vs mutations (LUAD)
res <- tcga_correlation(
  var1 = c("Age", "Gender", "Stage"), var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = c("KRAS", "EGFR", "TP53"), var2_modal = "Mutation", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 6: Clinical vs ImmuneCell
# ==============================================================================

# Case 28: Age vs T cells (BRCA)
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = c("T_cells_CD8", "T_cells_CD4_memory_resting"), var2_modal = "ImmuneCell",
  var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 29: Stage vs ALL_IMMUNE_CELLS (LUAD, heatmap)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "cibersort",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 30: Gender vs macrophages (KIRC)
res <- tcga_correlation(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "KIRC",
  var2 = c("Macrophages_M1", "Macrophages_M2"), var2_modal = "ImmuneCell",
  var2_cancers = "KIRC",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 31: Grade vs immune cells (BRCA)
res <- tcga_correlation(
  var1 = "Grade", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = c("T_cells_CD8", "Macrophages_M1", "NK_cells_activated"),
  var2_modal = "ImmuneCell", var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 32: ER status vs immune cells (BRCA, heatmap)
res <- tcga_correlation(
  var1 = "ER", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell", var2_cancers = "BRCA",
  immune_algorithm = "timer",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 7: Clinical vs Clinical
# ==============================================================================

# Case 33: Age vs Gender (BRCA)
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "Gender", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 34: Stage vs Age (LUAD)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "Age", var2_modal = "Clinical", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 35: Grade vs Stage (BRCA)
res <- tcga_correlation(
  var1 = "Grade", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "Stage", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 36: Gender vs Race (HNSC)
res <- tcga_correlation(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "HNSC",
  var2 = "Race", var2_modal = "Clinical", var2_cancers = "HNSC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 37: ER vs PR status (BRCA)
res <- tcga_correlation(
  var1 = "ER", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "PR", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 38: ER vs HER2 status (BRCA)
res <- tcga_correlation(
  var1 = "ER", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "HER2", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 8: Clinical vs CNV
# ==============================================================================

# Case 39: Stage vs MYC CNV (BRCA)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "MYC", var2_modal = "CNV", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 40: Grade vs ERBB2 CNV (BRCA)
res <- tcga_correlation(
  var1 = "Grade", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "ERBB2", var2_modal = "CNV", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 41: Age vs EGFR CNV (LUAD)
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "EGFR", var2_modal = "CNV", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 9: Clinical vs Methylation
# ==============================================================================

# Case 42: Stage vs TP53 methylation (BRCA)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "Methylation", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 43: Age vs BRCA1 methylation (BRCA)
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "BRCA1", var2_modal = "Methylation", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 44: Gender vs MLH1 methylation (COAD)
res <- tcga_correlation(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "COAD",
  var2 = "MLH1", var2_modal = "Methylation", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 10: Clinical vs miRNA
# ==============================================================================

# Case 45: Age vs hsa-let-7a-1 (BRCA)
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "hsa-let-7a-1", var2_modal = "miRNA", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 46: Stage vs hsa-mir-21 (LUAD)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "hsa-mir-21", var2_modal = "miRNA", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 47: Gender vs multiple miRNAs (BRCA)
res <- tcga_correlation(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = c("hsa-let-7a-1", "hsa-mir-21"), var2_modal = "miRNA", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 11: Multiple Clinical Variables
# ==============================================================================

# Case 48: Multiple demographics vs gene (BRCA)
res <- tcga_correlation(
  var1 = c("Age", "Gender", "Race"), var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 49: Multiple stage variables vs gene (LUAD)
res <- tcga_correlation(
  var1 = c("Stage", "T", "N"), var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "EGFR", var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 50: Multiple receptor status vs genes (BRCA)
res <- tcga_correlation(
  var1 = c("ER", "PR", "HER2"), var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = c("ESR1", "PGR", "ERBB2"), var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 51: Multiple clinical vs signatures (LUAD)
res <- tcga_correlation(
  var1 = c("Age", "Gender", "Stage"), var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = c("TMB", "TIL_Score", "Purity"), var2_modal = "Signature", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 12: Clinical Survival Analysis
# ==============================================================================

# Case 52: Age survival (BRCA)
res <- tcga_survival(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 53: Stage survival (LUAD) - auto forest plot
res <- tcga_survival(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 54: Gender survival (KIRC)
res <- tcga_survival(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "KIRC",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 55: Grade survival (BRCA)
res <- tcga_survival(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 56: Race survival (BRCA) - auto forest plot
res <- tcga_survival(
  var1 = "Race", var1_modal = "Clinical", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 57: Multiple clinical survival (BRCA) - forest plot
res <- tcga_survival(
  var1 = c("Age", "Gender", "Race"), var1_modal = "Clinical", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 58: Stage + Grade survival (BRCA)
res <- tcga_survival(
  var1 = c("Stage", "Age"), var1_modal = "Clinical", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 59: ER + PR + HER2 survival (BRCA)
res <- tcga_survival(
  var1 = c("ER", "PR", "HER2"), var1_modal = "Clinical", var1_cancers = "BRCA",
  surv_type = "PFS"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 13: Multi-Cancer Clinical Analysis
# ==============================================================================

# Case 60: Age across pan-cancer
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical",
  var1_cancers = c("BRCA", "LUAD", "COAD"),
  var2 = "TP53", var2_modal = "RNAseq",
  var2_cancers = c("BRCA", "LUAD", "COAD")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 61: Stage pan-cancer
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical",
  var1_cancers = c("LUAD", "COAD", "KIRC"),
  var2 = "TMB", var2_modal = "Signature",
  var2_cancers = c("LUAD", "COAD", "KIRC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 62: Gender pan-cancer
res <- tcga_correlation(
  var1 = "Gender", var1_modal = "Clinical",
  var1_cancers = c("BRCA", "LUAD", "HNSC"),
  var2 = "TIL_Score", var2_modal = "Signature",
  var2_cancers = c("BRCA", "LUAD", "HNSC")
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 14: Rare Clinical Variables
# ==============================================================================

# Case 63: Smoking vs KRAS (LUAD)
res <- tcga_correlation(
  var1 = "Smoking", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "KRAS", var2_modal = "Mutation", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 64: Alcohol vs gene (LIHC)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LIHC",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "LIHC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 65: VitalStatus vs signature (BRCA)
res <- tcga_correlation(
  var1 = "VitalStatus", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 15: Edge Cases - Clinical Variables
# ==============================================================================

# Case 66: All-NA variable test (BMI in some cancers)
res <- tcga_correlation(
  var1 = c("Age", "Gender", "Race", "BMI"), var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 67: Stage in different cancers
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "HNSC",
  var2 = "CD274", var2_modal = "RNAseq", var2_cancers = "HNSC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 68: Grade in GBM
res <- tcga_correlation(
  var1 = "Age", var1_modal = "Clinical", var1_cancers = "GBM",
  var2 = "IDH1", var2_modal = "Mutation", var2_cancers = "GBM"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 69: T stage numeric
res <- tcga_correlation(
  var1 = "TNum", var1_modal = "Clinical", var1_cancers = "BRCA",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 70: N stage numeric
res <- tcga_correlation(
  var1 = "NNum", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 71: M stage numeric
res <- tcga_correlation(
  var1 = "MNum", var1_modal = "Clinical", var1_cancers = "COAD",
  var2 = "TIL_Score", var2_modal = "Signature", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 72: StageNum vs multiple genes
res <- tcga_correlation(
  var1 = "StageNum", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = c("EGFR", "KRAS", "TP53"), var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)


