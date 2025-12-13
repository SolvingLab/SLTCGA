# SLTCGA Test Suite - RNAseq Modality
# Setup
library(devtools)
load_all()
Sys.setenv(SL_BULK_DATA = "/Users/liuzaoqu/Desktop/develop/DataMiner_Dev/bulk_data")

# ==============================================================================
# RNAseq vs Signature
# ==============================================================================

# Case 1: TP53 vs TMB (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 2: ESR1 vs Purity (BRCA)
res <- tcga_correlation(
  var1 = "ESR1", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 3: KRAS vs Stemness (LUAD)
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "Stemness", var2_modal = "Signature", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 4: BRCA1 vs HRD (OV)
res <- tcga_correlation(
  var1 = "BRCA1", var1_modal = "RNAseq", var1_cancers = "OV",
  var2 = "HRD", var2_modal = "Signature", var2_cancers = "OV"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 5: MYC vs Ploidy (COAD)
res <- tcga_correlation(
  var1 = "MYC", var1_modal = "RNAseq", var1_cancers = "COAD",
  var2 = "Ploidy", var2_modal = "Signature", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 6: CD274 vs TIL_Score (LUAD)
res <- tcga_correlation(
  var1 = "CD274", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "TIL_Score", var2_modal = "Signature", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 7: PDCD1 vs Leukocyte (SKCM)
res <- tcga_correlation(
  var1 = "PDCD1", var1_modal = "RNAseq", var1_cancers = "SKCM",
  var2 = "Leukocyte", var2_modal = "Signature", var2_cancers = "SKCM"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 8: CTLA4 vs Stromal (KIRC)
res <- tcga_correlation(
  var1 = "CTLA4", var1_modal = "RNAseq", var1_cancers = "KIRC",
  var2 = "Stromal", var2_modal = "Signature", var2_cancers = "KIRC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 9: IFNG vs IFN_Gamma (HNSC)
res <- tcga_correlation(
  var1 = "IFNG", var1_modal = "RNAseq", var1_cancers = "HNSC",
  var2 = "IFN_Gamma", var2_modal = "Signature", var2_cancers = "HNSC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 10: TGFB1 vs TGF_Beta (LIHC)
res <- tcga_correlation(
  var1 = "TGFB1", var1_modal = "RNAseq", var1_cancers = "LIHC",
  var2 = "TGF_Beta", var2_modal = "Signature", var2_cancers = "LIHC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 11: EGFR vs TMB (LUAD)
res <- tcga_correlation(
  var1 = "EGFR", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 12: PIK3CA vs Purity (UCEC)
res <- tcga_correlation(
  var1 = "PIK3CA", var1_modal = "RNAseq", var1_cancers = "UCEC",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "UCEC"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# RNAseq vs RNAseq
# ==============================================================================

# Case 13: TP53 vs MDM2 (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "MDM2", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 14: EGFR vs KRAS (LUAD)
res <- tcga_correlation(
  var1 = "EGFR", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "KRAS", var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 15: MYC vs CCND1 (HNSC)
res <- tcga_correlation(
  var1 = "MYC", var1_modal = "RNAseq", var1_cancers = "HNSC",
  var2 = "CCND1", var2_modal = "RNAseq", var2_cancers = "HNSC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 16: PIK3CA vs AKT1 (UCEC)
res <- tcga_correlation(
  var1 = "PIK3CA", var1_modal = "RNAseq", var1_cancers = "UCEC",
  var2 = "AKT1", var2_modal = "RNAseq", var2_cancers = "UCEC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 17: BRCA1 vs BRCA2 (OV)
res <- tcga_correlation(
  var1 = "BRCA1", var1_modal = "RNAseq", var1_cancers = "OV",
  var2 = "BRCA2", var2_modal = "RNAseq", var2_cancers = "OV"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 18: Multiple genes (BRCA)
res <- tcga_correlation(
  var1 = c("TP53", "ESR1", "ERBB2"), var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = c("PGR", "AR", "GATA3"), var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 19: Multiple genes (LUAD)
res <- tcga_correlation(
  var1 = c("EGFR", "KRAS", "ALK"), var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = c("STK11", "KEAP1", "NF1"), var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 20: Multi-cancer comparison
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = c("BRCA", "LUAD", "COAD"),
  var2 = "MYC", var2_modal = "RNAseq", var2_cancers = c("BRCA", "LUAD", "COAD")
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# RNAseq vs Methylation
# ==============================================================================

# Case 21: TP53 expression vs TP53 methylation (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "Methylation", var2_cancers = "BRCA",
  methylation_region = "Promoter_mean"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 22: BRCA1 vs BRCA1 methylation (OV)
res <- tcga_correlation(
  var1 = "BRCA1", var1_modal = "RNAseq", var1_cancers = "OV",
  var2 = "BRCA1", var2_modal = "Methylation", var2_cancers = "OV"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 23: MLH1 vs MLH1 methylation (COAD)
res <- tcga_correlation(
  var1 = "MLH1", var1_modal = "RNAseq", var1_cancers = "COAD",
  var2 = "MLH1", var2_modal = "Methylation", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 24: CDKN2A vs CDKN2A methylation (LUAD)
res <- tcga_correlation(
  var1 = "CDKN2A", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "CDKN2A", var2_modal = "Methylation", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# RNAseq vs miRNA
# ==============================================================================

# Case 25: MYC vs hsa-let-7a-1 (BRCA)
res <- tcga_correlation(
  var1 = "MYC", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "hsa-let-7a-1", var2_modal = "miRNA", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 26: KRAS vs hsa-mir-21 (LUAD)
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "hsa-mir-21", var2_modal = "miRNA", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 27: Multiple genes vs multiple miRNAs (BRCA)
res <- tcga_correlation(
  var1 = c("MYC", "KRAS"), var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = c("hsa-let-7a-1", "hsa-mir-21"), var2_modal = "miRNA", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# RNAseq vs Clinical
# ==============================================================================

# Case 28: TP53 vs Age (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "Age", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 29: ESR1 vs Gender (BRCA)
res <- tcga_correlation(
  var1 = "ESR1", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "Gender", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 30: EGFR vs Stage (LUAD)
res <- tcga_correlation(
  var1 = "EGFR", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "Stage", var2_modal = "Clinical", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 31: Multiple genes vs multiple clinical (BRCA)
res <- tcga_correlation(
  var1 = c("TP53", "ESR1"), var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = c("Age", "Gender", "Race"), var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# RNAseq vs ImmuneCell
# ==============================================================================

# Case 32: CD274 vs T_cells_CD8 (LUAD)
res <- tcga_correlation(
  var1 = "CD274", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "T_cells_CD8", var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 33: PDCD1 vs NK_cells (SKCM)
res <- tcga_correlation(
  var1 = "PDCD1", var1_modal = "RNAseq", var1_cancers = "SKCM",
  var2 = "NK_cells_activated", var2_modal = "ImmuneCell", var2_cancers = "SKCM",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 34: Multiple genes vs immune cells (BRCA)
res <- tcga_correlation(
  var1 = c("CD274", "PDCD1", "CTLA4"), var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = c("T_cells_CD8", "T_cells_CD4_memory_resting", "Macrophages_M1"),
  var2_modal = "ImmuneCell", var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# RNAseq vs Mutation
# ==============================================================================

# Case 35: TP53 expression vs TP53 mutation (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "Mutation", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 36: EGFR vs KRAS mutation (LUAD)
res <- tcga_correlation(
  var1 = "EGFR", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "KRAS", var2_modal = "Mutation", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 37: ESR1 vs PIK3CA mutation (BRCA)
res <- tcga_correlation(
  var1 = "ESR1", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "PIK3CA", var2_modal = "Mutation", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 38: Multiple genes vs multiple mutations (LUAD)
res <- tcga_correlation(
  var1 = c("EGFR", "ALK", "RET"), var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = c("KRAS", "STK11", "KEAP1"), var2_modal = "Mutation", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# RNAseq vs CNV
# ==============================================================================

# Case 39: MYC expression vs MYC CNV (BRCA)
res <- tcga_correlation(
  var1 = "MYC", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "MYC", var2_modal = "CNV", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 40: ERBB2 vs ERBB2 CNV (BRCA)
res <- tcga_correlation(
  var1 = "ERBB2", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "ERBB2", var2_modal = "CNV", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 41: EGFR vs EGFR CNV (LUAD)
res <- tcga_correlation(
  var1 = "EGFR", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "EGFR", var2_modal = "CNV", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 42: CCND1 vs CCND1 CNV (HNSC)
res <- tcga_correlation(
  var1 = "CCND1", var1_modal = "RNAseq", var1_cancers = "HNSC",
  var2 = "CCND1", var2_modal = "CNV", var2_cancers = "HNSC"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# RNAseq Enrichment Analysis
# ==============================================================================

# Case 43: KRAS RNAseq enrichment (LUAD)
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 44: TP53 RNAseq enrichment (BRCA)
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 45: EGFR RNAseq enrichment GO (LUAD)
res <- tcga_enrichment(
  var1 = "EGFR", var1_modal = "RNAseq", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "BP",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 46: MYC RNAseq genome scan (HNSC)
res <- tcga_enrichment(
  var1 = "MYC", var1_modal = "RNAseq", var1_cancers = "HNSC",
  analysis_type = "genome",
  top_n = 50
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 47: Multiple genes enrichment (BRCA)
res <- tcga_enrichment(
  var1 = c("TP53", "ESR1"), var1_modal = "RNAseq", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# RNAseq Survival Analysis
# ==============================================================================

# Case 48: TP53 survival (BRCA)
res <- tcga_survival(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  surv_type = "OS",
  cutoff_type = "optimal"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 49: KRAS survival (LUAD)
res <- tcga_survival(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = "LUAD",
  surv_type = "PFS",
  cutoff_type = "median"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 50: Multiple genes survival forest (BRCA)
res <- tcga_survival(
  var1 = c("TP53", "ESR1", "ERBB2"), var1_modal = "RNAseq", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 51: Multi-cancer survival (TP53)
res <- tcga_survival(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = c("BRCA", "LUAD"),
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# RNAseq Subtype Analysis
# ==============================================================================

# Case 52: ESR1 across BRCA subtypes
res <- tcga_correlation(
  var1 = "ESR1", var1_modal = "RNAseq",
  var1_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC"),
  var2 = "Purity", var2_modal = "Signature",
  var2_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 53: PGR across BRCA subtypes
res <- tcga_correlation(
  var1 = "PGR", var1_modal = "RNAseq",
  var1_cancers = c("BRCA_IDC", "BRCA_ILC"),
  var2 = "TMB", var2_modal = "Signature",
  var2_cancers = c("BRCA_IDC", "BRCA_ILC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 54: ERBB2 in BRCA subtypes
res <- tcga_correlation(
  var1 = "ERBB2", var1_modal = "RNAseq",
  var1_cancers = c("BRCA_IDC", "BRCA_TNBC"),
  var2 = "Stromal", var2_modal = "Signature",
  var2_cancers = c("BRCA_IDC", "BRCA_TNBC")
)
res$plot
head(res$stats)
head(res$raw_data)


