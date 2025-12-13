# SLTCGA Test Suite - Mutation Modality
# Setup
library(devtools)
load_all()
Sys.setenv(SL_BULK_DATA = "/Users/liuzaoqu/Desktop/develop/DataMiner_Dev/bulk_data")

# ==============================================================================
# Mutation vs RNAseq
# ==============================================================================

# Case 1: TP53 mutation vs TP53 expression (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 2: KRAS mutation vs EGFR expression (LUAD)
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
  var2 = "EGFR", var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 3: PIK3CA mutation vs AKT1 expression (BRCA)
res <- tcga_correlation(
  var1 = "PIK3CA", var1_modal = "Mutation", var1_cancers = "BRCA",
  var2 = "AKT1", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 4: BRAF mutation vs downstream targets (SKCM)
res <- tcga_correlation(
  var1 = "BRAF", var1_modal = "Mutation", var1_cancers = "SKCM",
  var2 = c("THBS1", "THBS2", "THBS3"), var2_modal = "RNAseq", var2_cancers = "SKCM"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 5: Multiple mutations vs gene expression (LUAD)
res <- tcga_correlation(
  var1 = c("KRAS", "EGFR", "ALK"), var1_modal = "Mutation", var1_cancers = "LUAD",
  var2 = c("STK11", "KEAP1"), var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Mutation vs Signature
# ==============================================================================

# Case 6: TP53 mutation vs TMB (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 7: KRAS mutation vs TIL_Score (LUAD)
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
  var2 = "TIL_Score", var2_modal = "Signature", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 8: PIK3CA mutation vs Purity (UCEC)
res <- tcga_correlation(
  var1 = "PIK3CA", var1_modal = "Mutation", var1_cancers = "UCEC",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "UCEC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 9: BRCA1 mutation vs HRD (OV)
res <- tcga_correlation(
  var1 = "BRCA1", var1_modal = "Mutation", var1_cancers = "OV",
  var2 = "HRD", var2_modal = "Signature", var2_cancers = "OV"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 10: MLH1 mutation vs MSI (COAD)
res <- tcga_correlation(
  var1 = "MLH1", var1_modal = "Mutation", var1_cancers = "COAD",
  var2 = "MSI", var2_modal = "Signature", var2_cancers = "COAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Mutation vs ImmuneCell
# ==============================================================================

# Case 11: TP53 mutation vs immune cells (LUAD, heatmap)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "LUAD",
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "cibersort",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 12: KRAS mutation vs immune cells (LUAD)
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 13: PIK3CA mutation vs T cells (BRCA)
res <- tcga_correlation(
  var1 = "PIK3CA", var1_modal = "Mutation", var1_cancers = "BRCA",
  var2 = c("T_cells_CD8", "T_cells_CD4_memory_resting"), var2_modal = "ImmuneCell",
  var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 14: BRAF mutation vs macrophages (SKCM)
res <- tcga_correlation(
  var1 = "BRAF", var1_modal = "Mutation", var1_cancers = "SKCM",
  var2 = c("Macrophages_M1", "Macrophages_M2"), var2_modal = "ImmuneCell",
  var2_cancers = "SKCM",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Mutation vs Mutation
# ==============================================================================

# Case 15: TP53 vs PIK3CA mutation (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA",
  var2 = "PIK3CA", var2_modal = "Mutation", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 16: KRAS vs EGFR mutation (LUAD)
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
  var2 = "EGFR", var2_modal = "Mutation", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 17: Multiple mutations (LUAD)
res <- tcga_correlation(
  var1 = c("KRAS", "EGFR", "ALK"), var1_modal = "Mutation", var1_cancers = "LUAD",
  var2 = c("TP53", "STK11", "KEAP1"), var2_modal = "Mutation", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Mutation vs Clinical
# ==============================================================================

# Case 18: TP53 mutation vs Age (BRCA)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA",
  var2 = "Age", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 19: KRAS mutation vs Stage (LUAD)
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
  var2 = "Stage", var2_modal = "Clinical", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 20: PIK3CA mutation vs Gender (BRCA)
res <- tcga_correlation(
  var1 = "PIK3CA", var1_modal = "Mutation", var1_cancers = "BRCA",
  var2 = "Gender", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Mutation Enrichment Analysis
# ==============================================================================

# Case 21: KRAS mutation enrichment (LUAD)
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 22: TP53 mutation enrichment (BRCA)
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 23: PIK3CA mutation enrichment GO (UCEC)
res <- tcga_enrichment(
  var1 = "PIK3CA", var1_modal = "Mutation", var1_cancers = "UCEC",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "BP",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 24: BRAF mutation enrichment (SKCM)
res <- tcga_enrichment(
  var1 = "BRAF", var1_modal = "Mutation", var1_cancers = "SKCM",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 25: Multiple mutations enrichment (LUAD)
res <- tcga_enrichment(
  var1 = c("KRAS", "EGFR"), var1_modal = "Mutation", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Mutation Survival Analysis
# ==============================================================================

# Case 26: TP53 mutation survival (BRCA)
res <- tcga_survival(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 27: KRAS mutation survival (LUAD)
res <- tcga_survival(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
  surv_type = "PFS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 28: Multiple mutations survival (LUAD)
res <- tcga_survival(
  var1 = c("KRAS", "EGFR", "TP53"), var1_modal = "Mutation", var1_cancers = "LUAD",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 29: PIK3CA mutation survival (BRCA)
res <- tcga_survival(
  var1 = "PIK3CA", var1_modal = "Mutation", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 30: Multi-cancer mutation survival
res <- tcga_survival(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = c("BRCA", "LUAD"),
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)


