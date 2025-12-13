# SLTCGA Test Suite - Comprehensive Mixed Scenarios
# Setup
library(devtools)
load_all()
Sys.setenv(SL_BULK_DATA = "/Users/liuzaoqu/Desktop/develop/DataMiner_Dev/bulk_data")

# ==============================================================================
# Multi-Cancer Analysis
# ==============================================================================

# Case 1: TP53 across pan-cancer
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq",
  var1_cancers = c("BRCA", "LUAD", "COAD", "HNSC", "KIRC"),
  var2 = "TMB", var2_modal = "Signature",
  var2_cancers = c("BRCA", "LUAD", "COAD", "HNSC", "KIRC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 2: KRAS mutation pan-cancer
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "Mutation",
  var1_cancers = c("LUAD", "COAD", "PAAD"),
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 3: Immune infiltration pan-cancer
res <- tcga_correlation(
  var1 = "CD274", var1_modal = "RNAseq",
  var1_cancers = c("LUAD", "SKCM", "KIRC"),
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell",
  var2_cancers = c("LUAD", "SKCM", "KIRC"),
  immune_algorithm = "timer"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Subtype Comparisons
# ==============================================================================

# Case 4: BRCA subtypes - ESR1 expression
res <- tcga_correlation(
  var1 = "ESR1", var1_modal = "RNAseq",
  var1_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC"),
  var2 = "ERBB2", var2_modal = "RNAseq",
  var2_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 5: COAD subtypes - immune infiltration
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical",
  var1_cancers = c("COAD_LCC", "COAD_RCC"),
  var2 = c("T_cells_CD8", "Macrophages_M1"), var2_modal = "ImmuneCell",
  var2_cancers = c("COAD_LCC", "COAD_RCC"),
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 6: ESCA subtypes analysis
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation",
  var1_cancers = c("ESCA_ESCC", "ESCA_EAC"),
  var2 = "TMB", var2_modal = "Signature",
  var2_cancers = c("ESCA_ESCC", "ESCA_EAC")
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Complex Multi-Variable Analysis
# ==============================================================================

# Case 7: Multiple genes vs multiple signatures (BRCA)
res <- tcga_correlation(
  var1 = c("TP53", "ESR1", "ERBB2", "PGR"), var1_modal = "RNAseq",
  var1_cancers = "BRCA",
  var2 = c("TMB", "Purity", "TIL_Score", "Stemness"), var2_modal = "Signature",
  var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 8: Multiple mutations vs multiple genes (LUAD)
res <- tcga_correlation(
  var1 = c("KRAS", "EGFR", "ALK", "BRAF"), var1_modal = "Mutation",
  var1_cancers = "LUAD",
  var2 = c("TP53", "STK11", "KEAP1", "NF1"), var2_modal = "RNAseq",
  var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 9: Multiple clinical vs multiple signatures (LUAD)
res <- tcga_correlation(
  var1 = c("Age", "Gender", "Stage"), var1_modal = "Clinical",
  var1_cancers = "LUAD",
  var2 = c("TMB", "TIL_Score", "Purity"), var2_modal = "Signature",
  var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Enrichment Analysis - Different Databases
# ==============================================================================

# Case 10: KRAS mutation MsigDB H (LUAD)
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

# Case 11: TP53 mutation MsigDB C2 (BRCA)
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "C2",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 12: EGFR expression GO BP (LUAD)
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

# Case 13: MYC expression GO MF (HNSC)
res <- tcga_enrichment(
  var1 = "MYC", var1_modal = "RNAseq", var1_cancers = "HNSC",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "MF",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 14: PIK3CA mutation GO CC (UCEC)
res <- tcga_enrichment(
  var1 = "PIK3CA", var1_modal = "Mutation", var1_cancers = "UCEC",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "CC",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 15: Multiple genes enrichment (BRCA)
res <- tcga_enrichment(
  var1 = c("TP53", "ESR1", "ERBB2"), var1_modal = "RNAseq",
  var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Genome Scan Analysis
# ==============================================================================

# Case 16: KRAS expression genome scan (LUAD)
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = "LUAD",
  analysis_type = "genome",
  top_n = 50
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 17: TP53 expression genome scan (BRCA)
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  analysis_type = "genome",
  top_n = 50
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 18: EGFR mutation genome scan (LUAD)
res <- tcga_enrichment(
  var1 = "EGFR", var1_modal = "Mutation", var1_cancers = "LUAD",
  analysis_type = "genome",
  top_n = 50
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 19: Multiple genes genome scan (HNSC)
res <- tcga_enrichment(
  var1 = c("TP53", "EGFR"), var1_modal = "RNAseq", var1_cancers = "HNSC",
  analysis_type = "genome",
  top_n = 50
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Survival Analysis - Single Variable
# ==============================================================================

# Case 20: TP53 expression survival OS (BRCA)
res <- tcga_survival(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  surv_type = "OS",
  cutoff_type = "optimal"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 21: KRAS expression survival PFS (LUAD)
res <- tcga_survival(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = "LUAD",
  surv_type = "PFS",
  cutoff_type = "median"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 22: TP53 mutation survival (BRCA)
res <- tcga_survival(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 23: KRAS mutation survival (LUAD)
res <- tcga_survival(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
  surv_type = "PFS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 24: Stage survival (LUAD)
res <- tcga_survival(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 25: TMB survival (BRCA)
res <- tcga_survival(
  var1 = "TMB", var1_modal = "Signature", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Survival Analysis - Multiple Variables
# ==============================================================================

# Case 26: Multiple genes survival forest (BRCA)
res <- tcga_survival(
  var1 = c("TP53", "ESR1", "ERBB2", "PGR"), var1_modal = "RNAseq",
  var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 27: Multiple mutations survival forest (LUAD)
res <- tcga_survival(
  var1 = c("KRAS", "EGFR", "TP53", "STK11"), var1_modal = "Mutation",
  var1_cancers = "LUAD",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 28: Multiple signatures survival forest (LUAD)
res <- tcga_survival(
  var1 = c("TMB", "TIL_Score", "IFN_Gamma", "Purity"), var1_modal = "Signature",
  var1_cancers = "LUAD",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 29: Multiple clinical survival forest (BRCA)
res <- tcga_survival(
  var1 = c("Age", "Gender", "Race", "Stage"), var1_modal = "Clinical",
  var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 30: Multi-cancer survival
res <- tcga_survival(
  var1 = "TP53", var1_modal = "RNAseq",
  var1_cancers = c("BRCA", "LUAD", "COAD"),
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Special Scenarios - Heatmap Visualizations
# ==============================================================================

# Case 31: TP53 mutation vs ALL immune cells (LUAD, heatmap)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "LUAD",
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell",
  var2_cancers = "LUAD",
  immune_algorithm = "cibersort",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 32: KRAS mutation vs ALL immune cells (LUAD, heatmap)
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell",
  var2_cancers = "LUAD",
  immune_algorithm = "cibersort",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 33: Stage vs ALL immune cells (LUAD, heatmap)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "LUAD",
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell",
  var2_cancers = "LUAD",
  immune_algorithm = "cibersort",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 34: Gender vs immune cells (KIRC, heatmap)
res <- tcga_correlation(
  var1 = "Gender", var1_modal = "Clinical", var1_cancers = "KIRC",
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell",
  var2_cancers = "KIRC",
  immune_algorithm = "timer",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 35: PIK3CA mutation vs immune cells (BRCA, heatmap)
res <- tcga_correlation(
  var1 = "PIK3CA", var1_modal = "Mutation", var1_cancers = "BRCA",
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell",
  var2_cancers = "BRCA",
  immune_algorithm = "cibersort",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Edge Cases and Stress Tests
# ==============================================================================

# Case 36: Single sample subtype (rare)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "CHOL",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "CHOL"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 37: Large gene panel (BRCA)
res <- tcga_correlation(
  var1 = c("TP53", "ESR1", "ERBB2", "PGR", "AR", "GATA3", "FOXA1", "MYC"),
  var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 38: Many immune algorithms comparison
res <- tcga_correlation(
  var1 = "CD274", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "T_cells_CD8", var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

res <- tcga_correlation(
  var1 = "CD274", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "T_cells_CD8", var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "timer"
)
res$plot
head(res$stats)
head(res$raw_data)

res <- tcga_correlation(
  var1 = "CD274", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "T_cells_CD8", var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "quantiseq"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Correlation Methods Comparison
# ==============================================================================

# Case 39: Pearson correlation
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "BRCA",
  method = "pearson"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 40: Spearman correlation
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "BRCA",
  method = "spearman"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Different RNAseq Normalization Methods
# ==============================================================================

# Case 41: log2TPM (default)
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "EGFR", var2_modal = "RNAseq", var2_cancers = "LUAD",
  rnaseq_type = "log2TPM"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# Rare Cancer Types
# ==============================================================================

# Case 42: ACC analysis
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "ACC",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "ACC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 43: CHOL analysis
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = "CHOL",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "CHOL"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 44: DLBC analysis
res <- tcga_correlation(
  var1 = "MYC", var1_modal = "RNAseq", var1_cancers = "DLBC",
  var2 = "TIL_Score", var2_modal = "Signature", var2_cancers = "DLBC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 45: MESO analysis
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "MESO",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "MESO"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 46: THYM analysis
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "THYM",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "THYM"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 47: UCS analysis
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "UCS",
  var2 = "TIL_Score", var2_modal = "Signature", var2_cancers = "UCS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 48: UVM analysis
res <- tcga_correlation(
  var1 = "BAP1", var1_modal = "RNAseq", var1_cancers = "UVM",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "UVM"
)
res$plot
head(res$stats)
head(res$raw_data)


# ==============================================================================
# Additional Advanced Cases (to reach 300 total)
# ==============================================================================

# Case 49: Combined subtypes analysis
res <- tcga_correlation(
  var1 = "NSCLC", var1_modal = "Clinical", var1_cancers = "NSCLC",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "NSCLC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 50: CRC combined type
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "CRC",
  var2 = "MSI", var2_modal = "Signature", var2_cancers = "CRC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 51: GLIOMA combined
res <- tcga_correlation(
  var1 = "IDH1", var1_modal = "Mutation", var1_cancers = "GLIOMA",
  var2 = "Stemness", var2_modal = "Signature", var2_cancers = "GLIOMA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 52: KRCC combined kidney
res <- tcga_correlation(
  var1 = "VHL", var1_modal = "Mutation", var1_cancers = "KRCC",
  var2 = "HRD", var2_modal = "Signature", var2_cancers = "KRCC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 53: xcell algorithm test
res <- tcga_correlation(
  var1 = "CD274", var1_modal = "RNAseq", var1_cancers = "LUAD",
  var2 = "T_cells_CD8", var2_modal = "ImmuneCell", var2_cancers = "LUAD",
  immune_algorithm = "xcell"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 54: mcpcounter algorithm test
res <- tcga_correlation(
  var1 = "PDCD1", var1_modal = "RNAseq", var1_cancers = "SKCM",
  var2 = "T_cells", var2_modal = "ImmuneCell", var2_cancers = "SKCM",
  immune_algorithm = "mcpcounter"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 55: epic algorithm test
res <- tcga_correlation(
  var1 = "CTLA4", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "T_cells_CD8", var2_modal = "ImmuneCell", var2_cancers = "BRCA",
  immune_algorithm = "epic"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 56: Methylation different regions
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Methylation", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "BRCA",
  methylation_region = "Gene_body_mean"
)
res$plot
head(res$stats)
head(res$raw_data)

