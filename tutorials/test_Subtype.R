# SLTCGA Test Suite - Cancer Subtype Analysis
# Testing all 32 TCGA subtypes + 4 combined types
# Setup
library(devtools)
load_all()
Sys.setenv(SL_BULK_DATA = "/Users/liuzaoqu/Desktop/develop/DataMiner_Dev/bulk_data")

# ==============================================================================
# SECTION 1: BRCA Subtypes (3 subtypes)
# ==============================================================================

# Case 1: ESR1 across BRCA subtypes
res <- tcga_correlation(
  var1 = "ESR1", var1_modal = "RNAseq",
  var1_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC"),
  var2 = "PGR", var2_modal = "RNAseq",
  var2_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 2: TP53 mutation across BRCA subtypes
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation",
  var1_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC"),
  var2 = "TMB", var2_modal = "Signature",
  var2_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 3: PIK3CA mutation across BRCA subtypes
res <- tcga_correlation(
  var1 = "PIK3CA", var1_modal = "Mutation",
  var1_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC"),
  var2 = "Purity", var2_modal = "Signature",
  var2_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 4: ERBB2 CNV across BRCA subtypes
res <- tcga_correlation(
  var1 = "ERBB2", var1_modal = "CNV",
  var1_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC"),
  var2 = "TIL_Score", var2_modal = "Signature",
  var2_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 5: Immune infiltration across BRCA subtypes
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical",
  var1_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC"),
  var2 = c("T_cells_CD8", "Macrophages_M1"), var2_modal = "ImmuneCell",
  var2_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC"),
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 6: BRCA_IDC vs BRCA_ILC only
res <- tcga_correlation(
  var1 = "ESR1", var1_modal = "RNAseq",
  var1_cancers = c("BRCA_IDC", "BRCA_ILC"),
  var2 = "ERBB2", var2_modal = "RNAseq",
  var2_cancers = c("BRCA_IDC", "BRCA_ILC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 7: BRCA_IDC survival
res <- tcga_survival(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA_IDC",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 8: BRCA_TNBC enrichment
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA_TNBC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 2: COAD Subtypes (3 subtypes)
# ==============================================================================

# Case 9: KRAS across COAD subtypes
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "Mutation",
  var1_cancers = c("COAD_LCC", "COAD_RCC", "COAD_MAC"),
  var2 = "MSI", var2_modal = "Signature",
  var2_cancers = c("COAD_LCC", "COAD_RCC", "COAD_MAC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 10: APC mutation across COAD subtypes
res <- tcga_correlation(
  var1 = "APC", var1_modal = "Mutation",
  var1_cancers = c("COAD_LCC", "COAD_RCC"),
  var2 = "Stemness", var2_modal = "Signature",
  var2_cancers = c("COAD_LCC", "COAD_RCC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 11: Stage across COAD subtypes
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical",
  var1_cancers = c("COAD_LCC", "COAD_MAC"),
  var2 = "TIL_Score", var2_modal = "Signature",
  var2_cancers = c("COAD_LCC", "COAD_MAC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 12: COAD_LCC enrichment
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "COAD_LCC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 3: ESCA Subtypes (2 subtypes)
# ==============================================================================

# Case 13: TP53 across ESCA subtypes
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation",
  var1_cancers = c("ESCA_ESCC", "ESCA_EAC"),
  var2 = "TMB", var2_modal = "Signature",
  var2_cancers = c("ESCA_ESCC", "ESCA_EAC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 14: CDKN2A across ESCA subtypes
res <- tcga_correlation(
  var1 = "CDKN2A", var1_modal = "RNAseq",
  var1_cancers = c("ESCA_ESCC", "ESCA_EAC"),
  var2 = "Purity", var2_modal = "Signature",
  var2_cancers = c("ESCA_ESCC", "ESCA_EAC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 15: ESCA_ESCC enrichment
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "ESCA_ESCC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 4: HNSC Subtypes (2 subtypes)
# ==============================================================================

# Case 16: TP53 across HNSC subtypes
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation",
  var1_cancers = c("HNSC_OSCC", "HNSC_LSCC"),
  var2 = "TIL_Score", var2_modal = "Signature",
  var2_cancers = c("HNSC_OSCC", "HNSC_LSCC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 17: EGFR across HNSC subtypes
res <- tcga_correlation(
  var1 = "EGFR", var1_modal = "RNAseq",
  var1_cancers = c("HNSC_OSCC", "HNSC_LSCC"),
  var2 = "Purity", var2_modal = "Signature",
  var2_cancers = c("HNSC_OSCC", "HNSC_LSCC")
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 5: LGG Subtypes (3 subtypes)
# ==============================================================================

# Case 18: IDH1 across LGG subtypes
res <- tcga_correlation(
  var1 = "IDH1", var1_modal = "Mutation",
  var1_cancers = c("LGG_ASTROCYTOMA", "LGG_OLIGODENDROGLIOMA"),
  var2 = "Stemness", var2_modal = "Signature",
  var2_cancers = c("LGG_ASTROCYTOMA", "LGG_OLIGODENDROGLIOMA")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 19: TP53 across LGG subtypes
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq",
  var1_cancers = c("LGG_ASTROCYTOMA", "LGG_OLIGOASTROCYTOMA"),
  var2 = "TIL_Score", var2_modal = "Signature",
  var2_cancers = c("LGG_ASTROCYTOMA", "LGG_OLIGOASTROCYTOMA")
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 6: SKCM Subtypes (2 subtypes)
# ==============================================================================

# Case 20: BRAF across SKCM subtypes
res <- tcga_correlation(
  var1 = "BRAF", var1_modal = "Mutation",
  var1_cancers = c("SKCM_MCM", "SKCM_PCM"),
  var2 = "TMB", var2_modal = "Signature",
  var2_cancers = c("SKCM_MCM", "SKCM_PCM")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 21: CD274 across SKCM subtypes
res <- tcga_correlation(
  var1 = "CD274", var1_modal = "RNAseq",
  var1_cancers = c("SKCM_MCM", "SKCM_PCM"),
  var2 = "TIL_Score", var2_modal = "Signature",
  var2_cancers = c("SKCM_MCM", "SKCM_PCM")
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 7: STAD Subtypes (2 subtypes)
# ==============================================================================

# Case 22: TP53 across STAD subtypes
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation",
  var1_cancers = c("STAD_DGA", "STAD_SRCC"),
  var2 = "MSI", var2_modal = "Signature",
  var2_cancers = c("STAD_DGA", "STAD_SRCC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 23: CDH1 across STAD subtypes
res <- tcga_correlation(
  var1 = "CDH1", var1_modal = "RNAseq",
  var1_cancers = c("STAD_DGA", "STAD_SRCC"),
  var2 = "Purity", var2_modal = "Signature",
  var2_cancers = c("STAD_DGA", "STAD_SRCC")
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 8: UCEC Subtypes (2 subtypes)
# ==============================================================================

# Case 24: PIK3CA across UCEC subtypes
res <- tcga_correlation(
  var1 = "PIK3CA", var1_modal = "Mutation",
  var1_cancers = c("UCEC_EEA", "UCEC_SEA"),
  var2 = "TMB", var2_modal = "Signature",
  var2_cancers = c("UCEC_EEA", "UCEC_SEA")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 25: TP53 across UCEC subtypes
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq",
  var1_cancers = c("UCEC_EEA", "UCEC_SEA"),
  var2 = "Purity", var2_modal = "Signature",
  var2_cancers = c("UCEC_EEA", "UCEC_SEA")
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 9: CESC Subtypes (1 subtype)
# ==============================================================================

# Case 26: CESC_CSCC analysis
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "CESC_CSCC",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "CESC_CSCC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 27: CESC_CSCC enrichment
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "CESC_CSCC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 10: PAAD Subtypes (1 subtype)
# ==============================================================================

# Case 28: PAAD_PADC analysis
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "PAAD_PADC",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "PAAD_PADC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 29: PAAD_PADC enrichment
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "PAAD_PADC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 11: PCPG Subtypes (2 subtypes)
# ==============================================================================

# Case 30: PCPG subtypes comparison
res <- tcga_correlation(
  var1 = "VHL", var1_modal = "Mutation",
  var1_cancers = c("PCPG_PHEOCHROMOCYTOMA", "PCPG_PARAGANGLIOMA"),
  var2 = "Purity", var2_modal = "Signature",
  var2_cancers = c("PCPG_PHEOCHROMOCYTOMA", "PCPG_PARAGANGLIOMA")
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 12: SARC Subtypes (2 subtypes)
# ==============================================================================

# Case 31: TP53 across SARC subtypes
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation",
  var1_cancers = c("SARC_DDLPS", "SARC_LMS"),
  var2 = "TMB", var2_modal = "Signature",
  var2_cancers = c("SARC_DDLPS", "SARC_LMS")
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 13: TGCT Subtypes (1 subtype)
# ==============================================================================

# Case 32: TGCT_SEMINOMA analysis
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "TGCT_SEMINOMA",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "TGCT_SEMINOMA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 14: THCA Subtypes (1 subtype)
# ==============================================================================

# Case 33: THCA_CPTC analysis
res <- tcga_correlation(
  var1 = "BRAF", var1_modal = "Mutation", var1_cancers = "THCA_CPTC",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "THCA_CPTC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 34: THCA_CPTC RNAseq
res <- tcga_correlation(
  var1 = "BRAF", var1_modal = "RNAseq", var1_cancers = "THCA_CPTC",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "THCA_CPTC"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 15: Combined Cancer Types (4 types)
# ==============================================================================

# Case 35: NSCLC (LUAD + LUSC combined)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "NSCLC",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "NSCLC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 36: CRC (COAD + READ combined)
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "CRC",
  var2 = "MSI", var2_modal = "Signature", var2_cancers = "CRC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 37: GLIOMA (GBM + LGG combined)
res <- tcga_correlation(
  var1 = "IDH1", var1_modal = "Mutation", var1_cancers = "GLIOMA",
  var2 = "Stemness", var2_modal = "Signature", var2_cancers = "GLIOMA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 38: KRCC (KIRC + KIRP + KICH combined)
res <- tcga_correlation(
  var1 = "VHL", var1_modal = "Mutation", var1_cancers = "KRCC",
  var2 = "HRD", var2_modal = "Signature", var2_cancers = "KRCC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 39: NSCLC enrichment
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "NSCLC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 40: CRC enrichment
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "CRC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 16: Cross-Subtype Comparisons (Main vs Subtype)
# ==============================================================================

# Case 41: BRCA (all) vs BRCA_IDC
res <- tcga_correlation(
  var1 = "ESR1", var1_modal = "RNAseq",
  var1_cancers = c("BRCA", "BRCA_IDC"),
  var2 = "PGR", var2_modal = "RNAseq",
  var2_cancers = c("BRCA", "BRCA_IDC")
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 42: COAD (all) vs COAD_LCC
res <- tcga_correlation(
  var1 = "KRAS", var1_modal = "Mutation",
  var1_cancers = c("COAD", "COAD_LCC"),
  var2 = "MSI", var2_modal = "Signature",
  var2_cancers = c("COAD", "COAD_LCC")
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 17: Subtype-Specific Survival Analysis
# ==============================================================================

# Case 43: BRCA_IDC survival
res <- tcga_survival(
  var1 = "ESR1", var1_modal = "RNAseq", var1_cancers = "BRCA_IDC",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 44: BRCA_TNBC survival
res <- tcga_survival(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA_TNBC",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 45: COAD_LCC survival
res <- tcga_survival(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "COAD_LCC",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 46: ESCA_ESCC survival
res <- tcga_survival(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "ESCA_ESCC",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 18: Multi-Subtype Enrichment
# ==============================================================================

# Case 47: BRCA subtypes enrichment comparison
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation",
  var1_cancers = c("BRCA_IDC", "BRCA_TNBC"),
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 48: COAD subtypes enrichment
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "Mutation",
  var1_cancers = c("COAD_LCC", "COAD_RCC"),
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 19: Subtype Immune Profiling
# ==============================================================================

# Case 49: Immune infiltration across BRCA subtypes (heatmap)
res <- tcga_correlation(
  var1 = "TP53", var1_modal = "Mutation",
  var1_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC"),
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell",
  var2_cancers = c("BRCA_IDC", "BRCA_ILC", "BRCA_TNBC"),
  immune_algorithm = "cibersort",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 50: Stage vs immune cells in BRCA_IDC (heatmap)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "BRCA_IDC",
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell", var2_cancers = "BRCA_IDC",
  immune_algorithm = "cibersort",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 51: Stage vs immune cells in COAD_LCC (heatmap)
res <- tcga_correlation(
  var1 = "Stage", var1_modal = "Clinical", var1_cancers = "COAD_LCC",
  var2 = "ALL_IMMUNE_CELLS", var2_modal = "ImmuneCell", var2_cancers = "COAD_LCC",
  immune_algorithm = "cibersort",
  plot_type = "heatmap"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 20: All Subtypes Summary Test
# ==============================================================================

# Case 52: Test data availability for each subtype
subtypes <- c(
  "BRCA_IDC", "BRCA_ILC", "BRCA_TNBC",
  "CESC_CSCC",
  "COAD_LCC", "COAD_RCC", "COAD_MAC",
  "ESCA_ESCC", "ESCA_EAC",
  "HNSC_OSCC", "HNSC_LSCC",
  "LGG_ASTROCYTOMA", "LGG_OLIGOASTROCYTOMA", "LGG_OLIGODENDROGLIOMA",
  "PAAD_PADC",
  "PCPG_PHEOCHROMOCYTOMA", "PCPG_PARAGANGLIOMA",
  "SARC_DDLPS", "SARC_LMS",
  "SKCM_MCM", "SKCM_PCM",
  "STAD_DGA", "STAD_SRCC",
  "TGCT_SEMINOMA",
  "THCA_CPTC",
  "UCEC_EEA", "UCEC_SEA"
)

for (subtype in subtypes) {
  res <- tcga_correlation(
    var1 = "TP53", var1_modal = "RNAseq",
    var1_cancers = subtype,
    var2 = "TMB", var2_modal = "Signature",
    var2_cancers = subtype
  )
  res$plot
  head(res$stats)
}


