# SLTCGA Test Suite - Enrichment Analysis (tcga_enrichment)
# Setup
library(devtools)
load_all()
Sys.setenv(SL_BULK_DATA = "/Users/liuzaoqu/Desktop/develop/DataMiner_Dev/bulk_data")

# ==============================================================================
# SECTION 1: Mutation-based Enrichment (DEA -> GSEA)
# ==============================================================================

# Case 1: KRAS mutation MsigDB Hallmark (LUAD)
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

# Case 2: TP53 mutation MsigDB Hallmark (BRCA)
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

# Case 3: PIK3CA mutation MsigDB Hallmark (UCEC)
res <- tcga_enrichment(
  var1 = "PIK3CA", var1_modal = "Mutation", var1_cancers = "UCEC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 4: EGFR mutation MsigDB Hallmark (LUAD)
res <- tcga_enrichment(
  var1 = "EGFR", var1_modal = "Mutation", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 5: BRAF mutation MsigDB Hallmark (SKCM)
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

# Case 6: IDH1 mutation MsigDB Hallmark (GBM)
res <- tcga_enrichment(
  var1 = "IDH1", var1_modal = "Mutation", var1_cancers = "GBM",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 7: VHL mutation MsigDB Hallmark (KIRC)
res <- tcga_enrichment(
  var1 = "VHL", var1_modal = "Mutation", var1_cancers = "KIRC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 2: Mutation-based Enrichment - Different MsigDB Categories
# ==============================================================================

# Case 8: KRAS mutation MsigDB C2 (canonical pathways) (LUAD)
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 9: TP53 mutation MsigDB C6 (oncogenic signatures) (BRCA)
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "C6",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 3: Mutation-based Enrichment - GO Database
# ==============================================================================

# Case 10: KRAS mutation GO BP (LUAD)
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "BP",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 11: TP53 mutation GO MF (BRCA)
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "MF",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 12: PIK3CA mutation GO CC (UCEC)
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

# ==============================================================================
# SECTION 4: Mutation-based Genome Scan (DEA -> Top genes)
# ==============================================================================

# Case 13: KRAS mutation genome scan (LUAD)
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
  analysis_type = "genome",
  top_n = 50
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 14: TP53 mutation genome scan (BRCA)
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA",
  analysis_type = "genome",
  top_n = 50
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 15: BRAF mutation genome scan (SKCM)
res <- tcga_enrichment(
  var1 = "BRAF", var1_modal = "Mutation", var1_cancers = "SKCM",
  analysis_type = "genome",
  top_n = 50
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 5: RNAseq-based Enrichment (Correlation -> GSEA)
# ==============================================================================

# Case 16: KRAS expression MsigDB Hallmark (LUAD)
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

# Case 17: TP53 expression MsigDB Hallmark (BRCA)
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

# Case 18: EGFR expression MsigDB Hallmark (LUAD)
res <- tcga_enrichment(
  var1 = "EGFR", var1_modal = "RNAseq", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 19: MYC expression MsigDB Hallmark (HNSC)
res <- tcga_enrichment(
  var1 = "MYC", var1_modal = "RNAseq", var1_cancers = "HNSC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 20: CD274 expression MsigDB Hallmark (SKCM)
res <- tcga_enrichment(
  var1 = "CD274", var1_modal = "RNAseq", var1_cancers = "SKCM",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 6: RNAseq-based Enrichment - Different MsigDB Categories
# ==============================================================================

# Case 21: KRAS expression MsigDB C2-CP-KEGG (LUAD)
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "C2-CP-KEGG",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 22: TP53 expression MsigDB C6 (oncogenic) (BRCA)
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "C6",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 22b: EGFR expression MsigDB C2-CP-REACTOME (LUAD)
res <- tcga_enrichment(
  var1 = "EGFR", var1_modal = "RNAseq", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "C2-CP-REACTOME",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 7: RNAseq-based Enrichment - GO Database
# ==============================================================================

# Case 23: KRAS expression GO BP (LUAD)
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "BP",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 24: TP53 expression GO MF (BRCA)
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "MF",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 25: EGFR expression GO CC (LUAD)
res <- tcga_enrichment(
  var1 = "EGFR", var1_modal = "RNAseq", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "CC",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 8: RNAseq-based Genome Scan (Correlation -> Top genes)
# ==============================================================================

# Case 26: KRAS expression genome scan (LUAD)
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = "LUAD",
  analysis_type = "genome",
  top_n = 50
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 27: TP53 expression genome scan (BRCA)
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  analysis_type = "genome",
  top_n = 50
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 28: EGFR expression genome scan (LUAD)
res <- tcga_enrichment(
  var1 = "EGFR", var1_modal = "RNAseq", var1_cancers = "LUAD",
  analysis_type = "genome",
  top_n = 50
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 29: MYC expression genome scan (HNSC)
res <- tcga_enrichment(
  var1 = "MYC", var1_modal = "RNAseq", var1_cancers = "HNSC",
  analysis_type = "genome",
  top_n = 50
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 9: Multiple Variables Enrichment (Mutation)
# ==============================================================================

# Case 30: Multiple mutations MsigDB Hallmark (LUAD)
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

# Case 31: Multiple mutations GO BP (LUAD)
res <- tcga_enrichment(
  var1 = c("KRAS", "TP53"), var1_modal = "Mutation", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "BP",
  top_n = 10
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 32: Multiple mutations genome scan (LUAD)
res <- tcga_enrichment(
  var1 = c("KRAS", "EGFR", "TP53"), var1_modal = "Mutation", var1_cancers = "LUAD",
  analysis_type = "genome",
  top_n = 50
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 10: Multiple Variables Enrichment (RNAseq)
# ==============================================================================

# Case 33: Multiple genes MsigDB Hallmark (BRCA)
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

# Case 34: Multiple genes GO BP (BRCA)
res <- tcga_enrichment(
  var1 = c("TP53", "ESR1", "ERBB2"), var1_modal = "RNAseq", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "BP",
  top_n = 10
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 35: Multiple genes genome scan (LUAD)
res <- tcga_enrichment(
  var1 = c("EGFR", "KRAS"), var1_modal = "RNAseq", var1_cancers = "LUAD",
  analysis_type = "genome",
  top_n = 50
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 11: Multi-Cancer Enrichment
# ==============================================================================

# Case 36: KRAS mutation pan-cancer enrichment
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = c("LUAD", "COAD", "PAAD"),
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 37: TP53 mutation pan-cancer enrichment
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = c("BRCA", "LUAD", "COAD"),
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 38: KRAS expression pan-cancer enrichment
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = c("LUAD", "COAD"),
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 12: Different Cancer Types
# ==============================================================================

# Case 39: Mutation enrichment in COAD
res <- tcga_enrichment(
  var1 = "APC", var1_modal = "Mutation", var1_cancers = "COAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 40: Mutation enrichment in HNSC
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "HNSC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 41: Mutation enrichment in LIHC
res <- tcga_enrichment(
  var1 = "CTNNB1", var1_modal = "Mutation", var1_cancers = "LIHC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 42: RNAseq enrichment in KIRC
res <- tcga_enrichment(
  var1 = "VHL", var1_modal = "RNAseq", var1_cancers = "KIRC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 43: RNAseq enrichment in PRAD
res <- tcga_enrichment(
  var1 = "AR", var1_modal = "RNAseq", var1_cancers = "PRAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 44: RNAseq enrichment in THCA
res <- tcga_enrichment(
  var1 = "BRAF", var1_modal = "RNAseq", var1_cancers = "THCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 13: Different Enrichment Databases (Advanced)
# ==============================================================================

# Case 45: KRAS mutation KEGG pathway (LUAD)
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "Mutation", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "KEGG",
  kegg_category = "pathway",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 46: TP53 mutation WikiPathways (BRCA)
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "Wiki",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 47: PIK3CA mutation Reactome (UCEC)
res <- tcga_enrichment(
  var1 = "PIK3CA", var1_modal = "Mutation", var1_cancers = "UCEC",
  analysis_type = "enrichment",
  enrich_database = "Reactome",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 14: Signature-based Enrichment
# ==============================================================================

# Case 48: TMB signature enrichment (BRCA)
res <- tcga_enrichment(
  var1 = "TMB", var1_modal = "Signature", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 49: TIL_Score signature enrichment (LUAD)
res <- tcga_enrichment(
  var1 = "TIL_Score", var1_modal = "Signature", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 50: Purity signature enrichment (COAD)
res <- tcga_enrichment(
  var1 = "Purity", var1_modal = "Signature", var1_cancers = "COAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 51: Stemness signature enrichment (GBM)
res <- tcga_enrichment(
  var1 = "Stemness", var1_modal = "Signature", var1_cancers = "GBM",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 15: ImmuneCell-based Enrichment
# ==============================================================================

# Case 52: T_cells_CD8 enrichment (LUAD)
res <- tcga_enrichment(
  var1 = "T_cells_CD8", var1_modal = "ImmuneCell", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  immune_algorithm = "cibersort",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 53: Macrophages_M1 enrichment (BRCA)
res <- tcga_enrichment(
  var1 = "Macrophages_M1", var1_modal = "ImmuneCell", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  immune_algorithm = "cibersort",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 54: NK_cells enrichment (SKCM)
res <- tcga_enrichment(
  var1 = "NK_cells_activated", var1_modal = "ImmuneCell", var1_cancers = "SKCM",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  immune_algorithm = "cibersort",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 16: CNV-based Enrichment
# ==============================================================================

# Case 55: MYC CNV enrichment (BRCA)
res <- tcga_enrichment(
  var1 = "MYC", var1_modal = "CNV", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 56: ERBB2 CNV enrichment (BRCA)
res <- tcga_enrichment(
  var1 = "ERBB2", var1_modal = "CNV", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 57: EGFR CNV enrichment (LUAD)
res <- tcga_enrichment(
  var1 = "EGFR", var1_modal = "CNV", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 17: Methylation-based Enrichment
# ==============================================================================

# Case 58: TP53 methylation enrichment (BRCA)
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Methylation", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 59: BRCA1 methylation enrichment (BRCA)
res <- tcga_enrichment(
  var1 = "BRCA1", var1_modal = "Methylation", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 60: MLH1 methylation enrichment (COAD)
res <- tcga_enrichment(
  var1 = "MLH1", var1_modal = "Methylation", var1_cancers = "COAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 18: miRNA-based Enrichment
# ==============================================================================

# Case 61: hsa-let-7a-1 enrichment (BRCA)
res <- tcga_enrichment(
  var1 = "hsa-let-7a-1", var1_modal = "miRNA", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 62: hsa-mir-21 enrichment (LUAD)
res <- tcga_enrichment(
  var1 = "hsa-mir-21", var1_modal = "miRNA", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 63: Multiple miRNAs enrichment (BRCA)
res <- tcga_enrichment(
  var1 = c("hsa-let-7a-1", "hsa-mir-21"), var1_modal = "miRNA", var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 15
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 19: Different Correlation Methods
# ==============================================================================

# Case 64: KRAS expression Pearson correlation (LUAD)
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  method = "pearson",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 65: KRAS expression Spearman correlation (LUAD)
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  method = "spearman",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 20: Edge Cases and Stress Tests
# ==============================================================================

# Case 66: Rare cancer (CHOL)
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "CHOL",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 67: Small cancer (DLBC)
res <- tcga_enrichment(
  var1 = "MYC", var1_modal = "RNAseq", var1_cancers = "DLBC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 68: High mutation burden cancer (STAD)
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "STAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# SECTION 21: Subtype-specific Enrichment
# ==============================================================================

# Case 69: BRCA_IDC enrichment
res <- tcga_enrichment(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "BRCA_IDC",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  msigdb_category = "H",
  top_n = 20
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 70: BRCA_TNBC enrichment
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

# Case 71: COAD_LCC enrichment
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

# Case 72: LUAD with KRAS RNAseq genome scan
res <- tcga_enrichment(
  var1 = "KRAS", var1_modal = "RNAseq", var1_cancers = "LUAD",
  analysis_type = "genome",
  top_n = 100
)
res$plot
head(res$stats, 20)
head(res$raw_data)


