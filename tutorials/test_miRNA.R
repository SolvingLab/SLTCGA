# SLTCGA Test Suite - miRNA Modality
# Setup
library(devtools)
load_all()
Sys.setenv(SL_BULK_DATA = "/Users/liuzaoqu/Desktop/develop/DataMiner_Dev/bulk_data")

# ==============================================================================
# miRNA vs RNAseq
# ==============================================================================

# Case 1: hsa-let-7a-1 vs MYC (BRCA)
res <- tcga_correlation(
  var1 = "hsa-let-7a-1", var1_modal = "miRNA", var1_cancers = "BRCA",
  var2 = "MYC", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 2: hsa-mir-21 vs KRAS (LUAD)
res <- tcga_correlation(
  var1 = "hsa-mir-21", var1_modal = "miRNA", var1_cancers = "LUAD",
  var2 = "KRAS", var2_modal = "RNAseq", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 3: hsa-mir-200c vs CDH1 (BRCA)
res <- tcga_correlation(
  var1 = "hsa-mir-200c", var1_modal = "miRNA", var1_cancers = "BRCA",
  var2 = "CDH1", var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 4: hsa-mir-155 vs TP53 (DLBC)
res <- tcga_correlation(
  var1 = "hsa-mir-155", var1_modal = "miRNA", var1_cancers = "DLBC",
  var2 = "TP53", var2_modal = "RNAseq", var2_cancers = "DLBC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 5: Multiple miRNAs vs multiple genes (BRCA)
res <- tcga_correlation(
  var1 = c("hsa-let-7a-1", "hsa-mir-21", "hsa-mir-200c"),
  var1_modal = "miRNA", var1_cancers = "BRCA",
  var2 = c("MYC", "KRAS", "CDH1"), var2_modal = "RNAseq", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# miRNA vs Signature
# ==============================================================================

# Case 6: hsa-let-7a-1 vs TMB (BRCA)
res <- tcga_correlation(
  var1 = "hsa-let-7a-1", var1_modal = "miRNA", var1_cancers = "BRCA",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 7: hsa-mir-21 vs Purity (LUAD)
res <- tcga_correlation(
  var1 = "hsa-mir-21", var1_modal = "miRNA", var1_cancers = "LUAD",
  var2 = "Purity", var2_modal = "Signature", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 8: hsa-mir-155 vs TIL_Score (DLBC)
res <- tcga_correlation(
  var1 = "hsa-mir-155", var1_modal = "miRNA", var1_cancers = "DLBC",
  var2 = "TIL_Score", var2_modal = "Signature", var2_cancers = "DLBC"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 9: Multiple miRNAs vs signatures (BRCA)
res <- tcga_correlation(
  var1 = c("hsa-let-7a-1", "hsa-mir-21"), var1_modal = "miRNA", var1_cancers = "BRCA",
  var2 = c("TMB", "Purity", "TIL_Score"), var2_modal = "Signature", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# miRNA vs Mutation
# ==============================================================================

# Case 10: hsa-let-7a-1 vs KRAS mutation (LUAD)
res <- tcga_correlation(
  var1 = "hsa-let-7a-1", var1_modal = "miRNA", var1_cancers = "LUAD",
  var2 = "KRAS", var2_modal = "Mutation", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 11: hsa-mir-21 vs TP53 mutation (BRCA)
res <- tcga_correlation(
  var1 = "hsa-mir-21", var1_modal = "miRNA", var1_cancers = "BRCA",
  var2 = "TP53", var2_modal = "Mutation", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 12: hsa-mir-200c vs CDH1 mutation (BRCA)
res <- tcga_correlation(
  var1 = "hsa-mir-200c", var1_modal = "miRNA", var1_cancers = "BRCA",
  var2 = "CDH1", var2_modal = "Mutation", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# miRNA vs ImmuneCell
# ==============================================================================

# Case 13: hsa-let-7a-1 vs T cells (LUAD)
res <- tcga_correlation(
  var1 = "hsa-let-7a-1", var1_modal = "miRNA", var1_cancers = "LUAD",
  var2 = c("T_cells_CD8", "T_cells_CD4_memory_resting"), var2_modal = "ImmuneCell",
  var2_cancers = "LUAD",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 14: hsa-mir-21 vs macrophages (BRCA)
res <- tcga_correlation(
  var1 = "hsa-mir-21", var1_modal = "miRNA", var1_cancers = "BRCA",
  var2 = c("Macrophages_M1", "Macrophages_M2"), var2_modal = "ImmuneCell",
  var2_cancers = "BRCA",
  immune_algorithm = "cibersort"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# miRNA vs miRNA
# ==============================================================================

# Case 15: hsa-let-7a-1 vs hsa-mir-21 (BRCA)
res <- tcga_correlation(
  var1 = "hsa-let-7a-1", var1_modal = "miRNA", var1_cancers = "BRCA",
  var2 = "hsa-mir-21", var2_modal = "miRNA", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 16: Multiple miRNAs correlation (LUAD)
res <- tcga_correlation(
  var1 = c("hsa-let-7a-1", "hsa-mir-21"), var1_modal = "miRNA", var1_cancers = "LUAD",
  var2 = c("hsa-mir-200c", "hsa-mir-155"), var2_modal = "miRNA", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# miRNA vs Clinical
# ==============================================================================

# Case 17: hsa-let-7a-1 vs Age (BRCA)
res <- tcga_correlation(
  var1 = "hsa-let-7a-1", var1_modal = "miRNA", var1_cancers = "BRCA",
  var2 = "Age", var2_modal = "Clinical", var2_cancers = "BRCA"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 18: hsa-mir-21 vs Stage (LUAD)
res <- tcga_correlation(
  var1 = "hsa-mir-21", var1_modal = "miRNA", var1_cancers = "LUAD",
  var2 = "Stage", var2_modal = "Clinical", var2_cancers = "LUAD"
)
res$plot
head(res$stats)
head(res$raw_data)

# ==============================================================================
# miRNA Survival Analysis
# ==============================================================================

# Case 19: hsa-let-7a-1 survival (BRCA)
res <- tcga_survival(
  var1 = "hsa-let-7a-1", var1_modal = "miRNA", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 20: hsa-mir-21 survival (LUAD)
res <- tcga_survival(
  var1 = "hsa-mir-21", var1_modal = "miRNA", var1_cancers = "LUAD",
  surv_type = "PFS"
)
res$plot
head(res$stats)
head(res$raw_data)

# Case 21: Multiple miRNAs survival (BRCA)
res <- tcga_survival(
  var1 = c("hsa-let-7a-1", "hsa-mir-21", "hsa-mir-200c"),
  var1_modal = "miRNA", var1_cancers = "BRCA",
  surv_type = "OS"
)
res$plot
head(res$stats)
head(res$raw_data)


