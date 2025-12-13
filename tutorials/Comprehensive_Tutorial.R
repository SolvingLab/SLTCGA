################################################################################
# SLCPTAC 完整教程
# 覆盖所有模态、所有场景、所有组合
# 特别关注：磷酸化、突变、临床变量
################################################################################

library(SLCPTAC)

################################################################################
# 第一部分：数据加载测试（cptac_load_modality）
################################################################################

# ==================== 单基因单癌种 ====================

# 1.1 RNAseq - BRCA
load_1_1 <- cptac_load_modality(
  var1 = "TP53",
  var1_modal = "RNAseq",
  var1_cancers = "BRCA"
)
load_1_1$var1_features
load_1_1$var1_types
head(load_1_1$data)

# 1.2 Protein - LUAD
load_1_2 <- cptac_load_modality(
  var1 = "EGFR",
  var1_modal = "Protein",
  var1_cancers = "LUAD"
)
load_1_2$var1_features
dim(load_1_2$data)

# 1.3 Phospho - BRCA（重点）
load_1_3 <- cptac_load_modality(
  var1 = "AKT1",
  var1_modal = "Phospho",
  var1_cancers = "BRCA"
)
load_1_3$var1_features
length(load_1_3$var1_features)

# 1.4 Mutation - LUAD（重点）
load_1_4 <- cptac_load_modality(
  var1 = "KRAS",
  var1_modal = "Mutation",
  var1_cancers = "LUAD"
)
load_1_4$var1_features
load_1_4$var1_types
table(load_1_4$data[,1])

# 1.5 Clinical - BRCA（重点）
load_1_5 <- cptac_load_modality(
  var1 = "Tumor_Stage",
  var1_modal = "Clinical",
  var1_cancers = "BRCA"
)
load_1_5$var1_features
load_1_5$var1_types
table(load_1_5$data[,1])

# 1.6 logCNA - PDAC
load_1_6 <- cptac_load_modality(
  var1 = "MYC",
  var1_modal = "logCNA",
  var1_cancers = "PDAC"
)
load_1_6$var1_features
summary(load_1_6$data[,1])

# ==================== 多基因单癌种 ====================

# 2.1 多RNAseq
load_2_1 <- cptac_load_modality(
  var1 = c("TP53", "EGFR", "KRAS"),
  var1_modal = "RNAseq",
  var1_cancers = "BRCA"
)
load_2_1$var1_features
length(load_2_1$var1_features)

# 2.2 多Phospho（重点）
load_2_2 <- cptac_load_modality(
  var1 = c("AKT1", "MTOR", "RPS6"),
  var1_modal = "Phospho",
  var1_cancers = "BRCA"
)
load_2_2$var1_features
length(load_2_2$var1_features)

# 2.3 多Mutation（重点）
load_2_3 <- cptac_load_modality(
  var1 = c("KRAS", "EGFR", "TP53", "ALK"),
  var1_modal = "Mutation",
  var1_cancers = "LUAD"
)
load_2_3$var1_features
sapply(load_2_3$data[,1:4], table)

# 2.4 多Clinical（重点）
load_2_4 <- cptac_load_modality(
  var1 = c("Age", "Tumor_Stage"),
  var1_modal = "Clinical",
  var1_cancers = "BRCA"
)
load_2_4$var1_features
load_2_4$var1_types

# ==================== 单基因多癌种 ====================

# 3.1 RNAseq跨癌种
load_3_1 <- cptac_load_modality(
  var1 = "TP53",
  var1_modal = "RNAseq",
  var1_cancers = c("BRCA", "LUAD", "COAD")
)
load_3_1$var1_features
length(load_3_1$var1_features)

# 3.2 Phospho跨癌种（重点）
load_3_2 <- cptac_load_modality(
  var1 = "AKT1",
  var1_modal = "Phospho",
  var1_cancers = c("BRCA", "LUAD", "CCRCC")
)
load_3_2$var1_features
length(load_3_2$var1_features)

# 3.3 Mutation跨癌种（重点）
load_3_3 <- cptac_load_modality(
  var1 = "TP53",
  var1_modal = "Mutation",
  var1_cancers = c("BRCA", "LUAD", "COAD")
)
load_3_3$var1_features
sapply(load_3_3$data[,1:3], table)

# ==================== 双变量组合 ====================

# 4.1 RNAseq vs Protein
load_4_1 <- cptac_load_modality(
  var1 = "TP53",
  var1_modal = "RNAseq",
  var1_cancers = "BRCA",
  var2 = "TP53",
  var2_modal = "Protein",
  var2_cancers = "BRCA"
)
load_4_1$var1_features
load_4_1$var2_features
dim(load_4_1$data)

# 4.2 RNAseq vs Phospho（重点）
load_4_2 <- cptac_load_modality(
  var1 = "AKT1",
  var1_modal = "RNAseq",
  var1_cancers = "BRCA",
  var2 = "AKT1",
  var2_modal = "Phospho",
  var2_cancers = "BRCA"
)
load_4_2$var1_features
load_4_2$var2_features
ncol(load_4_2$data)

# 4.3 Mutation vs Phospho（重点）
load_4_3 <- cptac_load_modality(
  var1 = "PIK3CA",
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  var2 = "AKT1",
  var2_modal = "Phospho",
  var2_cancers = "BRCA"
)
load_4_3$var1_types
load_4_3$var2_types

# 4.4 Clinical vs Phospho（重点）
load_4_4 <- cptac_load_modality(
  var1 = "Tumor_Stage",
  var1_modal = "Clinical",
  var1_cancers = "LUAD",
  var2 = "AKT1",
  var2_modal = "Phospho",
  var2_cancers = "LUAD"
)
load_4_4$var1_types
load_4_4$var2_types

################################################################################
# 第二部分：相关性分析（cptac_correlation）- 场景1-7
################################################################################

# ==================== 场景1: 1 continuous vs 1 continuous ====================

# 场景1.1: RNAseq vs Protein - 单癌种（CorPlot）
corr_1_1 <- cptac_correlation(
  var1 = "TP53",
  var1_modal = "RNAseq",
  var1_cancers = "BRCA",
  var2 = "TP53",
  var2_modal = "Protein",
  var2_cancers = "BRCA"
)
head(corr_1_1$stats)
head(corr_1_1$raw_data)
print(corr_1_1$plot)

# 场景1.2: RNAseq vs Phospho - 单癌种（重点）
corr_1_2 <- cptac_correlation(
  var1 = "AKT1",
  var1_modal = "RNAseq",
  var1_cancers = "BRCA",
  var2 = "AKT1",
  var2_modal = "Phospho",
  var2_cancers = "BRCA"
)
head(corr_1_2$stats)
head(corr_1_2$raw_data)
print(corr_1_2$plot)

# 场景1.3: Protein vs Phospho - 单癌种（重点）
corr_1_3 <- cptac_correlation(
  var1 = "MTOR",
  var1_modal = "Protein",
  var1_cancers = "LUAD",
  var2 = "MTOR",
  var2_modal = "Phospho",
  var2_cancers = "LUAD"
)
head(corr_1_3$stats)
head(corr_1_3$raw_data)
print(corr_1_3$plot)

# 场景1.4: RNAseq vs Protein - 多癌种（LollipopPlot）
corr_1_4 <- cptac_correlation(
  var1 = "TP53",
  var1_modal = "RNAseq",
  var1_cancers = c("BRCA", "LUAD", "COAD"),
  var2 = "TP53",
  var2_modal = "Protein",
  var2_cancers = c("BRCA", "LUAD", "COAD")
)
head(corr_1_4$stats)
head(corr_1_4$raw_data)
print(corr_1_4$plot)

# 场景1.5: Phospho vs Phospho - 多癌种（重点）
corr_1_5 <- cptac_correlation(
  var1 = "AKT1",
  var1_modal = "Phospho",
  var1_cancers = c("BRCA", "LUAD"),
  var2 = "MTOR",
  var2_modal = "Phospho",
  var2_cancers = c("BRCA", "LUAD")
)
head(corr_1_5$stats)
head(corr_1_5$raw_data)
print(corr_1_5$plot)

# ==================== 场景2: 1 continuous vs multiple continuous ====================

# 场景2.1: 1 RNAseq vs 多 Protein - 单癌种（LollipopPlot）
corr_2_1 <- cptac_correlation(
  var1 = "TP53",
  var1_modal = "RNAseq",
  var1_cancers = "BRCA",
  var2 = c("AKT1", "MTOR", "PTEN"),
  var2_modal = "Protein",
  var2_cancers = "BRCA"
)
head(corr_2_1$stats)
head(corr_2_1$raw_data)
print(corr_2_1$plot)

# 场景2.2: 1 Protein vs 多 Phospho - 单癌种（重点）
corr_2_2 <- cptac_correlation(
  var1 = "AKT1",
  var1_modal = "Protein",
  var1_cancers = "BRCA",
  var2 = c("AKT1", "MTOR", "RPS6"),
  var2_modal = "Phospho",
  var2_cancers = "BRCA"
)
head(corr_2_2$stats)
head(corr_2_2$raw_data)
print(corr_2_2$plot)

# 场景2.3: 1 RNAseq vs 多 Phospho - 单癌种（重点）
corr_2_3 <- cptac_correlation(
  var1 = "AKT1",
  var1_modal = "RNAseq",
  var1_cancers = "LUAD",
  var2 = c("AKT1", "EGFR", "SRC"),
  var2_modal = "Phospho",
  var2_cancers = "LUAD"
)
head(corr_2_3$stats)
head(corr_2_3$raw_data)
print(corr_2_3$plot)

# 场景2.4: 1 RNAseq vs 多 Protein - 多癌种（DotPlot）
corr_2_4 <- cptac_correlation(
  var1 = "TP53",
  var1_modal = "RNAseq",
  var1_cancers = c("BRCA", "LUAD"),
  var2 = c("AKT1", "MTOR"),
  var2_modal = "Protein",
  var2_cancers = c("BRCA", "LUAD")
)
head(corr_2_4$stats)
head(corr_2_4$raw_data)
print(corr_2_4$plot)

# ==================== 场景3: multiple continuous vs multiple continuous ====================

# 场景3.1: 多 RNAseq vs 多 Protein - 单癌种（DotPlot）
corr_3_1 <- cptac_correlation(
  var1 = c("TP53", "EGFR"),
  var1_modal = "RNAseq",
  var1_cancers = "BRCA",
  var2 = c("AKT1", "MTOR"),
  var2_modal = "Protein",
  var2_cancers = "BRCA"
)
head(corr_3_1$stats)
head(corr_3_1$raw_data)
print(corr_3_1$plot)

# 场景3.2: 多 Protein vs 多 Phospho - 单癌种（重点）
corr_3_2 <- cptac_correlation(
  var1 = c("AKT1", "MTOR", "PTEN"),
  var1_modal = "Protein",
  var1_cancers = "BRCA",
  var2 = c("AKT1", "MTOR", "RPS6"),
  var2_modal = "Phospho",
  var2_cancers = "BRCA"
)
head(corr_3_2$stats)
head(corr_3_2$raw_data)
print(corr_3_2$plot)

# 场景3.3: 多 RNAseq vs 多 Phospho - 单癌种（重点）
corr_3_3 <- cptac_correlation(
  var1 = c("AKT1", "EGFR"),
  var1_modal = "RNAseq",
  var1_cancers = "LUAD",
  var2 = c("AKT1", "EGFR", "SRC"),
  var2_modal = "Phospho",
  var2_cancers = "LUAD"
)
head(corr_3_3$stats)
head(corr_3_3$raw_data)
print(corr_3_3$plot)

# 场景3.4: 多 RNAseq vs 多 Protein - 多癌种（DotPlot）
corr_3_4 <- cptac_correlation(
  var1 = c("TP53", "EGFR"),
  var1_modal = "RNAseq",
  var1_cancers = c("BRCA", "LUAD"),
  var2 = c("AKT1", "MTOR"),
  var2_modal = "Protein",
  var2_cancers = c("BRCA", "LUAD")
)
head(corr_3_4$stats)
head(corr_3_4$raw_data)
print(corr_3_4$plot)

# 场景3.5: 同基因Phospho自相关 - 去掉对角线
corr_3_5 <- cptac_correlation(
  var1 = "AKT1",
  var1_modal = "Phospho",
  var1_cancers = "BRCA",
  var2 = "AKT1",
  var2_modal = "Phospho",
  var2_cancers = "BRCA"
)
head(corr_3_5$stats)
head(corr_3_5$raw_data)
print(corr_3_5$plot)

# ==================== 场景4: 1 categorical vs 1 continuous (BoxPlot) ====================

# 场景4.1: Mutation vs RNAseq（重点）
corr_4_1 <- cptac_correlation(
  var1 = "KRAS",
  var1_modal = "Mutation",
  var1_cancers = "LUAD",
  var2 = "EGFR",
  var2_modal = "RNAseq",
  var2_cancers = "LUAD"
)
head(corr_4_1$stats)
head(corr_4_1$raw_data)
print(corr_4_1$plot)

# 场景4.2: Mutation vs Protein（重点）
corr_4_2 <- cptac_correlation(
  var1 = "TP53",
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  var2 = "AKT1",
  var2_modal = "Protein",
  var2_cancers = "BRCA"
)
head(corr_4_2$stats)
head(corr_4_2$raw_data)
print(corr_4_2$plot)

# 场景4.3: Mutation vs Phospho（重点）
corr_4_3 <- cptac_correlation(
  var1 = "PIK3CA",
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  var2 = "AKT1",
  var2_modal = "Phospho",
  var2_cancers = "BRCA"
)
head(corr_4_3$stats)
head(corr_4_3$raw_data)
print(corr_4_3$plot)

# 场景4.4: Clinical vs RNAseq（重点）
corr_4_4 <- cptac_correlation(
  var1 = "Tumor_Stage",
  var1_modal = "Clinical",
  var1_cancers = "BRCA",
  var2 = "TP53",
  var2_modal = "RNAseq",
  var2_cancers = "BRCA"
)
head(corr_4_4$stats)
head(corr_4_4$raw_data)
print(corr_4_4$plot)

# 场景4.5: Clinical vs Phospho（重点）
corr_4_5 <- cptac_correlation(
  var1 = "Tumor_Stage",
  var1_modal = "Clinical",
  var1_cancers = "LUAD",
  var2 = "AKT1",
  var2_modal = "Phospho",
  var2_cancers = "LUAD"
)
head(corr_4_5$stats)
head(corr_4_5$raw_data)
print(corr_4_5$plot)

# 场景4.6: Clinical vs Protein
corr_4_6 <- cptac_correlation(
  var1 = "Age",
  var1_modal = "Clinical",
  var1_cancers = "BRCA",
  var2 = "MTOR",
  var2_modal = "Protein",
  var2_cancers = "BRCA"
)
head(corr_4_6$stats)
head(corr_4_6$raw_data)
print(corr_4_6$plot)

# 场景4.7: Mutation vs RNAseq - 多癌种
corr_4_7 <- cptac_correlation(
  var1 = "KRAS",
  var1_modal = "Mutation",
  var1_cancers = c("LUAD", "COAD"),
  var2 = "EGFR",
  var2_modal = "RNAseq",
  var2_cancers = c("LUAD", "COAD")
)
head(corr_4_7$stats)
head(corr_4_7$raw_data)
print(corr_4_7$plot)

# ==================== 场景5: 1 continuous vs multiple categorical ====================

# 场景5.1: 1 RNAseq vs 多 Mutation（重点）
corr_5_1 <- cptac_correlation(
  var1 = "EGFR",
  var1_modal = "RNAseq",
  var1_cancers = "LUAD",
  var2 = c("KRAS", "EGFR", "TP53"),
  var2_modal = "Mutation",
  var2_cancers = "LUAD"
)
head(corr_5_1$stats)
head(corr_5_1$raw_data)
print(corr_5_1$plot)

# 场景5.2: 1 Protein vs 多 Mutation（重点）
corr_5_2 <- cptac_correlation(
  var1 = "AKT1",
  var1_modal = "Protein",
  var1_cancers = "BRCA",
  var2 = c("PIK3CA", "TP53", "GATA3"),
  var2_modal = "Mutation",
  var2_cancers = "BRCA"
)
head(corr_5_2$stats)
head(corr_5_2$raw_data)
print(corr_5_2$plot)

# 场景5.3: 1 Phospho vs 多 Mutation（重点）
corr_5_3 <- cptac_correlation(
  var1 = "AKT1",
  var1_modal = "Phospho",
  var1_cancers = "BRCA",
  var2 = c("PIK3CA", "TP53"),
  var2_modal = "Mutation",
  var2_cancers = "BRCA"
)
head(corr_5_3$stats)
head(corr_5_3$raw_data)
print(corr_5_3$plot)

# 场景5.4: 1 RNAseq vs 多 Clinical（重点）
corr_5_4 <- cptac_correlation(
  var1 = "TP53",
  var1_modal = "RNAseq",
  var1_cancers = "LUAD",
  var2 = c("Tumor_Stage", "Age"),
  var2_modal = "Clinical",
  var2_cancers = "LUAD"
)
head(corr_5_4$stats)
head(corr_5_4$raw_data)
print(corr_5_4$plot)

# ==================== 场景6: multiple continuous vs 1 categorical ====================

# 场景6.1: 多 RNAseq vs 1 Mutation（重点）
corr_6_1 <- cptac_correlation(
  var1 = c("EGFR", "AKT1", "MTOR"),
  var1_modal = "RNAseq",
  var1_cancers = "LUAD",
  var2 = "KRAS",
  var2_modal = "Mutation",
  var2_cancers = "LUAD"
)
head(corr_6_1$stats)
head(corr_6_1$raw_data)
print(corr_6_1$plot)

# 场景6.2: 多 Protein vs 1 Mutation（重点）
corr_6_2 <- cptac_correlation(
  var1 = c("AKT1", "MTOR", "PTEN"),
  var1_modal = "Protein",
  var1_cancers = "BRCA",
  var2 = "PIK3CA",
  var2_modal = "Mutation",
  var2_cancers = "BRCA"
)
head(corr_6_2$stats)
head(corr_6_2$raw_data)
print(corr_6_2$plot)

# 场景6.3: 多 Phospho vs 1 Mutation（重点）
corr_6_3 <- cptac_correlation(
  var1 = c("AKT1", "MTOR", "RPS6"),
  var1_modal = "Phospho",
  var1_cancers = "BRCA",
  var2 = "PIK3CA",
  var2_modal = "Mutation",
  var2_cancers = "BRCA"
)
head(corr_6_3$stats)
head(corr_6_3$raw_data)
print(corr_6_3$plot)

# 场景6.4: 多 RNAseq vs 1 Clinical（重点）
corr_6_4 <- cptac_correlation(
  var1 = c("TP53", "EGFR", "KRAS"),
  var1_modal = "RNAseq",
  var1_cancers = "LUAD",
  var2 = "Tumor_Stage",
  var2_modal = "Clinical",
  var2_cancers = "LUAD"
)
head(corr_6_4$stats)
head(corr_6_4$raw_data)
print(corr_6_4$plot)

# 场景6.5: 多 Phospho vs 1 Clinical（重点）
corr_6_5 <- cptac_correlation(
  var1 = c("AKT1", "EGFR", "SRC"),
  var1_modal = "Phospho",
  var1_cancers = "LUAD",
  var2 = "Tumor_Stage",
  var2_modal = "Clinical",
  var2_cancers = "LUAD"
)
head(corr_6_5$stats)
head(corr_6_5$raw_data)
print(corr_6_5$plot)

# ==================== 场景7: categorical vs categorical ====================

# 场景7.1: 1 Mutation vs 1 Mutation - 百分比柱状图（重点）
corr_7_1 <- cptac_correlation(
  var1 = "KRAS",
  var1_modal = "Mutation",
  var1_cancers = "LUAD",
  var2 = "EGFR",
  var2_modal = "Mutation",
  var2_cancers = "LUAD"
)
head(corr_7_1$stats)
head(corr_7_1$raw_data)
print(corr_7_1$plot)

# 场景7.2: 1 Mutation vs 多 Mutation - 多个百分比柱状图（重点）
corr_7_2 <- cptac_correlation(
  var1 = c("KRAS", "EGFR"),
  var1_modal = "Mutation",
  var1_cancers = "LUAD",
  var2 = "TP53",
  var2_modal = "Mutation",
  var2_cancers = "LUAD"
)
head(corr_7_2$stats)
head(corr_7_2$raw_data)
print(corr_7_2$plot)

# 场景7.3: 多 Mutation vs 多 Mutation - Heatmap（重点：共突变和互斥）
corr_7_3 <- cptac_correlation(
  var1 = c("KRAS", "EGFR", "ALK", "BRAF"),
  var1_modal = "Mutation",
  var1_cancers = "LUAD",
  var2 = c("TP53", "STK11", "KEAP1"),
  var2_modal = "Mutation",
  var2_cancers = "LUAD"
)
head(corr_7_3$stats)
head(corr_7_3$raw_data)
print(corr_7_3$plot)

# 场景7.4: 1 Clinical vs 1 Mutation - 百分比柱状图（重点）
corr_7_4 <- cptac_correlation(
  var1 = "Tumor_Stage",
  var1_modal = "Clinical",
  var1_cancers = "BRCA",
  var2 = "PIK3CA",
  var2_modal = "Mutation",
  var2_cancers = "BRCA"
)
head(corr_7_4$stats)
head(corr_7_4$raw_data)
print(corr_7_4$plot)

# 场景7.5: 多 Clinical vs 多 Mutation - Heatmap
corr_7_5 <- cptac_correlation(
  var1 = c("Tumor_Stage", "Age"),
  var1_modal = "Clinical",
  var1_cancers = "LUAD",
  var2 = c("KRAS", "EGFR", "TP53"),
  var2_modal = "Mutation",
  var2_cancers = "LUAD"
)
head(corr_7_5$stats)
head(corr_7_5$raw_data)
print(corr_7_5$plot)

################################################################################
# 第三部分：富集分析（cptac_enrichment）- 场景8-15
################################################################################

# ==================== 场景8: 1 categorical vs genome-wide ====================

# 场景8.1: Mutation vs Protein genome - NetworkPlot（重点）
enrich_8_1 <- cptac_enrichment(
  var1 = "KRAS",
  var1_modal = "Mutation",
  var1_cancers = "LUAD",
  analysis_type = "genome",
  genome_modal = "Protein",
  top_n = 30
)
head(enrich_8_1$stats)
head(enrich_8_1$raw_data)
print(enrich_8_1$plot)

# 场景8.2: Mutation vs RNAseq genome（重点）
enrich_8_2 <- cptac_enrichment(
  var1 = "PIK3CA",
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  analysis_type = "genome",
  genome_modal = "RNAseq",
  top_n = 30
)
head(enrich_8_2$stats)
head(enrich_8_2$raw_data)
print(enrich_8_2$plot)

# 场景8.3: Mutation vs Phospho genome（重点）
enrich_8_3 <- cptac_enrichment(
  var1 = "TP53",
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  analysis_type = "genome",
  genome_modal = "Phospho",
  top_n = 30
)
head(enrich_8_3$stats)
head(enrich_8_3$raw_data)
print(enrich_8_3$plot)

# ==================== 场景9: 1 categorical vs enrichment ====================

# 场景9.1: Mutation vs MsigDB Hallmark（默认，重点）
enrich_9_1 <- cptac_enrichment(
  var1 = "KRAS",
  var1_modal = "Mutation",
  var1_cancers = "LUAD",
  analysis_type = "enrichment",
  top_n = 20
)
head(enrich_9_1$stats)
head(enrich_9_1$raw_data)
print(enrich_9_1$plot)

# 场景9.2: Mutation vs GO BP（重点）
enrich_9_2 <- cptac_enrichment(
  var1 = "PIK3CA",
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "BP",
  top_n = 20
)
head(enrich_9_2$stats)
head(enrich_9_2$raw_data)
print(enrich_9_2$plot)

# 场景9.3: Mutation vs KEGG（重点）
enrich_9_3 <- cptac_enrichment(
  var1 = "TP53",
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "KEGG",
  top_n = 20
)
head(enrich_9_3$stats)
head(enrich_9_3$raw_data)
print(enrich_9_3$plot)

# 场景9.4: Mutation vs Reactome（重点）
enrich_9_4 <- cptac_enrichment(
  var1 = "PIK3CA",
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "Reactome",
  top_n = 20
)
head(enrich_9_4$stats)
head(enrich_9_4$raw_data)
print(enrich_9_4$plot)

# ==================== 场景10: multiple categorical vs genome-wide ====================

# 场景10.1: 多 Mutation vs Protein genome - DotPlot Paired（重点）
enrich_10_1 <- cptac_enrichment(
  var1 = c("KRAS", "EGFR", "TP53"),
  var1_modal = "Mutation",
  var1_cancers = "LUAD",
  analysis_type = "genome",
  genome_modal = "Protein",
  top_n = 30
)
head(enrich_10_1$stats)
head(enrich_10_1$raw_data[[1]])
print(enrich_10_1$plot)

# 场景10.2: 多 Mutation vs RNAseq genome（重点）
enrich_10_2 <- cptac_enrichment(
  var1 = c("PIK3CA", "TP53"),
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  analysis_type = "genome",
  genome_modal = "RNAseq",
  top_n = 30
)
head(enrich_10_2$stats)
head(enrich_10_2$raw_data[[1]])
print(enrich_10_2$plot)

# 场景10.3: 多 Mutation vs Phospho genome（重点）
enrich_10_3 <- cptac_enrichment(
  var1 = c("PIK3CA", "TP53", "GATA3"),
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  analysis_type = "genome",
  genome_modal = "Phospho",
  top_n = 30
)
head(enrich_10_3$stats)
head(enrich_10_3$raw_data[[1]])
print(enrich_10_3$plot)

# ==================== 场景11: multiple categorical vs enrichment ====================

# 场景11.1: 多 Mutation vs MsigDB（重点）
enrich_11_1 <- cptac_enrichment(
  var1 = c("KRAS", "EGFR"),
  var1_modal = "Mutation",
  var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  top_n = 15
)
head(enrich_11_1$stats)
head(enrich_11_1$raw_data[[1]])
print(enrich_11_1$plot)

# 场景11.2: 多 Mutation vs GO BP（重点）
enrich_11_2 <- cptac_enrichment(
  var1 = c("PIK3CA", "TP53"),
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "BP",
  top_n = 15
)
head(enrich_11_2$stats)
head(enrich_11_2$raw_data[[1]])
print(enrich_11_2$plot)

# 场景11.3: 多 Mutation vs KEGG（重点）
enrich_11_3 <- cptac_enrichment(
  var1 = c("KRAS", "EGFR", "TP53"),
  var1_modal = "Mutation",
  var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "KEGG",
  top_n = 15
)
head(enrich_11_3$stats)
head(enrich_11_3$raw_data[[1]])
print(enrich_11_3$plot)

# 场景11.4: 多 Mutation vs Reactome（重点）
enrich_11_4 <- cptac_enrichment(
  var1 = c("PIK3CA", "TP53"),
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "Reactome",
  top_n = 15
)
head(enrich_11_4$stats)
head(enrich_11_4$raw_data[[1]])
print(enrich_11_4$plot)

# ==================== 场景12: 1 continuous vs genome-wide ====================

# 场景12.1: 1 RNAseq vs Protein genome - NetworkPlot
enrich_12_1 <- cptac_enrichment(
  var1 = "TP53",
  var1_modal = "RNAseq",
  var1_cancers = "BRCA",
  analysis_type = "genome",
  genome_modal = "Protein",
  method = "pearson",
  top_n = 30
)
head(enrich_12_1$stats)
head(enrich_12_1$raw_data)
print(enrich_12_1$plot)

# 场景12.2: 1 Protein vs Phospho genome（重点）
enrich_12_2 <- cptac_enrichment(
  var1 = "AKT1",
  var1_modal = "Protein",
  var1_cancers = "BRCA",
  analysis_type = "genome",
  genome_modal = "Phospho",
  method = "pearson",
  top_n = 30
)
head(enrich_12_2$stats)
head(enrich_12_2$raw_data)
print(enrich_12_2$plot)

# 场景12.3: 1 Phospho vs Protein genome（重点）
enrich_12_3 <- cptac_enrichment(
  var1 = "MTOR",
  var1_modal = "Phospho",
  var1_cancers = "LUAD",
  analysis_type = "genome",
  genome_modal = "Protein",
  method = "spearman",
  top_n = 30
)
head(enrich_12_3$stats)
head(enrich_12_3$raw_data)
print(enrich_12_3$plot)

# 场景12.4: 1 Phospho vs RNAseq genome（重点）
enrich_12_4 <- cptac_enrichment(
  var1 = "AKT1",
  var1_modal = "Phospho",
  var1_cancers = "BRCA",
  analysis_type = "genome",
  genome_modal = "RNAseq",
  method = "pearson",
  top_n = 30
)
head(enrich_12_4$stats)
head(enrich_12_4$raw_data)
print(enrich_12_4$plot)

# ==================== 场景13: 1 continuous vs enrichment ====================

# 场景13.1: 1 RNAseq vs MsigDB
enrich_13_1 <- cptac_enrichment(
  var1 = "TP53",
  var1_modal = "RNAseq",
  var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  method = "pearson",
  top_n = 15
)
head(enrich_13_1$stats)
head(enrich_13_1$raw_data)
print(enrich_13_1$plot)

# 场景13.2: 1 Protein vs GO BP
enrich_13_2 <- cptac_enrichment(
  var1 = "AKT1",
  var1_modal = "Protein",
  var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "BP",
  method = "spearman",
  top_n = 15
)
head(enrich_13_2$stats)
head(enrich_13_2$raw_data)
print(enrich_13_2$plot)

# 场景13.3: 1 Phospho vs KEGG（重点）
enrich_13_3 <- cptac_enrichment(
  var1 = "AKT1",
  var1_modal = "Phospho",
  var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "KEGG",
  method = "pearson",
  top_n = 15
)
head(enrich_13_3$stats)
head(enrich_13_3$raw_data)
print(enrich_13_3$plot)

# 场景13.4: 1 Phospho vs Reactome（重点）
enrich_13_4 <- cptac_enrichment(
  var1 = "MTOR",
  var1_modal = "Phospho",
  var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "Reactome",
  method = "pearson",
  top_n = 15
)
head(enrich_13_4$stats)
head(enrich_13_4$raw_data)
print(enrich_13_4$plot)

# 场景13.5: 1 Protein vs Reactome
enrich_13_5 <- cptac_enrichment(
  var1 = "TP53",
  var1_modal = "Protein",
  var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "Reactome",
  method = "pearson",
  top_n = 15
)
head(enrich_13_5$stats)
head(enrich_13_5$raw_data)
print(enrich_13_5$plot)

# ==================== 场景14: multiple continuous vs genome-wide ====================

# 场景14.1: 多 RNAseq vs Protein genome - DotPlot Paired
enrich_14_1 <- cptac_enrichment(
  var1 = c("TP53", "EGFR"),
  var1_modal = "RNAseq",
  var1_cancers = "BRCA",
  analysis_type = "genome",
  genome_modal = "Protein",
  method = "pearson",
  top_n = 30
)
head(enrich_14_1$stats)
head(enrich_14_1$raw_data[[1]])
print(enrich_14_1$plot)

# 场景14.2: 多 Protein vs Phospho genome（重点）
enrich_14_2 <- cptac_enrichment(
  var1 = c("AKT1", "MTOR"),
  var1_modal = "Protein",
  var1_cancers = "BRCA",
  analysis_type = "genome",
  genome_modal = "Phospho",
  method = "pearson",
  top_n = 30
)
head(enrich_14_2$stats)
head(enrich_14_2$raw_data[[1]])
print(enrich_14_2$plot)

# 场景14.3: 多 Phospho vs Protein genome（重点）
enrich_14_3 <- cptac_enrichment(
  var1 = c("AKT1", "EGFR", "SRC"),
  var1_modal = "Phospho",
  var1_cancers = "LUAD",
  analysis_type = "genome",
  genome_modal = "Protein",
  method = "spearman",
  top_n = 30
)
head(enrich_14_3$stats)
head(enrich_14_3$raw_data[[1]])
print(enrich_14_3$plot)

# 场景14.4: 多 Phospho vs RNAseq genome（重点）
enrich_14_4 <- cptac_enrichment(
  var1 = c("AKT1", "MTOR"),
  var1_modal = "Phospho",
  var1_cancers = "BRCA",
  analysis_type = "genome",
  genome_modal = "RNAseq",
  method = "pearson",
  top_n = 30
)
head(enrich_14_4$stats)
head(enrich_14_4$raw_data[[1]])
print(enrich_14_4$plot)

# ==================== 场景15: multiple continuous vs enrichment ====================

# 场景15.1: 多 RNAseq vs MsigDB - GSEA Matrix
enrich_15_1 <- cptac_enrichment(
  var1 = c("TP53", "EGFR"),
  var1_modal = "RNAseq",
  var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  method = "pearson",
  top_n = 15
)
head(enrich_15_1$stats)
head(enrich_15_1$raw_data[[1]])
print(enrich_15_1$plot)

# 场景15.2: 多 Protein vs GO BP
enrich_15_2 <- cptac_enrichment(
  var1 = c("AKT1", "MTOR", "PTEN"),
  var1_modal = "Protein",
  var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "GO",
  enrich_ont = "BP",
  method = "spearman",
  top_n = 15
)
head(enrich_15_2$stats)
head(enrich_15_2$raw_data[[1]])
print(enrich_15_2$plot)

# 场景15.3: 多 Phospho vs KEGG（重点）
enrich_15_3 <- cptac_enrichment(
  var1 = c("AKT1", "MTOR"),
  var1_modal = "Phospho",
  var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "KEGG",
  method = "pearson",
  top_n = 15
)
head(enrich_15_3$stats)
head(enrich_15_3$raw_data[[1]])
print(enrich_15_3$plot)

# 场景15.4: 多 Phospho vs Reactome（重点）
enrich_15_4 <- cptac_enrichment(
  var1 = c("AKT1", "EGFR", "SRC"),
  var1_modal = "Phospho",
  var1_cancers = "LUAD",
  analysis_type = "enrichment",
  enrich_database = "Reactome",
  method = "pearson",
  top_n = 15
)
head(enrich_15_4$stats)
head(enrich_15_4$raw_data[[1]])
print(enrich_15_4$plot)

################################################################################
# 第四部分：生存分析（cptac_survival）- 场景16-17
################################################################################

# ==================== 场景16: 1 feature vs survival (KM + Cox) ====================

# 场景16.1: RNAseq vs OS - 单癌种
surv_16_1 <- cptac_survival(
  var1 = "TP53",
  var1_modal = "RNAseq",
  var1_cancers = "BRCA",
  surv_type = "OS",
  cutoff_type = "optimal"
)
head(surv_16_1$stats)
head(surv_16_1$raw_data)
print(surv_16_1$plot)

# 场景16.2: Protein vs OS - 单癌种
surv_16_2 <- cptac_survival(
  var1 = "AKT1",
  var1_modal = "Protein",
  var1_cancers = "LUAD",
  surv_type = "OS",
  cutoff_type = "median"
)
head(surv_16_2$stats)
head(surv_16_2$raw_data)
print(surv_16_2$plot)

# 场景16.3: Phospho vs OS - 单癌种（重点）
surv_16_3 <- cptac_survival(
  var1 = "AKT1",
  var1_modal = "Phospho",
  var1_cancers = "BRCA",
  surv_type = "OS",
  cutoff_type = "optimal"
)
head(surv_16_3$stats)
head(surv_16_3$raw_data)
print(surv_16_3$plot)

# 场景16.4: Phospho vs PFS - 单癌种（重点）
surv_16_4 <- cptac_survival(
  var1 = "MTOR",
  var1_modal = "Phospho",
  var1_cancers = "LUAD",
  surv_type = "PFS",
  cutoff_type = "median"
)
head(surv_16_4$stats)
head(surv_16_4$raw_data)
print(surv_16_4$plot)

# 场景16.5: Mutation vs OS - 单癌种（重点）
surv_16_5 <- cptac_survival(
  var1 = "KRAS",
  var1_modal = "Mutation",
  var1_cancers = "LUAD",
  surv_type = "OS"
)
head(surv_16_5$stats)
head(surv_16_5$raw_data)
print(surv_16_5$plot)

# 场景16.6: Mutation vs OS - 单癌种（重点）
surv_16_6 <- cptac_survival(
  var1 = "PIK3CA",
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  surv_type = "OS"
)
head(surv_16_6$stats)
head(surv_16_6$raw_data)
print(surv_16_6$plot)

# 场景16.7: Clinical vs OS - 单癌种（重点）
surv_16_7 <- cptac_survival(
  var1 = "Tumor_Stage",
  var1_modal = "Clinical",
  var1_cancers = "BRCA",
  surv_type = "OS"
)
head(surv_16_7$stats)
head(surv_16_7$raw_data)
print(surv_16_7$plot)

# 场景16.8: Clinical vs PFS - 单癌种（重点）
surv_16_8 <- cptac_survival(
  var1 = "Age",
  var1_modal = "Clinical",
  var1_cancers = "LUAD",
  surv_type = "PFS",
  cutoff_type = "median"
)
head(surv_16_8$stats)
head(surv_16_8$raw_data)
print(surv_16_8$plot)

# ==================== 场景17: multiple features vs survival (Forest Plot) ====================

# 场景17.1: 多 RNAseq vs OS - 单癌种
surv_17_1 <- cptac_survival(
  var1 = c("TP53", "EGFR", "KRAS"),
  var1_modal = "RNAseq",
  var1_cancers = "LUAD",
  surv_type = "OS",
  cutoff_type = "optimal"
)
head(surv_17_1$stats)
head(surv_17_1$raw_data)
print(surv_17_1$plot)

# 场景17.2: 多 Protein vs OS - 单癌种
surv_17_2 <- cptac_survival(
  var1 = c("AKT1", "MTOR", "PTEN"),
  var1_modal = "Protein",
  var1_cancers = "BRCA",
  surv_type = "OS",
  cutoff_type = "median"
)
head(surv_17_2$stats)
head(surv_17_2$raw_data)
print(surv_17_2$plot)

# 场景17.3: 多 Phospho vs OS - 单癌种（重点：每个phospho site独立）
surv_17_3 <- cptac_survival(
  var1 = c("AKT1", "MTOR"),
  var1_modal = "Phospho",
  var1_cancers = "BRCA",
  surv_type = "OS",
  cutoff_type = "optimal"
)
head(surv_17_3$stats)
head(surv_17_3$raw_data)
print(surv_17_3$plot)

# 场景17.4: 多 Mutation vs OS - 单癌种（重点）
surv_17_4 <- cptac_survival(
  var1 = c("KRAS", "EGFR", "TP53"),
  var1_modal = "Mutation",
  var1_cancers = "LUAD",
  surv_type = "OS"
)
head(surv_17_4$stats)
head(surv_17_4$raw_data)
print(surv_17_4$plot)

# 场景17.5: 多 Mutation vs PFS - 单癌种（重点）
surv_17_5 <- cptac_survival(
  var1 = c("PIK3CA", "TP53"),
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  surv_type = "PFS"
)
head(surv_17_5$stats)
head(surv_17_5$raw_data)
print(surv_17_5$plot)

# 场景17.6: 多 Clinical vs OS - 单癌种（重点）
surv_17_6 <- cptac_survival(
  var1 = c("Age", "BMI"),
  var1_modal = "Clinical",
  var1_cancers = "BRCA",
  surv_type = "OS"
)
head(surv_17_6$stats)
head(surv_17_6$raw_data)
print(surv_17_6$plot)

# 场景17.7: 1基因多癌种 - Forest Plot（新：多癌种→场景17）
surv_17_7 <- cptac_survival(
  var1 = "TP53",
  var1_modal = "RNAseq",
  var1_cancers = c("BRCA", "LUAD"),
  surv_type = "OS",
  cutoff_type = "optimal"
)
head(surv_17_7$stats)
head(surv_17_7$raw_data)
print(surv_17_7$plot)

# 场景17.8: Phospho多癌种 - Forest Plot（重点：每个site×每个癌种）
surv_17_8 <- cptac_survival(
  var1 = "AKT1",
  var1_modal = "Phospho",
  var1_cancers = c("BRCA", "LUAD"),
  surv_type = "OS",
  cutoff_type = "optimal"
)
head(surv_17_8$stats)
head(surv_17_8$raw_data)
print(surv_17_8$plot)

# 场景17.9: 多基因多癌种 - Forest Plot
surv_17_9 <- cptac_survival(
  var1 = c("TP53", "EGFR"),
  var1_modal = "RNAseq",
  var1_cancers = c("BRCA", "LUAD"),
  surv_type = "OS",
  cutoff_type = "optimal"
)
head(surv_17_9$stats)
head(surv_17_9$raw_data)
print(surv_17_9$plot)

################################################################################
# 第五部分：高级专题测试
################################################################################

# ==================== 磷酸化深度测试 ====================

# 专题1: Phospho vs Phospho - 跨5个癌种
phospho_adv_1 <- cptac_correlation(
  var1 = "AKT1",
  var1_modal = "Phospho",
  var1_cancers = c("BRCA", "LUAD", "CCRCC", "UCEC", "PDAC"),
  var2 = "MTOR",
  var2_modal = "Phospho",
  var2_cancers = c("BRCA", "LUAD", "CCRCC", "UCEC", "PDAC")
)
head(phospho_adv_1$stats)
head(phospho_adv_1$raw_data)
print(phospho_adv_1$plot)

# 专题2: 多Phospho基因组扫描
phospho_adv_2 <- cptac_enrichment(
  var1 = c("AKT1", "MTOR", "RPS6"),
  var1_modal = "Phospho",
  var1_cancers = "BRCA",
  analysis_type = "genome",
  genome_modal = "Protein",
  top_n = 50
)
head(phospho_adv_2$stats)
head(phospho_adv_2$raw_data[[1]])
print(phospho_adv_2$plot)

# 专题3: Phospho多癌种生存分析
phospho_adv_3 <- cptac_survival(
  var1 = "AKT1",
  var1_modal = "Phospho",
  var1_cancers = c("BRCA", "LUAD", "CCRCC"),
  surv_type = "OS",
  cutoff_type = "optimal"
)
head(phospho_adv_3$stats)
head(phospho_adv_3$raw_data)
print(phospho_adv_3$plot)

# ==================== 突变深度测试 ====================

# 专题4: 多突变共突变/互斥分析 - Heatmap
mutation_adv_1 <- cptac_correlation(
  var1 = c("KRAS", "EGFR", "ALK", "BRAF"),
  var1_modal = "Mutation",
  var1_cancers = "LUAD",
  var2 = c("TP53", "STK11", "KEAP1"),
  var2_modal = "Mutation",
  var2_cancers = "LUAD"
)
head(mutation_adv_1$stats)
head(mutation_adv_1$raw_data)
print(mutation_adv_1$plot)

# 专题5: 突变对磷酸化的影响
mutation_adv_2 <- cptac_correlation(
  var1 = c("PIK3CA", "TP53"),
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  var2 = c("AKT1", "MTOR", "RPS6"),
  var2_modal = "Phospho",
  var2_cancers = "BRCA"
)
head(mutation_adv_2$stats)
head(mutation_adv_2$raw_data)
print(mutation_adv_2$plot)

# 专题6: 突变对磷酸化组的影响 - DotPlot Paired
mutation_adv_3 <- cptac_enrichment(
  var1 = c("PIK3CA", "TP53"),
  var1_modal = "Mutation",
  var1_cancers = "BRCA",
  analysis_type = "genome",
  genome_modal = "Phospho",
  top_n = 50
)
head(mutation_adv_3$stats)
head(mutation_adv_3$raw_data[[1]])
print(mutation_adv_3$plot)

# 专题7: 多突变多癌种生存分析 - Forest Plot
mutation_adv_4 <- cptac_survival(
  var1 = c("TP53", "KRAS"),
  var1_modal = "Mutation",
  var1_cancers = c("BRCA", "LUAD", "COAD"),
  surv_type = "OS"
)
head(mutation_adv_4$stats)
head(mutation_adv_4$raw_data)
print(mutation_adv_4$plot)

# ==================== 临床深度测试 ====================

# 专题8: 临床变量 vs 磷酸化 - 多个BoxPlot
clinical_adv_1 <- cptac_correlation(
  var1 = c("Age", "Tumor_Stage"),
  var1_modal = "Clinical",
  var1_cancers = "BRCA",
  var2 = c("AKT1", "MTOR"),
  var2_modal = "Phospho",
  var2_cancers = "BRCA"
)
head(clinical_adv_1$stats)
head(clinical_adv_1$raw_data)
print(clinical_adv_1$plot)

# 专题9: 临床变量 vs 突变 - Heatmap
clinical_adv_2 <- cptac_correlation(
  var1 = c("Tumor_Stage", "Age"),
  var1_modal = "Clinical",
  var1_cancers = "LUAD",
  var2 = c("KRAS", "EGFR", "TP53"),
  var2_modal = "Mutation",
  var2_cancers = "LUAD"
)
head(clinical_adv_2$stats)
head(clinical_adv_2$raw_data)
print(clinical_adv_2$plot)

# 专题10: 多临床变量生存分析 - Forest Plot
clinical_adv_3 <- cptac_survival(
  var1 = c("Age", "BMI"),
  var1_modal = "Clinical",
  var1_cancers = c("BRCA", "LUAD"),
  surv_type = "OS"
)
head(clinical_adv_3$stats)
head(clinical_adv_3$raw_data)
print(clinical_adv_3$plot)

# ==================== 跨模态综合测试 ====================

# 专题11: Protein vs Phospho跨多癌种
cross_modal_1 <- cptac_correlation(
  var1 = "AKT1",
  var1_modal = "Protein",
  var1_cancers = c("BRCA", "LUAD", "CCRCC"),
  var2 = "AKT1",
  var2_modal = "Phospho",
  var2_cancers = c("BRCA", "LUAD", "CCRCC")
)
head(cross_modal_1$stats)
head(cross_modal_1$raw_data)
print(cross_modal_1$plot)

# 专题12: RNAseq vs Mutation - 多基因多癌种
cross_modal_2 <- cptac_correlation(
  var1 = c("TP53", "KRAS"),
  var1_modal = "RNAseq",
  var1_cancers = c("BRCA", "LUAD"),
  var2 = c("TP53", "KRAS"),
  var2_modal = "Mutation",
  var2_cancers = c("BRCA", "LUAD")
)
head(cross_modal_2$stats)
head(cross_modal_2$raw_data)
print(cross_modal_2$plot)

# 专题13: 多Protein vs Phospho enrichment - GSEA Matrix
cross_modal_3 <- cptac_enrichment(
  var1 = c("AKT1", "MTOR", "PTEN"),
  var1_modal = "Protein",
  var1_cancers = "BRCA",
  analysis_type = "enrichment",
  enrich_database = "MsigDB",
  method = "pearson",
  top_n = 15
)
head(cross_modal_3$stats)
head(cross_modal_3$raw_data[[1]])
print(cross_modal_3$plot)

# ==================== 大规模测试 ====================

# 专题14: 5个基因 × 5个突变 - Heatmap（测试大矩阵）
large_test_1 <- cptac_correlation(
  var1 = c("KRAS", "EGFR", "ALK", "BRAF", "RET"),
  var1_modal = "Mutation",
  var1_cancers = "LUAD",
  var2 = c("TP53", "STK11", "KEAP1", "NF1", "RBM10"),
  var2_modal = "Mutation",
  var2_cancers = "LUAD"
)
head(large_test_1$stats)
head(large_test_1$raw_data)
print(large_test_1$plot)

# 专题15: 多Phospho vs Protein genome（测试大图）
large_test_2 <- cptac_enrichment(
  var1 = c("AKT1", "MTOR", "RPS6"),
  var1_modal = "Phospho",
  var1_cancers = "BRCA",
  analysis_type = "genome",
  genome_modal = "Protein",
  method = "pearson",
  top_n = 50
)
head(large_test_2$stats)
head(large_test_2$raw_data[[1]])
print(large_test_2$plot)

################################################################################
# 教程完成
################################################################################
