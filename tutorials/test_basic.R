# ==============================================================================
# SLTCGA Basic Test Script
# ==============================================================================
# Quick test to verify core functionality
# ==============================================================================

library(devtools)
load_all()

# Set data path
Sys.setenv(SL_BULK_DATA = "/Users/liuzaoqu/Desktop/develop/DataMiner_Dev/bulk_data")

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║              SLTCGA Basic Functionality Test              ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

# ==============================================================================
# Test 1: Helper Functions
# ==============================================================================

cat("Test 1: Helper Functions\n")
cat("─────────────────────────────────────────────────────────────\n")

cat("\n1.1 list_modalities()\n")
list_modalities()

cat("\n1.2 list_cancer_types()\n")
list_cancer_types()

cat("\n1.3 list_variables(modal='Signature')\n")
list_variables(modal = "Signature")

cat("\n1.4 search_variables('TMB')\n")
search_variables("TMB")

cat("\n✓ Test 1 passed: All helper functions work\n\n")

# ==============================================================================
# Test 2: Data Loading
# ==============================================================================

cat("Test 2: Data Loading\n")
cat("─────────────────────────────────────────────────────────────\n")

cat("\n2.1 Loading RNAseq data\n")
test_rnaseq <- tcga_load_modality(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA"
)
cat(sprintf("  ✓ Loaded: %d samples, %d features\n", 
           nrow(test_rnaseq$data), length(test_rnaseq$var1_features)))

cat("\n2.2 Loading Signature data\n")
test_sig <- tcga_load_modality(
  var1 = "TMB", var1_modal = "Signature", var1_cancers = "BRCA"
)
cat(sprintf("  ✓ Loaded: %d samples, %d features\n", 
           nrow(test_sig$data), length(test_sig$var1_features)))

cat("\n2.3 Loading Mutation data\n")
test_mut <- tcga_load_modality(
  var1 = "TP53", var1_modal = "Mutation", var1_cancers = "LUAD"
)
cat(sprintf("  ✓ Loaded: %d samples, %d features\n", 
           nrow(test_mut$data), length(test_mut$var1_features)))

cat("\n✓ Test 2 passed: Data loading works\n\n")

# ==============================================================================
# Test 3: Correlation Analysis (Scenario 1)
# ==============================================================================

cat("Test 3: Correlation Analysis (Scenario 1)\n")
cat("─────────────────────────────────────────────────────────────\n")

cat("\nRunning: TP53 (RNAseq) vs TMB (Signature) in BRCA\n")

result1 <- tcga_correlation(
  var1 = "TP53", var1_modal = "RNAseq", var1_cancers = "BRCA",
  var2 = "TMB", var2_modal = "Signature", var2_cancers = "BRCA"
)

cat("\nResults:\n")
print(result1$stats)
cat(sprintf("\nPlot dimensions: %.1f × %.1f inches\n", 
           attr(result1$plot, "width"), attr(result1$plot, "height")))

cat("\n✓ Test 3 passed: Correlation analysis works\n\n")

# ==============================================================================
# Summary
# ==============================================================================

cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║                  All Tests Passed! ✓                      ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

cat("SLTCGA is ready to use!\n\n")
cat("Next steps:\n")
cat("  1. Run tutorials/Quick_Start.R for more examples\n")
cat("  2. Try different scenarios (mutation, enrichment, survival)\n")
cat("  3. Explore immune cell analysis\n\n")

