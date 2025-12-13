# ==============================================================================
# Package Startup
# ==============================================================================

.onLoad <- function(libname, pkgname) {
  # Check SL_BULK_DATA environment variable
  bulk_data <- Sys.getenv("SL_BULK_DATA")

  if (bulk_data == "") {
    packageStartupMessage(
      "Note: SL_BULK_DATA environment variable not set.\n",
      "Please set it using:\n",
      "  Sys.setenv(SL_BULK_DATA = '/path/to/bulk_data')\n",
      "Required data files should be in:\n",
      "  - TCGA_Omics_Split/\n",
      "  - TCGA_RNAseq_log2TPM_SplitCancer/\n",
      "  - TCGA_Clinical_Final_SovingLab.qs\n",
      "  - TCGA_miRNA_pancancer_scores.qs\n",
      "  - TCGA_DeconvCell_All_scores.qs"
    )
  } else {
    packageStartupMessage("SLTCGA loaded successfully!\n")
  }
}

#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "\nUse list_modalities() to see available data types.\n\n",
    "╔════════════════════════════════════════════════════════════╗\n",
    "║  8 Modalities | 33 Main + 27 Subtypes + 5 Combined | 17 Scenarios║\n",
    "╚════════════════════════════════════════════════════════════╝\n"
  )
}
