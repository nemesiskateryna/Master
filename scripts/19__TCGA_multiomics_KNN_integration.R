
###############################################################################
# SCRIPT 16. TCGA-SKCM KNN MULTI-OMICS INTEGRATION
#
# Purpose:
# Integrate KNN-imputed DNA methylation, RNA-seq expression
# and clinical data for the TCGA-SKCM discovery cohort.
#
# Methylation: 409,556 CpGs x 460 samples
# RNA-seq:      25,042 genes x 460 samples
# Clinical:        460 samples x 15 variables
#
# Methylation is stored as an HDF5-backed DelayedMatrix.
# The original integrated object is preserved.
###############################################################################

# 1. Load packages
library(HDF5Array)
library(rhdf5)

# 2. Define input and output paths
project_dir <- getwd()

imputed_h5 <- file.path(
  project_dir,
  "data/TCGA/processed",
  "TCGA_SKCM_beta_matrix_460_StandardCpG_20pct_KNN.h5"
)

original_rds <- file.path(
  project_dir,
  "data/TCGA/processed",
  "TCGA_SKCM_integrated_multiomics.rds"
)

output_rds <- file.path(
  project_dir,
  "data/TCGA/processed",
  "TCGA_SKCM_integrated_multiomics_KNN.rds"
)

# 3. Check input files
stopifnot(
  file.exists(imputed_h5),
  file.exists(original_rds)
)

if (file.exists(output_rds)) {
  stop("Output file already exists. No file was overwritten.")
}

# 4. Load existing RNA-seq and clinical data
integrated <- readRDS(original_rds)

rna_matrix <- integrated$expression_matrix
clinical_data <- integrated$clinical_data

rna_samples <- as.character(colnames(rna_matrix))
clinical_samples <- as.character(rownames(clinical_data))

# 5. Read methylation identifiers
knn_cpg_ids <- as.character(
  rhdf5::h5read(imputed_h5, "CpG_ID")
)

knn_sample_ids <- as.character(
  rhdf5::h5read(imputed_h5, "Sample_Barcode")
)

# 6. Verify sample matching
rna_id_16 <- substr(rna_samples, 1, 16)

stopifnot(
  length(knn_sample_ids) == length(rna_samples),
  !anyDuplicated(knn_sample_ids),
  !anyDuplicated(rna_id_16),
  all(substr(knn_sample_ids, 1, 16) == rna_id_16),
  identical(clinical_samples, rna_samples)
)

# 7. Create HDF5-backed methylation matrix
beta_knn <- HDF5Array::HDF5Array(
  imputed_h5,
  "beta"
)

# 8. Verify dimensions and assign identifiers
stopifnot(
  identical(dim(beta_knn), c(409556L, 460L)),
  length(knn_cpg_ids) == nrow(beta_knn),
  length(rna_samples) == ncol(beta_knn),
  !anyDuplicated(knn_cpg_ids)
)

rownames(beta_knn) <- knn_cpg_ids
colnames(beta_knn) <- rna_samples

# 9. Update metadata
knn_metadata <- integrated$metadata

knn_metadata$methylation_imputation <- list(
  method = "Blockwise KNN imputation",
  package = "impute",
  k = 10,
  block_size = 1000,
  imputed_values = 1137342,
  benchmark = "Five repeated stratified masking experiments",
  validation = "Final HDF5 QC: PASS"
)

knn_metadata$methylation_source <- basename(imputed_h5)

knn_metadata$matching_method <-
  "Matching by first 16 characters of TCGA sample barcode"

knn_metadata$clinical_alignment <-
  "Exact full barcode order; verified against RNA-seq"

# 10. Create new multi-omics object
integrated_knn <- list(
  beta_matrix = beta_knn,
  expression_matrix = rna_matrix,
  clinical_data = clinical_data,
  metadata = knn_metadata
)

# 11. Final integration QC
stopifnot(
  identical(dim(integrated_knn$beta_matrix), c(409556L, 460L)),
  identical(dim(integrated_knn$expression_matrix), c(25042L, 460L)),
  identical(dim(integrated_knn$clinical_data), c(460L, 15L)),
  identical(
    colnames(integrated_knn$beta_matrix),
    colnames(integrated_knn$expression_matrix)
  ),
  identical(
    colnames(integrated_knn$expression_matrix),
    rownames(integrated_knn$clinical_data)
  ),
  !anyDuplicated(rownames(integrated_knn$beta_matrix))
)

# 12. Save new object
saveRDS(
  integrated_knn,
  file = output_rds,
  compress = FALSE
)

# 13. Report
cat("\n========================================\n")
cat("TCGA-SKCM MULTI-OMICS INTEGRATION\n")
cat("========================================\n")

cat("Methylation:", paste(dim(beta_knn), collapse = " x "), "\n")
cat("RNA-seq:", paste(dim(rna_matrix), collapse = " x "), "\n")
cat("Clinical:", paste(dim(clinical_data), collapse = " x "), "\n")

cat("Matched samples:", length(rna_samples), "\n")
cat("CpG features:", length(knn_cpg_ids), "\n")

cat("Output file:", output_rds, "\n")
cat("File exists:", file.exists(output_rds), "\n")

cat("\nINTEGRATION QC: PASS\n")
