
###############################################################################
# SCRIPT 17
# TCGA-SKCM: INTEGRATION WITH IMMUNE SCORES
#
# Purpose:
# Integrate transcriptional immune scores with the existing
# TCGA-SKCM KNN-imputed multi-omics dataset.
#
# Input:
#   data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN.rds
#   results/tables/TCGA_SKCM_immune_scores.csv
#
# Output:
#   data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune_portable.rds
#   results/tables/TCGA_SKCM_immune_integration_QC.csv
###############################################################################

# 0. PROJECT DIRECTORY --------------------------------------------------------

setwd("~/master/Master")

# 1. LOAD INPUT DATA -----------------------------------------------------------

integrated_knn <- readRDS(
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN.rds"
)

immune_scores <- read.csv(
  "results/tables/TCGA_SKCM_immune_scores.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n--- INPUT DIMENSIONS ---\n")

cat(
  "Methylation:",
  paste(dim(integrated_knn$beta_matrix), collapse = " x "),
  "\n"
)

cat(
  "Expression:",
  paste(dim(integrated_knn$expression_matrix), collapse = " x "),
  "\n"
)

cat(
  "Clinical:",
  paste(dim(integrated_knn$clinical_data), collapse = " x "),
  "\n"
)

cat(
  "Immune scores:",
  paste(dim(immune_scores), collapse = " x "),
  "\n"
)

# 2. VALIDATE IMMUNE SCORE DATA -----------------------------------------------

required_columns <- c(
  "Sample_ID",
  "TIS_mean_zscore",
  "CYT_mean_zscore",
  "CD8_mean_zscore",
  "IFNG_mean_zscore",
  "immune_group",
  "sample_type"
)

stopifnot(
  all(required_columns %in% colnames(immune_scores))
)

stopifnot(
  !anyDuplicated(immune_scores$Sample_ID)
)

stopifnot(
  !anyNA(immune_scores[, required_columns])
)

# 3. VALIDATE SAMPLE IDENTIFIERS -----------------------------------------------

rna_ids <- colnames(integrated_knn$expression_matrix)
immune_ids <- immune_scores$Sample_ID

stopifnot(!is.null(rna_ids))
stopifnot(!anyDuplicated(rna_ids))

# Every immune-score sample must have an RNA-seq counterpart.
immune_match <- match(rna_ids, immune_ids)

stopifnot(!anyNA(immune_match))

# Reorder immune scores to exactly match RNA-seq sample order.
immune_scores_aligned <- immune_scores[immune_match, , drop = FALSE]

stopifnot(
  identical(
    as.character(immune_scores_aligned$Sample_ID),
    as.character(rna_ids)
  )
)

# 4. VALIDATE METHYLATION SAMPLE ALIGNMENT -------------------------------------

beta_ids <- colnames(integrated_knn$beta_matrix)

stopifnot(!is.null(beta_ids))
stopifnot(length(beta_ids) == length(rna_ids))

# HDF5 sample IDs may be 16-character TCGA barcodes.
# Verify the first 16 characters before assigning full RNA barcodes.
stopifnot(
  identical(
    substr(beta_ids, 1, 16),
    substr(rna_ids, 1, 16)
  )
)

# Assign full RNA sample IDs only after verified positional matching.
colnames(integrated_knn$beta_matrix) <- rna_ids

# 5. VALIDATE CLINICAL ALIGNMENT -----------------------------------------------

clinical_ids <- rownames(integrated_knn$clinical_data)

stopifnot(!is.null(clinical_ids))

stopifnot(
  identical(
    as.character(clinical_ids),
    as.character(rna_ids)
  )
)

# 6. CHECK IMMUNE PHENOTYPE -----------------------------------------------------

immune_scores_aligned$immune_group <- factor(
  immune_scores_aligned$immune_group,
  levels = c("Immune-low", "Immune-high")
)

group_distribution <- table(
  immune_scores_aligned$immune_group
)

cat("\n--- IMMUNE GROUP DISTRIBUTION ---\n")
print(group_distribution)

stopifnot(
  all(group_distribution > 0)
)

# Confirm that the labels correspond to the cohort-median TIS definition.
tis_threshold <- median(
  immune_scores_aligned$TIS_mean_zscore
)

expected_group <- ifelse(
  immune_scores_aligned$TIS_mean_zscore > tis_threshold,
  "Immune-high",
  "Immune-low"
)

stopifnot(
  identical(
    as.character(immune_scores_aligned$immune_group),
    expected_group
  )
)

# 7. INTEGRATE IMMUNE SCORES ----------------------------------------------------

integrated_knn_immune <- integrated_knn

integrated_knn_immune$immune_scores <- immune_scores_aligned

integrated_knn_immune$metadata$immune_phenotype <- list(
  primary_score = "TIS_mean_zscore",
  threshold_method = "Cohort median",
  threshold_value = tis_threshold,
  group_definition = paste(
    "Immune-high: TIS above cohort median;",
    "Immune-low: TIS at or below cohort median"
  ),
  group_counts = as.list(group_distribution)
)

# 8. FINAL ALIGNMENT QC --------------------------------------------------------

qc_results <- data.frame(
  Check = c(
    "RNA sample count",
    "Immune score sample count",
    "Unique RNA sample IDs",
    "Unique immune sample IDs",
    "Immune/RNA alignment",
    "Methylation/RNA barcode alignment",
    "Clinical/RNA alignment",
    "Missing immune scores",
    "Immune group balance",
    "TIS label definition"
  ),
  Result = c(
    ncol(integrated_knn_immune$expression_matrix),
    nrow(integrated_knn_immune$immune_scores),
    length(unique(rna_ids)),
    length(unique(immune_ids)),
    identical(
      integrated_knn_immune$immune_scores$Sample_ID,
      rna_ids
    ),
    identical(
      substr(colnames(integrated_knn_immune$beta_matrix), 1, 16),
      substr(rna_ids, 1, 16)
    ),
    identical(
      rownames(integrated_knn_immune$clinical_data),
      rna_ids
    ),
    sum(is.na(integrated_knn_immune$immune_scores)),
    paste(group_distribution, collapse = " / "),
    identical(
      as.character(integrated_knn_immune$immune_scores$immune_group),
      expected_group
    )
  ),
  stringsAsFactors = FALSE
)

cat("\n--- FINAL INTEGRATION QC ---\n")
print(qc_results, row.names = FALSE)

# 9. STOP IF ALIGNMENT FAILS ----------------------------------------------------

alignment_pass <- (
  identical(
    integrated_knn_immune$immune_scores$Sample_ID,
    rna_ids
  ) &&
    identical(
      substr(colnames(integrated_knn_immune$beta_matrix), 1, 16),
      substr(rna_ids, 1, 16)
    ) &&
    identical(
      rownames(integrated_knn_immune$clinical_data),
      rna_ids
    ) &&
    !anyNA(integrated_knn_immune$immune_scores)
)

stopifnot(alignment_pass)

# 10. SAVE QC REPORT ------------------------------------------------------------

write.csv(
  qc_results,
  "results/tables/TCGA_SKCM_immune_integration_QC.csv",
  row.names = FALSE
)

# 11. SAVE PORTABLE OBJECT ------------------------------------------------------

# Keep the large HDF5 matrix in its separate HDF5 file.
# Store its path in the portable RDS instead of serializing the
# out-of-memory matrix object.

integrated_knn_immune_portable <- integrated_knn_immune

integrated_knn_immune_portable$beta_matrix <- NULL

integrated_knn_immune_portable$beta_matrix_reference <-
  "data/TCGA/processed/TCGA_SKCM_beta_matrix_460_StandardCpG_20pct_KNN.h5"

saveRDS(
  integrated_knn_immune_portable,
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune_portable.rds"
)

cat("\n--- FINAL STATUS ---\n")
cat("INTEGRATION QC:", if (alignment_pass) "PASS" else "FAIL", "\n")
cat("Portable object saved successfully.\n")
cat("HDF5 methylation matrix remains a separate file.\n")
