
# ============================================================
# 22_TCGA_SKCM_sensitivity_analysis.R
#
# PROJECT:
# Integrative Multi-Omics Analysis of Epigenetic Signatures
# Associated with the Immune Phenotype of Melanoma
#
# ANALYSIS:
# Sensitivity Analysis Using a Training-Only TIS Threshold
#
# PURPOSE:
# 1. Recalculate the TIS threshold using training samples only.
# 2. Compare original and alternative test labels.
# 3. Re-evaluate the final Elastic Net model using unchanged
#    holdout predictions.
# 4. Assess the robustness of holdout performance.
#
# IMPORTANT:
# - The final model is NOT retrained.
# - The original holdout predictions are NOT modified.
# - No feature selection or hyperparameter tuning is performed.
#
# INPUT FILES:
# - TCGA_SKCM_integrated_multiomics_KNN_immune.rds
# - TCGA_SKCM_ML_train_test_split.rds
# - TCGA_SKCM_final_holdout_evaluation.rds
#
# OUTPUT FILES:
# - results/models/TCGA_SKCM_sensitivity_analysis.rds
# - results/tables/TCGA_SKCM_sensitivity_analysis_metrics.csv
# - results/tables/TCGA_SKCM_sensitivity_analysis_predictions.csv
# - results/tables/TCGA_SKCM_sensitivity_analysis_label_comparison.csv
# - results/tables/TCGA_SKCM_sensitivity_analysis_QC.csv
#
# RANDOM SEED:
# Not required: no random sampling or model fitting.
#
# ============================================================


# ============================================================
# 0. INITIALIZATION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("TCGA-SKCM SENSITIVITY ANALYSIS\n")
cat("============================================================\n")

# Set project working directory before running this script.
# Expected working directory: /Users/katie/master/Master

project_dir <- getwd()

cat("\nWorking directory:", project_dir, "\n")

# Create output directories if necessary.

dir.create(
  "results/models",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. REQUIRED INPUT FILES
# ============================================================

input_files <- c(
  integrated_data =
    "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune.rds",
  
  split_manifest =
    "results/models/TCGA_SKCM_ML_train_test_split.rds",
  
  holdout_evaluation =
    "results/models/TCGA_SKCM_final_holdout_evaluation.rds"
)

file_qc <- file.exists(input_files)

cat("\nINPUT FILE QC:\n")

print(
  data.frame(
    File = names(input_files),
    Path = unname(input_files),
    Exists = unname(file_qc)
  ),
  row.names = FALSE
)

if (!all(file_qc)) {
  stop("ERROR: One or more required input files are missing.")
}


# ============================================================
# 2. LOAD DATA
# ============================================================

cat("\n--- Loading input data ---\n")

ml_source <- readRDS(input_files["integrated_data"])

split_manifest <- readRDS(input_files["split_manifest"])

holdout_eval <- readRDS(input_files["holdout_evaluation"])

immune_scores <- ml_source$immune_scores

holdout_predictions <- holdout_eval$predictions

cat("Immune scores dimensions:",
    paste(dim(immune_scores), collapse = " x "), "\n")

cat("Training samples:",
    length(split_manifest$train_sample_ids), "\n")

cat("Test samples:",
    length(split_manifest$test_sample_ids), "\n")

cat("Holdout predictions:",
    nrow(holdout_predictions), "\n")


# ============================================================
# 3. INITIAL DATA QC
# ============================================================

required_immune_columns <- c(
  "Sample_ID",
  "TIS_mean_zscore",
  "immune_group"
)

required_prediction_columns <- c(
  "Sample_ID",
  "Actual",
  "Predicted",
  "Probability_Immune_high",
  "Classification_threshold"
)

initial_qc <- c(
  correct_immune_dimensions =
    nrow(immune_scores) == 460,
  
  required_immune_columns_present =
    all(required_immune_columns %in% colnames(immune_scores)),
  
  unique_immune_sample_IDs =
    !anyDuplicated(immune_scores$Sample_ID),
  
  correct_training_sample_count =
    length(split_manifest$train_sample_ids) == 368,
  
  correct_test_sample_count =
    length(split_manifest$test_sample_ids) == 92,
  
  required_prediction_columns_present =
    all(required_prediction_columns %in%
          colnames(holdout_predictions)),
  
  correct_holdout_prediction_count =
    nrow(holdout_predictions) == 92
)

cat("\nINITIAL QC:\n")
print(initial_qc)

if (!all(initial_qc)) {
  stop("ERROR: Initial data QC failed.")
}


# ============================================================
# 4. SEPARATE TRAINING AND TEST SAMPLES
# ============================================================

training_scores <- immune_scores[
  immune_scores$Sample_ID %in%
    split_manifest$train_sample_ids,
  ,
  drop = FALSE
]

test_scores <- immune_scores[
  immune_scores$Sample_ID %in%
    split_manifest$test_sample_ids,
  ,
  drop = FALSE
]

# Restore the exact order defined in the split manifest.

training_scores <- training_scores[
  match(
    split_manifest$train_sample_ids,
    training_scores$Sample_ID
  ),
  ,
  drop = FALSE
]

test_scores <- test_scores[
  match(
    split_manifest$test_sample_ids,
    test_scores$Sample_ID
  ),
  ,
  drop = FALSE
]

stopifnot(
  nrow(training_scores) == 368,
  nrow(test_scores) == 92,
  identical(
    training_scores$Sample_ID,
    split_manifest$train_sample_ids
  ),
  identical(
    test_scores$Sample_ID,
    split_manifest$test_sample_ids
  )
)

cat("\nSAMPLE DISTRIBUTION:\n")

cat("Training:", nrow(training_scores), "\n")
cat("Test:", nrow(test_scores), "\n")


# ============================================================
# 5. CALCULATE TRAINING-ONLY TIS THRESHOLD
# ============================================================

training_tis_median <- median(
  training_scores$TIS_mean_zscore,
  na.rm = TRUE
)

if (!is.finite(training_tis_median)) {
  stop("ERROR: Training TIS median is not finite.")
}

cat("\nTRAINING-ONLY TIS MEDIAN:\n")
print(training_tis_median)


# ============================================================
# 6. GENERATE ALTERNATIVE TEST LABELS
# ============================================================

test_scores$immune_group_training_median <- ifelse(
  test_scores$TIS_mean_zscore > training_tis_median,
  "Immune-high",
  "Immune-low"
)

test_scores$immune_group_training_median <- factor(
  test_scores$immune_group_training_median,
  levels = c("Immune-low", "Immune-high")
)

if (anyNA(test_scores$immune_group_training_median)) {
  stop("ERROR: Missing alternative test labels.")
}

cat("\nNEW TEST LABEL DISTRIBUTION:\n")

print(
  table(test_scores$immune_group_training_median)
)

cat("\nORIGINAL TEST LABEL DISTRIBUTION:\n")

print(
  table(test_scores$immune_group)
)


# ============================================================
# 7. COMPARE ORIGINAL AND ALTERNATIVE LABELS
# ============================================================

label_comparison <- table(
  Original = factor(
    test_scores$immune_group,
    levels = c("Immune-low", "Immune-high")
  ),
  Training_only_median =
    test_scores$immune_group_training_median
)

n_changed <- sum(
  as.character(test_scores$immune_group) !=
    as.character(test_scores$immune_group_training_median)
)

label_agreement <- mean(
  as.character(test_scores$immune_group) ==
    as.character(test_scores$immune_group_training_median)
)

cat("\nLABEL COMPARISON:\n")
print(label_comparison)

cat("\nChanged labels:", n_changed, "of 92\n")

cat(
  "Label agreement:",
  round(label_agreement * 100, 2),
  "%\n"
)


# ============================================================
# 8. ALIGN ALTERNATIVE LABELS WITH HOLDOUT PREDICTIONS
# ============================================================

prediction_indices <- match(
  holdout_predictions$Sample_ID,
  test_scores$Sample_ID
)

if (anyNA(prediction_indices)) {
  stop("ERROR: Holdout prediction IDs could not be matched.")
}

new_labels <- test_scores$immune_group_training_median[
  prediction_indices
]

if (!identical(
  holdout_predictions$Sample_ID,
  test_scores$Sample_ID[prediction_indices]
)) {
  stop("ERROR: Sample ID alignment failed.")
}

actual_new <- factor(
  new_labels,
  levels = c("Immune-low", "Immune-high")
)

predicted_existing <- factor(
  holdout_predictions$Predicted,
  levels = c("Immune-low", "Immune-high")
)

probability_high <-
  holdout_predictions$Probability_Immune_high

if (
  anyNA(actual_new) ||
  anyNA(predicted_existing) ||
  any(!is.finite(probability_high)) ||
  any(probability_high < 0 | probability_high > 1)
) {
  stop("ERROR: Invalid labels, predictions or probabilities.")
}


# ============================================================
# 9. CONFUSION MATRIX
# ============================================================

sensitivity_confusion <- table(
  Actual = actual_new,
  Predicted = predicted_existing
)

cat("\nSENSITIVITY CONFUSION MATRIX:\n")

print(sensitivity_confusion)

TN <- unname(
  sensitivity_confusion["Immune-low", "Immune-low"]
)

FP <- unname(
  sensitivity_confusion["Immune-low", "Immune-high"]
)

FN <- unname(
  sensitivity_confusion["Immune-high", "Immune-low"]
)

TP <- unname(
  sensitivity_confusion["Immune-high", "Immune-high"]
)


# ============================================================
# 10. CLASSIFICATION METRICS
# ============================================================

accuracy <- (TP + TN) / (TP + TN + FP + FN)

sensitivity <- TP / (TP + FN)

specificity <- TN / (TN + FP)

precision <- if ((TP + FP) > 0) {
  TP / (TP + FP)
} else {
  NA_real_
}

f1_score <- if (
  is.finite(precision) &&
  (precision + sensitivity) > 0
) {
  2 * precision * sensitivity /
    (precision + sensitivity)
} else {
  NA_real_
}

balanced_accuracy <- (
  sensitivity + specificity
) / 2


# ============================================================
# 11. ROC-AUC AND PROBABILITY METRICS
# ============================================================

if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Package pROC is required. Install it using BiocManager or install.packages().")
}

sensitivity_roc <- pROC::roc(
  response = actual_new,
  predictor = probability_high,
  levels = c("Immune-low", "Immune-high"),
  direction = "<",
  quiet = TRUE
)

sensitivity_auc <- as.numeric(
  pROC::auc(sensitivity_roc)
)

sensitivity_auc_ci <- as.numeric(
  pROC::ci.auc(sensitivity_roc)
)

actual_numeric <- as.numeric(
  actual_new == "Immune-high"
)

epsilon <- 1e-15

probability_clipped <- pmin(
  pmax(probability_high, epsilon),
  1 - epsilon
)

brier_score <- mean(
  (probability_high - actual_numeric)^2
)

log_loss <- -mean(
  actual_numeric * log(probability_clipped) +
    (1 - actual_numeric) *
    log(1 - probability_clipped)
)


# ============================================================
# 12. SENSITIVITY METRICS TABLE
# ============================================================

sensitivity_metrics <- data.frame(
  Metric = c(
    "Accuracy",
    "Sensitivity",
    "Specificity",
    "Precision",
    "F1_score",
    "Balanced_Accuracy",
    "ROC_AUC",
    "ROC_AUC_CI_95_Lower",
    "ROC_AUC_CI_95_Upper",
    "Brier_Score",
    "Log_Loss"
  ),
  
  Value = c(
    accuracy,
    sensitivity,
    specificity,
    precision,
    f1_score,
    balanced_accuracy,
    sensitivity_auc,
    sensitivity_auc_ci[1],
    sensitivity_auc_ci[3],
    brier_score,
    log_loss
  )
)

cat("\nSENSITIVITY ANALYSIS METRICS:\n")

print(
  sensitivity_metrics,
  row.names = FALSE
)


# ============================================================
# 13. ORIGINAL VS SENSITIVITY METRICS
# ============================================================

original_metrics <- holdout_eval$metrics

if (!all(c("Metric", "Value") %in% colnames(original_metrics))) {
  stop("ERROR: Unexpected structure of original holdout metrics.")
}

comparison_metrics <- merge(
  original_metrics,
  sensitivity_metrics,
  by = "Metric",
  suffixes = c("_Original", "_Sensitivity")
)

comparison_metrics$Difference <-
  comparison_metrics$Value_Sensitivity -
  comparison_metrics$Value_Original

cat("\nORIGINAL VS SENSITIVITY METRICS:\n")

print(
  comparison_metrics,
  row.names = FALSE
)


# ============================================================
# 14. SAVE PREDICTIONS WITH ALTERNATIVE LABELS
# ============================================================

sensitivity_predictions_final <- data.frame(
  Sample_ID = holdout_predictions$Sample_ID,
  
  Original_Actual =
    holdout_predictions$Actual,
  
  Training_only_Actual =
    as.character(actual_new),
  
  Predicted =
    holdout_predictions$Predicted,
  
  Probability_Immune_high =
    probability_high,
  
  Classification_threshold =
    holdout_predictions$Classification_threshold
)


# ============================================================
# 15. QUALITY CONTROL
# ============================================================

sensitivity_qc <- c(
  
  input_files_exist =
    all(file.exists(input_files)),
  
  correct_training_count =
    nrow(training_scores) == 368,
  
  correct_test_count =
    nrow(test_scores) == 92,
  
  unique_test_sample_IDs =
    !anyDuplicated(test_scores$Sample_ID),
  
  exact_test_ID_alignment =
    identical(
      holdout_predictions$Sample_ID,
      test_scores$Sample_ID[prediction_indices]
    ),
  
  finite_training_median =
    is.finite(training_tis_median),
  
  no_missing_alternative_labels =
    !anyNA(actual_new),
  
  correct_confusion_matrix_count =
    sum(sensitivity_confusion) == 92,
  
  valid_probabilities =
    all(
      is.finite(probability_high) &
        probability_high >= 0 &
        probability_high <= 1
    ),
  
  finite_metrics =
    all(is.finite(sensitivity_metrics$Value)),
  
  original_predictions_preserved =
    isTRUE(all.equal(
      holdout_predictions$Probability_Immune_high,
      probability_high
    )),
  
  model_not_retrained =
    TRUE
)

cat("\nSENSITIVITY ANALYSIS QC:\n")

print(sensitivity_qc)

if (!all(sensitivity_qc)) {
  stop("ERROR: Sensitivity analysis QC failed.")
}


# ============================================================
# 16. CREATE FINAL RESULTS OBJECT
# ============================================================

sensitivity_results <- list(
  
  analysis =
    "Training-only TIS threshold sensitivity analysis",
  
  training_TIS_median =
    training_tis_median,
  
  training_sample_count =
    nrow(training_scores),
  
  test_sample_count =
    nrow(test_scores),
  
  changed_labels =
    n_changed,
  
  label_agreement =
    label_agreement,
  
  label_comparison =
    label_comparison,
  
  confusion_matrix =
    sensitivity_confusion,
  
  metrics =
    sensitivity_metrics,
  
  original_metrics =
    original_metrics,
  
  comparison_metrics =
    comparison_metrics,
  
  predictions =
    sensitivity_predictions_final,
  
  roc =
    sensitivity_roc,
  
  roc_auc_ci =
    sensitivity_auc_ci,
  
  qc =
    sensitivity_qc,
  
  model_retrained =
    FALSE,
  
  original_model_predictions_preserved =
    TRUE
)


# ============================================================
# 17. SAVE RESULTS
# ============================================================

saveRDS(
  sensitivity_results,
  "results/models/TCGA_SKCM_sensitivity_analysis.rds"
)

write.csv(
  sensitivity_metrics,
  "results/tables/TCGA_SKCM_sensitivity_analysis_metrics.csv",
  row.names = FALSE
)

write.csv(
  sensitivity_predictions_final,
  "results/tables/TCGA_SKCM_sensitivity_analysis_predictions.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame.matrix(label_comparison),
  "results/tables/TCGA_SKCM_sensitivity_analysis_label_comparison.csv"
)

write.csv(
  comparison_metrics,
  "results/tables/TCGA_SKCM_sensitivity_analysis_metric_comparison.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    QC_Check = names(sensitivity_qc),
    Passed = as.logical(sensitivity_qc)
  ),
  "results/tables/TCGA_SKCM_sensitivity_analysis_QC.csv",
  row.names = FALSE
)


# ============================================================
# 18. FINAL FILE VERIFICATION
# ============================================================

output_files <- c(
  "results/models/TCGA_SKCM_sensitivity_analysis.rds",
  
  "results/tables/TCGA_SKCM_sensitivity_analysis_metrics.csv",
  
  "results/tables/TCGA_SKCM_sensitivity_analysis_predictions.csv",
  
  "results/tables/TCGA_SKCM_sensitivity_analysis_label_comparison.csv",
  
  "results/tables/TCGA_SKCM_sensitivity_analysis_metric_comparison.csv",
  
  "results/tables/TCGA_SKCM_sensitivity_analysis_QC.csv"
)

output_qc <- file.exists(output_files)

cat("\nOUTPUT FILE QC:\n")

print(
  data.frame(
    File = output_files,
    Exists = output_qc
  ),
  row.names = FALSE
)

if (!all(output_qc)) {
  stop("ERROR: One or more output files were not saved.")
}


# ============================================================
# 19. FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SENSITIVITY ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat(
  "Training-only TIS median:",
  training_tis_median,
  "\n"
)

cat(
  "Changed test labels:",
  n_changed,
  "of",
  nrow(test_scores),
  "\n"
)

cat(
  "Label agreement:",
  round(label_agreement * 100, 2),
  "%\n"
)

cat(
  "Sensitivity ROC-AUC:",
  round(sensitivity_auc, 4),
  "\n"
)

cat(
  "Sensitivity Accuracy:",
  round(accuracy, 4),
  "\n"
)

cat(
  "Original model predictions preserved: YES\n"
)

cat(
  "Model retrained: NO\n"
)

cat(
  "All QC checks passed:",
  all(sensitivity_qc),
  "\n"
)

cat(
  "All output files saved:",
  all(output_qc),
  "\n"
)

cat("\nFINAL SENSITIVITY ANALYSIS: PASS\n")

# ============================================================
# END OF SCRIPT
# ============================================================
