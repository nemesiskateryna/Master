
###############################################################################
# SCRIPT 18
# TCGA-SKCM: ELASTIC NET PILOT MODEL - AUDIT AND RESULTS EXPORT
#
# Purpose:
# Validate and document the previously trained Elastic Net pilot model
# evaluated on outer validation fold 1.
#
# Input:
#   results/models/TCGA_SKCM_ElasticNet_fold1_pilot.rds
#
# Outputs:
#   results/tables/TCGA_SKCM_ElasticNet_fold1_pilot_metrics.csv
#   results/tables/TCGA_SKCM_ElasticNet_fold1_pilot_predictions.csv
#   results/tables/TCGA_SKCM_ElasticNet_fold1_pilot_selected_CpG.csv
#   results/tables/TCGA_SKCM_ElasticNet_fold1_pilot_QC.csv
###############################################################################

# 0. PROJECT DIRECTORY --------------------------------------------------------

setwd("~/master/Master")

# 1. LOAD SAVED PILOT MODEL ----------------------------------------------------

pilot_check <- readRDS(
  "results/models/TCGA_SKCM_ElasticNet_fold1_pilot.rds"
)

# 2. CHECK OBJECT STRUCTURE ----------------------------------------------------

required_components <- c(
  "model",
  "metrics",
  "predictions",
  "selected_CpG_lambda_min",
  "selected_CpG_lambda_1se",
  "top_5000_CpG",
  "confusion_matrix",
  "fold_number",
  "seed"
)

stopifnot(
  all(required_components %in% names(pilot_check))
)

cat("\n--- PILOT MODEL INFORMATION ---\n")

cat("Fold:", pilot_check$fold_number, "\n")
cat("Seed:", pilot_check$seed, "\n")

cat(
  "Features at lambda.min:",
  length(pilot_check$selected_CpG_lambda_min),
  "\n"
)

cat(
  "Features at lambda.1se:",
  length(pilot_check$selected_CpG_lambda_1se),
  "\n"
)

# 3. VALIDATE SELECTED FEATURES ------------------------------------------------

selected_cpg <- pilot_check$selected_CpG_lambda_1se

stopifnot(
  length(selected_cpg) > 0,
  !anyDuplicated(selected_cpg),
  all(grepl("^cg[0-9]+$", selected_cpg))
)

selected_cpg_table <- data.frame(
  CpG_ID = selected_cpg,
  Selection_method = "Elastic Net lambda.1se",
  Fold = pilot_check$fold_number,
  stringsAsFactors = FALSE
)

# 4. VALIDATE PREDICTIONS -------------------------------------------------------

predictions <- pilot_check$predictions

required_prediction_columns <- c(
  "Probability_Immune_high"
)

stopifnot(
  all(required_prediction_columns %in% colnames(predictions))
)

probabilities <- predictions$Probability_Immune_high

stopifnot(
  all(is.finite(probabilities)),
  all(probabilities >= 0 & probabilities <= 1)
)

stopifnot(
  nrow(predictions) == 74
)

# 5. VALIDATE METRICS -----------------------------------------------------------

metrics <- pilot_check$metrics

required_metric_columns <- c(
  "Fold",
  "Model",
  "Alpha",
  "Lambda_1se",
  "Selected_CpG",
  "Accuracy",
  "Balanced_accuracy",
  "Sensitivity",
  "Specificity",
  "Precision",
  "F1",
  "AUC"
)

stopifnot(
  all(required_metric_columns %in% colnames(metrics))
)

stopifnot(
  nrow(metrics) == 1
)

# 6. CONFUSION MATRIX -----------------------------------------------------------

cat("\n--- CONFUSION MATRIX ---\n")
print(pilot_check$confusion_matrix)

# 7. FINAL PILOT QC -------------------------------------------------------------

pilot_checks <- c(
  required_object_components =
    all(required_components %in% names(pilot_check)),
  
  selected_CpG_count =
    length(selected_cpg) == 89,
  
  selected_CpG_unique =
    !anyDuplicated(selected_cpg),
  
  standard_CpG_identifiers =
    all(grepl("^cg[0-9]+$", selected_cpg)),
  
  validation_prediction_count =
    nrow(predictions) == 74,
  
  finite_prediction_probabilities =
    all(is.finite(probabilities)),
  
  probabilities_within_range =
    all(probabilities >= 0 & probabilities <= 1),
  
  metrics_available =
    all(required_metric_columns %in% colnames(metrics)),
  
  exactly_one_metrics_row =
    nrow(metrics) == 1
)

pilot_qc <- data.frame(
  Check = names(pilot_checks),
  Result = ifelse(pilot_checks, "PASS", "FAIL"),
  stringsAsFactors = FALSE
)

cat("\n--- CORRECTED PILOT MODEL QC ---\n")
print(pilot_qc, row.names = FALSE)

pilot_qc_pass <- all(pilot_checks)

cat(
  "\nPILOT MODEL QC:",
  if (pilot_qc_pass) "PASS" else "FAIL",
  "\n"
)

stopifnot(pilot_qc_pass)

###############################################################################
# STEP 8. EXPORT QC AND RESULTS
###############################################################################

write.csv(
  metrics,
  "results/tables/TCGA_SKCM_ElasticNet_fold1_pilot_metrics.csv",
  row.names = FALSE
)

write.csv(
  predictions,
  "results/tables/TCGA_SKCM_ElasticNet_fold1_pilot_predictions.csv",
  row.names = FALSE
)

write.csv(
  selected_cpg_table,
  "results/tables/TCGA_SKCM_ElasticNet_fold1_pilot_selected_CpG.csv",
  row.names = FALSE
)

write.csv(
  pilot_qc,
  "results/tables/TCGA_SKCM_ElasticNet_fold1_pilot_QC.csv",
  row.names = FALSE
)

cat("\n--- EXPORTED FILES ---\n")
cat("Pilot metrics: PASS\n")
cat("Pilot predictions: PASS\n")
cat("Selected CpG list: PASS\n")
cat("Pilot QC report: PASS\n")