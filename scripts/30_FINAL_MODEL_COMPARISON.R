###############################################################################
# FINAL MODEL COMPARISON — STEP 1
# Read-only audit of saved nested CV and holdout results
###############################################################################

files_to_check <- c(
  nested_cv = "results/models/TCGA_SKCM_nested_CV_OOF_results.rds",
  final_mlp = "results/models/TCGA_SKCM_MLP_final_model.rds",
  final_en = "results/models/TCGA_SKCM_final_ElasticNet_model.rds",
  holdout = "results/models/TCGA_SKCM_MLP_holdout_evaluation.rds"
)

cat("\n========== FILE AVAILABILITY ==========\n")
print(data.frame(
  Object = names(files_to_check),
  File = unname(files_to_check),
  Exists = file.exists(files_to_check)
))

stopifnot(all(file.exists(files_to_check)))

nested_cv_object <- readRDS(files_to_check["nested_cv"])
final_mlp_object <- readRDS(files_to_check["final_mlp"])
final_en_object <- readRDS(files_to_check["final_en"])
holdout_object <- readRDS(files_to_check["holdout"])

cat("\n========== NESTED CV OBJECT ==========\n")
print(names(nested_cv_object))
str(nested_cv_object, max.level = 2)

cat("\n========== MLP NESTED CV METRICS ==========\n")
print(final_mlp_object$nested_cv_oof_metrics)

cat("\n========== ELASTIC NET NESTED CV OBJECT ==========\n")
print(names(nested_cv_object))

cat("\n========== HOLDOUT COMPARISON METRICS ==========\n")
print(holdout_object$comparison_metrics)

cat("\n========== AUDIT COMPLETED ==========\n")


###############################################################################
# FIX — STANDARDIZE AND COMBINE MODEL METRICS
###############################################################################

# Common metrics for all four evaluations
common_columns <- c(
  "N",
  "Accuracy",
  "Sensitivity",
  "Specificity",
  "Precision",
  "F1",
  "Balanced_Accuracy",
  "ROC_AUC",
  "Brier_Score",
  "Log_Loss"
)

# Function to standardize one evaluation row
standardize_metrics <- function(data, model, evaluation) {
  
  data <- as.data.frame(data)
  
  missing_columns <- setdiff(common_columns, names(data))
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  result <- data[, common_columns, drop = FALSE]
  
  result$Model <- model
  result$Evaluation <- evaluation
  
  result <- result[
    ,
    c("Model", "Evaluation", common_columns),
    drop = FALSE
  ]
  
  return(result)
}

# Four standardized rows
mlp_nested_row <- standardize_metrics(
  mlp_nested_metrics,
  "MLP",
  "Nested_CV_OOF"
)

en_nested_row <- standardize_metrics(
  en_nested_metrics,
  "Elastic_Net",
  "Nested_CV_OOF"
)

mlp_holdout_row <- standardize_metrics(
  holdout_object$MLP$metrics,
  "MLP",
  "Independent_holdout"
)

en_holdout_row <- standardize_metrics(
  holdout_object$Elastic_Net$metrics,
  "Elastic_Net",
  "Independent_holdout"
)

# Combine
final_model_comparison <- rbind(
  mlp_nested_row,
  en_nested_row,
  mlp_holdout_row,
  en_holdout_row
)

rownames(final_model_comparison) <- NULL

cat("\n========== FINAL MODEL COMPARISON ==========\n")
print(final_model_comparison, row.names = FALSE)

cat("\n========== DIMENSIONS ==========\n")
print(dim(final_model_comparison))

cat("\n========== STRUCTURE CHECK ==========\n")
print(names(final_model_comparison))



###############################################################################
# FINAL MODEL COMPARISON — STEP 3
# Final quality control before saving
# Read-only
###############################################################################

metric_columns <- c(
  "Accuracy", "Sensitivity", "Specificity",
  "Precision", "F1", "Balanced_Accuracy",
  "ROC_AUC", "Brier_Score", "Log_Loss"
)

comparison_qc <- c(
  correct_dimensions =
    identical(dim(final_model_comparison), c(4L, 12L)),
  
  correct_model_names =
    setequal(unique(final_model_comparison$Model),
             c("MLP", "Elastic_Net")),
  
  correct_evaluation_types =
    setequal(unique(final_model_comparison$Evaluation),
             c("Nested_CV_OOF", "Independent_holdout")),
  
  correct_sample_counts =
    all(final_model_comparison$N[
      final_model_comparison$Evaluation == "Nested_CV_OOF"
    ] == 368) &&
    all(final_model_comparison$N[
      final_model_comparison$Evaluation == "Independent_holdout"
    ] == 92),
  
  all_metrics_finite =
    all(is.finite(as.matrix(
      final_model_comparison[, metric_columns]
    ))),
  
  classification_metrics_in_range =
    all(as.matrix(
      final_model_comparison[
        ,
        c("Accuracy", "Sensitivity", "Specificity",
          "Precision", "F1", "Balanced_Accuracy", "ROC_AUC")
      ]
    ) >= 0) &&
    all(as.matrix(
      final_model_comparison[
        ,
        c("Accuracy", "Sensitivity", "Specificity",
          "Precision", "F1", "Balanced_Accuracy", "ROC_AUC")
      ]
    ) <= 1),
  
  probability_metrics_nonnegative =
    all(final_model_comparison$Brier_Score >= 0) &&
    all(final_model_comparison$Log_Loss >= 0),
  
  elastic_net_higher_accuracy_both_evaluations =
    all(
      final_model_comparison$Accuracy[
        final_model_comparison$Model == "Elastic_Net"
      ] >
        final_model_comparison$Accuracy[
          final_model_comparison$Model == "MLP"
        ]
    ),
  
  elastic_net_higher_auc_both_evaluations =
    all(
      final_model_comparison$ROC_AUC[
        final_model_comparison$Model == "Elastic_Net"
      ] >
        final_model_comparison$ROC_AUC[
          final_model_comparison$Model == "MLP"
        ]
    )
)

cat("\n========== FINAL COMPARISON QC ==========\n")
print(comparison_qc)

cat("\nAll checks passed:\n")
print(all(comparison_qc))

stopifnot(all(comparison_qc))




###############################################################################
# FINAL MODEL COMPARISON — STEP 4
# Save final comparison table
# Never overwrite existing files
###############################################################################

output_files <- c(
  csv = "results/tables/TCGA_SKCM_MLP_vs_ElasticNet_final_comparison.csv",
  rds = "results/models/TCGA_SKCM_MLP_vs_ElasticNet_final_comparison.rds"
)

# Create directories if necessary
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/models", recursive = TRUE, showWarnings = FALSE)

# Protect existing files
if (any(file.exists(output_files))) {
  stop(
    "Output file already exists. Nothing was overwritten. ",
    "Check: ",
    paste(output_files[file.exists(output_files)], collapse = "; ")
  )
}

# Save
write.csv(
  final_model_comparison,
  output_files["csv"],
  row.names = FALSE
)

saveRDS(
  list(
    comparison_table = final_model_comparison,
    quality_control = comparison_qc,
    evaluation_design = list(
      nested_cv_samples = 368,
      holdout_samples = 92,
      positive_class = "Immune-high",
      classification_threshold = 0.5
    ),
    models = c("MLP", "Elastic_Net"),
    source_files = files_to_check,
    script = "Final model comparison",
    evaluation_date = Sys.Date()
  ),
  output_files["rds"]
)

# Verify saved files
saved_check <- data.frame(
  File = unname(output_files),
  Exists = file.exists(output_files),
  Size_MB = round(
    file.info(output_files)$size / 1024^2,
    4
  )
)

cat("\n========== SAVED FILE CHECK ==========\n")
print(saved_check, row.names = FALSE)

stopifnot(
  all(saved_check$Exists),
  all(saved_check$Size_MB > 0)
)

cat("\nFINAL MODEL COMPARISON SAVED SUCCESSFULLY\n")





