###############################################################################
# 
# TCGA-SKCM: FINAL MLP HOLDOUT EVALUATION
#
# Purpose:
#   Independent evaluation of the locked final MLP model
#   on the untouched 92-sample holdout set.
#
# IMPORTANT:
#   - No model retraining
#   - No hyperparameter tuning
#   - No feature reselection
#   - No preprocessing refitting
#   - No overwriting existing files
###############################################################################

cat("\n========== SCRIPT 25: INITIAL CHECK ==========\n")

# 1. Define input files
input_files <- c(
  integrated_data =
    "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune.rds",
  
  split_manifest =
    "results/models/TCGA_SKCM_ML_train_test_split.rds",
  
  final_mlp_model =
    "results/models/TCGA_SKCM_MLP_final_model.rds",
  
  final_elastic_net_model =
    "results/models/TCGA_SKCM_final_ElasticNet_model.rds"
)

# 2. Check that all inputs exist
input_check <- data.frame(
  Object = names(input_files),
  File = unname(input_files),
  Exists = file.exists(input_files)
)

print(input_check, row.names = FALSE)

if (!all(input_check$Exists)) {
  stop("One or more required input files are missing.")
}

# 3. Define new output files
output_files <- c(
  "results/models/TCGA_SKCM_MLP_holdout_evaluation.rds",
  "results/models/TCGA_SKCM_MLP_vs_ElasticNet_holdout.rds",
  "results/tables/TCGA_SKCM_MLP_holdout_metrics.csv",
  "results/tables/TCGA_SKCM_MLP_vs_ElasticNet_holdout_comparison.csv"
)





###############################################################################
# STEP 2. INSPECT SAVED MODEL STRUCTURES
###############################################################################

cat("\n========== LOADING SAVED OBJECTS ==========\n")

integrated_data <- readRDS(input_files["integrated_data"])

split_manifest <- readRDS(input_files["split_manifest"])

final_mlp_object <- readRDS(input_files["final_mlp_model"])

final_en_object <- readRDS(input_files["final_elastic_net_model"])


cat("\n========== MLP MODEL STRUCTURE ==========\n")

cat("\nTop-level elements:\n")
print(names(final_mlp_object))

cat("\nMLP model class:\n")
print(class(final_mlp_object$model))

cat("\nMLP model convergence:\n")
print(final_mlp_object$model$convergence)

cat("\nMLP model parameters:\n")
print(final_mlp_object$hyperparameter_selection)


cat("\n========== ELASTIC NET MODEL STRUCTURE ==========\n")

cat("\nTop-level elements:\n")
print(names(final_en_object))

cat("\nElastic Net model class:\n")
print(class(final_en_object$model))


cat("\n========== SPLIT MANIFEST ==========\n")

cat("\nTraining samples:\n")
print(length(split_manifest$train_sample_ids))

cat("\nHoldout samples:\n")
print(length(split_manifest$test_sample_ids))

cat("\nHoldout label distribution:\n")
print(table(split_manifest$test_labels))


cat("\n========== INTEGRATED DATA ==========\n")

cat("\nIntegrated object elements:\n")
print(names(integrated_data))

cat("\nBeta matrix dimensions:\n")
print(dim(integrated_data$beta_matrix))

cat("\n============================================\n")

###############################################################################
# STEP 3. MODEL FEATURE AND PREPROCESSING AUDIT
###############################################################################

cat("\n========== MLP FEATURE AUDIT ==========\n")

cat("\nSelected CpGs:\n")
print(length(final_mlp_object$selected_cpgs))

cat("\nNumber of PCA components:\n")
print(final_mlp_object$n_pcs)

cat("\nNumber of input features:\n")
print(final_mlp_object$n_features)

cat("\nPCA object class:\n")
print(class(final_mlp_object$pca))

cat("\nPCA rotation dimensions:\n")
print(dim(final_mlp_object$pca$rotation))

cat("\nScore means length:\n")
print(length(final_mlp_object$score_means))

cat("\nScore SDs length:\n")
print(length(final_mlp_object$score_sds))


cat("\n========== ELASTIC NET FEATURE AUDIT ==========\n")

cat("\nSelected CpGs:\n")
print(length(final_en_object$selected_cpg_ids))

cat("\nNumber of features:\n")
print(final_en_object$n_features)

cat("\nAlpha:\n")
print(final_en_object$alpha)

cat("\nLambda:\n")
print(final_en_object$lambda)

cat("\nLambda index:\n")
print(final_en_object$lambda_index)

cat("\nTraining sample count:\n")
print(final_en_object$n_training_samples)


cat("\n========== ELASTIC NET MODEL DETAILS ==========\n")

cat("\nModel dimensions:\n")
print(dim(final_en_object$model$beta))

cat("\nModel coefficients at selected lambda:\n")

en_coef <- as.matrix(
  coef(
    final_en_object$model,
    s = final_en_object$lambda
  )
)

print(dim(en_coef))

cat("\nFirst coefficient names:\n")
print(head(rownames(en_coef), 10))


cat("\n========== SAMPLE ID ALIGNMENT ==========\n")

beta_matrix <- integrated_data$beta_matrix

cat("\nHoldout IDs present in beta matrix:\n")
print(all(split_manifest$test_sample_ids %in% colnames(beta_matrix)))

cat("\nMLP training IDs match split manifest:\n")
print(setequal(
  final_mlp_object$training_sample_ids,
  split_manifest$train_sample_ids
))

cat("\nElastic Net training IDs match split manifest:\n")
print(setequal(
  final_en_object$training_sample_ids,
  split_manifest$train_sample_ids
))


cat("\n========== FEATURE ID ALIGNMENT ==========\n")

cat("\nMLP selected CpGs present in beta matrix:\n")
print(all(final_mlp_object$selected_cpgs %in% rownames(beta_matrix)))

cat("\nElastic Net selected CpGs present in beta matrix:\n")
print(all(final_en_object$selected_cpg_ids %in% rownames(beta_matrix)))

cat("\n============================================\n")




###############################################################################
# STEP 4. ELASTIC NET CLASS AND FEATURE ORDER AUDIT
###############################################################################

cat("\n========== ELASTIC NET CLASS AUDIT ==========\n")

cat("\nModel class names:\n")
print(final_en_object$model$classnames)

cat("\nTraining label levels:\n")
print(levels(final_en_object$training_labels))

cat("\nSplit manifest label levels:\n")
print(levels(split_manifest$train_labels))


cat("\n========== FEATURE ORDER AUDIT ==========\n")

en_feature_names <- rownames(en_coef)[-1]

cat("\nSelected CpG IDs equal coefficient names:\n")
print(setequal(
  final_en_object$selected_cpg_ids,
  en_feature_names
))

cat("\nSelected CpG order identical to coefficient order:\n")
print(identical(
  final_en_object$selected_cpg_ids,
  en_feature_names
))

cat("\nNumber of coefficient features:\n")
print(length(en_feature_names))


cat("\n========== MLP CLASS AUDIT ==========\n")

cat("\nMLP class levels:\n")
print(final_mlp_object$class_levels)

cat("\nMLP output names:\n")
print(colnames(final_mlp_object$model$fitted.values))


cat("\n========== HOLDOUT LABEL ALIGNMENT ==========\n")

cat("\nHoldout label IDs match holdout sample IDs:\n")
print(setequal(
  names(split_manifest$test_labels),
  split_manifest$test_sample_ids
))

cat("\n============================================\n")


###############################################################################
# STEP 5. ELASTIC NET CLASS MAPPING
# Training samples only — holdout is not accessed
###############################################################################

cat("\n========== ELASTIC NET CLASS MAPPING ==========\n")

# Prepare training matrix in the exact coefficient order
en_train_matrix <- as.matrix(
  t(
    beta_matrix[
      final_en_object$selected_cpg_ids,
      final_en_object$training_sample_ids
    ]
  )
)

# Align training labels to matrix rows
en_train_labels <- split_manifest$train_labels[
  rownames(en_train_matrix)
]

stopifnot(
  identical(
    rownames(en_train_matrix),
    names(en_train_labels)
  ),
  identical(
    colnames(en_train_matrix),
    final_en_object$selected_cpg_ids
  ),
  all(is.finite(en_train_matrix)),
  !anyNA(en_train_labels)
)

# Predict probabilities for numeric class 1
en_train_prob_class1 <- as.numeric(
  predict(
    final_en_object$model,
    newx = en_train_matrix,
    s = final_en_object$lambda,
    type = "response"
  )
)

# Compare mean probability by known biological class
class_mapping_audit <- data.frame(
  Actual_class = levels(en_train_labels),
  Mean_probability_class1 = vapply(
    levels(en_train_labels),
    function(cl) {
      mean(en_train_prob_class1[en_train_labels == cl])
    },
    numeric(1)
  )
)

cat("\nMean probability of numeric class 1 by biological class:\n")
print(class_mapping_audit, row.names = FALSE)

cat("\nTraining prediction dimensions:\n")
print(length(en_train_prob_class1))

cat("\nProbability range:\n")
print(range(en_train_prob_class1))

cat("\nClass mapping interpretation:\n")

if (
  class_mapping_audit$Mean_probability_class1[
    class_mapping_audit$Actual_class == "Immune-high"
  ] >
  class_mapping_audit$Mean_probability_class1[
    class_mapping_audit$Actual_class == "Immune-low"
  ]
) {
  cat("Numeric class 1 is consistent with Immune-high.\n")
} else {
  cat("Numeric class 1 is NOT consistent with Immune-high.\n")
}

cat("\n============================================\n")



###############################################################################
# STEP 6. GENERATE FINAL HOLDOUT PREDICTIONS
# Apply locked models; do not retrain or refit preprocessing
###############################################################################

cat("\n========== PREPARING HOLDOUT DATA ==========\n")

holdout_ids <- split_manifest$test_sample_ids

holdout_labels <- split_manifest$test_labels[holdout_ids]

stopifnot(
  length(holdout_ids) == 92,
  length(holdout_labels) == 92,
  identical(names(holdout_labels), holdout_ids),
  !anyNA(holdout_labels)
)

###############################################################################
# 1. MLP HOLDOUT PREDICTIONS
###############################################################################

cat("\n========== MLP HOLDOUT PREDICTIONS ==========\n")

# Extract exactly the CpGs used by the final MLP
mlp_holdout_matrix <- as.matrix(
  t(
    beta_matrix[
      final_mlp_object$selected_cpgs,
      holdout_ids
    ]
  )
)

stopifnot(
  identical(
    colnames(mlp_holdout_matrix),
    final_mlp_object$selected_cpgs
  ),
  identical(
    rownames(mlp_holdout_matrix),
    holdout_ids
  ),
  all(is.finite(mlp_holdout_matrix))
)

# Apply the saved PCA transformation
mlp_holdout_pca <- predict(
  final_mlp_object$pca,
  newdata = mlp_holdout_matrix
)

# Retain the selected number of components
mlp_holdout_pca <- mlp_holdout_pca[
  ,
  seq_len(final_mlp_object$n_pcs),
  drop = FALSE
]

# Apply training-derived scaling parameters
mlp_holdout_scaled <- sweep(
  mlp_holdout_pca,
  2,
  final_mlp_object$score_means,
  FUN = "-"
)

mlp_holdout_scaled <- sweep(
  mlp_holdout_scaled,
  2,
  final_mlp_object$score_sds,
  FUN = "/"
)

stopifnot(
  nrow(mlp_holdout_scaled) == 92,
  ncol(mlp_holdout_scaled) == final_mlp_object$n_pcs,
  all(is.finite(mlp_holdout_scaled))
)

# Predict using the locked final MLP
mlp_holdout_probabilities <- predict(
  final_mlp_object$model,
  newdata = mlp_holdout_scaled,
  type = "raw"
)

mlp_holdout_probabilities <- as.matrix(
  mlp_holdout_probabilities
)

stopifnot(
  nrow(mlp_holdout_probabilities) == 92,
  all(c("Immune-low", "Immune-high") %in%
        colnames(mlp_holdout_probabilities))
)

mlp_holdout_predictions <- data.frame(
  Sample_ID = holdout_ids,
  Actual = as.character(holdout_labels),
  Probability_Immune_low =
    mlp_holdout_probabilities[, "Immune-low"],
  Probability_Immune_high =
    mlp_holdout_probabilities[, "Immune-high"],
  Predicted = ifelse(
    mlp_holdout_probabilities[, "Immune-high"] >= 0.5,
    "Immune-high",
    "Immune-low"
  ),
  stringsAsFactors = FALSE
)

###############################################################################
# 2. ELASTIC NET HOLDOUT PREDICTIONS
###############################################################################

cat("\n========== ELASTIC NET HOLDOUT PREDICTIONS ==========\n")

# Extract CpGs in the exact coefficient order
en_holdout_matrix <- as.matrix(
  t(
    beta_matrix[
      final_en_object$selected_cpg_ids,
      holdout_ids
    ]
  )
)

stopifnot(
  identical(
    colnames(en_holdout_matrix),
    final_en_object$selected_cpg_ids
  ),
  identical(
    rownames(en_holdout_matrix),
    holdout_ids
  ),
  all(is.finite(en_holdout_matrix))
)

# Numeric class 1 was verified as Immune-high
en_holdout_prob_high <- as.numeric(
  predict(
    final_en_object$model,
    newx = en_holdout_matrix,
    s = final_en_object$lambda,
    type = "response"
  )
)

en_holdout_predictions <- data.frame(
  Sample_ID = holdout_ids,
  Actual = as.character(holdout_labels),
  Probability_Immune_high = en_holdout_prob_high,
  Predicted = ifelse(
    en_holdout_prob_high >= 0.5,
    "Immune-high",
    "Immune-low"
  ),
  stringsAsFactors = FALSE
)

###############################################################################
# 3. PREDICTION QUALITY CONTROL
###############################################################################

cat("\n========== HOLDOUT PREDICTION QC ==========\n")

prediction_qc <- c(
  
  holdout_sample_count =
    length(holdout_ids) == 92,
  
  holdout_IDs_unique =
    !anyDuplicated(holdout_ids),
  
  holdout_labels_complete =
    !anyNA(holdout_labels),
  
  MLP_predictions_cover_92_samples =
    nrow(mlp_holdout_predictions) == 92,
  
  Elastic_Net_predictions_cover_92_samples =
    nrow(en_holdout_predictions) == 92,
  
  MLP_IDs_preserved =
    identical(mlp_holdout_predictions$Sample_ID, holdout_ids),
  
  Elastic_Net_IDs_preserved =
    identical(en_holdout_predictions$Sample_ID, holdout_ids),
  
  MLP_probabilities_finite =
    all(is.finite(mlp_holdout_predictions$Probability_Immune_high)),
  
  Elastic_Net_probabilities_finite =
    all(is.finite(en_holdout_predictions$Probability_Immune_high)),
  
  MLP_probabilities_valid =
    all(
      mlp_holdout_predictions$Probability_Immune_high >= 0 &
        mlp_holdout_predictions$Probability_Immune_high <= 1
    ),
  
  Elastic_Net_probabilities_valid =
    all(
      en_holdout_predictions$Probability_Immune_high >= 0 &
        en_holdout_predictions$Probability_Immune_high <= 1
    ),
  
  Actual_labels_identical =
    identical(
      mlp_holdout_predictions$Actual,
      en_holdout_predictions$Actual
    )
)

print(prediction_qc)

stopifnot(all(prediction_qc))

cat("\nMLP holdout predictions:\n")
print(head(mlp_holdout_predictions))

cat("\nElastic Net holdout predictions:\n")
print(head(en_holdout_predictions))

cat("\nMLP probability range:\n")
print(range(mlp_holdout_predictions$Probability_Immune_high))

cat("\nElastic Net probability range:\n")
print(range(en_holdout_predictions$Probability_Immune_high))

cat("\nHoldout prediction generation completed.\n")
cat("No model was retrained.\n")
cat("============================================\n")



###############################################################################
# STEP 7. HOLDOUT PERFORMANCE METRICS
###############################################################################

cat("\n========== HOLDOUT PERFORMANCE METRICS ==========\n")

# Metric calculation function
calculate_metrics <- function(prediction_data, model_name) {
  
  actual <- factor(
    prediction_data$Actual,
    levels = c("Immune-low", "Immune-high")
  )
  
  predicted <- factor(
    prediction_data$Predicted,
    levels = c("Immune-low", "Immune-high")
  )
  
  probability <- prediction_data$Probability_Immune_high
  
  stopifnot(
    length(actual) == 92,
    length(predicted) == 92,
    length(probability) == 92,
    !anyNA(actual),
    !anyNA(predicted),
    all(is.finite(probability)),
    all(probability >= 0 & probability <= 1)
  )
  
  # Confusion matrix
  cm <- table(
    Actual = actual,
    Predicted = predicted
  )
  
  TN <- cm["Immune-low", "Immune-low"]
  FP <- cm["Immune-low", "Immune-high"]
  FN <- cm["Immune-high", "Immune-low"]
  TP <- cm["Immune-high", "Immune-high"]
  
  # Classification metrics
  accuracy <- (TP + TN) / sum(cm)
  
  sensitivity <- TP / (TP + FN)
  
  specificity <- TN / (TN + FP)
  
  precision <- if ((TP + FP) > 0) {
    TP / (TP + FP)
  } else {
    NA_real_
  }
  
  f1 <- if (is.finite(precision + sensitivity) &&
            (precision + sensitivity) > 0) {
    2 * precision * sensitivity /
      (precision + sensitivity)
  } else {
    NA_real_
  }
  
  balanced_accuracy <- mean(
    c(sensitivity, specificity)
  )
  
  # ROC-AUC using the Mann-Whitney rank statistic
  actual_binary <- as.numeric(actual == "Immune-high")
  
  n_positive <- sum(actual_binary == 1)
  n_negative <- sum(actual_binary == 0)
  
  auc <- (
    sum(rank(probability)[actual_binary == 1]) -
      n_positive * (n_positive + 1) / 2
  ) / (n_positive * n_negative)
  
  # Brier score
  brier <- mean(
    (probability - actual_binary)^2
  )
  
  # Log-loss
  epsilon <- 1e-15
  
  probability_clipped <- pmin(
    pmax(probability, epsilon),
    1 - epsilon
  )
  
  log_loss <- -mean(
    actual_binary * log(probability_clipped) +
      (1 - actual_binary) * log(1 - probability_clipped)
  )
  
  metrics <- data.frame(
    Model = model_name,
    N = length(actual),
    TN = as.numeric(TN),
    FP = as.numeric(FP),
    FN = as.numeric(FN),
    TP = as.numeric(TP),
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1 = f1,
    Balanced_Accuracy = balanced_accuracy,
    ROC_AUC = auc,
    Brier_Score = brier,
    Log_Loss = log_loss
  )
  
  list(
    confusion_matrix = cm,
    metrics = metrics
  )
}

###############################################################################
# Calculate metrics for both models
###############################################################################

mlp_holdout_evaluation <- calculate_metrics(
  mlp_holdout_predictions,
  "MLP"
)

en_holdout_evaluation <- calculate_metrics(
  en_holdout_predictions,
  "Elastic_Net"
)

###############################################################################
# Display confusion matrices
###############################################################################

cat("\n========== MLP CONFUSION MATRIX ==========\n")
print(mlp_holdout_evaluation$confusion_matrix)

cat("\n========== ELASTIC NET CONFUSION MATRIX ==========\n")
print(en_holdout_evaluation$confusion_matrix)

###############################################################################
# Combined metrics
###############################################################################

holdout_metrics <- rbind(
  mlp_holdout_evaluation$metrics,
  en_holdout_evaluation$metrics
)

cat("\n========== HOLDOUT METRICS COMPARISON ==========\n")

print(
  holdout_metrics,
  row.names = FALSE,
  digits = 4
)

cat("\n========== METRIC QC ==========\n")

metric_columns <- c(
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

metric_qc <- c(
  both_models_have_92_samples =
    all(holdout_metrics$N == 92),
  
  all_metrics_finite =
    all(is.finite(
      as.matrix(holdout_metrics[, metric_columns])
    )),
  
  classification_metrics_in_range =
    all(
      as.matrix(
        holdout_metrics[
          ,
          c(
            "Accuracy",
            "Sensitivity",
            "Specificity",
            "Precision",
            "F1",
            "Balanced_Accuracy",
            "ROC_AUC"
          )
        ]
      ) >= 0 &
        as.matrix(
          holdout_metrics[
            ,
            c(
              "Accuracy",
              "Sensitivity",
              "Specificity",
              "Precision",
              "F1",
              "Balanced_Accuracy",
              "ROC_AUC"
            )
          ]
        ) <= 1
    ),
  
  probability_metrics_nonnegative =
    all(
      holdout_metrics[
        ,
        c("Brier_Score", "Log_Loss")
      ] >= 0
    )
)

print(metric_qc)

stopifnot(all(metric_qc))

cat("\nHoldout metrics calculated successfully.\n")
cat("============================================\n")

###############################################################################
# STEP 8. SAVE FINAL HOLDOUT EVALUATION ARTIFACTS
###############################################################################

cat("\n========== SAVING HOLDOUT RESULTS ==========\n")

# Extend evaluation objects with their sample-level predictions
mlp_holdout_evaluation$predictions <- mlp_holdout_predictions
en_holdout_evaluation$predictions <- en_holdout_predictions

# Complete evaluation object
holdout_evaluation_object <- list(
  evaluation_type = "Independent holdout evaluation",
  cohort = "TCGA-SKCM",
  n_holdout = 92,
  positive_class = "Immune-high",
  classification_threshold = 0.5,
  
  MLP = mlp_holdout_evaluation,
  Elastic_Net = en_holdout_evaluation,
  
  comparison_metrics = holdout_metrics,
  
  prediction_qc = prediction_qc,
  metric_qc = metric_qc,
  
  holdout_sample_ids = holdout_ids,
  holdout_labels = holdout_labels,
  
  model_training_samples = 368,
  holdout_used_for_training = FALSE,
  
  script = "29_TCGA_SKCM_MLP_holdout_evaluation.R",
  evaluation_date = Sys.Date()
)

# Dedicated model comparison object
holdout_comparison_object <- list(
  models = c("MLP", "Elastic_Net"),
  sample_ids = holdout_ids,
  actual_labels = as.character(holdout_labels),
  
  MLP_predictions = mlp_holdout_predictions,
  Elastic_Net_predictions = en_holdout_predictions,
  
  metrics = holdout_metrics,
  
  confusion_matrices = list(
    MLP = mlp_holdout_evaluation$confusion_matrix,
    Elastic_Net = en_holdout_evaluation$confusion_matrix
  ),
  
  threshold = 0.5,
  positive_class = "Immune-high",
  same_holdout_samples = TRUE
)

# Output paths
output_holdout_rds <-
  "results/models/TCGA_SKCM_MLP_holdout_evaluation.rds"

output_comparison_rds <-
  "results/models/TCGA_SKCM_MLP_vs_ElasticNet_holdout.rds"

output_metrics_csv <-
  "results/tables/TCGA_SKCM_MLP_holdout_metrics.csv"

output_comparison_csv <-
  "results/tables/TCGA_SKCM_MLP_vs_ElasticNet_holdout_comparison.csv"

final_output_paths <- c(
  output_holdout_rds,
  output_comparison_rds,
  output_metrics_csv,
  output_comparison_csv
)

# Create directories if needed
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

# Never overwrite existing files
if (any(file.exists(final_output_paths))) {
  stop(
    "At least one output file already exists. ",
    "Nothing was overwritten."
  )
}

# Save RDS objects
saveRDS(
  holdout_evaluation_object,
  output_holdout_rds
)

saveRDS(
  holdout_comparison_object,
  output_comparison_rds
)

# Save CSV tables
write.csv(
  holdout_metrics,
  output_metrics_csv,
  row.names = FALSE
)

write.csv(
  holdout_metrics,
  output_comparison_csv,
  row.names = FALSE
)

###############################################################################
# VERIFY SAVED FILES
###############################################################################

cat("\n========== SAVED FILE VERIFICATION ==========\n")

saved_file_check <- data.frame(
  File = final_output_paths,
  Exists = file.exists(final_output_paths),
  Size_MB = round(
    file.info(final_output_paths)$size / 1024^2,
    4
  )
)

print(saved_file_check, row.names = FALSE)

stopifnot(
  all(saved_file_check$Exists),
  all(saved_file_check$Size_MB > 0)
)

cat("\nSCRIPT 29 HOLDOUT EVALUATION COMPLETED.\n")
cat("Both locked models were evaluated on the same 92 samples.\n")
cat("No model was retrained or tuned using holdout data.\n")
cat("All results were saved without overwriting existing files.\n")




