###############################################################################
# ExplainMyMelanoma AI
# SCRIPT 08 — PREDICTION CONSISTENCY AUDIT
# READ-ONLY
###############################################################################

project_root <- "/Users/katie/master/Master"

tcga_file <- file.path(
  project_root,
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune.rds"
)

mlp_file <- file.path(
  project_root,
  "results/models/TCGA_SKCM_MLP_final_model.rds"
)

tcga_object <- readRDS(tcga_file)
mlp <- readRDS(mlp_file)

# Reconstruct the exact training input
beta <- tcga_object$beta_matrix

training_beta <- as.matrix(
  beta[
    mlp$selected_cpgs,
    mlp$training_sample_ids,
    drop = FALSE
  ]
)

X_train <- t(training_beta)

X_scaled <- scale(
  X_train,
  center = mlp$pca$center,
  scale = mlp$pca$scale
)

pca_scores <- X_scaled %*% mlp$pca$rotation

pca_scores_scaled <- scale(
  pca_scores,
  center = mlp$score_means,
  scale = mlp$score_sds
)

# Obtain class probabilities
mlp_output <- predict(
  mlp$model,
  newdata = pca_scores_scaled,
  type = "raw"
)

cat("\n========== 1. CLASS ORDER ==========\n")

cat("Expected class order:\n")
print(mlp$class_levels)

cat("Prediction column order:\n")
print(colnames(mlp_output))

cat(
  "Column order matches model class levels:",
  identical(colnames(mlp_output), mlp$class_levels),
  "\n"
)

cat("\n========== 2. PROBABILITY VALIDATION ==========\n")

row_probability_sums <- rowSums(mlp_output)

cat(
  "Minimum row probability sum:",
  min(row_probability_sums),
  "\n"
)

cat(
  "Maximum row probability sum:",
  max(row_probability_sums),
  "\n"
)

cat(
  "Maximum deviation from 1:",
  max(abs(row_probability_sums - 1)),
  "\n"
)

cat(
  "All probabilities within [0, 1]:",
  all(mlp_output >= 0 & mlp_output <= 1),
  "\n"
)

cat(
  "All row sums approximately 1:",
  all(abs(row_probability_sums - 1) < 1e-10),
  "\n"
)

cat("\n========== 3. SAMPLE ORDER ==========\n")

cat(
  "Prediction sample order matches training IDs:",
  identical(rownames(mlp_output), mlp$training_sample_ids),
  "\n"
)

cat(
  "Training labels are in the same sample order:",
  identical(names(mlp$training_labels), mlp$training_sample_ids),
  "\n"
)

cat("\n========== 4. PREDICTION VS ACTUAL LABELS ==========\n")

predicted_class <- factor(
  colnames(mlp_output)[max.col(mlp_output)],
  levels = mlp$class_levels
)

actual_class <- factor(
  as.character(mlp$training_labels),
  levels = mlp$class_levels
)

confusion_matrix <- table(
  Actual = actual_class,
  Predicted = predicted_class
)

print(confusion_matrix)

accuracy <- mean(predicted_class == actual_class)

cat("\nTraining-set prediction accuracy:", accuracy, "\n")

cat("\n========== 5. PREDICTION SUMMARY ==========\n")

prediction_summary <- data.frame(
  Sample_ID = rownames(mlp_output),
  Actual_Class = as.character(actual_class),
  Predicted_Class = as.character(predicted_class),
  Probability_Immune_Low = mlp_output[, "Immune-low"],
  Probability_Immune_High = mlp_output[, "Immune-high"],
  Correct = predicted_class == actual_class
)

cat("Prediction summary dimensions:\n")
print(dim(prediction_summary))

cat("Correct predictions:\n")
print(sum(prediction_summary$Correct))

cat("Incorrect predictions:\n")
print(sum(!prediction_summary$Correct))

cat("\nSCRIPT 08 — COMPLETED\n")
###############################################################################
