###############################################################################
# ExplainMyMelanoma AI
# SCRIPT 07 — MLP PREDICTION OUTPUT AUDIT
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

# Extract training samples and reconstruct PCA
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

# Apply stored training score scaling
pca_scores_scaled <- scale(
  pca_scores,
  center = mlp$score_means,
  scale = mlp$score_sds
)

# Predict with the saved neural network
cat("\n========== 1. MODEL OUTPUT ==========\n")

mlp_output <- predict(
  mlp$model,
  newdata = pca_scores_scaled,
  type = "raw"
)

cat("Output class:", paste(class(mlp_output), collapse = ", "), "\n")
cat("Output dimensions:", paste(dim(mlp_output), collapse = " x "), "\n")

cat("First 10 predictions:\n")
print(head(mlp_output, 10))

cat("\n========== 2. MODEL CONFIGURATION ==========\n")

print(mlp$selected_configuration)
cat("Class levels:\n")
print(mlp$class_levels)

cat("\n========== 3. OUTPUT RANGE ==========\n")

cat("Minimum output:", min(mlp_output), "\n")
cat("Maximum output:", max(mlp_output), "\n")
cat("Missing outputs:", sum(is.na(mlp_output)), "\n")

cat("\n========== 4. TRAINING LABEL COMPARISON ==========\n")

predicted_class <- if (is.matrix(mlp_output)) {
  colnames(mlp_output)[max.col(mlp_output)]
} else {
  ifelse(
    mlp_output >= 0.5,
    "Immune-high",
    "Immune-low"
  )
}

cat("Predicted class counts:\n")
print(table(predicted_class))

cat("Actual training class counts:\n")
print(table(mlp$training_labels))

cat("\nSCRIPT 07 — COMPLETED\n")
###############################################################################
