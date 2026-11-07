###############################################################################
# ExplainMyMelanoma AI
# SCRIPT 09 — MLP HOLDOUT PIPELINE AUDIT
# READ-ONLY — NO MODEL RETRAINING
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

split_file <- file.path(
  project_root,
  "results/models/TCGA_SKCM_ML_train_test_split.rds"
)

tcga_object <- readRDS(tcga_file)
mlp <- readRDS(mlp_file)
split <- readRDS(split_file)

beta <- tcga_object$beta_matrix

# 1. Verify train/test separation
cat("\n========== 1. TRAIN / TEST SEPARATION ==========\n")

train_ids <- mlp$training_sample_ids
test_ids <- split$test_sample_ids

cat("Training samples:", length(train_ids), "\n")
cat("Test samples:", length(test_ids), "\n")

cat(
  "No patient/sample overlap:",
  length(intersect(train_ids, test_ids)) == 0,
  "\n"
)

cat(
  "All test samples present in beta matrix:",
  all(test_ids %in% colnames(beta)),
  "\n"
)

cat(
  "Test labels correspond to test IDs:",
  identical(names(split$test_labels), test_ids),
  "\n"
)

cat("\n========== 2. HOLDOUT INPUT ==========\n")

test_beta <- as.matrix(
  beta[
    mlp$selected_cpgs,
    test_ids,
    drop = FALSE
  ]
)

X_test <- t(test_beta)

cat(
  "Test input dimensions:",
  paste(dim(X_test), collapse = " x "),
  "\n"
)

cat("Missing values:", sum(is.na(X_test)), "\n")
cat("Non-finite values:", sum(!is.finite(X_test)), "\n")

# 3. Apply the exact training preprocessing
cat("\n========== 3. APPLY TRAINING PREPROCESSING ==========\n")

X_test_scaled <- scale(
  X_test,
  center = mlp$pca$center,
  scale = mlp$pca$scale
)

test_pca_scores <- X_test_scaled %*% mlp$pca$rotation

test_pca_scaled <- scale(
  test_pca_scores,
  center = mlp$score_means,
  scale = mlp$score_sds
)

cat(
  "PCA dimensions:",
  paste(dim(test_pca_scores), collapse = " x "),
  "\n"
)

cat(
  "Scaled PCA dimensions:",
  paste(dim(test_pca_scaled), collapse = " x "),
  "\n"
)

# 4. Predict using the saved MLP
cat("\n========== 4. HOLDOUT PREDICTIONS ==========\n")

holdout_probabilities <- predict(
  mlp$model,
  newdata = test_pca_scaled,
  type = "raw"
)

cat(
  "Probability dimensions:",
  paste(dim(holdout_probabilities), collapse = " x "),
  "\n"
)

cat("Prediction columns:\n")
print(colnames(holdout_probabilities))

cat(
  "Prediction row order matches test IDs:",
  identical(rownames(holdout_probabilities), test_ids),
  "\n"
)

cat(
  "Missing probabilities:",
  sum(is.na(holdout_probabilities)),
  "\n"
)

cat(
  "Maximum probability sum deviation:",
  max(abs(rowSums(holdout_probabilities) - 1)),
  "\n"
)

# 5. Class predictions and confusion matrix
cat("\n========== 5. HOLDOUT CLASSIFICATION ==========\n")

actual <- factor(
  as.character(split$test_labels),
  levels = mlp$class_levels
)

predicted <- factor(
  colnames(holdout_probabilities)[
    max.col(holdout_probabilities)
  ],
  levels = mlp$class_levels
)

confusion <- table(
  Actual = actual,
  Predicted = predicted
)

print(confusion)

# Positive class: Immune-high
TP <- confusion["Immune-high", "Immune-high"]
TN <- confusion["Immune-low", "Immune-low"]
FP <- confusion["Immune-low", "Immune-high"]
FN <- confusion["Immune-high", "Immune-low"]

accuracy <- (TP + TN) / sum(confusion)
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
precision <- TP / (TP + FP)
F1 <- 2 * precision * sensitivity / (precision + sensitivity)
balanced_accuracy <- (sensitivity + specificity) / 2

# 6. Probability metrics
cat("\n========== 6. PROBABILITY METRICS ==========\n")

p_high <- holdout_probabilities[, "Immune-high"]

y <- as.integer(actual == "Immune-high")

brier_score <- mean((p_high - y)^2)

epsilon <- 1e-15
p_clamped <- pmin(pmax(p_high, epsilon), 1 - epsilon)

log_loss <- -mean(
  y * log(p_clamped) +
    (1 - y) * log(1 - p_clamped)
)

# AUC using the rank-sum method
positive_scores <- p_high[y == 1]
negative_scores <- p_high[y == 0]

n_positive <- length(positive_scores)
n_negative <- length(negative_scores)

auc <- (
  sum(rank(p_high)[y == 1]) -
    n_positive * (n_positive + 1) / 2
) / (n_positive * n_negative)

metrics <- data.frame(
  Metric = c(
    "Accuracy",
    "Sensitivity",
    "Specificity",
    "Precision",
    "F1",
    "Balanced accuracy",
    "AUC",
    "Brier score",
    "Log loss"
  ),
  Value = c(
    accuracy,
    sensitivity,
    specificity,
    precision,
    F1,
    balanced_accuracy,
    auc,
    brier_score,
    log_loss
  )
)

print(metrics, row.names = FALSE)

# 7. Compare with previously established holdout results
cat("\n========== 7. REFERENCE METRIC COMPARISON ==========\n")

reference_metrics <- data.frame(
  Metric = c(
    "Accuracy",
    "Sensitivity",
    "Specificity",
    "Precision",
    "F1",
    "Balanced accuracy",
    "AUC",
    "Brier score",
    "Log loss"
  ),
  Reference = c(
    0.83696,
    0.84783,
    0.82609,
    0.82979,
    0.83871,
    0.83696,
    0.91635,
    0.13897,
    0.58235
  )
)

comparison <- merge(
  reference_metrics,
  metrics,
  by = "Metric",
  sort = FALSE
)

comparison$Difference <- comparison$Value - comparison$Reference

print(comparison, row.names = FALSE)

cat("\nSCRIPT 09 — COMPLETED\n")
###############################################################################
