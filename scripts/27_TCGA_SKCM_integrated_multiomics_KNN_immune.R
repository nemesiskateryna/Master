
###############################################################################
# SCRIPT 23
# TCGA-SKCM MULTILAYER PERCEPTRON (MLP)
# Nested cross-validation preparation
#
# Objective:
# Predict Immune-high / Immune-low using DNA methylation CpG beta-values.
#
# Important:
# - Reuse existing patient-level train/test split.
# - Reuse existing nested CV manifests.
# - No holdout-based feature selection or tuning.
# - No modification of existing Elastic Net results.
###############################################################################

# 1. Packages ---------------------------------------------------------------

if (!requireNamespace("nnet", quietly = TRUE)) {
  install.packages("nnet")
}

library(nnet)

# 2. Project directory -------------------------------------------------------

project_dir <- "/Users/katie/master/Master"

setwd(project_dir)

# 3. Load existing methylation data -----------------------------------------

ml_source <- readRDS(
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune.rds"
)

# 4. Load existing train/test split ------------------------------------------

split_manifest <- readRDS(
  "results/models/TCGA_SKCM_ML_train_test_split.rds"
)

# 5. Load existing nested CV manifest ---------------------------------------

nested_manifest <- readRDS(
  "results/models/TCGA_SKCM_nested_CV_manifest.rds"
)

# 6. Extract methylation matrix and phenotype -------------------------------

beta_matrix <- ml_source$beta_matrix

immune_group <- ml_source$immune_group

# Explicit factor order
immune_group <- factor(
  as.character(immune_group),
  levels = c("Immune-low", "Immune-high")
)

# 7. Extract sample IDs from existing split ---------------------------------

train_ids <- split_manifest$train_ids
test_ids  <- split_manifest$test_ids

# 8. Basic object initialization --------------------------------------------

cat("\n========== MLP INITIALIZATION ==========\n")

cat(
  "Methylation matrix dimensions:",
  paste(dim(beta_matrix), collapse = " x "),
  "\n"
)

cat("Training samples:", length(train_ids), "\n")
cat("Holdout samples:", length(test_ids), "\n")

cat("\nTraining phenotype distribution:\n")
print(table(immune_group[train_ids]))

cat("\nHoldout phenotype distribution:\n")
print(table(immune_group[test_ids]))

cat("\nNested CV outer folds:",
    length(nested_manifest$outer_folds),
    "\n"
)

cat("\nNested CV inner folds:",
    length(nested_manifest$inner_folds),
    "\n"
)

cat("\nMLP package version:\n")
print(packageVersion("nnet"))

cat("\nMLP initialization completed.\n")



###############################################################################
# STEP 2 — READ-ONLY STRUCTURE AUDIT FOR MLP
###############################################################################

cat("\n========== SPLIT MANIFEST STRUCTURE ==========\n")

cat("\nTop-level names:\n")
print(names(split_manifest))

cat("\nStructure:\n")
str(split_manifest, max.level = 3)

cat("\n========== NESTED MANIFEST STRUCTURE ==========\n")

cat("\nTop-level names:\n")
print(names(nested_manifest))

cat("\nStructure:\n")
str(nested_manifest, max.level = 3)


cat("\n========== SPLIT ID CANDIDATES ==========\n")

if (is.list(split_manifest)) {
  
  for (nm in names(split_manifest)) {
    
    obj <- split_manifest[[nm]]
    
    cat(
      "\nField:", nm,
      "| Class:", paste(class(obj), collapse = ","),
      "| Length:", length(obj),
      "\n"
    )
    
    if (is.character(obj)) {
      cat("First IDs:\n")
      print(head(obj))
    }
    
    if (is.data.frame(obj)) {
      cat("Dimensions:", paste(dim(obj), collapse = " x "), "\n")
      print(head(obj))
    }
  }
}


cat("\n========== NESTED FOLD CANDIDATES ==========\n")

if (is.list(nested_manifest)) {
  
  for (nm in names(nested_manifest)) {
    
    obj <- nested_manifest[[nm]]
    
    cat(
      "\nField:", nm,
      "| Class:", paste(class(obj), collapse = ","),
      "| Length:", length(obj),
      "\n"
    )
    
    if (is.list(obj)) {
      cat("Child names:\n")
      print(names(obj))
    }
    
    if (is.data.frame(obj)) {
      cat("Dimensions:", paste(dim(obj), collapse = " x "), "\n")
      print(head(obj))
    }
  }
}

cat("\n========== PHENOTYPE ALIGNMENT ==========\n")

cat("Phenotype length:", length(immune_group), "\n")
cat("Phenotype names present:", !is.null(names(immune_group)), "\n")
cat("First phenotype IDs:\n")
print(head(names(immune_group)))

cat("\n========== AUDIT COMPLETE ==========\n")


###############################################################################
# STEP 3 — CORRECT DATA ALIGNMENT FOR MLP
# Read-only: no files are created or overwritten
###############################################################################

# 1. Correct sample ID fields ------------------------------------------------

train_ids <- split_manifest$train_sample_ids
test_ids  <- split_manifest$test_sample_ids

# 2. Correct phenotype labels ------------------------------------------------

y_train <- split_manifest$train_labels
y_test  <- split_manifest$test_labels

# Preserve explicit class order
y_train <- factor(
  as.character(y_train),
  levels = c("Immune-low", "Immune-high")
)

y_test <- factor(
  as.character(y_test),
  levels = c("Immune-low", "Immune-high")
)

names(y_train) <- train_ids
names(y_test)  <- test_ids

# 3. Check methylation matrix sample IDs -------------------------------------

matrix_sample_ids <- colnames(beta_matrix)

cat("\n========== MATRIX ALIGNMENT ==========\n")

cat("Matrix dimensions:",
    paste(dim(beta_matrix), collapse = " x "), "\n")

cat("Matrix sample IDs available:",
    !is.null(matrix_sample_ids), "\n")

cat("Training IDs:", length(train_ids), "\n")
cat("Holdout IDs:", length(test_ids), "\n")

cat("Training IDs found in matrix:",
    sum(train_ids %in% matrix_sample_ids), "\n")

cat("Holdout IDs found in matrix:",
    sum(test_ids %in% matrix_sample_ids), "\n")

# 4. Construct matrices in the exact manifest order -------------------------

train_positions <- match(train_ids, matrix_sample_ids)
test_positions  <- match(test_ids, matrix_sample_ids)

stopifnot(
  !is.null(matrix_sample_ids),
  !anyNA(train_positions),
  !anyNA(test_positions),
  !anyDuplicated(train_ids),
  !anyDuplicated(test_ids)
)

X_train <- beta_matrix[
  ,
  train_positions,
  drop = FALSE
]

X_test <- beta_matrix[
  ,
  test_positions,
  drop = FALSE
]

# 5. Verify exact alignment ---------------------------------------------------

cat("\n========== TRAIN / HOLDOUT OBJECTS ==========\n")

cat("X_train:",
    paste(dim(X_train), collapse = " x "), "\n")

cat("X_test:",
    paste(dim(X_test), collapse = " x "), "\n")

cat("Training column order matches manifest:",
    identical(colnames(X_train), train_ids), "\n")

cat("Holdout column order matches manifest:",
    identical(colnames(X_test), test_ids), "\n")

cat("\nTraining labels:\n")
print(table(y_train))

cat("\nHoldout labels:\n")
print(table(y_test))

cat("\nTraining labels named correctly:",
    identical(names(y_train), train_ids), "\n")

cat("Holdout labels named correctly:",
    identical(names(y_test), test_ids), "\n")

# 6. Verify nested outer folds against training sample order -----------------

cat("\n========== NESTED CV ALIGNMENT ==========\n")

outer_alignment <- lapply(
  seq_along(nested_manifest$nested_cv),
  function(i) {
    
    fold <- nested_manifest$nested_cv[[i]]
    
    expected_train_ids <- train_ids[fold$outer_train_idx]
    expected_valid_ids <- train_ids[fold$outer_valid_idx]
    
    data.frame(
      Outer_Fold = fold$outer_fold,
      Train_Size = length(fold$outer_train_idx),
      Validation_Size = length(fold$outer_valid_idx),
      Train_IDs_Match = identical(
        expected_train_ids,
        fold$outer_train_sample_ids
      ),
      Validation_IDs_Match = identical(
        expected_valid_ids,
        fold$outer_valid_sample_ids
      ),
      Inner_Folds = length(fold$inner_folds)
    )
  }
)

outer_alignment <- do.call(rbind, outer_alignment)

print(outer_alignment)

# 7. Inspect inner fold structure --------------------------------------------

cat("\n========== FIRST INNER FOLD STRUCTURE ==========\n")

first_inner_fold <- nested_manifest$nested_cv[[1]]$inner_folds[[1]]

cat("Field names:\n")
print(names(first_inner_fold))

cat("\nStructure:\n")
str(first_inner_fold, max.level = 2)

# 8. Final QC ----------------------------------------------------------------

mlp_alignment_qc <- c(
  correct_training_size = ncol(X_train) == 368,
  correct_holdout_size = ncol(X_test) == 92,
  correct_feature_count = nrow(X_train) == 409556,
  training_ID_order = identical(colnames(X_train), train_ids),
  holdout_ID_order = identical(colnames(X_test), test_ids),
  training_label_order = identical(names(y_train), train_ids),
  holdout_label_order = identical(names(y_test), test_ids),
  all_outer_train_IDs_match = all(outer_alignment$Train_IDs_Match),
  all_outer_validation_IDs_match = all(outer_alignment$Validation_IDs_Match),
  four_inner_folds_each = all(outer_alignment$Inner_Folds == 4)
)

cat("\n========== MLP ALIGNMENT QC ==========\n")
print(mlp_alignment_qc)

cat("\nAll QC checks passed:",
    all(mlp_alignment_qc), "\n")





###############################################################################
# STEP 6 — MLP COMPUTATION ENVIRONMENT AUDIT
# Read-only
###############################################################################

cat("\n========== MATRIX CLASS ==========\n")
print(class(beta_matrix))

cat("\nIs DelayedMatrix:\n")
print(methods::is(beta_matrix, "DelayedMatrix"))

cat("\n========== REQUIRED PACKAGES ==========\n")

packages <- c(
  "nnet",
  "DelayedArray",
  "DelayedMatrixStats",
  "matrixStats",
  "BiocParallel"
)

package_status <- data.frame(
  Package = packages,
  Installed = vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
)

print(package_status, row.names = FALSE)

cat("\n========== PACKAGE VERSIONS ==========\n")

for (pkg in packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(
      pkg, ":",
      as.character(packageVersion(pkg)),
      "\n"
    )
  }
}

cat("\n========== MATRIX DIMENSIONS ==========\n")

cat("Total CpGs:", nrow(beta_matrix), "\n")
cat("Total samples:", ncol(beta_matrix), "\n")
cat("Training samples:", ncol(X_train), "\n")

cat("\n========== MLP ENVIRONMENT AUDIT COMPLETE ==========\n")



###############################################################################
# STEP 7 — TEST DELAYEDMATRIX VARIANCE COMPUTATION
# Read-only: no files are created or modified
###############################################################################

cat("\n========== TEST: CpG VARIANCE COMPUTATION ==========\n")

# Беремо лише перші 1000 CpG для технічного тесту
test_matrix <- X_train[
  seq_len(1000),
  ,
  drop = FALSE
]

# Обчислюємо дисперсію між training-зразками
test_variance <- DelayedMatrixStats::rowVars(
  test_matrix,
  na.rm = TRUE
)

cat("Input class:", class(test_matrix), "\n")
cat("Number of variance values:", length(test_variance), "\n")

cat(
  "Finite variance values:",
  sum(is.finite(test_variance)),
  "\n"
)

cat(
  "NA variance values:",
  sum(is.na(test_variance)),
  "\n"
)

cat(
  "Zero variance CpGs:",
  sum(test_variance == 0, na.rm = TRUE),
  "\n"
)

cat(
  "Positive variance CpGs:",
  sum(test_variance > 0, na.rm = TRUE),
  "\n"
)

cat("\nVariance summary:\n")
print(summary(test_variance))

cat("\n========== STEP 7 QC ==========\n")

variance_qc <- c(
  correct_length = length(test_variance) == 1000,
  all_finite = all(is.finite(test_variance)),
  at_least_one_variable_cpg =
    any(test_variance > 0, na.rm = TRUE)
)

print(variance_qc)

stopifnot(all(variance_qc))


###############################################################################
# STEP 8 — FULL CpG VARIANCE RANKING IN ONE INNER FOLD
# Read-only: no files are created or modified
###############################################################################

cat("\n========== STEP 8: INNER FOLD 1 FEATURE SELECTION ==========\n")

# 1. Отримуємо перший outer fold
outer1 <- nested_manifest$nested_cv[[1]]

outer_train_idx <- outer1$outer_train_idx
inner_valid_idx <- outer1$inner_folds[[1]]

# 2. Формуємо inner-training індекси
inner_train_idx <- setdiff(
  seq_along(outer_train_idx),
  inner_valid_idx
)

# 3. Відбираємо лише inner-training зразки
inner_train_ids <- train_ids[
  outer_train_idx[inner_train_idx]
]

X_inner_train <- beta_matrix[
  ,
  match(inner_train_ids, colnames(beta_matrix)),
  drop = FALSE
]

cat("Inner-training samples:", ncol(X_inner_train), "\n")
cat("CpGs:", nrow(X_inner_train), "\n")

stopifnot(
  ncol(X_inner_train) == length(inner_train_ids),
  identical(colnames(X_inner_train), inner_train_ids)
)

# 4. Обчислюємо дисперсію CpG
variance_time <- system.time({
  
  cpg_variance <- DelayedMatrixStats::rowVars(
    X_inner_train,
    na.rm = TRUE
  )
  
})

cat("\nVariance computation time:\n")
print(variance_time)

# 5. Перевіряємо результати
variance_qc <- c(
  correct_length = length(cpg_variance) == nrow(beta_matrix),
  all_finite = all(is.finite(cpg_variance)),
  no_negative_values = all(cpg_variance >= 0),
  enough_variable_cpgs =
    sum(cpg_variance > 0) >= 5000
)

cat("\nVariance QC:\n")
print(variance_qc)

stopifnot(all(variance_qc))

# 6. Відбираємо 5000 CpG з найбільшою дисперсією
top5000_idx <- order(
  cpg_variance,
  decreasing = TRUE
)[seq_len(5000)]

top5000_cpg <- rownames(beta_matrix)[top5000_idx]

cat("\nSelected CpGs:", length(top5000_idx), "\n")
cat("Unique CpG IDs:", length(unique(top5000_cpg)), "\n")

cat("\nTop 10 CpGs:\n")
print(data.frame(
  CpG_ID = top5000_cpg[1:10],
  Variance = cpg_variance[top5000_idx[1:10]]
))

# 7. Фінальний QC
feature_selection_qc <- c(
  exactly_5000 = length(top5000_idx) == 5000,
  unique_features = length(unique(top5000_idx)) == 5000,
  valid_indices = all(top5000_idx >= 1 &
                        top5000_idx <= nrow(beta_matrix)),
  correct_training_samples =
    identical(colnames(X_inner_train), inner_train_ids)
)

cat("\n========== STEP 8 FINAL QC ==========\n")
print(feature_selection_qc)

stopifnot(all(feature_selection_qc))

###############################################################################
# DEFINE INNER FOLD SAMPLE IDs
###############################################################################

# Індекси поточного outer fold
outer_train_idx <- nested_manifest$nested_cv[[outer_fold]]$outer_train_idx

# Inner folds поточного outer fold
inner_folds <- nested_manifest$nested_cv[[outer_fold]]$inner_folds

# ID зовнішньої навчальної частини
outer_ids <- train_ids[outer_train_idx]

# Валідаційні позиції поточного inner fold
valid_idx <- inner_folds[[inner_fold]]

# ID inner-training та inner-validation
inner_train_idx <- setdiff(
  seq_along(outer_ids),
  valid_idx
)

inner_train_ids <- outer_ids[inner_train_idx]
inner_valid_ids <- outer_ids[valid_idx]

cat("Outer fold:", outer_fold, "\n")
cat("Inner fold:", inner_fold, "\n")
cat("Inner-training N:", length(inner_train_ids), "\n")
cat("Inner-validation N:", length(inner_valid_ids), "\n")

cat(
  "No ID overlap:",
  !any(inner_train_ids %in% inner_valid_ids),
  "\n"
)

stopifnot(
  length(inner_train_ids) > 0,
  length(inner_valid_ids) > 0,
  !any(inner_train_ids %in% inner_valid_ids)
)

###############################################################################
# STEP 9 — RESELECT CpGs AND FIT PCA
# Outer fold 5 / Inner fold 4
# Read-only
###############################################################################

# 1. Обчислюємо дисперсію лише на 222 training-зразках
cpg_variance <- DelayedMatrixStats::rowVars(
  X_inner_train,
  na.rm = TRUE
)

stopifnot(
  length(cpg_variance) == nrow(beta_matrix),
  all(is.finite(cpg_variance)),
  sum(cpg_variance > 0) >= 5000
)

# 2. Відбираємо 5000 найбільш варіабельних CpG
top5000_idx <- order(
  cpg_variance,
  decreasing = TRUE
)[seq_len(5000)]

top5000_cpg <- rownames(beta_matrix)[top5000_idx]

stopifnot(
  length(top5000_idx) == 5000,
  length(unique(top5000_idx)) == 5000
)

# 3. Витягуємо training і validation
train_5000 <- as.matrix(
  X_inner_train[
    top5000_idx,
    ,
    drop = FALSE
  ]
)

valid_5000 <- as.matrix(
  beta_matrix[
    top5000_idx,
    inner_valid_positions,
    drop = FALSE
  ]
)

# 4. Транспонуємо
train_pca_input <- t(train_5000)
valid_pca_input <- t(valid_5000)

# 5. Перевірка матриць
stopifnot(
  nrow(train_pca_input) == 222,
  nrow(valid_pca_input) == 74,
  ncol(train_pca_input) == 5000,
  ncol(valid_pca_input) == 5000,
  all(is.finite(train_pca_input)),
  all(is.finite(valid_pca_input))
)

# 6. Перевірка нульової дисперсії
train_sd <- apply(
  train_pca_input,
  2,
  sd
)

stopifnot(
  all(is.finite(train_sd)),
  all(train_sd > 0)
)

# 7. PCA — навчання тільки на training
pca_time <- system.time({
  
  pca_fit <- prcomp(
    train_pca_input,
    center = TRUE,
    scale. = TRUE,
    rank. = 30
  )
  
})

cat("\nPCA computation time:\n")
print(pca_time)

# 8. Training та validation scores
train_pca_scores <- pca_fit$x

valid_pca_scores <- predict(
  pca_fit,
  newdata = valid_pca_input
)

# 9. Перевірка відповідності ID
rownames(train_pca_scores) <- inner_train_ids
rownames(valid_pca_scores) <- inner_valid_ids

# 10. QC
pca_qc <- c(
  correct_training_N =
    nrow(train_pca_scores) == 222,
  
  correct_validation_N =
    nrow(valid_pca_scores) == 74,
  
  exactly_30_components =
    ncol(train_pca_scores) == 30,
  
  finite_training_scores =
    all(is.finite(train_pca_scores)),
  
  finite_validation_scores =
    all(is.finite(valid_pca_scores)),
  
  correct_training_ID_order =
    identical(rownames(train_pca_scores), inner_train_ids),
  
  correct_validation_ID_order =
    identical(rownames(valid_pca_scores), inner_valid_ids)
)

cat("\nPCA dimensions:\n")
cat("Training:", paste(dim(train_pca_scores), collapse = " x "), "\n")
cat("Validation:", paste(dim(valid_pca_scores), collapse = " x "), "\n")

cat("\nPCA QC:\n")
print(pca_qc)

stopifnot(all(pca_qc))

cat("\nSTEP 9 COMPLETED SUCCESSFULLY\n")



###############################################################################
# STEP 10 — VERIFY PHENOTYPE LABEL ALIGNMENT
# Outer fold 5 / Inner fold 4
# Read-only
###############################################################################

# 1. Отримуємо мітки за ID зразків
y_inner_train <- factor(
  as.character(y_train[inner_train_ids]),
  levels = c("Immune-low", "Immune-high")
)

y_inner_valid <- factor(
  as.character(y_train[inner_valid_ids]),
  levels = c("Immune-low", "Immune-high")
)

# 2. Призначаємо імена
names(y_inner_train) <- inner_train_ids
names(y_inner_valid) <- inner_valid_ids

# 3. Перевіряємо відповідність ID та PCA scores
label_alignment_qc <- c(
  training_label_count =
    length(y_inner_train) == nrow(train_pca_scores),
  
  validation_label_count =
    length(y_inner_valid) == nrow(valid_pca_scores),
  
  training_ID_order =
    identical(names(y_inner_train), rownames(train_pca_scores)),
  
  validation_ID_order =
    identical(names(y_inner_valid), rownames(valid_pca_scores)),
  
  no_missing_training_labels =
    !anyNA(y_inner_train),
  
  no_missing_validation_labels =
    !anyNA(y_inner_valid),
  
  both_training_classes =
    length(unique(y_inner_train)) == 2,
  
  both_validation_classes =
    length(unique(y_inner_valid)) == 2
)

cat("\n========== TRAINING LABEL DISTRIBUTION ==========\n")
print(table(y_inner_train))

cat("\n========== VALIDATION LABEL DISTRIBUTION ==========\n")
print(table(y_inner_valid))

cat("\n========== LABEL ALIGNMENT QC ==========\n")
print(label_alignment_qc)

stopifnot(all(label_alignment_qc))

cat("\nSTEP 10 COMPLETED SUCCESSFULLY\n")




###############################################################################
# STEP 11 — PILOT MLP TRAINING (CORRECTED ONE-HOT RESPONSE)
# Outer fold 5 / Inner fold 4
# No file writing
###############################################################################

set.seed(20261106)

# 1. Перевіряємо вхідні дані
stopifnot(
  nrow(train_pca_scores) == length(y_inner_train),
  nrow(valid_pca_scores) == length(y_inner_valid),
  identical(rownames(train_pca_scores), names(y_inner_train)),
  identical(rownames(valid_pca_scores), names(y_inner_valid))
)

# 2. Перетворюємо класи на one-hot матрицю
y_train_onehot <- nnet::class.ind(y_inner_train)

colnames(y_train_onehot) <- levels(y_inner_train)

# 3. Перевіряємо кодування
cat("\nOne-hot response dimensions:\n")
print(dim(y_train_onehot))

cat("\nOne-hot response columns:\n")
print(colnames(y_train_onehot))

stopifnot(
  nrow(y_train_onehot) == nrow(train_pca_scores),
  ncol(y_train_onehot) == 2,
  all(rowSums(y_train_onehot) == 1)
)

# 4. Навчаємо MLP
mlp_time <- system.time({
  
  mlp_pilot <- nnet::nnet(
    x = train_pca_scores,
    y = y_train_onehot,
    size = 5,
    decay = 0.1,
    maxit = 500,
    MaxNWts = 5000,
    trace = FALSE,
    softmax = TRUE
  )
  
})

cat("\nMLP training time:\n")
print(mlp_time)

# 5. Прогнозуємо ймовірності
train_prob <- predict(
  mlp_pilot,
  newdata = train_pca_scores,
  type = "raw"
)

valid_prob <- predict(
  mlp_pilot,
  newdata = valid_pca_scores,
  type = "raw"
)

# 6. Перевіряємо назви та розмірності
cat("\nTraining probability dimensions:\n")
print(dim(train_prob))

cat("\nValidation probability dimensions:\n")
print(dim(valid_prob))

cat("\nOutput class names:\n")
print(colnames(valid_prob))

# 7. QC прогнозованих імовірностей
mlp_probability_qc <- c(
  correct_training_rows =
    nrow(train_prob) == 222,
  
  correct_validation_rows =
    nrow(valid_prob) == 74,
  
  two_output_classes =
    ncol(valid_prob) == 2,
  
  correct_class_names =
    identical(colnames(valid_prob), levels(y_inner_train)),
  
  finite_training_probabilities =
    all(is.finite(train_prob)),
  
  finite_validation_probabilities =
    all(is.finite(valid_prob)),
  
  probabilities_in_range =
    all(valid_prob >= 0 & valid_prob <= 1),
  
  rows_sum_to_one =
    all(abs(rowSums(valid_prob) - 1) < 1e-6)
)

cat("\n========== MLP PROBABILITY QC ==========\n")
print(mlp_probability_qc)

stopifnot(all(mlp_probability_qc))

# 8. Класифікація validation-зразків
valid_predicted_class <- colnames(valid_prob)[
  max.col(valid_prob, ties.method = "first")
]

valid_predicted_class <- factor(
  valid_predicted_class,
  levels = levels(y_inner_valid)
)

cat("\nValidation confusion matrix:\n")
print(
  table(
    Observed = y_inner_valid,
    Predicted = valid_predicted_class
  )
)

# 9. Інформація про навчання
cat("\nMLP convergence information:\n")
cat("Iterations:", mlp_pilot$niter, "\n")
cat("Convergence code:", mlp_pilot$convergence, "\n")
cat("Final error:", mlp_pilot$value, "\n")

cat("\nSTEP 11 PILOT COMPLETED\n")


###############################################################################
# STEP 12 — STANDARDIZE PCA SCORES
# Training statistics only; no file writing
###############################################################################

# 1. Обчислюємо параметри масштабування на training
pca_means <- colMeans(train_pca_scores)
pca_sds <- apply(train_pca_scores, 2, sd)

stopifnot(
  all(is.finite(pca_means)),
  all(is.finite(pca_sds)),
  all(pca_sds > 0)
)

# 2. Масштабуємо training
train_pca_scaled <- scale(
  train_pca_scores,
  center = pca_means,
  scale = pca_sds
)

# 3. Масштабуємо validation тими самими параметрами
valid_pca_scaled <- scale(
  valid_pca_scores,
  center = pca_means,
  scale = pca_sds
)

# 4. Перевірка
scaling_qc <- c(
  training_dimensions_preserved =
    identical(dim(train_pca_scaled), dim(train_pca_scores)),
  
  validation_dimensions_preserved =
    identical(dim(valid_pca_scaled), dim(valid_pca_scores)),
  
  training_means_approximately_zero =
    all(abs(colMeans(train_pca_scaled)) < 1e-8),
  
  training_sds_approximately_one =
    all(abs(apply(train_pca_scaled, 2, sd) - 1) < 1e-8),
  
  validation_values_finite =
    all(is.finite(valid_pca_scaled)),
  
  training_IDs_preserved =
    identical(rownames(train_pca_scaled), inner_train_ids),
  
  validation_IDs_preserved =
    identical(rownames(valid_pca_scaled), inner_valid_ids)
)

cat("\n========== PCA SCALING QC ==========\n")
print(scaling_qc)

stopifnot(all(scaling_qc))

cat("\nTraining scaled dimensions:",
    paste(dim(train_pca_scaled), collapse = " x "), "\n")

cat("Validation scaled dimensions:",
    paste(dim(valid_pca_scaled), collapse = " x "), "\n")

cat("\nSTEP 12 COMPLETED\n")


###############################################################################
# STEP 13 — REUSABLE MLP FOLD PREPROCESSING FUNCTION
# No files are created or modified
###############################################################################

prepare_mlp_fold <- function(
    outer_fold,
    inner_fold,
    n_features = 5000,
    n_pcs = 30
) {
  
  # 1. Отримуємо структуру поточного fold
  fold_object <- nested_manifest$nested_cv[[outer_fold]]
  
  outer_train_idx <- fold_object$outer_train_idx
  inner_folds <- fold_object$inner_folds
  
  stopifnot(
    inner_fold >= 1,
    inner_fold <= length(inner_folds)
  )
  
  outer_ids <- train_ids[outer_train_idx]
  
  valid_idx <- inner_folds[[inner_fold]]
  
  inner_train_idx <- setdiff(
    seq_along(outer_ids),
    valid_idx
  )
  
  inner_train_ids <- outer_ids[inner_train_idx]
  inner_valid_ids <- outer_ids[valid_idx]
  
  # 2. Зіставляємо ID з beta_matrix
  train_pos <- match(
    inner_train_ids,
    colnames(beta_matrix)
  )
  
  valid_pos <- match(
    inner_valid_ids,
    colnames(beta_matrix)
  )
  
  stopifnot(
    !anyNA(train_pos),
    !anyNA(valid_pos)
  )
  
  # 3. Витягуємо training matrix
  X_train_fold <- beta_matrix[
    ,
    train_pos,
    drop = FALSE
  ]
  
  stopifnot(
    identical(colnames(X_train_fold), inner_train_ids)
  )
  
  # 4. Feature selection: training only
  cpg_variance <- DelayedMatrixStats::rowVars(
    X_train_fold,
    na.rm = TRUE
  )
  
  stopifnot(
    all(is.finite(cpg_variance)),
    sum(cpg_variance > 0) >= n_features
  )
  
  selected_idx <- order(
    cpg_variance,
    decreasing = TRUE
  )[seq_len(n_features)]
  
  # 5. Витягуємо вибрані CpG
  X_train_selected <- as.matrix(
    X_train_fold[
      selected_idx,
      ,
      drop = FALSE
    ]
  )
  
  X_valid_selected <- as.matrix(
    beta_matrix[
      selected_idx,
      valid_pos,
      drop = FALSE
    ]
  )
  
  # 6. Транспонуємо: samples x CpGs
  X_train_selected <- t(X_train_selected)
  X_valid_selected <- t(X_valid_selected)
  
  stopifnot(
    all(is.finite(X_train_selected)),
    all(is.finite(X_valid_selected))
  )
  
  # 7. PCA: fit тільки на training
  pca_fit <- prcomp(
    X_train_selected,
    center = TRUE,
    scale. = TRUE,
    rank. = n_pcs
  )
  
  train_scores <- pca_fit$x
  
  valid_scores <- predict(
    pca_fit,
    newdata = X_valid_selected
  )
  
  # 8. Стандартизація PCA scores: training only
  pca_means <- colMeans(train_scores)
  pca_sds <- apply(train_scores, 2, sd)
  
  stopifnot(
    all(is.finite(pca_means)),
    all(is.finite(pca_sds)),
    all(pca_sds > 0)
  )
  
  train_scaled <- scale(
    train_scores,
    center = pca_means,
    scale = pca_sds
  )
  
  valid_scaled <- scale(
    valid_scores,
    center = pca_means,
    scale = pca_sds
  )
  
  rownames(train_scaled) <- inner_train_ids
  rownames(valid_scaled) <- inner_valid_ids
  
  # 9. Отримуємо phenotype labels
  y_train_fold <- factor(
    as.character(y_train[inner_train_ids]),
    levels = c("Immune-low", "Immune-high")
  )
  
  y_valid_fold <- factor(
    as.character(y_train[inner_valid_ids]),
    levels = c("Immune-low", "Immune-high")
  )
  
  names(y_train_fold) <- inner_train_ids
  names(y_valid_fold) <- inner_valid_ids
  
  # 10. Final QC
  stopifnot(
    nrow(train_scaled) == length(y_train_fold),
    nrow(valid_scaled) == length(y_valid_fold),
    ncol(train_scaled) == n_pcs,
    ncol(valid_scaled) == n_pcs,
    identical(rownames(train_scaled), names(y_train_fold)),
    identical(rownames(valid_scaled), names(y_valid_fold)),
    !anyNA(y_train_fold),
    !anyNA(y_valid_fold),
    all(is.finite(train_scaled)),
    all(is.finite(valid_scaled))
  )
  
  return(list(
    train_x = train_scaled,
    valid_x = valid_scaled,
    train_y = y_train_fold,
    valid_y = y_valid_fold,
    selected_cpgs = rownames(beta_matrix)[selected_idx],
    pca = pca_fit,
    pca_means = pca_means,
    pca_sds = pca_sds,
    train_ids = inner_train_ids,
    valid_ids = inner_valid_ids
  ))
}

cat("Preprocessing function created successfully.\n")




test_fold <- prepare_mlp_fold(
  outer_fold = 5,
  inner_fold = 4,
  n_features = 5000,
  n_pcs = 30
)

cat("\nTraining dimensions:\n")
print(dim(test_fold$train_x))

cat("\nValidation dimensions:\n")
print(dim(test_fold$valid_x))

cat("\nTraining labels:\n")
print(table(test_fold$train_y))

cat("\nValidation labels:\n")
print(table(test_fold$valid_y))

cat("\nSelected CpGs:", length(test_fold$selected_cpgs), "\n")

cat("\n========== REUSABLE FUNCTION QC ==========\n")

function_qc <- c(
  training_N = nrow(test_fold$train_x) == 222,
  validation_N = nrow(test_fold$valid_x) == 74,
  training_PCs = ncol(test_fold$train_x) == 30,
  validation_PCs = ncol(test_fold$valid_x) == 30,
  selected_5000_CpGs = length(test_fold$selected_cpgs) == 5000,
  training_ID_alignment =
    identical(rownames(test_fold$train_x), test_fold$train_ids),
  validation_ID_alignment =
    identical(rownames(test_fold$valid_x), test_fold$valid_ids),
  training_label_alignment =
    identical(names(test_fold$train_y), test_fold$train_ids),
  validation_label_alignment =
    identical(names(test_fold$valid_y), test_fold$valid_ids)
)

print(function_qc)

stopifnot(all(function_qc))

cat("\nSTEP 13 COMPLETED SUCCESSFULLY\n")



###############################################################################
# STEP 14 — PILOT INNER CV FOR MLP
# Outer fold 5; one hyperparameter configuration
# No files are created or modified
###############################################################################

pilot_outer_fold <- 5

pilot_config <- list(
  n_pcs = 20,
  size = 5,
  decay = 0.1
)

pilot_results <- vector("list", 4)

for (j in seq_len(4)) {
  
  cat("\nProcessing Outer fold:",
      pilot_outer_fold,
      "| Inner fold:", j, "\n")
  
  # 1. Fold-specific preprocessing
  fold_data <- prepare_mlp_fold(
    outer_fold = pilot_outer_fold,
    inner_fold = j,
    n_features = 5000,
    n_pcs = 30
  )
  
  # 2. Use selected number of PCA components
  x_train <- fold_data$train_x[
    ,
    seq_len(pilot_config$n_pcs),
    drop = FALSE
  ]
  
  x_valid <- fold_data$valid_x[
    ,
    seq_len(pilot_config$n_pcs),
    drop = FALSE
  ]
  
  # 3. Convert training labels to one-hot matrix
  y_onehot <- nnet::class.ind(fold_data$train_y)
  colnames(y_onehot) <- levels(fold_data$train_y)
  
  # 4. Train MLP
  set.seed(20261105 + j)
  
  mlp_fit <- nnet::nnet(
    x = x_train,
    y = y_onehot,
    size = pilot_config$size,
    decay = pilot_config$decay,
    maxit = 500,
    MaxNWts = 5000,
    trace = FALSE,
    softmax = TRUE
  )
  
  # 5. Predict validation probabilities
  probabilities <- predict(
    mlp_fit,
    newdata = x_valid,
    type = "raw"
  )
  
  probabilities <- as.matrix(probabilities)
  
  # 6. Align probability columns with class levels
  probabilities <- probabilities[
    ,
    levels(fold_data$valid_y),
    drop = FALSE
  ]
  
  # 7. Calculate validation log-loss
  actual_class_index <- match(
    as.character(fold_data$valid_y),
    colnames(probabilities)
  )
  
  actual_probabilities <- probabilities[
    cbind(
      seq_len(nrow(probabilities)),
      actual_class_index
    )
  ]
  
  eps <- 1e-15
  
  fold_logloss <- -mean(
    log(pmax(actual_probabilities, eps))
  )
  
  # 8. Store results
  pilot_results[[j]] <- data.frame(
    Outer_Fold = pilot_outer_fold,
    Inner_Fold = j,
    N_train = nrow(x_train),
    N_valid = nrow(x_valid),
    PCA_components = pilot_config$n_pcs,
    Hidden_neurons = pilot_config$size,
    Decay = pilot_config$decay,
    LogLoss = fold_logloss,
    Convergence = mlp_fit$convergence,
    Probabilities_finite = all(is.finite(probabilities)),
    Probabilities_in_range =
      all(probabilities >= 0 & probabilities <= 1),
    Probability_rows_sum_to_one =
      all(abs(rowSums(probabilities) - 1) < 1e-6)
  )
  
  cat("Validation log-loss:",
      round(fold_logloss, 4), "\n")
  
  cat("Convergence code:",
      mlp_fit$convergence, "\n")
}

###############################################################################
# 9. Summarize pilot
###############################################################################

pilot_results_df <- do.call(
  rbind,
  pilot_results
)

cat("\n========== PILOT INNER CV RESULTS ==========\n")
print(pilot_results_df, row.names = FALSE)

cat("\nMean log-loss:",
    mean(pilot_results_df$LogLoss), "\n")

cat("SD log-loss:",
    sd(pilot_results_df$LogLoss), "\n")

###############################################################################
# 10. Quality control
###############################################################################

pilot_qc <- c(
  four_inner_folds = nrow(pilot_results_df) == 4,
  all_logloss_finite = all(is.finite(pilot_results_df$LogLoss)),
  all_converged = all(pilot_results_df$Convergence == 0),
  probabilities_finite =
    all(pilot_results_df$Probabilities_finite),
  probabilities_in_range =
    all(pilot_results_df$Probabilities_in_range),
  probability_rows_sum_to_one =
    all(pilot_results_df$Probability_rows_sum_to_one)
)

cat("\n========== PILOT QC ==========\n")
print(pilot_qc)

stopifnot(all(pilot_qc))

cat("\nSTEP 14 COMPLETED SUCCESSFULLY\n")


###############################################################################
# STEP 15 — FULL NESTED INNER GRID SEARCH FOR MLP
# 18 configurations x 20 inner folds = 360 model fits
# No holdout data used
# No files are created or modified
###############################################################################

mlp_grid <- expand.grid(
  n_pcs = c(10, 20, 30),
  size = c(3, 5, 8),
  decay = c(0.01, 0.1),
  KEEP.OUT.ATTRS = FALSE
)

mlp_grid$Config_ID <- seq_len(nrow(mlp_grid))

cat("Number of configurations:", nrow(mlp_grid), "\n")
print(mlp_grid)

mlp_grid_results <- list()
result_counter <- 1

###############################################################################
# OUTER / INNER FOLD LOOPS
###############################################################################

for (outer in seq_len(5)) {
  
  cat("\n\n========================================\n")
  cat("OUTER FOLD:", outer, "OF 5\n")
  cat("========================================\n")
  
  for (inner in seq_len(4)) {
    
    cat("\nPreparing Outer", outer,
        "| Inner", inner, "\n")
    
    # Fold-specific preprocessing:
    # feature selection, PCA and scaling use inner-training only
    fold_data <- prepare_mlp_fold(
      outer_fold = outer,
      inner_fold = inner,
      n_features = 5000,
      n_pcs = 30
    )
    
    y_onehot <- nnet::class.ind(fold_data$train_y)
    colnames(y_onehot) <- levels(fold_data$train_y)
    
    for (config_idx in seq_len(nrow(mlp_grid))) {
      
      config <- mlp_grid[config_idx, ]
      
      x_train <- fold_data$train_x[
        ,
        seq_len(config$n_pcs),
        drop = FALSE
      ]
      
      x_valid <- fold_data$valid_x[
        ,
        seq_len(config$n_pcs),
        drop = FALSE
      ]
      
      # Deterministic seed for each fold/configuration
      set.seed(
        20261105 +
          outer * 10000 +
          inner * 1000 +
          config$Config_ID
      )
      
      fit_result <- tryCatch({
        
        fit <- nnet::nnet(
          x = x_train,
          y = y_onehot,
          size = config$size,
          decay = config$decay,
          maxit = 1000,
          MaxNWts = 5000,
          trace = FALSE,
          softmax = TRUE
        )
        
        probabilities <- predict(
          fit,
          newdata = x_valid,
          type = "raw"
        )
        
        probabilities <- as.matrix(probabilities)
        
        probabilities <- probabilities[
          ,
          levels(fold_data$valid_y),
          drop = FALSE
        ]
        
        actual_class_index <- match(
          as.character(fold_data$valid_y),
          colnames(probabilities)
        )
        
        actual_probabilities <- probabilities[
          cbind(
            seq_len(nrow(probabilities)),
            actual_class_index
          )
        ]
        
        logloss <- -mean(
          log(pmax(actual_probabilities, 1e-15))
        )
        
        list(
          LogLoss = logloss,
          Convergence = fit$convergence,
          Probabilities_finite =
            all(is.finite(probabilities)),
          Probabilities_in_range =
            all(probabilities >= 0 &
                  probabilities <= 1),
          Probability_rows_sum_to_one =
            all(abs(rowSums(probabilities) - 1) < 1e-6),
          Error = NA_character_
        )
        
      }, error = function(e) {
        
        list(
          LogLoss = NA_real_,
          Convergence = NA_integer_,
          Probabilities_finite = FALSE,
          Probabilities_in_range = FALSE,
          Probability_rows_sum_to_one = FALSE,
          Error = conditionMessage(e)
        )
      })
      
      mlp_grid_results[[result_counter]] <- data.frame(
        Outer_Fold = outer,
        Inner_Fold = inner,
        Config_ID = config$Config_ID,
        PCA_components = config$n_pcs,
        Hidden_neurons = config$size,
        Decay = config$decay,
        N_train = nrow(x_train),
        N_valid = nrow(x_valid),
        LogLoss = fit_result$LogLoss,
        Convergence = fit_result$Convergence,
        Probabilities_finite =
          fit_result$Probabilities_finite,
        Probabilities_in_range =
          fit_result$Probabilities_in_range,
        Probability_rows_sum_to_one =
          fit_result$Probability_rows_sum_to_one,
        Error = fit_result$Error,
        stringsAsFactors = FALSE
      )
      
      result_counter <- result_counter + 1
    }
    
    cat("Completed Outer", outer,
        "| Inner", inner,
        "| 18 configurations\n")
  }
}

###############################################################################
# COMBINE RESULTS
###############################################################################

mlp_grid_results_df <- do.call(
  rbind,
  mlp_grid_results
)

rownames(mlp_grid_results_df) <- NULL

cat("\n========== GRID SEARCH SUMMARY ==========\n")

cat("Total model fits:",
    nrow(mlp_grid_results_df), "\n")

cat("Finite log-loss values:",
    sum(is.finite(mlp_grid_results_df$LogLoss)), "\n")

cat("Converged models:",
    sum(mlp_grid_results_df$Convergence == 0, na.rm = TRUE),
    "\n")

cat("Failed fits:",
    sum(!is.na(mlp_grid_results_df$Error)), "\n")

###############################################################################
# SUMMARIZE EACH CONFIGURATION
###############################################################################

config_summaries <- lapply(
  split(mlp_grid_results_df, mlp_grid_results_df$Config_ID),
  function(df) {
    
    data.frame(
      Config_ID = df$Config_ID[1],
      PCA_components = df$PCA_components[1],
      Hidden_neurons = df$Hidden_neurons[1],
      Decay = df$Decay[1],
      Folds_completed = sum(is.finite(df$LogLoss)),
      Converged_folds = sum(df$Convergence == 0, na.rm = TRUE),
      Mean_logloss = if (all(is.finite(df$LogLoss))) {
        mean(df$LogLoss)
      } else {
        NA_real_
      },
      SD_logloss = if (all(is.finite(df$LogLoss))) {
        sd(df$LogLoss)
      } else {
        NA_real_
      },
      All_folds_converged =
        all(df$Convergence == 0),
      All_probabilities_valid =
        all(df$Probabilities_finite) &&
        all(df$Probabilities_in_range) &&
        all(df$Probability_rows_sum_to_one)
    )
  }
)

mlp_config_summary <- do.call(
  rbind,
  config_summaries
)

mlp_config_summary <- mlp_config_summary[
  order(mlp_config_summary$Mean_logloss),
]

rownames(mlp_config_summary) <- NULL

cat("\n========== CONFIGURATION PERFORMANCE ==========\n")
print(mlp_config_summary, row.names = FALSE)

###############################################################################
# QUALITY CONTROL
###############################################################################

grid_qc <- c(
  expected_360_fits = nrow(mlp_grid_results_df) == 360,
  all_20_inner_folds_present =
    length(unique(
      paste(
        mlp_grid_results_df$Outer_Fold,
        mlp_grid_results_df$Inner_Fold
      )
    )) == 20,
  all_18_configs_per_fold =
    all(table(
      mlp_grid_results_df$Outer_Fold,
      mlp_grid_results_df$Inner_Fold
    ) == 18),
  all_logloss_finite =
    all(is.finite(mlp_grid_results_df$LogLoss)),
  all_models_converged =
    all(mlp_grid_results_df$Convergence == 0),
  all_probabilities_valid =
    all(mlp_grid_results_df$Probabilities_finite) &&
    all(mlp_grid_results_df$Probabilities_in_range) &&
    all(mlp_grid_results_df$Probability_rows_sum_to_one),
  all_configurations_have_20_folds =
    all(mlp_config_summary$Folds_completed == 20)
)

cat("\n========== GRID SEARCH QC ==========\n")
print(grid_qc)

cat("\nSTEP 15 FINISHED\n")


###############################################################################
# STEP 16 — OUTER-SPECIFIC HYPERPARAMETER SELECTION
# Each outer fold uses ONLY its own 4 inner folds
# No files are created or modified
###############################################################################

outer_config_summaries <- list()
selected_configs_by_outer <- list()

for (outer in seq_len(5)) {
  
  outer_results <- subset(
    mlp_grid_results_df,
    Outer_Fold == outer
  )
  
  outer_summary <- do.call(
    rbind,
    lapply(
      split(outer_results, outer_results$Config_ID),
      function(df) {
        
        data.frame(
          Outer_Fold = outer,
          Config_ID = df$Config_ID[1],
          PCA_components = df$PCA_components[1],
          Hidden_neurons = df$Hidden_neurons[1],
          Decay = df$Decay[1],
          Inner_folds = sum(is.finite(df$LogLoss)),
          Mean_logloss = mean(df$LogLoss),
          SD_logloss = sd(df$LogLoss),
          All_converged =
            all(df$Convergence == 0),
          All_probabilities_valid =
            all(df$Probabilities_finite) &&
            all(df$Probabilities_in_range) &&
            all(df$Probability_rows_sum_to_one)
        )
      }
    )
  )
  
  outer_summary <- outer_summary[
    order(outer_summary$Mean_logloss),
  ]
  
  rownames(outer_summary) <- NULL
  
  outer_config_summaries[[outer]] <- outer_summary
  
  selected_configs_by_outer[[outer]] <- outer_summary[1, ]
  
  cat("\n====================================\n")
  cat("OUTER FOLD:", outer, "\n")
  cat("====================================\n")
  
  print(head(outer_summary, 5), row.names = FALSE)
  
  cat("\nSelected configuration:\n")
  print(outer_summary[1, ], row.names = FALSE)
}

###############################################################################
# COMBINE SELECTED CONFIGURATIONS
###############################################################################

selected_outer_configs <- do.call(
  rbind,
  selected_configs_by_outer
)

rownames(selected_outer_configs) <- NULL

cat("\n========== SELECTED CONFIGURATIONS ==========\n")
print(selected_outer_configs, row.names = FALSE)

###############################################################################
# QC
###############################################################################

outer_selection_qc <- c(
  five_outer_folds =
    nrow(selected_outer_configs) == 5,
  
  four_inner_folds_each =
    all(vapply(
      outer_config_summaries,
      function(x) all(x$Inner_folds == 4),
      logical(1)
    )),
  
  all_selected_logloss_finite =
    all(is.finite(selected_outer_configs$Mean_logloss)),
  
  all_selected_models_converged =
    all(selected_outer_configs$All_converged),
  
  all_selected_probabilities_valid =
    all(selected_outer_configs$All_probabilities_valid),
  
  one_configuration_per_outer =
    length(unique(selected_outer_configs$Outer_Fold)) == 5
)

cat("\n========== OUTER SELECTION QC ==========\n")
print(outer_selection_qc)

stopifnot(all(outer_selection_qc))

cat("\nSTEP 16 COMPLETED SUCCESSFULLY\n")





###############################################################################
# STEP 17 — OUTER-FOLD MLP TRAINING AND OOF PREDICTIONS
# No holdout usage
# No files are created or modified
###############################################################################

# 1. Ініціалізація OOF prediction matrix
class_levels <- c("Immune-low", "Immune-high")

outer_oof_probabilities <- matrix(
  NA_real_,
  nrow = length(train_ids),
  ncol = length(class_levels),
  dimnames = list(train_ids, class_levels)
)

outer_oof_predictions <- rep(
  NA_character_,
  length(train_ids)
)

outer_fold_assignment <- rep(
  NA_integer_,
  length(train_ids)
)

outer_model_diagnostics <- vector("list", 5)

###############################################################################
# 2. Outer-fold loop
###############################################################################

for (outer in seq_len(5)) {
  
  cat("\n========================================\n")
  cat("TRAINING OUTER FOLD:", outer, "OF 5\n")
  cat("========================================\n")
  
  fold_object <- nested_manifest$nested_cv[[outer]]
  
  outer_train_idx <- fold_object$outer_train_idx
  outer_valid_idx <- fold_object$outer_valid_idx
  
  outer_train_ids <- train_ids[outer_train_idx]
  outer_valid_ids <- train_ids[outer_valid_idx]
  
  # Selected hyperparameters for this outer fold
  config <- selected_outer_configs[
    selected_outer_configs$Outer_Fold == outer,
  ]
  
  n_pcs <- config$PCA_components
  hidden_neurons <- config$Hidden_neurons
  decay_value <- config$Decay
  
  cat(
    "Selected configuration:",
    n_pcs, "PCs |",
    hidden_neurons, "neurons | decay",
    decay_value, "\n"
  )
  
  # 3. Match sample IDs to beta matrix
  train_pos <- match(
    outer_train_ids,
    colnames(beta_matrix)
  )
  
  valid_pos <- match(
    outer_valid_ids,
    colnames(beta_matrix)
  )
  
  stopifnot(
    !anyNA(train_pos),
    !anyNA(valid_pos)
  )
  
  # 4. Extract outer-training matrix
  X_outer_train <- beta_matrix[
    ,
    train_pos,
    drop = FALSE
  ]
  
  # 5. Feature selection on outer-training only
  outer_variance <- DelayedMatrixStats::rowVars(
    X_outer_train,
    na.rm = TRUE
  )
  
  stopifnot(
    all(is.finite(outer_variance)),
    sum(outer_variance > 0) >= 5000
  )
  
  selected_cpg_idx <- order(
    outer_variance,
    decreasing = TRUE
  )[seq_len(5000)]
  
  # 6. Extract selected CpGs
  X_train_selected <- as.matrix(
    X_outer_train[
      selected_cpg_idx,
      ,
      drop = FALSE
    ]
  )
  
  X_valid_selected <- as.matrix(
    beta_matrix[
      selected_cpg_idx,
      valid_pos,
      drop = FALSE
    ]
  )
  
  X_train_selected <- t(X_train_selected)
  X_valid_selected <- t(X_valid_selected)
  
  stopifnot(
    all(is.finite(X_train_selected)),
    all(is.finite(X_valid_selected))
  )
  
  # 7. PCA fit on outer-training only
  pca_fit <- prcomp(
    X_train_selected,
    center = TRUE,
    scale. = TRUE,
    rank. = n_pcs
  )
  
  train_scores <- pca_fit$x
  
  valid_scores <- predict(
    pca_fit,
    newdata = X_valid_selected
  )
  
  # 8. Scale PCA scores using outer-training parameters
  pca_means <- colMeans(train_scores)
  pca_sds <- apply(train_scores, 2, sd)
  
  stopifnot(
    all(is.finite(pca_means)),
    all(is.finite(pca_sds)),
    all(pca_sds > 0)
  )
  
  train_scaled <- scale(
    train_scores,
    center = pca_means,
    scale = pca_sds
  )
  
  valid_scaled <- scale(
    valid_scores,
    center = pca_means,
    scale = pca_sds
  )
  
  rownames(train_scaled) <- outer_train_ids
  rownames(valid_scaled) <- outer_valid_ids
  
  # 9. Training labels
  y_outer_train <- factor(
    as.character(y_train[outer_train_ids]),
    levels = class_levels
  )
  
  y_outer_valid <- factor(
    as.character(y_train[outer_valid_ids]),
    levels = class_levels
  )
  
  names(y_outer_train) <- outer_train_ids
  names(y_outer_valid) <- outer_valid_ids
  
  stopifnot(
    identical(rownames(train_scaled), names(y_outer_train)),
    identical(rownames(valid_scaled), names(y_outer_valid)),
    !anyNA(y_outer_train),
    !anyNA(y_outer_valid)
  )
  
  y_onehot <- nnet::class.ind(y_outer_train)
  colnames(y_onehot) <- class_levels
  
  # 10. Fit final model for this outer fold
  set.seed(20261105 + outer * 100000)
  
  mlp_outer_fit <- nnet::nnet(
    x = train_scaled[
      ,
      seq_len(n_pcs),
      drop = FALSE
    ],
    y = y_onehot,
    size = hidden_neurons,
    decay = decay_value,
    maxit = 1000,
    MaxNWts = 5000,
    trace = FALSE,
    softmax = TRUE
  )
  
  # 11. Predict outer-validation probabilities
  outer_prob <- predict(
    mlp_outer_fit,
    newdata = valid_scaled[
      ,
      seq_len(n_pcs),
      drop = FALSE
    ],
    type = "raw"
  )
  
  outer_prob <- as.matrix(outer_prob)
  
  outer_prob <- outer_prob[
    ,
    class_levels,
    drop = FALSE
  ]
  
  # 12. Predicted class
  outer_pred <- class_levels[
    max.col(outer_prob, ties.method = "first")
  ]
  
  # 13. Store OOF predictions in original training order
  outer_oof_probabilities[
    outer_valid_ids,
  ] <- outer_prob
  
  outer_oof_predictions[
    match(outer_valid_ids, train_ids)
  ] <- outer_pred
  
  outer_fold_assignment[
    outer_valid_idx
  ] <- outer
  
  # 14. Diagnostics
  outer_model_diagnostics[[outer]] <- data.frame(
    Outer_Fold = outer,
    N_train = length(outer_train_ids),
    N_valid = length(outer_valid_ids),
    PCA_components = n_pcs,
    Hidden_neurons = hidden_neurons,
    Decay = decay_value,
    Convergence = mlp_outer_fit$convergence,
    Training_IDs_aligned =
      identical(rownames(train_scaled), outer_train_ids),
    Validation_IDs_aligned =
      identical(rownames(valid_scaled), outer_valid_ids),
    Probabilities_finite =
      all(is.finite(outer_prob)),
    Probabilities_in_range =
      all(outer_prob >= 0 & outer_prob <= 1),
    Probability_rows_sum_to_one =
      all(abs(rowSums(outer_prob) - 1) < 1e-6)
  )
  
  cat(
    "Outer fold", outer,
    "completed | convergence:",
    mlp_outer_fit$convergence,
    "| validation samples:",
    length(outer_valid_ids), "\n"
  )
}

###############################################################################
# 15. Combine diagnostics
###############################################################################

outer_model_diagnostics_df <- do.call(
  rbind,
  outer_model_diagnostics
)

cat("\n========== OUTER MODEL DIAGNOSTICS ==========\n")
print(outer_model_diagnostics_df, row.names = FALSE)

###############################################################################
# 16. OOF prediction QC
###############################################################################

outer_oof_qc <- c(
  all_368_samples_predicted =
    all(!is.na(outer_oof_predictions)),
  
  all_probabilities_present =
    all(is.finite(outer_oof_probabilities)),
  
  exactly_368_predictions =
    length(outer_oof_predictions) == 368,
  
  all_sample_IDs_preserved =
    identical(
      rownames(outer_oof_probabilities),
      train_ids
    ),
  
  all_samples_assigned_to_outer_fold =
    all(!is.na(outer_fold_assignment)),
  
  each_sample_predicted_once =
    sum(!is.na(outer_oof_predictions)) == 368,
  
  all_models_converged =
    all(outer_model_diagnostics_df$Convergence == 0),
  
  all_probabilities_valid =
    all(outer_model_diagnostics_df$Probabilities_finite) &&
    all(outer_model_diagnostics_df$Probabilities_in_range) &&
    all(outer_model_diagnostics_df$Probability_rows_sum_to_one)
)

cat("\n========== OUTER OOF QC ==========\n")
print(outer_oof_qc)

stopifnot(all(outer_oof_qc))

###############################################################################
# 17. OOF prediction overview
###############################################################################

cat("\nOOF prediction dimensions:\n")
print(dim(outer_oof_probabilities))

cat("\nPredicted class distribution:\n")
print(table(outer_oof_predictions))

cat("\nActual class distribution:\n")
print(table(y_train))

cat("\nSTEP 17 COMPLETED SUCCESSFULLY\n")




###############################################################################
# STEP 18 — MLP OUTER OOF PERFORMANCE METRICS
# Evaluation on 368 out-of-fold predictions
# No holdout usage
# No files are created or modified
###############################################################################

# 1. Align actual and predicted classes
actual <- factor(
  as.character(y_train),
  levels = c("Immune-low", "Immune-high")
)

predicted <- factor(
  outer_oof_predictions,
  levels = c("Immune-low", "Immune-high")
)

prob_high <- outer_oof_probabilities[
  train_ids,
  "Immune-high"
]

stopifnot(
  length(actual) == 368,
  length(predicted) == 368,
  length(prob_high) == 368,
  !anyNA(actual),
  !anyNA(predicted),
  all(is.finite(prob_high))
)

# 2. Confusion matrix
confusion <- table(
  Observed = actual,
  Predicted = predicted
)

cat("\n========== MLP OOF CONFUSION MATRIX ==========\n")
print(confusion)

TN <- confusion["Immune-low", "Immune-low"]
FP <- confusion["Immune-low", "Immune-high"]
FN <- confusion["Immune-high", "Immune-low"]
TP <- confusion["Immune-high", "Immune-high"]

# 3. Classification metrics
accuracy <- (TP + TN) / sum(confusion)

sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)

precision <- if ((TP + FP) > 0) {
  TP / (TP + FP)
} else {
  NA_real_
}

f1_score <- if (is.finite(precision + sensitivity) &&
                (precision + sensitivity) > 0) {
  2 * precision * sensitivity /
    (precision + sensitivity)
} else {
  NA_real_
}

balanced_accuracy <- mean(
  c(sensitivity, specificity)
)

# 4. ROC-AUC using rank statistic
actual_binary <- as.integer(
  actual == "Immune-high"
)

auc <- (
  sum(rank(prob_high, ties.method = "average")[actual_binary == 1]) -
    sum(actual_binary == 1) *
    (sum(actual_binary == 1) + 1) / 2
) / (
  sum(actual_binary == 1) *
    sum(actual_binary == 0)
)

# 5. Brier score
brier_score <- mean(
  (prob_high - actual_binary)^2
)

# 6. Log-loss
eps <- 1e-15

logloss <- -mean(
  actual_binary * log(pmax(prob_high, eps)) +
    (1 - actual_binary) *
    log(pmax(1 - prob_high, eps))
)

# 7. Metrics table
mlp_oof_metrics <- data.frame(
  Metric = c(
    "Accuracy",
    "Sensitivity_Immune_high",
    "Specificity_Immune_low",
    "Precision_Immune_high",
    "F1_score",
    "Balanced_accuracy",
    "ROC_AUC",
    "Brier_score",
    "Log_loss"
  ),
  Value = c(
    accuracy,
    sensitivity,
    specificity,
    precision,
    f1_score,
    balanced_accuracy,
    auc,
    brier_score,
    logloss
  )
)

cat("\n========== MLP OUTER OOF METRICS ==========\n")
print(mlp_oof_metrics, row.names = FALSE)

# 8. Fold-specific accuracy and AUC
fold_metrics <- lapply(
  sort(unique(outer_fold_assignment)),
  function(fold) {
    
    idx <- which(outer_fold_assignment == fold)
    
    fold_actual <- actual[idx]
    fold_pred <- predicted[idx]
    fold_prob <- prob_high[idx]
    
    fold_binary <- as.integer(
      fold_actual == "Immune-high"
    )
    
    fold_auc <- (
      sum(
        rank(fold_prob, ties.method = "average")[
          fold_binary == 1
        ]
      ) -
        sum(fold_binary == 1) *
        (sum(fold_binary == 1) + 1) / 2
    ) / (
      sum(fold_binary == 1) *
        sum(fold_binary == 0)
    )
    
    data.frame(
      Outer_Fold = fold,
      N = length(idx),
      Accuracy = mean(fold_pred == fold_actual),
      AUC = fold_auc,
      LogLoss = -mean(
        fold_binary * log(pmax(fold_prob, eps)) +
          (1 - fold_binary) *
          log(pmax(1 - fold_prob, eps))
      )
    )
  }
)

mlp_fold_metrics <- do.call(rbind, fold_metrics)

cat("\n========== FOLD-SPECIFIC METRICS ==========\n")
print(mlp_fold_metrics, row.names = FALSE)

# 9. QC
mlp_metrics_qc <- c(
  confusion_matrix_total_368 = sum(confusion) == 368,
  all_metrics_finite = all(is.finite(mlp_oof_metrics$Value)),
  auc_in_valid_range = auc >= 0 && auc <= 1,
  probabilities_in_valid_range =
    all(prob_high >= 0 & prob_high <= 1),
  five_outer_folds_reported =
    nrow(mlp_fold_metrics) == 5
)

cat("\n========== METRICS QC ==========\n")
print(mlp_metrics_qc)

stopifnot(all(mlp_metrics_qc))

cat("\nSTEP 18 COMPLETED SUCCESSFULLY\n")


###############################################################################
# STEP 19 — LOCATE ELASTIC NET OOF PREDICTIONS
# Read-only inventory
###############################################################################

cat("\n========== ELASTIC NET OBJECTS IN R SESSION ==========\n")

en_objects <- ls(envir = .GlobalEnv)

en_objects <- en_objects[
  grepl(
    "elastic|glmnet|nested|oof|prediction",
    en_objects,
    ignore.case = TRUE
  )
]

print(en_objects)

cat("\n========== RELEVANT MODEL FILES ==========\n")

model_files <- list.files(
  "results/models",
  pattern = "\\.(rds|RDS)$",
  full.names = TRUE,
  recursive = TRUE
)

relevant_files <- model_files[
  grepl(
    "elastic|glmnet|nested|oof|prediction",
    basename(model_files),
    ignore.case = TRUE
  )
]

print(data.frame(
  File = relevant_files,
  Size_MB = round(
    file.info(relevant_files)$size / 1024^2,
    3
  ),
  row.names = NULL
))

cat("\nSTEP 19 INVENTORY COMPLETED\n")

###############################################################################
# STEP 21 — FAIR OOF COMPARISON: MLP VS ELASTIC NET
# Same 368 samples, same labels, same outer folds
# Read-only
###############################################################################

# 1. Prepare Elastic Net OOF data
en_oof <- oof_pred

required_columns <- c(
  "Sample_ID",
  "Actual",
  "Predicted",
  "Probability_Immune_high",
  "Outer_Fold"
)

stopifnot(
  all(required_columns %in% colnames(en_oof)),
  nrow(en_oof) == 368,
  !anyDuplicated(en_oof$Sample_ID),
  all(is.finite(en_oof$Probability_Immune_high))
)

# 2. Align Elastic Net to MLP training ID order
en_match <- match(
  train_ids,
  en_oof$Sample_ID
)

stopifnot(
  !anyNA(en_match),
  length(unique(en_match)) == 368
)

en_oof_aligned <- en_oof[en_match, ]

# 3. Verify identical samples and labels
mlp_actual <- as.character(y_train)

en_actual <- as.character(
  en_oof_aligned$Actual
)

mlp_fold <- as.integer(
  outer_fold_assignment
)

en_fold <- as.integer(
  en_oof_aligned$Outer_Fold
)

comparison_qc <- c(
  same_sample_IDs = identical(
    train_ids,
    en_oof_aligned$Sample_ID
  ),
  same_actual_labels = identical(
    mlp_actual,
    en_actual
  ),
  same_outer_folds = identical(
    mlp_fold,
    en_fold
  ),
  MLP_predictions_complete =
    length(outer_oof_predictions) == 368,
  EN_predictions_complete =
    length(en_oof_aligned$Predicted) == 368,
  EN_probabilities_valid = all(
    en_oof_aligned$Probability_Immune_high >= 0 &
      en_oof_aligned$Probability_Immune_high <= 1
  )
)

cat("\n========== FAIR COMPARISON QC ==========\n")
print(comparison_qc)

stopifnot(all(comparison_qc))

# 4. Common metric calculation function
calculate_metrics <- function(
    actual,
    predicted,
    probability_high
) {
  
  actual <- factor(
    actual,
    levels = c("Immune-low", "Immune-high")
  )
  
  predicted <- factor(
    predicted,
    levels = c("Immune-low", "Immune-high")
  )
  
  binary <- as.integer(
    actual == "Immune-high"
  )
  
  cm <- table(
    Observed = actual,
    Predicted = predicted
  )
  
  TN <- cm["Immune-low", "Immune-low"]
  FP <- cm["Immune-low", "Immune-high"]
  FN <- cm["Immune-high", "Immune-low"]
  TP <- cm["Immune-high", "Immune-high"]
  
  accuracy <- (TP + TN) / sum(cm)
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  
  precision <- if ((TP + FP) > 0) {
    TP / (TP + FP)
  } else {
    NA_real_
  }
  
  f1 <- if (
    is.finite(precision + sensitivity) &&
    (precision + sensitivity) > 0
  ) {
    2 * precision * sensitivity /
      (precision + sensitivity)
  } else {
    NA_real_
  }
  
  ranks <- rank(
    probability_high,
    ties.method = "average"
  )
  
  auc <- (
    sum(ranks[binary == 1]) -
      sum(binary == 1) *
      (sum(binary == 1) + 1) / 2
  ) / (
    sum(binary == 1) *
      sum(binary == 0)
  )
  
  eps <- 1e-15
  
  brier <- mean(
    (probability_high - binary)^2
  )
  
  logloss <- -mean(
    binary * log(pmax(probability_high, eps)) +
      (1 - binary) *
      log(pmax(1 - probability_high, eps))
  )
  
  data.frame(
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1 = f1,
    Balanced_Accuracy = mean(
      c(sensitivity, specificity)
    ),
    ROC_AUC = auc,
    Brier_Score = brier,
    Log_Loss = logloss
  )
}

# 5. Calculate pooled OOF metrics
mlp_comparison_metrics <- calculate_metrics(
  actual = mlp_actual,
  predicted = as.character(outer_oof_predictions),
  probability_high = as.numeric(
    outer_oof_probabilities[
      train_ids,
      "Immune-high"
    ]
  )
)

en_comparison_metrics <- calculate_metrics(
  actual = en_actual,
  predicted = as.character(
    en_oof_aligned$Predicted
  ),
  probability_high = en_oof_aligned$Probability_Immune_high
)

comparison_table <- rbind(
  MLP = mlp_comparison_metrics,
  Elastic_Net = en_comparison_metrics
)

cat("\n========== POOLED OOF MODEL COMPARISON ==========\n")
print(round(comparison_table, 4))

# 6. Calculate metrics separately for each outer fold
fold_comparison <- do.call(
  rbind,
  lapply(sort(unique(mlp_fold)), function(fold) {
    
    idx <- which(mlp_fold == fold)
    
    mlp_fold_result <- calculate_metrics(
      actual = mlp_actual[idx],
      predicted = as.character(
        outer_oof_predictions[idx]
      ),
      probability_high = as.numeric(
        outer_oof_probabilities[
          train_ids[idx],
          "Immune-high"
        ]
      )
    )
    
    en_fold_result <- calculate_metrics(
      actual = en_actual[idx],
      predicted = as.character(
        en_oof_aligned$Predicted[idx]
      ),
      probability_high =
        en_oof_aligned$Probability_Immune_high[idx]
    )
    
    rbind(
      data.frame(
        Outer_Fold = fold,
        Model = "MLP",
        N = length(idx),
        mlp_fold_result
      ),
      data.frame(
        Outer_Fold = fold,
        Model = "Elastic_Net",
        N = length(idx),
        en_fold_result
      )
    )
  })
)

cat("\n========== OUTER FOLD COMPARISON ==========\n")

numeric_cols <- vapply(
  fold_comparison,
  is.numeric,
  logical(1)
)

fold_comparison_display <- fold_comparison

fold_comparison_display[numeric_cols] <- lapply(
  fold_comparison_display[numeric_cols],
  round,
  digits = 4
)

print(
  fold_comparison_display,
  row.names = FALSE
)

cat("\nSTEP 21 DISPLAY COMPLETED\n")




