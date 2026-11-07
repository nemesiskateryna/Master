
###############################################################################
# SCRIPT 21
# TCGA-SKCM: NESTED CROSS-VALIDATION — INNER MODEL TUNING
#
# Objective:
# Tune Elastic Net regularization using inner cross-validation.
#
# Pipeline:
#   1. Load existing training data and nested CV manifest.
#   2. Load checkpoint with previously completed folds.
#   3. Skip folds that have already passed QC.
#   4. Perform variance-based CpG selection within each inner-training fold.
#   5. Fit Elastic Net models (alpha = 0.5).
#   6. Evaluate 50 lambda positions using validation log-loss.
#   7. Save checkpoint and QC after every completed fold.
#
# IMPORTANT:
#   - Outer validation samples are never used in inner tuning.
#   - Feature selection is performed separately within each inner fold.
#   - The independent test set remains untouched.
###############################################################################


###############################################################################
# STEP 1. INITIALIZATION
###############################################################################


# 1. Робоча директорія
project_dir <- "/Users/katie/master/Master"

setwd(project_dir)

cat("Working directory:", getwd(), "\n")

# 2. Необхідні пакети
required_packages <- c(
  "HDF5Array",
  "DelayedArray",
  "DelayedMatrixStats",
  "glmnet",
  "pROC"
)

package_status <- vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

cat("\nPACKAGE STATUS:\n")
print(package_status)

if (!all(package_status)) {
  stop("Some required packages are missing.")
}

# 3. Завантаження пакетів
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
library(glmnet)
library(pROC)

# 4. Основні шляхи
integrated_file <- file.path(
  "data/TCGA/processed",
  "TCGA_SKCM_integrated_multiomics_KNN_immune.rds"
)

split_file <- file.path(
  "results/models",
  "TCGA_SKCM_ML_train_test_split.rds"
)

manifest_file <- file.path(
  "results/models",
  "TCGA_SKCM_nested_CV_manifest.rds"
)

checkpoint_file <- file.path(
  "results/models/nested_CV_inner",
  "TCGA_SKCM_inner_CV_partial_results.rds"
)

progress_file <- file.path(
  "results/tables",
  "TCGA_SKCM_inner_CV_progress_QC.csv"
)

# 5. Перевірка файлів
required_files <- c(
  integrated_file,
  split_file,
  manifest_file,
  checkpoint_file
)

file_status <- file.exists(required_files)

names(file_status) <- basename(required_files)

cat("\nFILE STATUS:\n")
print(file_status)

if (!all(file_status)) {
  stop("Required input or checkpoint file is missing.")
}

# 6. Кількість CpG та фолдів
n_features <- 5000
n_lambda <- 50
alpha_value <- 0.5

cat("\nINITIALIZATION PARAMETERS:\n")
cat("Selected CpG:", n_features, "\n")
cat("Lambda path:", n_lambda, "\n")
cat("Alpha:", alpha_value, "\n")

cat("\nSTAGE 1 INITIALIZATION: PASS\n")


###############################################################################
# ЕТАП 2. ЗАВАНТАЖЕННЯ ДАНИХ ТА ПЕРЕВІРКА СТРУКТУРИ
###############################################################################

# 1. Завантаження інтегрованого об'єкта
integrated <- readRDS(integrated_file)

# 2. Завантаження train/test split
split_manifest <- readRDS(split_file)

# 3. Завантаження nested CV manifest
manifest_check <- readRDS(manifest_file)

# 4. Перевірка компонентів інтегрованого об'єкта
cat("\nINTEGRATED OBJECT COMPONENTS:\n")
print(names(integrated))

# 5. Отримання матриці метилювання
X <- integrated$beta_matrix

# 6. Перевірка розмірності
cat("\nBETA MATRIX DIMENSIONS:\n")
print(dim(X))

# 7. Перевірка immune scores
cat("\nIMMUNE SCORES DIMENSIONS:\n")
print(dim(integrated$immune_scores))

# 8. Перевірка train/test manifest
cat("\nSPLIT MANIFEST COMPONENTS:\n")
print(names(split_manifest))

# 9. Перевірка nested manifest
cat("\nNESTED MANIFEST COMPONENTS:\n")
print(names(manifest_check))

cat("\nNUMBER OF OUTER FOLDS:\n")
print(length(manifest_check$nested_cv))

cat("\nINNER FOLDS PER OUTER FOLD:\n")
print(vapply(
  manifest_check$nested_cv,
  function(x) length(x$inner_folds),
  integer(1)
))

# 10. Базові перевірки
stopifnot(
  nrow(X) == 409556,
  ncol(X) == 460,
  length(manifest_check$nested_cv) == 5,
  all(vapply(
    manifest_check$nested_cv,
    function(x) length(x$inner_folds) == 4,
    logical(1)
  ))
)

cat("\nSTAGE 2 DATA LOADING: PASS\n")


###############################################################################
# ЕТАП 3. ФОРМУВАННЯ TRAIN / TEST DATA
###############################################################################

# 1. Отримання індексів
train_idx <- split_manifest$train_indices
test_idx  <- split_manifest$test_indices

# 2. Перевірка кількості зразків
cat("Training samples:", length(train_idx), "\n")
cat("Test samples:", length(test_idx), "\n")

# 3. Формування матриць
X_train <- X[, train_idx, drop = FALSE]
X_test  <- X[, test_idx, drop = FALSE]

# 4. Формування immune phenotype labels
immune_scores <- integrated$immune_scores

y <- factor(
  immune_scores$immune_group,
  levels = c("Immune-low", "Immune-high")
)

names(y) <- immune_scores$Sample_ID

# 5. Формування labels для train/test
y_train <- y[colnames(X_train)]
y_test  <- y[colnames(X_test)]

# 6. Перевірка розмірності
cat("\nTRAIN MATRIX:\n")
print(dim(X_train))

cat("\nTEST MATRIX:\n")
print(dim(X_test))

# 7. Перевірка розподілу класів
cat("\nTRAIN LABEL DISTRIBUTION:\n")
print(table(y_train))

cat("\nTEST LABEL DISTRIBUTION:\n")
print(table(y_test))

# 8. Перевірка відповідності ID
train_ids_match <- identical(
  colnames(X_train),
  split_manifest$train_sample_ids
)

test_ids_match <- identical(
  colnames(X_test),
  split_manifest$test_sample_ids
)

# 9. Перевірка labels
train_labels_match <- identical(
  as.character(y_train),
  as.character(split_manifest$train_labels)
)

test_labels_match <- identical(
  as.character(y_test),
  as.character(split_manifest$test_labels)
)

# 10. QC
stage3_qc <- c(
  correct_train_size = ncol(X_train) == 368,
  correct_test_size = ncol(X_test) == 92,
  train_IDs_match = train_ids_match,
  test_IDs_match = test_ids_match,
  train_labels_match = train_labels_match,
  test_labels_match = test_labels_match,
  no_missing_train_labels = !anyNA(y_train),
  no_missing_test_labels = !anyNA(y_test),
  correct_feature_count = nrow(X_train) == 409556
)

cat("\nSTAGE 3 QC:\n")
print(stage3_qc)

if (!all(stage3_qc)) {
  stop("STAGE 3 QC: FAIL")
}

cat("\nSTAGE 3 TRAIN/TEST PREPARATION: PASS\n")


###############################################################################
# ЕТАП 4. CHECKPOINT AUDIT
###############################################################################

# 1. Завантаження checkpoint
inner_cv_results <- readRDS(checkpoint_file)

# 2. Перевірка кількості записів
cat("Total checkpoint slots:",
    length(inner_cv_results), "\n")

cat("Completed entries:",
    sum(!vapply(inner_cv_results, is.null, logical(1))),
    "\n")

# 3. Аудит кожного збереженого фолду
checkpoint_qc <- list()

for (fold_name in names(inner_cv_results)) {
  
  result <- inner_cv_results[[fold_name]]
  
  if (is.null(result)) next
  
  outer_id <- as.integer(
    sub("Outer_([0-9]+)_Inner_.*", "\\1", fold_name)
  )
  
  inner_id <- as.integer(
    sub(".*_Inner_([0-9]+)", "\\1", fold_name)
  )
  
  outer_data <- manifest_check$nested_cv[[outer_id]]
  
  expected_train_idx <- outer_data$outer_train_idx[
    setdiff(
      seq_along(outer_data$outer_train_idx),
      outer_data$inner_folds[[inner_id]]
    )
  ]
  
  expected_valid_idx <- outer_data$outer_train_idx[
    outer_data$inner_folds[[inner_id]]
  ]
  
  fold_qc <- c(
    correct_outer_fold = result$outer_fold == outer_id,
    correct_inner_fold = result$inner_fold == inner_id,
    
    training_IDs_match = identical(
      result$training_sample_ids,
      names(y_train)[expected_train_idx]
    ),
    
    validation_IDs_match = identical(
      result$validation_sample_ids,
      names(y_train)[expected_valid_idx]
    ),
    
    correct_CpG_count =
      length(result$selected_cpg_ids) == 5000,
    
    unique_CpG_IDs =
      !anyDuplicated(result$selected_cpg_ids),
    
    correct_lambda_count =
      length(result$lambda_values) == 50,
    
    correct_logloss_count =
      length(result$logloss) == 50,
    
    finite_logloss =
      all(is.finite(result$logloss)),
    
    correct_best_lambda =
      result$best_lambda_index == which.min(result$logloss),
    
    original_QC_pass =
      isTRUE(all(result$qc))
  )
  
  checkpoint_qc[[fold_name]] <- fold_qc
  
  cat("\n", fold_name, "\n")
  print(fold_qc)
  
  if (!all(fold_qc)) {
    stop(paste("Checkpoint QC failed:", fold_name))
  }
}

cat("\nCHECKPOINT AUDIT: PASS\n")


###############################################################################
# ЕТАП 5. FUNCTION FOR PROCESSING ONE INNER FOLD
###############################################################################

process_inner_fold <- function(
    outer_id,
    inner_id,
    X_train,
    y_train,
    manifest,
    n_features = 5000,
    n_lambda = 50,
    alpha_value = 0.5
) {
  
  cat(
    "\nProcessing Outer", outer_id,
    "/ Inner", inner_id, "\n"
  )
  
  # 1. Отримання manifest
  outer_data <- manifest$nested_cv[[outer_id]]
  
  # 2. Індекси inner training / validation
  valid_rel <- outer_data$inner_folds[[inner_id]]
  
  train_rel <- setdiff(
    seq_along(outer_data$outer_train_idx),
    valid_rel
  )
  
  train_idx <- outer_data$outer_train_idx[train_rel]
  valid_idx <- outer_data$outer_train_idx[valid_rel]
  
  # 3. Перевірка індексів
  stopifnot(
    length(train_idx) > 0,
    length(valid_idx) > 0,
    length(intersect(train_idx, valid_idx)) == 0
  )
  
  # 4. Відбір CpG лише на inner training
  cat("Calculating CpG variance...\n")
  
  cpg_variance <- DelayedMatrixStats::rowVars(
    X_train[, train_idx, drop = FALSE],
    na.rm = TRUE
  )
  
  if (
    any(!is.finite(cpg_variance)) ||
    length(cpg_variance) != nrow(X_train)
  ) {
    stop("Invalid CpG variance values.")
  }
  
  top_idx <- order(
    cpg_variance,
    decreasing = TRUE
  )[seq_len(n_features)]
  
  selected_cpg_ids <- rownames(X_train)[top_idx]
  
  stopifnot(
    length(selected_cpg_ids) == n_features,
    !anyDuplicated(selected_cpg_ids),
    !anyNA(selected_cpg_ids)
  )
  
  # 5. Формування матриць
  X_fold_train <- t(as.matrix(
    X_train[top_idx, train_idx, drop = FALSE]
  ))
  
  X_fold_valid <- t(as.matrix(
    X_train[top_idx, valid_idx, drop = FALSE]
  ))
  
  y_fold_train <- as.integer(
    y_train[train_idx] == "Immune-high"
  )
  
  y_fold_valid <- as.integer(
    y_train[valid_idx] == "Immune-high"
  )
  
  # 6. Перевірка матриць
  stopifnot(
    all(dim(X_fold_train) == c(length(train_idx), n_features)),
    all(dim(X_fold_valid) == c(length(valid_idx), n_features)),
    all(is.finite(X_fold_train)),
    all(is.finite(X_fold_valid)),
    all(y_fold_train %in% c(0L, 1L)),
    all(y_fold_valid %in% c(0L, 1L))
  )
  
  # 7. Навчання Elastic Net
  set.seed(2026 + outer_id * 100 + inner_id)
  
  fit <- glmnet::glmnet(
    x = X_fold_train,
    y = y_fold_train,
    family = "binomial",
    alpha = alpha_value,
    nlambda = n_lambda,
    lambda.min.ratio = 0.01,
    standardize = TRUE
  )
  
  lambda_values <- fit$lambda
  
  stopifnot(length(lambda_values) == n_lambda)
  
  # 8. Прогнози
  predictions <- predict(
    fit,
    newx = X_fold_valid,
    type = "response"
  )
  
  predictions <- as.matrix(predictions)
  
  # 9. Log-loss
  predictions <- pmin(
    pmax(predictions, 1e-15),
    1 - 1e-15
  )
  
  logloss <- vapply(
    seq_len(ncol(predictions)),
    function(j) {
      
      p <- predictions[, j]
      
      -mean(
        y_fold_valid * log(p) +
          (1 - y_fold_valid) * log(1 - p)
      )
    },
    numeric(1)
  )
  
  best_lambda_index <- which.min(logloss)
  
  # 10. QC
  fold_qc <- c(
    correct_training_size =
      nrow(X_fold_train) == length(train_idx),
    
    correct_validation_size =
      nrow(X_fold_valid) == length(valid_idx),
    
    correct_feature_count =
      ncol(X_fold_train) == n_features,
    
    correct_lambda_count =
      length(lambda_values) == n_lambda,
    
    finite_predictions =
      all(is.finite(predictions)),
    
    finite_logloss =
      all(is.finite(logloss)),
    
    correct_best_lambda =
      best_lambda_index == which.min(logloss),
    
    no_missing_selected_CpG =
      !anyNA(selected_cpg_ids)
  )
  
  if (!all(fold_qc)) {
    stop("Inner fold QC failed.")
  }
  
  # 11. Результат
  fold_result <- list(
    outer_fold = outer_id,
    inner_fold = inner_id,
    fold_name = paste0(
      "Outer_", outer_id,
      "_Inner_", inner_id
    ),
    
    alpha = alpha_value,
    n_features = n_features,
    
    n_training = length(train_idx),
    n_validation = length(valid_idx),
    
    training_sample_ids = names(y_train)[train_idx],
    validation_sample_ids = names(y_train)[valid_idx],
    
    selected_cpg_ids = selected_cpg_ids,
    
    lambda_values = lambda_values,
    logloss = logloss,
    
    best_lambda_index = best_lambda_index,
    best_lambda = lambda_values[best_lambda_index],
    best_logloss = logloss[best_lambda_index],
    
    qc = fold_qc
  )
  
  cat("\nBEST LAMBDA INDEX:", best_lambda_index, "\n")
  cat("BEST LAMBDA:", fold_result$best_lambda, "\n")
  cat("BEST LOG-LOSS:", fold_result$best_logloss, "\n")
  
  cat("\nINNER FOLD QC:\n")
  print(fold_qc)
  
  return(fold_result)
}

cat("\nSTAGE 5 FUNCTION CREATED: PASS\n")

###############################################################################
# ЕТАП 6. FUNCTION AND CHECKPOINT READINESS
###############################################################################

# 1. Перевірка функції
stopifnot(
  exists("process_inner_fold"),
  is.function(process_inner_fold)
)

cat("Function exists: TRUE\n")

# 2. Перевірка аргументів
print(formals(process_inner_fold))

# 3. Перевірка матриці та labels
cat("\nTraining matrix:\n")
print(dim(X_train))

cat("\nTraining labels:\n")
print(table(y_train))

# 4. Перевірка manifest
cat("\nOuter folds:", length(manifest_check$nested_cv), "\n")

cat(
  "Inner folds per outer fold:",
  paste(
    vapply(
      manifest_check$nested_cv,
      function(x) length(x$inner_folds),
      integer(1)
    ),
    collapse = ", "
  ),
  "\n"
)

# 5. Перевірка checkpoint
inner_cv_results <- readRDS(checkpoint_file)

completed_names <- names(inner_cv_results)[
  !vapply(inner_cv_results, is.null, logical(1))
]

cat("\nCompleted folds:\n")
print(completed_names)

cat("\nRemaining folds:", 20 - length(completed_names), "\n")

# 6. Загальний QC
stage6_qc <- c(
  function_available = is.function(process_inner_fold),
  correct_training_dimensions =
    all(dim(X_train) == c(409556, 368)),
  correct_training_labels =
    length(y_train) == 368 && !anyNA(y_train),
  five_outer_folds =
    length(manifest_check$nested_cv) == 5,
  four_inner_folds_each =
    all(vapply(
      manifest_check$nested_cv,
      function(x) length(x$inner_folds) == 4,
      logical(1)
    )),
  checkpoint_has_20_slots =
    length(inner_cv_results) == 20,
  two_completed_folds =
    length(completed_names) == 2
)

cat("\nSTAGE 6 QC:\n")
print(stage6_qc)

if (!all(stage6_qc)) {
  stop("STAGE 6 QC: FAIL")
}

cat("\nSTAGE 6 READINESS: PASS\n")

###############################################################################
# ЕТАП 7. TEST ONE INNER FOLD — OUTER 1 / INNER 2
###############################################################################

test_fold <- process_inner_fold(
  outer_id = 1,
  inner_id = 2,
  X_train = X_train,
  y_train = y_train,
  manifest = manifest_check,
  n_features = 5000,
  n_lambda = 50,
  alpha_value = 0.5
)

cat("\nTEST FOLD RESULT:\n")

cat("Fold:", test_fold$fold_name, "\n")
cat("Training samples:", test_fold$n_training, "\n")
cat("Validation samples:", test_fold$n_validation, "\n")
cat("Selected CpG:", length(test_fold$selected_cpg_ids), "\n")
cat("Lambda values:", length(test_fold$lambda_values), "\n")
cat("Best lambda index:", test_fold$best_lambda_index, "\n")
cat("Best lambda:", test_fold$best_lambda, "\n")
cat("Best log-loss:", test_fold$best_logloss, "\n")

cat("\nFOLD QC:\n")
print(test_fold$qc)

stopifnot(
  test_fold$outer_fold == 1,
  test_fold$inner_fold == 2,
  test_fold$n_features == 5000,
  test_fold$n_training == 220,
  test_fold$n_validation == 74,
  length(test_fold$selected_cpg_ids) == 5000,
  length(test_fold$lambda_values) == 50,
  length(test_fold$logloss) == 50,
  all(is.finite(test_fold$logloss)),
  isTRUE(all(test_fold$qc))
)

cat("\nSTAGE 7 TEST FOLD: PASS\n")

###############################################################################
# ЕТАП 8. SAVE TEST FOLD TO CHECKPOINT
###############################################################################

# 1. Повторне завантаження актуального checkpoint
inner_cv_results <- readRDS(checkpoint_file)

# 2. Назва фолду
fold_name <- "Outer_1_Inner_2"

# 3. Перевірка, що слот ще не зайнятий
if (!is.null(inner_cv_results[[fold_name]])) {
  stop("This fold already exists. Checkpoint was not modified.")
}

# 4. Додавання тестового результату
inner_cv_results[[fold_name]] <- test_fold

# 5. Збереження checkpoint
saveRDS(
  inner_cv_results,
  checkpoint_file
)

# 6. Повторне завантаження
checkpoint_verify <- readRDS(checkpoint_file)

saved_fold <- checkpoint_verify[[fold_name]]

# 7. Перевірка збереження
save_qc <- c(
  fold_exists = !is.null(saved_fold),
  
  correct_fold_name =
    saved_fold$fold_name == fold_name,
  
  correct_training_size =
    saved_fold$n_training == 220,
  
  correct_validation_size =
    saved_fold$n_validation == 74,
  
  correct_CpG_count =
    length(saved_fold$selected_cpg_ids) == 5000,
  
  correct_lambda_count =
    length(saved_fold$lambda_values) == 50,
  
  logloss_preserved =
    identical(saved_fold$logloss, test_fold$logloss),
  
  lambda_preserved =
    identical(saved_fold$lambda_values, test_fold$lambda_values),
  
  QC_pass =
    isTRUE(all(saved_fold$qc))
)

cat("\nCHECKPOINT SAVE QC:\n")
print(save_qc)

if (!all(save_qc)) {
  stop("Checkpoint save QC: FAIL")
}

cat(
  "\nCompleted folds:",
  sum(!vapply(checkpoint_verify, is.null, logical(1))),
  "/20\n"
)

cat("\nSTAGE 8 CHECKPOINT SAVE: PASS\n")

###############################################################################
# ЕТАП 9. CHECKPOINT READINESS FOR BATCH PROCESSING
###############################################################################

inner_cv_results <- readRDS(checkpoint_file)

completed_names <- names(inner_cv_results)[
  !vapply(inner_cv_results, is.null, logical(1))
]

expected_completed <- c(
  "Outer_1_Inner_1",
  "Outer_1_Inner_2",
  "Outer_5_Inner_4"
)

cat("\nCOMPLETED FOLDS:\n")
print(completed_names)

# 1. Перевірка кількості
stopifnot(length(completed_names) == 3)

# 2. Перевірка очікуваних назв
stopifnot(setequal(completed_names, expected_completed))

# 3. Перевірка кожного завершеного результату
for (fold_name in completed_names) {
  
  result <- inner_cv_results[[fold_name]]
  
  stopifnot(
    length(result$selected_cpg_ids) == 5000,
    length(result$lambda_values) == 50,
    length(result$logloss) == 50,
    all(is.finite(result$logloss)),
    result$best_lambda_index == which.min(result$logloss),
    isTRUE(all(result$qc))
  )
}

# 4. Перевірка незавершених слотів
remaining_names <- names(inner_cv_results)[
  vapply(inner_cv_results, is.null, logical(1))
]

cat("\nRemaining folds:", length(remaining_names), "\n")

print(remaining_names)

stopifnot(length(remaining_names) == 17)

# 5. Підсумок
cat("\nBATCH READINESS QC: PASS\n")
cat("Completed:", length(completed_names), "/20\n")
cat("Remaining:", length(remaining_names), "/20\n")




###############################################################################
# ЕТАП 10. RUN REMAINING INNER FOLDS WITH CHECKPOINTING
###############################################################################

# 1. Завантажити актуальний checkpoint
inner_cv_results <- readRDS(checkpoint_file)

# 2. Перелік усіх фолдів
all_fold_names <- names(inner_cv_results)

# 3. Перелік завершених фолдів
completed_names <- all_fold_names[
  !vapply(inner_cv_results, is.null, logical(1))
]

# 4. Перелік незавершених фолдів
remaining_names <- all_fold_names[
  vapply(inner_cv_results, is.null, logical(1))
]

cat("Already completed:", length(completed_names), "\n")
cat("Remaining:", length(remaining_names), "\n")

stopifnot(
  length(all_fold_names) == 20,
  length(completed_names) == 3,
  length(remaining_names) == 17
)

# 5. Пакетний цикл
for (i in seq_along(remaining_names)) {
  
  fold_name <- remaining_names[i]
  
  outer_id <- as.integer(
    sub("Outer_([0-9]+)_Inner_.*", "\\1", fold_name)
  )
  
  inner_id <- as.integer(
    sub(".*_Inner_([0-9]+)", "\\1", fold_name)
  )
  
  cat("\n")
  cat("============================================\n")
  cat("FOLD", i, "OF", length(remaining_names), "\n")
  cat("Processing:", fold_name, "\n")
  cat("============================================\n")
  
  # Перевірка, що слот досі порожній
  if (!is.null(inner_cv_results[[fold_name]])) {
    stop(paste("Unexpected existing result:", fold_name))
  }
  
  # 6. Обробка одного фолду
  candidate <- process_inner_fold(
    outer_id = outer_id,
    inner_id = inner_id,
    X_train = X_train,
    y_train = y_train,
    manifest = manifest_check,
    n_features = 5000,
    n_lambda = 50,
    alpha_value = 0.5
  )
  
  # 7. Очікувані індекси з manifest
  outer_data <- manifest_check$nested_cv[[outer_id]]
  
  expected_train_idx <- outer_data$outer_train_idx[
    setdiff(
      seq_along(outer_data$outer_train_idx),
      outer_data$inner_folds[[inner_id]]
    )
  ]
  
  expected_valid_idx <- outer_data$outer_train_idx[
    outer_data$inner_folds[[inner_id]]
  ]
  
  # 8. Розширений QC
  batch_qc <- c(
    correct_fold_name =
      identical(candidate$fold_name, fold_name),
    
    training_IDs_match =
      identical(
        candidate$training_sample_ids,
        names(y_train)[expected_train_idx]
      ),
    
    validation_IDs_match =
      identical(
        candidate$validation_sample_ids,
        names(y_train)[expected_valid_idx]
      ),
    
    correct_CpG_count =
      length(candidate$selected_cpg_ids) == 5000,
    
    unique_CpG_IDs =
      !anyDuplicated(candidate$selected_cpg_ids),
    
    correct_lambda_count =
      length(candidate$lambda_values) == 50,
    
    correct_logloss_count =
      length(candidate$logloss) == 50,
    
    finite_logloss =
      all(is.finite(candidate$logloss)),
    
    correct_best_lambda =
      candidate$best_lambda_index ==
      which.min(candidate$logloss),
    
    original_QC_pass =
      isTRUE(all(candidate$qc))
  )
  
  cat("\nBATCH QC:\n")
  print(batch_qc)
  
  if (!all(batch_qc)) {
    stop(paste("Batch QC failed:", fold_name))
  }
  
  # 9. Додавання результату лише після успішного QC
  inner_cv_results[[fold_name]] <- candidate
  
  # 10. Збереження checkpoint після кожного фолду
  saveRDS(
    inner_cv_results,
    checkpoint_file
  )
  
  # 11. Перевірка повторного завантаження
  checkpoint_verify <- readRDS(checkpoint_file)
  
  saved_result <- checkpoint_verify[[fold_name]]
  
  stopifnot(
    !is.null(saved_result),
    identical(saved_result$logloss, candidate$logloss),
    identical(
      saved_result$selected_cpg_ids,
      candidate$selected_cpg_ids
    ),
    isTRUE(all(saved_result$qc))
  )
  
  # 12. Оновлення таблиці прогресу
  completed_now <- names(checkpoint_verify)[
    !vapply(checkpoint_verify, is.null, logical(1))
  ]
  
  progress_df <- data.frame(
    Fold = names(checkpoint_verify),
    Completed = !vapply(
      checkpoint_verify,
      is.null,
      logical(1)
    ),
    Best_Lambda_Index = vapply(
      checkpoint_verify,
      function(x) {
        if (is.null(x)) NA_integer_
        else as.integer(x$best_lambda_index)
      },
      integer(1)
    ),
    Best_Logloss = vapply(
      checkpoint_verify,
      function(x) {
        if (is.null(x)) NA_real_
        else as.numeric(x$best_logloss)
      },
      numeric(1)
    )
  )
  
  write.csv(
    progress_df,
    progress_file,
    row.names = FALSE
  )
  
  cat("\nCHECKPOINT SAVED:", fold_name, "\n")
  cat("Total completed:", length(completed_now), "/20\n")
  cat("Best log-loss:", saved_result$best_logloss, "\n")
}

cat("\n============================================\n")
cat("ALL REMAINING FOLDS PROCESSED\n")
cat("============================================\n")

cat(
  "Completed folds:",
  sum(!vapply(inner_cv_results, is.null, logical(1))),
  "/20\n"
)

###############################################################################
# ЕТАП 11. FULL AUDIT OF ALL 20 INNER FOLDS
###############################################################################

# 1. Завантаження checkpoint
inner_cv_results <- readRDS(checkpoint_file)

# 2. Перевірка кількості фолдів
stopifnot(
  length(inner_cv_results) == 20,
  all(!vapply(inner_cv_results, is.null, logical(1)))
)

# 3. Очікувані назви фолдів
expected_names <- unlist(
  lapply(1:5, function(o) {
    paste0("Outer_", o, "_Inner_", 1:4)
  })
)

stopifnot(setequal(names(inner_cv_results), expected_names))

# 4. Створення таблиці аудиту
audit_results <- vector("list", 20)

for (i in seq_along(expected_names)) {
  
  fold_name <- expected_names[i]
  result <- inner_cv_results[[fold_name]]
  
  outer_id <- as.integer(
    sub("Outer_([0-9]+)_Inner_.*", "\\1", fold_name)
  )
  
  inner_id <- as.integer(
    sub(".*_Inner_([0-9]+)", "\\1", fold_name)
  )
  
  outer_data <- manifest_check$nested_cv[[outer_id]]
  
  expected_valid_idx <- outer_data$outer_train_idx[
    outer_data$inner_folds[[inner_id]]
  ]
  
  expected_train_idx <- outer_data$outer_train_idx[
    setdiff(
      seq_along(outer_data$outer_train_idx),
      outer_data$inner_folds[[inner_id]]
    )
  ]
  
  # 5. Перевірки
  qc <- c(
    correct_fold_name =
      identical(result$fold_name, fold_name),
    
    correct_outer_fold =
      identical(as.integer(result$outer_fold), outer_id),
    
    correct_inner_fold =
      identical(as.integer(result$inner_fold), inner_id),
    
    training_IDs_match =
      identical(
        result$training_sample_ids,
        names(y_train)[expected_train_idx]
      ),
    
    validation_IDs_match =
      identical(
        result$validation_sample_ids,
        names(y_train)[expected_valid_idx]
      ),
    
    correct_training_size =
      length(result$training_sample_ids) ==
      length(expected_train_idx),
    
    correct_validation_size =
      length(result$validation_sample_ids) ==
      length(expected_valid_idx),
    
    correct_CpG_count =
      length(result$selected_cpg_ids) == 5000,
    
    unique_CpG_IDs =
      !anyDuplicated(result$selected_cpg_ids),
    
    CpG_IDs_exist_in_matrix =
      all(result$selected_cpg_ids %in% rownames(X_train)),
    
    correct_lambda_count =
      length(result$lambda_values) == 50,
    
    correct_logloss_count =
      length(result$logloss) == 50,
    
    finite_lambda =
      all(is.finite(result$lambda_values)),
    
    finite_logloss =
      all(is.finite(result$logloss)),
    
    correct_best_lambda_index =
      result$best_lambda_index ==
      which.min(result$logloss),
    
    correct_best_lambda_value =
      isTRUE(all.equal(
        as.numeric(result$best_lambda),
        as.numeric(result$lambda_values[
          result$best_lambda_index
        ])
      )),
    
    correct_best_logloss =
      isTRUE(all.equal(
        as.numeric(result$best_logloss),
        min(result$logloss)
      )),
    
    original_fold_QC =
      isTRUE(all(result$qc))
  )
  
  audit_results[[i]] <- data.frame(
    Fold = fold_name,
    Training_n = length(result$training_sample_ids),
    Validation_n = length(result$validation_sample_ids),
    CpG_n = length(result$selected_cpg_ids),
    Lambda_n = length(result$lambda_values),
    Best_lambda_index = result$best_lambda_index,
    Best_logloss = result$best_logloss,
    All_checks_pass = all(qc),
    Failed_checks = sum(!qc)
  )
  
  if (!all(qc)) {
    cat("\nFAILED:", fold_name, "\n")
    print(qc)
  }
}

# 6. Підсумкова таблиця
audit_table <- do.call(rbind, audit_results)

cat("\nFULL INNER CV AUDIT:\n")
print(audit_table, row.names = FALSE)

# 7. Загальний результат
cat("\nTOTAL FOLDS:", nrow(audit_table), "\n")
cat("PASSED:", sum(audit_table$All_checks_pass), "\n")
cat("FAILED:", sum(!audit_table$All_checks_pass), "\n")

stopifnot(
  nrow(audit_table) == 20,
  all(audit_table$All_checks_pass),
  all(audit_table$Failed_checks == 0)
)

# 8. Збереження аудиту
write.csv(
  audit_table,
  "results/tables/TCGA_SKCM_inner_CV_full_audit.csv",
  row.names = FALSE
)

cat("\nFULL INNER CV AUDIT: PASS\n")

###############################################################################
# ЕТАП 12. IDENTIFY FAILED INNER FOLD QC
###############################################################################

# 1. Назва фолду, який не пройшов аудит
failed_folds <- audit_table$Fold[
  !audit_table$All_checks_pass
]

cat("FAILED FOLD(S):\n")
print(failed_folds)

stopifnot(length(failed_folds) == 1)

# 2. Детальна перевірка проблемного фолду
fold_name <- failed_folds[1]

result <- inner_cv_results[[fold_name]]

outer_id <- as.integer(
  sub("Outer_([0-9]+)_Inner_.*", "\\1", fold_name)
)

inner_id <- as.integer(
  sub(".*_Inner_([0-9]+)", "\\1", fold_name)
)

outer_data <- manifest_check$nested_cv[[outer_id]]

expected_valid_idx <- outer_data$outer_train_idx[
  outer_data$inner_folds[[inner_id]]
]

expected_train_idx <- outer_data$outer_train_idx[
  setdiff(
    seq_along(outer_data$outer_train_idx),
    outer_data$inner_folds[[inner_id]]
  )
]

# 3. Повторне виконання всіх QC окремо
failed_qc <- c(
  correct_fold_name =
    identical(result$fold_name, fold_name),
  
  correct_outer_fold =
    identical(as.integer(result$outer_fold), outer_id),
  
  correct_inner_fold =
    identical(as.integer(result$inner_fold), inner_id),
  
  training_IDs_match =
    identical(
      result$training_sample_ids,
      names(y_train)[expected_train_idx]
    ),
  
  validation_IDs_match =
    identical(
      result$validation_sample_ids,
      names(y_train)[expected_valid_idx]
    ),
  
  correct_training_size =
    length(result$training_sample_ids) ==
    length(expected_train_idx),
  
  correct_validation_size =
    length(result$validation_sample_ids) ==
    length(expected_valid_idx),
  
  correct_CpG_count =
    length(result$selected_cpg_ids) == 5000,
  
  unique_CpG_IDs =
    !anyDuplicated(result$selected_cpg_ids),
  
  CpG_IDs_exist_in_matrix =
    all(result$selected_cpg_ids %in% rownames(X_train)),
  
  correct_lambda_count =
    length(result$lambda_values) == 50,
  
  correct_logloss_count =
    length(result$logloss) == 50,
  
  finite_lambda =
    all(is.finite(result$lambda_values)),
  
  finite_logloss =
    all(is.finite(result$logloss)),
  
  correct_best_lambda_index =
    result$best_lambda_index ==
    which.min(result$logloss),
  
  correct_best_lambda_value =
    isTRUE(all.equal(
      as.numeric(result$best_lambda),
      as.numeric(result$lambda_values[
        result$best_lambda_index
      ])
    )),
  
  correct_best_logloss =
    isTRUE(all.equal(
      as.numeric(result$best_logloss),
      min(result$logloss)
    )),
  
  original_fold_QC =
    isTRUE(all(result$qc))
)

cat("\nDETAILED QC FOR:", fold_name, "\n")
print(failed_qc)

cat("\nFAILED CHECKS:\n")
print(names(failed_qc)[!failed_qc])

cat("\nRESULT STRUCTURE:\n")
str(result[c(
  "outer_fold",
  "inner_fold",
  "n_training",
  "n_validation",
  "best_lambda_index",
  "best_lambda",
  "best_logloss"
)])


###############################################################################
# ЕТАП 13. REPAIR LEGACY PILOT FOLD NAME
###############################################################################

# 1. Завантаження checkpoint
inner_cv_results <- readRDS(checkpoint_file)

# 2. Отримання першого фолду
pilot <- inner_cv_results[["Outer_1_Inner_1"]]

# 3. Перевірка поточної назви
cat("Current fold_name:\n")
print(pilot$fold_name)

# 4. Встановлення правильної назви
pilot$fold_name <- "Outer_1_Inner_1"

# 5. Повернення в checkpoint
inner_cv_results[["Outer_1_Inner_1"]] <- pilot

# 6. Збереження
saveRDS(inner_cv_results, checkpoint_file)

# 7. Повторне завантаження
checkpoint_verify <- readRDS(checkpoint_file)

# 8. Перевірка
stopifnot(
  identical(
    checkpoint_verify[["Outer_1_Inner_1"]]$fold_name,
    "Outer_1_Inner_1"
  )
)

cat("\nFOLD NAME REPAIR: PASS\n")



###############################################################################
# ЕТАП 11. FULL AUDIT OF ALL 20 INNER FOLDS
###############################################################################

# 1. Завантаження checkpoint
inner_cv_results <- readRDS(checkpoint_file)

# 2. Перевірка кількості фолдів
stopifnot(
  length(inner_cv_results) == 20,
  all(!vapply(inner_cv_results, is.null, logical(1)))
)

# 3. Очікувані назви фолдів
expected_names <- unlist(
  lapply(1:5, function(o) {
    paste0("Outer_", o, "_Inner_", 1:4)
  })
)

stopifnot(setequal(names(inner_cv_results), expected_names))

# 4. Створення таблиці аудиту
audit_results <- vector("list", 20)

for (i in seq_along(expected_names)) {
  
  fold_name <- expected_names[i]
  result <- inner_cv_results[[fold_name]]
  
  outer_id <- as.integer(
    sub("Outer_([0-9]+)_Inner_.*", "\\1", fold_name)
  )
  
  inner_id <- as.integer(
    sub(".*_Inner_([0-9]+)", "\\1", fold_name)
  )
  
  outer_data <- manifest_check$nested_cv[[outer_id]]
  
  expected_valid_idx <- outer_data$outer_train_idx[
    outer_data$inner_folds[[inner_id]]
  ]
  
  expected_train_idx <- outer_data$outer_train_idx[
    setdiff(
      seq_along(outer_data$outer_train_idx),
      outer_data$inner_folds[[inner_id]]
    )
  ]
  
  # 5. Перевірки
  qc <- c(
    correct_fold_name =
      identical(result$fold_name, fold_name),
    
    correct_outer_fold =
      identical(as.integer(result$outer_fold), outer_id),
    
    correct_inner_fold =
      identical(as.integer(result$inner_fold), inner_id),
    
    training_IDs_match =
      identical(
        result$training_sample_ids,
        names(y_train)[expected_train_idx]
      ),
    
    validation_IDs_match =
      identical(
        result$validation_sample_ids,
        names(y_train)[expected_valid_idx]
      ),
    
    correct_training_size =
      length(result$training_sample_ids) ==
      length(expected_train_idx),
    
    correct_validation_size =
      length(result$validation_sample_ids) ==
      length(expected_valid_idx),
    
    correct_CpG_count =
      length(result$selected_cpg_ids) == 5000,
    
    unique_CpG_IDs =
      !anyDuplicated(result$selected_cpg_ids),
    
    CpG_IDs_exist_in_matrix =
      all(result$selected_cpg_ids %in% rownames(X_train)),
    
    correct_lambda_count =
      length(result$lambda_values) == 50,
    
    correct_logloss_count =
      length(result$logloss) == 50,
    
    finite_lambda =
      all(is.finite(result$lambda_values)),
    
    finite_logloss =
      all(is.finite(result$logloss)),
    
    correct_best_lambda_index =
      result$best_lambda_index ==
      which.min(result$logloss),
    
    correct_best_lambda_value =
      isTRUE(all.equal(
        as.numeric(result$best_lambda),
        as.numeric(result$lambda_values[
          result$best_lambda_index
        ])
      )),
    
    correct_best_logloss =
      isTRUE(all.equal(
        as.numeric(result$best_logloss),
        min(result$logloss)
      )),
    
    original_fold_QC =
      isTRUE(all(result$qc))
  )
  
  audit_results[[i]] <- data.frame(
    Fold = fold_name,
    Training_n = length(result$training_sample_ids),
    Validation_n = length(result$validation_sample_ids),
    CpG_n = length(result$selected_cpg_ids),
    Lambda_n = length(result$lambda_values),
    Best_lambda_index = result$best_lambda_index,
    Best_logloss = result$best_logloss,
    All_checks_pass = all(qc),
    Failed_checks = sum(!qc)
  )
  
  if (!all(qc)) {
    cat("\nFAILED:", fold_name, "\n")
    print(qc)
  }
}

# 6. Підсумкова таблиця
audit_table <- do.call(rbind, audit_results)

cat("\nFULL INNER CV AUDIT:\n")
print(audit_table, row.names = FALSE)

# 7. Загальний результат
cat("\nTOTAL FOLDS:", nrow(audit_table), "\n")
cat("PASSED:", sum(audit_table$All_checks_pass), "\n")
cat("FAILED:", sum(!audit_table$All_checks_pass), "\n")

stopifnot(
  nrow(audit_table) == 20,
  all(audit_table$All_checks_pass),
  all(audit_table$Failed_checks == 0)
)

# 8. Збереження аудиту
write.csv(
  audit_table,
  "results/tables/TCGA_SKCM_inner_CV_full_audit.csv",
  row.names = FALSE
)

cat("\nFULL INNER CV AUDIT: PASS\n")


###############################################################################
# ЕТАП 14. AGGREGATE INNER CV AND SELECT OPTIMAL LAMBDA INDEX
###############################################################################

# 1. Завантаження перевірених результатів
inner_cv_results <- readRDS(checkpoint_file)

stopifnot(
  length(inner_cv_results) == 20,
  all(!vapply(inner_cv_results, is.null, logical(1)))
)

# 2. Підготовка результатів для п'яти outer folds
outer_tuning_results <- vector("list", 5)

for (outer_id in 1:5) {
  
  cat("\n====================================\n")
  cat("OUTER FOLD:", outer_id, "\n")
  cat("====================================\n")
  
  # Назви чотирьох inner folds
  inner_names <- paste0(
    "Outer_", outer_id, "_Inner_", 1:4
  )
  
  # Отримання результатів
  fold_results <- inner_cv_results[inner_names]
  
  stopifnot(
    length(fold_results) == 4,
    all(!vapply(fold_results, is.null, logical(1)))
  )
  
  # Перевірка довжини lambda path
  stopifnot(
    all(vapply(
      fold_results,
      function(x) length(x$logloss) == 50,
      logical(1)
    ))
  )
  
  # 3. Матриця log-loss:
  # рядки = inner folds
  # стовпці = позиції lambda
  logloss_matrix <- do.call(
    rbind,
    lapply(fold_results, function(x) x$logloss)
  )
  
  stopifnot(
    all(dim(logloss_matrix) == c(4, 50)),
    all(is.finite(logloss_matrix))
  )
  
  # 4. Середній log-loss для кожної позиції
  mean_logloss <- colMeans(logloss_matrix)
  
  # Стандартне відхилення між inner folds
  sd_logloss <- apply(
    logloss_matrix,
    2,
    sd
  )
  
  # 5. Оптимальна позиція lambda
  best_index <- which.min(mean_logloss)
  
  # Lambda values у кожному inner fold
  selected_lambda_per_fold <- vapply(
    fold_results,
    function(x) x$lambda_values[best_index],
    numeric(1)
  )
  
  # 6. Збереження результатів outer fold
  outer_tuning_results[[outer_id]] <- list(
    outer_fold = outer_id,
    alpha = 0.5,
    n_features = 5000,
    n_inner_folds = 4,
    n_lambda = 50,
    mean_logloss = mean_logloss,
    sd_logloss = sd_logloss,
    logloss_matrix = logloss_matrix,
    best_lambda_index = best_index,
    best_mean_logloss = mean_logloss[best_index],
    selected_lambda_per_inner_fold =
      selected_lambda_per_fold
  )
  
  cat("Best lambda index:", best_index, "\n")
  cat("Mean log-loss:", mean_logloss[best_index], "\n")
  cat("SD log-loss:", sd_logloss[best_index], "\n")
  cat("Lambda values across inner folds:\n")
  print(selected_lambda_per_fold)
}

# 7. Загальна таблиця п'яти outer folds
outer_tuning_table <- do.call(
  rbind,
  lapply(outer_tuning_results, function(x) {
    
    data.frame(
      Outer_Fold = x$outer_fold,
      Alpha = x$alpha,
      Selected_CpG = x$n_features,
      Inner_Folds = x$n_inner_folds,
      Lambda_Path_Length = x$n_lambda,
      Best_Lambda_Index = x$best_lambda_index,
      Best_Mean_Logloss = x$best_mean_logloss,
      SD_Logloss = x$sd_logloss[
        x$best_lambda_index
      ]
    )
  })
)

cat("\nOUTER FOLD TUNING SUMMARY:\n")
print(outer_tuning_table, row.names = FALSE)

# 8. Фінальний QC
aggregation_qc <- c(
  five_outer_folds =
    length(outer_tuning_results) == 5,
  
  four_inner_folds_each =
    all(vapply(
      outer_tuning_results,
      function(x) x$n_inner_folds == 4,
      logical(1)
    )),
  
  fifty_lambda_positions =
    all(vapply(
      outer_tuning_results,
      function(x) length(x$mean_logloss) == 50,
      logical(1)
    )),
  
  finite_mean_logloss =
    all(vapply(
      outer_tuning_results,
      function(x) all(is.finite(x$mean_logloss)),
      logical(1)
    )),
  
  valid_best_indices =
    all(vapply(
      outer_tuning_results,
      function(x) x$best_lambda_index %in% 1:50,
      logical(1)
    )),
  
  correct_minimum_selection =
    all(vapply(
      outer_tuning_results,
      function(x) {
        x$best_lambda_index ==
          which.min(x$mean_logloss)
      },
      logical(1)
    ))
)

cat("\nAGGREGATION QC:\n")
print(aggregation_qc)

stopifnot(all(aggregation_qc))

# 9. Збереження результатів
saveRDS(
  outer_tuning_results,
  "results/models/TCGA_SKCM_outer_fold_tuning_results.rds"
)

write.csv(
  outer_tuning_table,
  "results/tables/TCGA_SKCM_outer_fold_tuning_summary.csv",
  row.names = FALSE
)

cat("\nINNER CV AGGREGATION: PASS\n")


###############################################################################
# ЕТАП 15. OUTER LOOP — OUTER FOLD 1
###############################################################################

outer_id <- 1

cat("\n====================================\n")
cat("OUTER FOLD:", outer_id, "\n")
cat("====================================\n")

# 1. Отримання індексів із manifest
outer_data <- manifest_check$nested_cv[[outer_id]]

outer_train_idx <- outer_data$outer_train_idx
outer_valid_idx <- outer_data$outer_valid_idx

stopifnot(
  length(intersect(outer_train_idx, outer_valid_idx)) == 0,
  length(outer_train_idx) + length(outer_valid_idx) == 368
)

# 2. Вибір оптимального lambda index із inner CV
best_index <- outer_tuning_results[[outer_id]]$best_lambda_index

cat("Selected lambda index:", best_index, "\n")

stopifnot(best_index %in% 1:50)

# 3. Відбір 5000 CpG виключно на outer-training samples
cat("\nCalculating CpG variance...\n")

outer_variance <- DelayedMatrixStats::rowVars(
  X_train[, outer_train_idx, drop = FALSE],
  na.rm = TRUE
)

stopifnot(
  all(is.finite(outer_variance)),
  length(outer_variance) == nrow(X_train)
)

outer_top_idx <- order(
  outer_variance,
  decreasing = TRUE
)[seq_len(5000)]

outer_selected_cpg <- rownames(X_train)[outer_top_idx]

stopifnot(
  length(outer_selected_cpg) == 5000,
  !anyDuplicated(outer_selected_cpg)
)

# 4. Формування матриць
X_outer_train <- t(as.matrix(
  X_train[outer_top_idx, outer_train_idx, drop = FALSE]
))

X_outer_valid <- t(as.matrix(
  X_train[outer_top_idx, outer_valid_idx, drop = FALSE]
))

y_outer_train <- as.integer(
  y_train[outer_train_idx] == "Immune-high"
)

y_outer_valid <- as.integer(
  y_train[outer_valid_idx] == "Immune-high"
)

stopifnot(
  all(dim(X_outer_train) == c(length(outer_train_idx), 5000)),
  all(dim(X_outer_valid) == c(length(outer_valid_idx), 5000)),
  all(is.finite(X_outer_train)),
  all(is.finite(X_outer_valid))
)

# 5. Навчання Elastic Net на всій outer-training вибірці
set.seed(2026 + outer_id)

outer_fit <- glmnet::glmnet(
  x = X_outer_train,
  y = y_outer_train,
  family = "binomial",
  alpha = 0.5,
  nlambda = 50,
  lambda.min.ratio = 0.01,
  standardize = TRUE
)

stopifnot(length(outer_fit$lambda) == 50)

# 6. Lambda для вибраної позиції шляху
selected_lambda <- outer_fit$lambda[best_index]

cat("Refitted lambda:", selected_lambda, "\n")

# 7. Прогноз на outer-validation
outer_prob <- as.numeric(
  predict(
    outer_fit,
    newx = X_outer_valid,
    s = selected_lambda,
    type = "response"
  )
)

stopifnot(
  length(outer_prob) == length(outer_valid_idx),
  all(is.finite(outer_prob)),
  all(outer_prob >= 0 & outer_prob <= 1)
)

# 8. Класифікація з порогом 0.5
outer_pred <- as.integer(outer_prob >= 0.5)

# 9. Confusion matrix
outer_confusion <- table(
  Actual = factor(
    y_outer_valid,
    levels = c(0, 1),
    labels = c("Immune-low", "Immune-high")
  ),
  Predicted = factor(
    outer_pred,
    levels = c(0, 1),
    labels = c("Immune-low", "Immune-high")
  )
)

cat("\nCONFUSION MATRIX:\n")
print(outer_confusion)

# 10. Метрики
TN <- outer_confusion[1, 1]
FP <- outer_confusion[1, 2]
FN <- outer_confusion[2, 1]
TP <- outer_confusion[2, 2]

accuracy <- (TP + TN) / sum(outer_confusion)

sensitivity <- if ((TP + FN) > 0) {
  TP / (TP + FN)
} else NA_real_

specificity <- if ((TN + FP) > 0) {
  TN / (TN + FP)
} else NA_real_

precision <- if ((TP + FP) > 0) {
  TP / (TP + FP)
} else NA_real_

f1_score <- if (is.finite(precision) &&
                is.finite(sensitivity) &&
                (precision + sensitivity) > 0) {
  2 * precision * sensitivity /
    (precision + sensitivity)
} else NA_real_

balanced_accuracy <- mean(
  c(sensitivity, specificity),
  na.rm = TRUE
)

# 11. AUC
outer_roc <- pROC::roc(
  response = factor(y_outer_valid, levels = c(0, 1)),
  predictor = outer_prob,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

outer_auc <- as.numeric(pROC::auc(outer_roc))

# 12. Таблиця метрик
outer_metrics <- data.frame(
  Outer_Fold = outer_id,
  Training_n = length(outer_train_idx),
  Validation_n = length(outer_valid_idx),
  Selected_CpG = 5000,
  Best_Lambda_Index = best_index,
  Lambda = selected_lambda,
  Accuracy = accuracy,
  Sensitivity = sensitivity,
  Specificity = specificity,
  Precision = precision,
  F1 = f1_score,
  Balanced_Accuracy = balanced_accuracy,
  AUC = outer_auc
)

cat("\nOUTER FOLD 1 METRICS:\n")
print(outer_metrics, row.names = FALSE)

# 13. Збереження прогнозів
outer_predictions <- data.frame(
  Sample_ID = names(y_train)[outer_valid_idx],
  Actual = ifelse(
    y_outer_valid == 1,
    "Immune-high",
    "Immune-low"
  ),
  Predicted = ifelse(
    outer_pred == 1,
    "Immune-high",
    "Immune-low"
  ),
  Probability_Immune_high = outer_prob
)

# 14. QC
outer_qc <- c(
  correct_training_size =
    nrow(X_outer_train) == length(outer_train_idx),
  
  correct_validation_size =
    nrow(X_outer_valid) == length(outer_valid_idx),
  
  correct_feature_count =
    ncol(X_outer_train) == 5000,
  
  correct_prediction_count =
    length(outer_prob) == length(outer_valid_idx),
  
  finite_predictions =
    all(is.finite(outer_prob)),
  
  valid_probability_range =
    all(outer_prob >= 0 & outer_prob <= 1),
  
  finite_metrics =
    all(is.finite(as.numeric(
      outer_metrics[, c(
        "Accuracy", "Sensitivity", "Specificity",
        "Precision", "F1", "Balanced_Accuracy", "AUC"
      )]
    )))
)

cat("\nOUTER FOLD QC:\n")
print(outer_qc)

stopifnot(all(outer_qc))

# 15. Створення checkpoint для outer evaluation
outer_output_dir <- "results/models/nested_CV_outer"

if (!dir.exists(outer_output_dir)) {
  dir.create(
    outer_output_dir,
    recursive = TRUE
  )
}

outer_evaluation_results <- vector("list", 5)
names(outer_evaluation_results) <- paste0(
  "Outer_", 1:5
)

outer_evaluation_results[[outer_id]] <- list(
  outer_fold = outer_id,
  selected_cpg_ids = outer_selected_cpg,
  best_lambda_index = best_index,
  lambda = selected_lambda,
  metrics = outer_metrics,
  confusion_matrix = outer_confusion,
  predictions = outer_predictions,
  qc = outer_qc
)

outer_checkpoint_file <- file.path(
  outer_output_dir,
  "TCGA_SKCM_outer_CV_partial_results.rds"
)

saveRDS(
  outer_evaluation_results,
  outer_checkpoint_file
)

cat("\nOUTER FOLD 1: PASS\n")
cat("Checkpoint saved:", outer_checkpoint_file, "\n")


###############################################################################
# ЕТАП 16. VERIFY OUTER FOLD 1 CHECKPOINT
###############################################################################

outer_checkpoint_check <- readRDS(outer_checkpoint_file)

saved_outer1 <- outer_checkpoint_check[["Outer_1"]]

outer1_save_qc <- c(
  fold_exists =
    !is.null(saved_outer1),
  
  correct_fold_id =
    identical(saved_outer1$outer_fold, 1),
  
  correct_CpG_count =
    length(saved_outer1$selected_cpg_ids) == 5000,
  
  correct_prediction_count =
    nrow(saved_outer1$predictions) == 74,
  
  correct_metric_count =
    nrow(saved_outer1$metrics) == 1,
  
  predictions_finite =
    all(is.finite(
      saved_outer1$predictions$Probability_Immune_high
    )),
  
  metrics_preserved =
    isTRUE(all.equal(
      saved_outer1$metrics$AUC,
      outer_metrics$AUC
    )),
  
  QC_pass =
    isTRUE(all(saved_outer1$qc))
)

cat("\nOUTER 1 CHECKPOINT QC:\n")
print(outer1_save_qc)

stopifnot(all(outer1_save_qc))

cat("\nOUTER 1 CHECKPOINT: PASS\n")



###############################################################################
# ЕТАП 17. OUTER FOLD 2 — TRAINING AND EVALUATION
###############################################################################

outer_id <- 2

cat("\n====================================\n")
cat("OUTER FOLD:", outer_id, "\n")
cat("====================================\n")

# 1. Індекси
outer_data <- manifest_check$nested_cv[[outer_id]]

outer_train_idx <- outer_data$outer_train_idx
outer_valid_idx <- outer_data$outer_valid_idx

stopifnot(
  length(intersect(outer_train_idx, outer_valid_idx)) == 0,
  length(outer_train_idx) + length(outer_valid_idx) == 368
)

# 2. Оптимальна позиція lambda
best_index <- outer_tuning_results[[outer_id]]$best_lambda_index

cat("Selected lambda index:", best_index, "\n")
stopifnot(best_index %in% 1:50)

# 3. Відбір CpG на outer-training
cat("\nCalculating CpG variance...\n")

outer_variance <- DelayedMatrixStats::rowVars(
  X_train[, outer_train_idx, drop = FALSE],
  na.rm = TRUE
)

stopifnot(all(is.finite(outer_variance)))

outer_top_idx <- order(
  outer_variance,
  decreasing = TRUE
)[seq_len(5000)]

outer_selected_cpg <- rownames(X_train)[outer_top_idx]

stopifnot(
  length(outer_selected_cpg) == 5000,
  !anyDuplicated(outer_selected_cpg)
)

# 4. Матриці та labels
X_outer_train <- t(as.matrix(
  X_train[outer_top_idx, outer_train_idx, drop = FALSE]
))

X_outer_valid <- t(as.matrix(
  X_train[outer_top_idx, outer_valid_idx, drop = FALSE]
))

y_outer_train <- as.integer(
  y_train[outer_train_idx] == "Immune-high"
)

y_outer_valid <- as.integer(
  y_train[outer_valid_idx] == "Immune-high"
)

stopifnot(
  all(dim(X_outer_train) == c(length(outer_train_idx), 5000)),
  all(dim(X_outer_valid) == c(length(outer_valid_idx), 5000)),
  all(is.finite(X_outer_train)),
  all(is.finite(X_outer_valid))
)

# 5. Навчання Elastic Net
set.seed(2026 + outer_id)

outer_fit <- glmnet::glmnet(
  x = X_outer_train,
  y = y_outer_train,
  family = "binomial",
  alpha = 0.5,
  nlambda = 50,
  lambda.min.ratio = 0.01,
  standardize = TRUE
)

stopifnot(length(outer_fit$lambda) == 50)

selected_lambda <- outer_fit$lambda[best_index]

cat("Refitted lambda:", selected_lambda, "\n")

# 6. Прогнозування
outer_prob <- as.numeric(
  predict(
    outer_fit,
    newx = X_outer_valid,
    s = selected_lambda,
    type = "response"
  )
)

stopifnot(
  length(outer_prob) == length(outer_valid_idx),
  all(is.finite(outer_prob)),
  all(outer_prob >= 0 & outer_prob <= 1)
)

outer_pred <- as.integer(outer_prob >= 0.5)

# 7. Confusion matrix
outer_confusion <- table(
  Actual = factor(
    y_outer_valid,
    levels = c(0, 1),
    labels = c("Immune-low", "Immune-high")
  ),
  Predicted = factor(
    outer_pred,
    levels = c(0, 1),
    labels = c("Immune-low", "Immune-high")
  )
)

cat("\nCONFUSION MATRIX:\n")
print(outer_confusion)

# 8. Метрики
TN <- outer_confusion[1, 1]
FP <- outer_confusion[1, 2]
FN <- outer_confusion[2, 1]
TP <- outer_confusion[2, 2]

accuracy <- (TP + TN) / sum(outer_confusion)
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_

f1_score <- if (
  is.finite(precision) &&
  is.finite(sensitivity) &&
  precision + sensitivity > 0
) {
  2 * precision * sensitivity / (precision + sensitivity)
} else {
  NA_real_
}

balanced_accuracy <- mean(c(sensitivity, specificity))

outer_roc <- pROC::roc(
  response = factor(y_outer_valid, levels = c(0, 1)),
  predictor = outer_prob,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

outer_auc <- as.numeric(pROC::auc(outer_roc))

outer_metrics <- data.frame(
  Outer_Fold = outer_id,
  Training_n = length(outer_train_idx),
  Validation_n = length(outer_valid_idx),
  Selected_CpG = 5000,
  Best_Lambda_Index = best_index,
  Lambda = selected_lambda,
  Accuracy = accuracy,
  Sensitivity = sensitivity,
  Specificity = specificity,
  Precision = precision,
  F1 = f1_score,
  Balanced_Accuracy = balanced_accuracy,
  AUC = outer_auc
)

cat("\nOUTER FOLD 2 METRICS:\n")
print(outer_metrics, row.names = FALSE)

# 9. Прогнози
outer_predictions <- data.frame(
  Sample_ID = names(y_train)[outer_valid_idx],
  Actual = ifelse(y_outer_valid == 1, "Immune-high", "Immune-low"),
  Predicted = ifelse(outer_pred == 1, "Immune-high", "Immune-low"),
  Probability_Immune_high = outer_prob
)

# 10. QC
outer_qc <- c(
  correct_training_size =
    nrow(X_outer_train) == length(outer_train_idx),
  
  correct_validation_size =
    nrow(X_outer_valid) == length(outer_valid_idx),
  
  correct_feature_count =
    ncol(X_outer_train) == 5000,
  
  correct_prediction_count =
    length(outer_prob) == length(outer_valid_idx),
  
  finite_predictions =
    all(is.finite(outer_prob)),
  
  valid_probability_range =
    all(outer_prob >= 0 & outer_prob <= 1),
  
  finite_metrics =
    all(is.finite(as.numeric(
      outer_metrics[, c(
        "Accuracy", "Sensitivity", "Specificity",
        "Precision", "F1", "Balanced_Accuracy", "AUC"
      )]
    )))
)

cat("\nOUTER FOLD QC:\n")
print(outer_qc)

stopifnot(all(outer_qc))

# 11. Оновлення checkpoint — без перезапису Outer 1
outer_checkpoint_check <- readRDS(outer_checkpoint_file)

stopifnot(
  is.null(outer_checkpoint_check[["Outer_2"]]),
  !is.null(outer_checkpoint_check[["Outer_1"]])
)

outer_checkpoint_check[["Outer_2"]] <- list(
  outer_fold = outer_id,
  selected_cpg_ids = outer_selected_cpg,
  best_lambda_index = best_index,
  lambda = selected_lambda,
  metrics = outer_metrics,
  confusion_matrix = outer_confusion,
  predictions = outer_predictions,
  qc = outer_qc
)

saveRDS(
  outer_checkpoint_check,
  outer_checkpoint_file
)

# 12. Повторна перевірка
outer_checkpoint_verify <- readRDS(outer_checkpoint_file)

stopifnot(
  !is.null(outer_checkpoint_verify[["Outer_1"]]),
  !is.null(outer_checkpoint_verify[["Outer_2"]]),
  isTRUE(all(outer_checkpoint_verify[["Outer_1"]]$qc)),
  isTRUE(all(outer_checkpoint_verify[["Outer_2"]]$qc)),
  identical(
    outer_checkpoint_verify[["Outer_2"]]$metrics$AUC,
    outer_metrics$AUC
  )
)

cat("\nOUTER FOLD 2 CHECKPOINT: PASS\n")






###############################################################################
# ЕТАП 18. SAVE AND VERIFY OUTER FOLD 2
###############################################################################

outer_checkpoint_check <- readRDS(outer_checkpoint_file)

# Переконуємося, що Outer 1 збережений
stopifnot(
  !is.null(outer_checkpoint_check[["Outer_1"]]),
  isTRUE(all(outer_checkpoint_check[["Outer_1"]]$qc))
)

# Якщо Outer 2 ще не збережений — додаємо його
if (is.null(outer_checkpoint_check[["Outer_2"]])) {
  
  outer_checkpoint_check[["Outer_2"]] <- list(
    outer_fold = outer_id,
    selected_cpg_ids = outer_selected_cpg,
    best_lambda_index = best_index,
    lambda = selected_lambda,
    metrics = outer_metrics,
    confusion_matrix = outer_confusion,
    predictions = outer_predictions,
    qc = outer_qc
  )
  
  saveRDS(
    outer_checkpoint_check,
    outer_checkpoint_file
  )
  
} else {
  
  cat("Outer 2 already exists. Checking saved result.\n")
}

# Повторне завантаження
outer_checkpoint_verify <- readRDS(outer_checkpoint_file)

saved_outer2 <- outer_checkpoint_verify[["Outer_2"]]

save_qc <- c(
  outer1_preserved =
    !is.null(outer_checkpoint_verify[["Outer_1"]]),
  
  outer2_exists =
    !is.null(saved_outer2),
  
  correct_fold_id =
    identical(as.integer(saved_outer2$outer_fold), 2L),
  
  correct_CpG_count =
    length(saved_outer2$selected_cpg_ids) == 5000,
  
  correct_prediction_count =
    nrow(saved_outer2$predictions) == 74,
  
  metrics_preserved =
    isTRUE(all.equal(
      saved_outer2$metrics$AUC,
      outer_metrics$AUC
    )),
  
  QC_pass =
    isTRUE(all(saved_outer2$qc))
)

cat("\nOUTER 2 CHECKPOINT QC:\n")
print(save_qc)

stopifnot(all(save_qc))

cat("\nOUTER FOLD 2 CHECKPOINT: PASS\n")


###############################################################################
# ЕТАП 19. OUTER FOLD 3 — TRAINING AND EVALUATION
###############################################################################

outer_id <- 3

cat("\n====================================\n")
cat("OUTER FOLD:", outer_id, "\n")
cat("====================================\n")

# 1. Індекси
outer_data <- manifest_check$nested_cv[[outer_id]]

outer_train_idx <- outer_data$outer_train_idx
outer_valid_idx <- outer_data$outer_valid_idx

stopifnot(
  length(intersect(outer_train_idx, outer_valid_idx)) == 0,
  length(outer_train_idx) + length(outer_valid_idx) == 368
)

# 2. Lambda index із inner CV
best_index <- outer_tuning_results[[outer_id]]$best_lambda_index

cat("Selected lambda index:", best_index, "\n")
stopifnot(best_index %in% 1:50)

# 3. Відбір CpG лише на outer-training
cat("\nCalculating CpG variance...\n")

outer_variance <- DelayedMatrixStats::rowVars(
  X_train[, outer_train_idx, drop = FALSE],
  na.rm = TRUE
)

stopifnot(all(is.finite(outer_variance)))

outer_top_idx <- order(
  outer_variance,
  decreasing = TRUE
)[seq_len(5000)]

outer_selected_cpg <- rownames(X_train)[outer_top_idx]

stopifnot(
  length(outer_selected_cpg) == 5000,
  !anyDuplicated(outer_selected_cpg)
)

# 4. Матриці та labels
X_outer_train <- t(as.matrix(
  X_train[outer_top_idx, outer_train_idx, drop = FALSE]
))

X_outer_valid <- t(as.matrix(
  X_train[outer_top_idx, outer_valid_idx, drop = FALSE]
))

y_outer_train <- as.integer(
  y_train[outer_train_idx] == "Immune-high"
)

y_outer_valid <- as.integer(
  y_train[outer_valid_idx] == "Immune-high"
)

stopifnot(
  all(dim(X_outer_train) == c(length(outer_train_idx), 5000)),
  all(dim(X_outer_valid) == c(length(outer_valid_idx), 5000)),
  all(is.finite(X_outer_train)),
  all(is.finite(X_outer_valid))
)

# 5. Навчання Elastic Net
set.seed(2026 + outer_id)

outer_fit <- glmnet::glmnet(
  x = X_outer_train,
  y = y_outer_train,
  family = "binomial",
  alpha = 0.5,
  nlambda = 50,
  lambda.min.ratio = 0.01,
  standardize = TRUE
)

stopifnot(length(outer_fit$lambda) == 50)

selected_lambda <- outer_fit$lambda[best_index]

cat("Refitted lambda:", selected_lambda, "\n")

# 6. Прогнозування
outer_prob <- as.numeric(
  predict(
    outer_fit,
    newx = X_outer_valid,
    s = selected_lambda,
    type = "response"
  )
)

stopifnot(
  length(outer_prob) == length(outer_valid_idx),
  all(is.finite(outer_prob)),
  all(outer_prob >= 0 & outer_prob <= 1)
)

outer_pred <- as.integer(outer_prob >= 0.5)

# 7. Confusion matrix
outer_confusion <- table(
  Actual = factor(
    y_outer_valid,
    levels = c(0, 1),
    labels = c("Immune-low", "Immune-high")
  ),
  Predicted = factor(
    outer_pred,
    levels = c(0, 1),
    labels = c("Immune-low", "Immune-high")
  )
)

cat("\nCONFUSION MATRIX:\n")
print(outer_confusion)

# 8. Метрики
TN <- outer_confusion[1, 1]
FP <- outer_confusion[1, 2]
FN <- outer_confusion[2, 1]
TP <- outer_confusion[2, 2]

accuracy <- (TP + TN) / sum(outer_confusion)
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_

f1_score <- if (
  is.finite(precision) &&
  is.finite(sensitivity) &&
  precision + sensitivity > 0
) {
  2 * precision * sensitivity / (precision + sensitivity)
} else {
  NA_real_
}

balanced_accuracy <- mean(c(sensitivity, specificity))

outer_roc <- pROC::roc(
  response = factor(y_outer_valid, levels = c(0, 1)),
  predictor = outer_prob,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

outer_auc <- as.numeric(pROC::auc(outer_roc))

outer_metrics <- data.frame(
  Outer_Fold = outer_id,
  Training_n = length(outer_train_idx),
  Validation_n = length(outer_valid_idx),
  Selected_CpG = 5000,
  Best_Lambda_Index = best_index,
  Lambda = selected_lambda,
  Accuracy = accuracy,
  Sensitivity = sensitivity,
  Specificity = specificity,
  Precision = precision,
  F1 = f1_score,
  Balanced_Accuracy = balanced_accuracy,
  AUC = outer_auc
)

cat("\nOUTER FOLD 3 METRICS:\n")
print(outer_metrics, row.names = FALSE)

# 9. Прогнози
outer_predictions <- data.frame(
  Sample_ID = names(y_train)[outer_valid_idx],
  Actual = ifelse(y_outer_valid == 1, "Immune-high", "Immune-low"),
  Predicted = ifelse(outer_pred == 1, "Immune-high", "Immune-low"),
  Probability_Immune_high = outer_prob
)

# 10. QC
outer_qc <- c(
  correct_training_size =
    nrow(X_outer_train) == length(outer_train_idx),
  
  correct_validation_size =
    nrow(X_outer_valid) == length(outer_valid_idx),
  
  correct_feature_count =
    ncol(X_outer_train) == 5000,
  
  correct_prediction_count =
    length(outer_prob) == length(outer_valid_idx),
  
  finite_predictions =
    all(is.finite(outer_prob)),
  
  valid_probability_range =
    all(outer_prob >= 0 & outer_prob <= 1),
  
  finite_metrics =
    all(is.finite(as.numeric(
      outer_metrics[, c(
        "Accuracy", "Sensitivity", "Specificity",
        "Precision", "F1", "Balanced_Accuracy", "AUC"
      )]
    )))
)

cat("\nOUTER FOLD QC:\n")
print(outer_qc)

stopifnot(all(outer_qc))

# 11. Оновлення checkpoint
outer_checkpoint_check <- readRDS(outer_checkpoint_file)

stopifnot(
  is.null(outer_checkpoint_check[["Outer_3"]]),
  !is.null(outer_checkpoint_check[["Outer_1"]]),
  !is.null(outer_checkpoint_check[["Outer_2"]])
)

outer_checkpoint_check[["Outer_3"]] <- list(
  outer_fold = outer_id,
  selected_cpg_ids = outer_selected_cpg,
  best_lambda_index = best_index,
  lambda = selected_lambda,
  metrics = outer_metrics,
  confusion_matrix = outer_confusion,
  predictions = outer_predictions,
  qc = outer_qc
)

saveRDS(
  outer_checkpoint_check,
  outer_checkpoint_file
)

# 12. Перевірка checkpoint
outer_checkpoint_verify <- readRDS(outer_checkpoint_file)

stopifnot(
  all(vapply(
    outer_checkpoint_verify[1:3],
    function(x) !is.null(x) && isTRUE(all(x$qc)),
    logical(1)
  ))
)

cat("\nOUTER FOLD 3 CHECKPOINT: PASS\n")

###############################################################################
# ЕТАП 20. OUTER FOLD 3 — CHECKPOINT VERIFICATION
###############################################################################

outer_checkpoint_verify <- readRDS(outer_checkpoint_file)

outer3_saved <- outer_checkpoint_verify[["Outer_3"]]

stopifnot(!is.null(outer3_saved))

save_qc <- c(
  outer1_preserved =
    !is.null(outer_checkpoint_verify[["Outer_1"]]),
  
  outer2_preserved =
    !is.null(outer_checkpoint_verify[["Outer_2"]]),
  
  outer3_exists =
    !is.null(outer3_saved),
  
  correct_fold_id =
    identical(outer3_saved$outer_fold, 3L) ||
    identical(outer3_saved$outer_fold, 3),
  
  correct_CpG_count =
    length(outer3_saved$selected_cpg_ids) == 5000,
  
  correct_prediction_count =
    nrow(outer3_saved$predictions) == 74,
  
  metrics_preserved =
    isTRUE(all.equal(
      outer3_saved$metrics$AUC,
      outer_metrics$AUC
    )),
  
  QC_pass =
    all(outer3_saved$qc)
)

cat("\nOUTER 3 CHECKPOINT QC:\n")
print(save_qc)

stopifnot(all(save_qc))

cat("\nOUTER FOLD 3 CHECKPOINT: PASS\n")



###############################################################################
# ЕТАП 21. OUTER FOLD 4 — TRAINING AND EVALUATION
###############################################################################

outer_id <- 4

cat("\n====================================\n")
cat("OUTER FOLD:", outer_id, "\n")
cat("====================================\n")

# 1. Індекси
outer_data <- manifest_check$nested_cv[[outer_id]]

outer_train_idx <- outer_data$outer_train_idx
outer_valid_idx <- outer_data$outer_valid_idx

stopifnot(
  length(intersect(outer_train_idx, outer_valid_idx)) == 0,
  length(outer_train_idx) + length(outer_valid_idx) == 368
)

cat("Training samples:", length(outer_train_idx), "\n")
cat("Validation samples:", length(outer_valid_idx), "\n")

# 2. Lambda index із inner CV
best_index <- outer_tuning_results[[outer_id]]$best_lambda_index

cat("Selected lambda index:", best_index, "\n")
stopifnot(best_index %in% 1:50)

# 3. Відбір CpG лише на outer-training
cat("\nCalculating CpG variance...\n")

outer_variance <- DelayedMatrixStats::rowVars(
  X_train[, outer_train_idx, drop = FALSE],
  na.rm = TRUE
)

stopifnot(all(is.finite(outer_variance)))

outer_top_idx <- order(
  outer_variance,
  decreasing = TRUE
)[seq_len(5000)]

outer_selected_cpg <- rownames(X_train)[outer_top_idx]

stopifnot(
  length(outer_selected_cpg) == 5000,
  !anyDuplicated(outer_selected_cpg)
)

# 4. Матриці та labels
X_outer_train <- t(as.matrix(
  X_train[outer_top_idx, outer_train_idx, drop = FALSE]
))

X_outer_valid <- t(as.matrix(
  X_train[outer_top_idx, outer_valid_idx, drop = FALSE]
))

y_outer_train <- as.integer(
  y_train[outer_train_idx] == "Immune-high"
)

y_outer_valid <- as.integer(
  y_train[outer_valid_idx] == "Immune-high"
)

stopifnot(
  all(dim(X_outer_train) == c(length(outer_train_idx), 5000)),
  all(dim(X_outer_valid) == c(length(outer_valid_idx), 5000)),
  all(is.finite(X_outer_train)),
  all(is.finite(X_outer_valid))
)

# 5. Навчання Elastic Net
set.seed(2026 + outer_id)

outer_fit <- glmnet::glmnet(
  x = X_outer_train,
  y = y_outer_train,
  family = "binomial",
  alpha = 0.5,
  nlambda = 50,
  lambda.min.ratio = 0.01,
  standardize = TRUE
)

stopifnot(length(outer_fit$lambda) == 50)

selected_lambda <- outer_fit$lambda[best_index]

cat("Refitted lambda:", selected_lambda, "\n")

# 6. Прогнозування
outer_prob <- as.numeric(
  predict(
    outer_fit,
    newx = X_outer_valid,
    s = selected_lambda,
    type = "response"
  )
)

stopifnot(
  length(outer_prob) == length(outer_valid_idx),
  all(is.finite(outer_prob)),
  all(outer_prob >= 0 & outer_prob <= 1)
)

outer_pred <- as.integer(outer_prob >= 0.5)

# 7. Confusion matrix
outer_confusion <- table(
  Actual = factor(
    y_outer_valid,
    levels = c(0, 1),
    labels = c("Immune-low", "Immune-high")
  ),
  Predicted = factor(
    outer_pred,
    levels = c(0, 1),
    labels = c("Immune-low", "Immune-high")
  )
)

cat("\nCONFUSION MATRIX:\n")
print(outer_confusion)

# 8. Метрики
TN <- outer_confusion[1, 1]
FP <- outer_confusion[1, 2]
FN <- outer_confusion[2, 1]
TP <- outer_confusion[2, 2]

accuracy <- (TP + TN) / sum(outer_confusion)
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_

f1_score <- if (
  is.finite(precision) &&
  is.finite(sensitivity) &&
  precision + sensitivity > 0
) {
  2 * precision * sensitivity / (precision + sensitivity)
} else {
  NA_real_
}

balanced_accuracy <- mean(c(sensitivity, specificity))

outer_roc <- pROC::roc(
  response = factor(y_outer_valid, levels = c(0, 1)),
  predictor = outer_prob,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

outer_auc <- as.numeric(pROC::auc(outer_roc))

outer_metrics <- data.frame(
  Outer_Fold = outer_id,
  Training_n = length(outer_train_idx),
  Validation_n = length(outer_valid_idx),
  Selected_CpG = 5000,
  Best_Lambda_Index = best_index,
  Lambda = selected_lambda,
  Accuracy = accuracy,
  Sensitivity = sensitivity,
  Specificity = specificity,
  Precision = precision,
  F1 = f1_score,
  Balanced_Accuracy = balanced_accuracy,
  AUC = outer_auc
)

cat("\nOUTER FOLD 4 METRICS:\n")
print(outer_metrics, row.names = FALSE)

# 9. Прогнози
outer_predictions <- data.frame(
  Sample_ID = names(y_train)[outer_valid_idx],
  Actual = ifelse(y_outer_valid == 1, "Immune-high", "Immune-low"),
  Predicted = ifelse(outer_pred == 1, "Immune-high", "Immune-low"),
  Probability_Immune_high = outer_prob
)

# 10. QC
outer_qc <- c(
  correct_training_size =
    nrow(X_outer_train) == length(outer_train_idx),
  
  correct_validation_size =
    nrow(X_outer_valid) == length(outer_valid_idx),
  
  correct_feature_count =
    ncol(X_outer_train) == 5000,
  
  correct_prediction_count =
    length(outer_prob) == length(outer_valid_idx),
  
  finite_predictions =
    all(is.finite(outer_prob)),
  
  valid_probability_range =
    all(outer_prob >= 0 & outer_prob <= 1),
  
  finite_metrics =
    all(is.finite(as.numeric(
      outer_metrics[, c(
        "Accuracy", "Sensitivity", "Specificity",
        "Precision", "F1", "Balanced_Accuracy", "AUC"
      )]
    )))
)

cat("\nOUTER FOLD QC:\n")
print(outer_qc)

stopifnot(all(outer_qc))


###############################################################################
# ЕТАП 22. OUTER FOLD 4 — CHECKPOINT SAVE AND VERIFICATION
###############################################################################

outer_checkpoint_check <- readRDS(outer_checkpoint_file)

# Переконуємося, що попередні folds збережені
stopifnot(
  !is.null(outer_checkpoint_check[["Outer_1"]]),
  !is.null(outer_checkpoint_check[["Outer_2"]]),
  !is.null(outer_checkpoint_check[["Outer_3"]])
)

# Не перезаписуємо вже збережений fold
if (is.null(outer_checkpoint_check[["Outer_4"]])) {
  
  outer_checkpoint_check[["Outer_4"]] <- list(
    outer_fold = outer_id,
    selected_cpg_ids = outer_selected_cpg,
    best_lambda_index = best_index,
    lambda = selected_lambda,
    metrics = outer_metrics,
    confusion_matrix = outer_confusion,
    predictions = outer_predictions,
    qc = outer_qc
  )
  
  saveRDS(
    outer_checkpoint_check,
    outer_checkpoint_file
  )
  
  cat("Outer 4 saved successfully.\n")
  
} else {
  stop("Outer_4 already exists. Check its contents before overwriting.")
}

# Перевірка після збереження
outer_checkpoint_verify <- readRDS(outer_checkpoint_file)

outer4_saved <- outer_checkpoint_verify[["Outer_4"]]

save_qc <- c(
  outer1_preserved =
    !is.null(outer_checkpoint_verify[["Outer_1"]]),
  
  outer2_preserved =
    !is.null(outer_checkpoint_verify[["Outer_2"]]),
  
  outer3_preserved =
    !is.null(outer_checkpoint_verify[["Outer_3"]]),
  
  outer4_exists =
    !is.null(outer4_saved),
  
  correct_fold_id =
    identical(as.integer(outer4_saved$outer_fold), 4L),
  
  correct_CpG_count =
    length(outer4_saved$selected_cpg_ids) == 5000,
  
  correct_prediction_count =
    nrow(outer4_saved$predictions) == length(outer_valid_idx),
  
  metrics_preserved =
    isTRUE(all.equal(
      outer4_saved$metrics$AUC,
      outer_metrics$AUC
    )),
  
  QC_pass =
    all(outer4_saved$qc)
)

cat("\nOUTER 4 CHECKPOINT QC:\n")
print(save_qc)

stopifnot(all(save_qc))

cat("\nOUTER FOLD 4 CHECKPOINT: PASS\n")



###############################################################################
# ЕТАП 23. OUTER FOLD 5 — TRAINING AND EVALUATION
###############################################################################

outer_id <- 5

cat("\n====================================\n")
cat("OUTER FOLD:", outer_id, "\n")
cat("====================================\n")

# 1. Індекси
outer_data <- manifest_check$nested_cv[[outer_id]]

outer_train_idx <- outer_data$outer_train_idx
outer_valid_idx <- outer_data$outer_valid_idx

stopifnot(
  length(intersect(outer_train_idx, outer_valid_idx)) == 0,
  length(outer_train_idx) + length(outer_valid_idx) == 368
)

cat("Training samples:", length(outer_train_idx), "\n")
cat("Validation samples:", length(outer_valid_idx), "\n")

# 2. Lambda index із inner CV
best_index <- outer_tuning_results[[outer_id]]$best_lambda_index

cat("Selected lambda index:", best_index, "\n")
stopifnot(best_index %in% 1:50)

# 3. Відбір CpG лише на outer-training
cat("\nCalculating CpG variance...\n")

outer_variance <- DelayedMatrixStats::rowVars(
  X_train[, outer_train_idx, drop = FALSE],
  na.rm = TRUE
)

stopifnot(all(is.finite(outer_variance)))

outer_top_idx <- order(
  outer_variance,
  decreasing = TRUE
)[seq_len(5000)]

outer_selected_cpg <- rownames(X_train)[outer_top_idx]

stopifnot(
  length(outer_selected_cpg) == 5000,
  !anyDuplicated(outer_selected_cpg)
)

# 4. Матриці та labels
X_outer_train <- t(as.matrix(
  X_train[outer_top_idx, outer_train_idx, drop = FALSE]
))

X_outer_valid <- t(as.matrix(
  X_train[outer_top_idx, outer_valid_idx, drop = FALSE]
))

y_outer_train <- as.integer(
  y_train[outer_train_idx] == "Immune-high"
)

y_outer_valid <- as.integer(
  y_train[outer_valid_idx] == "Immune-high"
)

stopifnot(
  all(dim(X_outer_train) == c(length(outer_train_idx), 5000)),
  all(dim(X_outer_valid) == c(length(outer_valid_idx), 5000)),
  all(is.finite(X_outer_train)),
  all(is.finite(X_outer_valid))
)

# 5. Навчання Elastic Net
set.seed(2026 + outer_id)

outer_fit <- glmnet::glmnet(
  x = X_outer_train,
  y = y_outer_train,
  family = "binomial",
  alpha = 0.5,
  nlambda = 50,
  lambda.min.ratio = 0.01,
  standardize = TRUE
)

stopifnot(length(outer_fit$lambda) == 50)

selected_lambda <- outer_fit$lambda[best_index]

cat("Refitted lambda:", selected_lambda, "\n")

# 6. Прогнозування
outer_prob <- as.numeric(
  predict(
    outer_fit,
    newx = X_outer_valid,
    s = selected_lambda,
    type = "response"
  )
)

stopifnot(
  length(outer_prob) == length(outer_valid_idx),
  all(is.finite(outer_prob)),
  all(outer_prob >= 0 & outer_prob <= 1)
)

outer_pred <- as.integer(outer_prob >= 0.5)

# 7. Confusion matrix
outer_confusion <- table(
  Actual = factor(
    y_outer_valid,
    levels = c(0, 1),
    labels = c("Immune-low", "Immune-high")
  ),
  Predicted = factor(
    outer_pred,
    levels = c(0, 1),
    labels = c("Immune-low", "Immune-high")
  )
)

cat("\nCONFUSION MATRIX:\n")
print(outer_confusion)

# 8. Метрики
TN <- outer_confusion[1, 1]
FP <- outer_confusion[1, 2]
FN <- outer_confusion[2, 1]
TP <- outer_confusion[2, 2]

accuracy <- (TP + TN) / sum(outer_confusion)
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_

f1_score <- if (
  is.finite(precision) &&
  is.finite(sensitivity) &&
  precision + sensitivity > 0
) {
  2 * precision * sensitivity / (precision + sensitivity)
} else {
  NA_real_
}

balanced_accuracy <- mean(c(sensitivity, specificity))

outer_roc <- pROC::roc(
  response = factor(y_outer_valid, levels = c(0, 1)),
  predictor = outer_prob,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

outer_auc <- as.numeric(pROC::auc(outer_roc))

outer_metrics <- data.frame(
  Outer_Fold = outer_id,
  Training_n = length(outer_train_idx),
  Validation_n = length(outer_valid_idx),
  Selected_CpG = 5000,
  Best_Lambda_Index = best_index,
  Lambda = selected_lambda,
  Accuracy = accuracy,
  Sensitivity = sensitivity,
  Specificity = specificity,
  Precision = precision,
  F1 = f1_score,
  Balanced_Accuracy = balanced_accuracy,
  AUC = outer_auc
)

cat("\nOUTER FOLD 5 METRICS:\n")
print(outer_metrics, row.names = FALSE)

# 9. Прогнози
outer_predictions <- data.frame(
  Sample_ID = names(y_train)[outer_valid_idx],
  Actual = ifelse(y_outer_valid == 1, "Immune-high", "Immune-low"),
  Predicted = ifelse(outer_pred == 1, "Immune-high", "Immune-low"),
  Probability_Immune_high = outer_prob
)

# 10. QC
outer_qc <- c(
  correct_training_size =
    nrow(X_outer_train) == length(outer_train_idx),
  
  correct_validation_size =
    nrow(X_outer_valid) == length(outer_valid_idx),
  
  correct_feature_count =
    ncol(X_outer_train) == 5000,
  
  correct_prediction_count =
    length(outer_prob) == length(outer_valid_idx),
  
  finite_predictions =
    all(is.finite(outer_prob)),
  
  valid_probability_range =
    all(outer_prob >= 0 & outer_prob <= 1),
  
  finite_metrics =
    all(is.finite(as.numeric(
      outer_metrics[, c(
        "Accuracy", "Sensitivity", "Specificity",
        "Precision", "F1", "Balanced_Accuracy", "AUC"
      )]
    )))
)

cat("\nOUTER FOLD QC:\n")
print(outer_qc)

stopifnot(all(outer_qc))


###############################################################################
# ЕТАП 24. OUTER FOLD 5 — CHECKPOINT SAVE AND VERIFICATION
###############################################################################

outer_checkpoint_check <- readRDS(outer_checkpoint_file)

# Перевіряємо попередні folds
stopifnot(
  !is.null(outer_checkpoint_check[["Outer_1"]]),
  !is.null(outer_checkpoint_check[["Outer_2"]]),
  !is.null(outer_checkpoint_check[["Outer_3"]]),
  !is.null(outer_checkpoint_check[["Outer_4"]])
)

# Переконуємося, що Outer 5 ще не записаний
if (is.null(outer_checkpoint_check[["Outer_5"]])) {
  
  outer_checkpoint_check[["Outer_5"]] <- list(
    outer_fold = outer_id,
    selected_cpg_ids = outer_selected_cpg,
    best_lambda_index = best_index,
    lambda = selected_lambda,
    metrics = outer_metrics,
    confusion_matrix = outer_confusion,
    predictions = outer_predictions,
    qc = outer_qc
  )
  
  saveRDS(
    outer_checkpoint_check,
    outer_checkpoint_file
  )
  
  cat("Outer 5 saved successfully.\n")
  
} else {
  stop("Outer_5 already exists. Check its contents before overwriting.")
}

# Перевірка збереженого checkpoint
outer_checkpoint_verify <- readRDS(outer_checkpoint_file)

outer5_saved <- outer_checkpoint_verify[["Outer_5"]]

save_qc <- c(
  outer1_preserved =
    !is.null(outer_checkpoint_verify[["Outer_1"]]),
  
  outer2_preserved =
    !is.null(outer_checkpoint_verify[["Outer_2"]]),
  
  outer3_preserved =
    !is.null(outer_checkpoint_verify[["Outer_3"]]),
  
  outer4_preserved =
    !is.null(outer_checkpoint_verify[["Outer_4"]]),
  
  outer5_exists =
    !is.null(outer5_saved),
  
  correct_fold_id =
    identical(as.integer(outer5_saved$outer_fold), 5L),
  
  correct_CpG_count =
    length(outer5_saved$selected_cpg_ids) == 5000,
  
  correct_prediction_count =
    nrow(outer5_saved$predictions) == length(outer_valid_idx),
  
  metrics_preserved =
    isTRUE(all.equal(
      outer5_saved$metrics$AUC,
      outer_metrics$AUC
    )),
  
  QC_pass =
    all(outer5_saved$qc)
)

cat("\nOUTER 5 CHECKPOINT QC:\n")
print(save_qc)

stopifnot(all(save_qc))

cat("\nOUTER FOLD 5 CHECKPOINT: PASS\n")



###############################################################################
# ЕТАП 25. NESTED CV — FULL OUTER AUDIT AND OOF EVALUATION
###############################################################################

cat("\n============================================\n")
cat("NESTED CV: FULL OUTER AUDIT\n")
cat("============================================\n")

# 1. Завантаження checkpoint
outer_checkpoint_verify <- readRDS(outer_checkpoint_file)

stopifnot(
  length(outer_checkpoint_verify) == 5,
  all(paste0("Outer_", 1:5) %in% names(outer_checkpoint_verify))
)

# 2. Перевірка кожного fold
fold_audit <- do.call(
  rbind,
  lapply(1:5, function(i) {
    
    fold_name <- paste0("Outer_", i)
    fold <- outer_checkpoint_verify[[fold_name]]
    
    data.frame(
      Outer_Fold = i,
      Fold_exists = !is.null(fold),
      Correct_fold_ID = identical(
        as.integer(fold$outer_fold),
        as.integer(i)
      ),
      Correct_CpG_count =
        length(fold$selected_cpg_ids) == 5000,
      Correct_prediction_count =
        nrow(fold$predictions) ==
        length(manifest_check$nested_cv[[i]]$outer_valid_idx),
      Unique_sample_IDs =
        !anyDuplicated(fold$predictions$Sample_ID),
      Finite_probabilities =
        all(is.finite(
          fold$predictions$Probability_Immune_high
        )),
      QC_pass =
        all(fold$qc)
    )
  })
)

cat("\nFOLD AUDIT:\n")
print(fold_audit, row.names = FALSE)

stopifnot(
  all(fold_audit[, -1])
)

# 3. Об'єднання всіх OOF-прогнозів
oof_predictions <- do.call(
  rbind,
  lapply(1:5, function(i) {
    
    fold_predictions <-
      outer_checkpoint_verify[[paste0("Outer_", i)]]$predictions
    
    fold_predictions$Outer_Fold <- i
    
    fold_predictions
  })
)

rownames(oof_predictions) <- NULL

# 4. Перевірка кількості та унікальності зразків
expected_training_ids <- names(y_train)

oof_sample_ids <- oof_predictions$Sample_ID

oof_audit <- c(
  total_predictions =
    nrow(oof_predictions) == 368,
  
  unique_predictions =
    !anyDuplicated(oof_sample_ids),
  
  all_training_samples_present =
    setequal(oof_sample_ids, expected_training_ids),
  
  no_missing_sample_IDs =
    !anyNA(oof_sample_ids),
  
  no_missing_probabilities =
    all(is.finite(oof_predictions$Probability_Immune_high)),
  
  probabilities_in_range =
    all(
      oof_predictions$Probability_Immune_high >= 0 &
        oof_predictions$Probability_Immune_high <= 1
    ),
  
  no_test_sample_overlap =
    length(intersect(
      oof_sample_ids,
      names(y_test)
    )) == 0
)

cat("\nOOF SAMPLE AUDIT:\n")
print(oof_audit)

stopifnot(all(oof_audit))

# 5. Перевірка фактичних labels
expected_actual <- ifelse(
  y_train[oof_sample_ids] == "Immune-high",
  "Immune-high",
  "Immune-low"
)

labels_match <- all(
  oof_predictions$Actual == expected_actual
)

cat("\nActual labels match original labels:", labels_match, "\n")

stopifnot(labels_match)

# 6. Сортування за початковим порядком training samples
oof_predictions <- oof_predictions[
  match(expected_training_ids, oof_predictions$Sample_ID),
]

rownames(oof_predictions) <- NULL

stopifnot(
  identical(
    oof_predictions$Sample_ID,
    expected_training_ids
  )
)

# 7. Загальна confusion matrix
oof_confusion <- table(
  Actual = factor(
    oof_predictions$Actual,
    levels = c("Immune-low", "Immune-high")
  ),
  Predicted = factor(
    oof_predictions$Predicted,
    levels = c("Immune-low", "Immune-high")
  )
)

cat("\nPOOLED OOF CONFUSION MATRIX:\n")
print(oof_confusion)

# 8. Загальні метрики
TN <- oof_confusion[1, 1]
FP <- oof_confusion[1, 2]
FN <- oof_confusion[2, 1]
TP <- oof_confusion[2, 2]

oof_accuracy <- (TP + TN) / sum(oof_confusion)
oof_sensitivity <- TP / (TP + FN)
oof_specificity <- TN / (TN + FP)
oof_precision <- TP / (TP + FP)

oof_f1 <- 2 * oof_precision * oof_sensitivity /
  (oof_precision + oof_sensitivity)

oof_balanced_accuracy <- mean(
  c(oof_sensitivity, oof_specificity)
)

oof_roc <- pROC::roc(
  response = factor(
    oof_predictions$Actual,
    levels = c("Immune-low", "Immune-high")
  ),
  predictor = oof_predictions$Probability_Immune_high,
  levels = c("Immune-low", "Immune-high"),
  direction = "<",
  quiet = TRUE
)

oof_auc <- as.numeric(pROC::auc(oof_roc))

oof_metrics <- data.frame(
  Metric = c(
    "Accuracy",
    "Sensitivity",
    "Specificity",
    "Precision",
    "F1",
    "Balanced_Accuracy",
    "ROC_AUC"
  ),
  Value = c(
    oof_accuracy,
    oof_sensitivity,
    oof_specificity,
    oof_precision,
    oof_f1,
    oof_balanced_accuracy,
    oof_auc
  )
)

cat("\nPOOLED NESTED CV METRICS:\n")
print(oof_metrics, row.names = FALSE)

###############################################################################
# ЕТАП 25.9. FINAL NESTED OOF QC
###############################################################################

nested_oof_qc <- c(
  five_folds_passed =
    all(fold_audit$QC_pass),
  
  unique_368_OOF_predictions =
    nrow(oof_predictions) == 368 &&
    !anyDuplicated(oof_predictions$Sample_ID),
  
  all_training_samples_predicted =
    setequal(
      oof_predictions$Sample_ID,
      expected_training_ids
    ),
  
  no_test_samples_used =
    length(intersect(
      oof_predictions$Sample_ID,
      names(y_test)
    )) == 0,
  
  actual_labels_match =
    labels_match,
  
  finite_pooled_metrics =
    all(is.finite(oof_metrics$Value))
)

cat("\nFINAL NESTED OOF QC:\n")
print(nested_oof_qc)

stopifnot(all(nested_oof_qc))

cat("\nNESTED CV OOF AUDIT: PASS\n")

# Збереження прогнозів
write.csv(
  oof_predictions,
  "results/tables/TCGA_SKCM_nested_CV_OOF_predictions.csv",
  row.names = FALSE
)

# Збереження метрик
write.csv(
  oof_metrics,
  "results/tables/TCGA_SKCM_nested_CV_OOF_metrics.csv",
  row.names = FALSE
)

# Збереження повного результату
saveRDS(
  list(
    predictions = oof_predictions,
    confusion_matrix = oof_confusion,
    metrics = oof_metrics,
    fold_audit = fold_audit,
    oof_qc = nested_oof_qc
  ),
  "results/models/TCGA_SKCM_nested_CV_OOF_results.rds"
)

cat("\nOOF RESULTS SAVED SUCCESSFULLY.\n")


###############################################################################
# ЕТАП 27. ПЕРЕВІРКА FINAL TRAINING CV MANIFEST
###############################################################################

cv_manifest_file <- "results/models/TCGA_SKCM_ML_5fold_CV_manifest.rds"

stopifnot(file.exists(cv_manifest_file))

cv_manifest <- readRDS(cv_manifest_file)

cat("\nFINAL TRAINING CV MANIFEST:\n")

cat("Number of folds:", length(cv_manifest$cv_folds), "\n")
cat("Training samples:", length(cv_manifest$training_sample_ids), "\n")

# 1. Перевірка структури
stopifnot(
  length(cv_manifest$cv_folds) == 5,
  length(cv_manifest$training_sample_ids) == 368
)

# 2. Перевірка відповідності sample IDs
manifest_ids_match <- identical(
  as.character(cv_manifest$training_sample_ids),
  as.character(names(y_train))
)

# 3. Перевірка розмірів folds
fold_sizes <- vapply(
  cv_manifest$cv_folds,
  length,
  integer(1)
)

cat("\nFOLD SIZES:\n")
print(fold_sizes)

# 4. Перевірка, що всі індекси коректні
all_indices <- unlist(
  cv_manifest$cv_folds,
  use.names = FALSE
)

manifest_qc <- c(
  five_folds = length(cv_manifest$cv_folds) == 5,
  
  correct_training_ID_count =
    length(cv_manifest$training_sample_ids) == 368,
  
  training_IDs_match =
    manifest_ids_match,
  
  all_indices_in_range =
    all(all_indices >= 1 & all_indices <= 368),
  
  all_training_samples_assigned_once =
    setequal(as.integer(all_indices), seq_len(368)) &&
    !anyDuplicated(all_indices),
  
  correct_fold_sizes =
    identical(as.integer(fold_sizes), c(74L, 74L, 74L, 74L, 72L)),
  
  balanced_folds =
    all(vapply(
      cv_manifest$cv_folds,
      function(idx) {
        fold_labels <- y_train[idx]
        sum(fold_labels == "Immune-low") ==
          sum(fold_labels == "Immune-high")
      },
      logical(1)
    ))
)

cat("\nFINAL CV MANIFEST QC:\n")
print(manifest_qc)

stopifnot(all(manifest_qc))

cat("\nFINAL TRAINING CV MANIFEST: PASS\n")

###############################################################################
# ЕТАП 28. FINAL TRAINING TUNING — FOLD 1
###############################################################################

set.seed(20261106)

# 1. Визначаємо training та validation індекси
fold_id <- 1

valid_idx <- cv_manifest$cv_folds[[fold_id]]
train_idx <- setdiff(seq_len(ncol(X_train)), valid_idx)

cat("FINAL TUNING — FOLD", fold_id, "\n")
cat("Training samples:", length(train_idx), "\n")
cat("Validation samples:", length(valid_idx), "\n")

# 2. Перевіряємо баланс класів
cat("\nTRAINING LABELS:\n")
print(table(y_train[train_idx]))

cat("\nVALIDATION LABELS:\n")
print(table(y_train[valid_idx]))

stopifnot(
  length(train_idx) == 294,
  length(valid_idx) == 74,
  !anyDuplicated(train_idx),
  !anyDuplicated(valid_idx),
  length(intersect(train_idx, valid_idx)) == 0
)

# 3. Відбір 5000 найбільш варіабельних CpG
# Варіабельність розраховується тільки на training-частині фолда

X_inner_train <- X_train[, train_idx, drop = FALSE]

cpg_variance <- DelayedMatrixStats::rowVars(
  X_inner_train,
  na.rm = TRUE
)

stopifnot(
  length(cpg_variance) == nrow(X_train),
  all(is.finite(cpg_variance))
)

top_cpg_idx <- order(
  cpg_variance,
  decreasing = TRUE
)[seq_len(5000)]

selected_cpg_ids <- rownames(X_train)[top_cpg_idx]

stopifnot(
  length(selected_cpg_ids) == 5000,
  !anyDuplicated(selected_cpg_ids)
)

cat("\nSelected CpG:", length(selected_cpg_ids), "\n")

# 4. Формуємо матриці для glmnet
X_fold_train <- t(
  as.matrix(X_train[top_cpg_idx, train_idx, drop = FALSE])
)

X_fold_valid <- t(
  as.matrix(X_train[top_cpg_idx, valid_idx, drop = FALSE])
)

y_fold_train <- as.integer(
  y_train[train_idx] == "Immune-high"
)

y_fold_valid <- as.integer(
  y_train[valid_idx] == "Immune-high"
)

stopifnot(
  identical(dim(X_fold_train), c(294L, 5000L)),
  identical(dim(X_fold_valid), c(74L, 5000L)),
  !anyNA(X_fold_train),
  !anyNA(X_fold_valid),
  all(is.finite(X_fold_train)),
  all(is.finite(X_fold_valid))
)

cat("\nTRAIN MATRIX:", paste(dim(X_fold_train), collapse = " x "), "\n")
cat("VALIDATION MATRIX:", paste(dim(X_fold_valid), collapse = " x "), "\n")

# 5. Навчання Elastic Net
enet_fit <- glmnet::glmnet(
  x = X_fold_train,
  y = y_fold_train,
  family = "binomial",
  alpha = 0.5,
  nlambda = 50,
  lambda.min.ratio = 0.01,
  standardize = TRUE
)

stopifnot(length(enet_fit$lambda) == 50)

# 6. Прогноз для всіх 50 lambda
valid_probabilities <- predict(
  enet_fit,
  newx = X_fold_valid,
  type = "response"
)

valid_probabilities <- as.matrix(valid_probabilities)

stopifnot(
  nrow(valid_probabilities) == 74,
  ncol(valid_probabilities) == 50,
  all(is.finite(valid_probabilities))
)

# 7. Розрахунок log-loss для кожної lambda
epsilon <- 1e-15

valid_probabilities <- pmin(
  pmax(valid_probabilities, epsilon),
  1 - epsilon
)

fold_logloss <- vapply(
  seq_len(ncol(valid_probabilities)),
  function(j) {
    p <- valid_probabilities[, j]
    
    -mean(
      y_fold_valid * log(p) +
        (1 - y_fold_valid) * log(1 - p)
    )
  },
  numeric(1)
)

stopifnot(
  length(fold_logloss) == 50,
  all(is.finite(fold_logloss))
)

# 8. Зберігаємо результат фолда
fold_result <- list(
  fold_id = fold_id,
  training_indices = train_idx,
  validation_indices = valid_idx,
  training_sample_ids = names(y_train)[train_idx],
  validation_sample_ids = names(y_train)[valid_idx],
  selected_cpg_ids = selected_cpg_ids,
  alpha = 0.5,
  n_features = 5000,
  lambda = enet_fit$lambda,
  logloss = fold_logloss,
  best_lambda_index = which.min(fold_logloss),
  best_logloss = min(fold_logloss),
  qc_pass = TRUE
)

# 9. Checkpoint — не втрачаємо результати попередніх фолдів
checkpoint_file <- paste0(
  "results/models/TCGA_SKCM_final_training_tuning_partial.rds"
)

dir.create(
  dirname(checkpoint_file),
  recursive = TRUE,
  showWarnings = FALSE
)

if (file.exists(checkpoint_file)) {
  final_tuning_fold_results <- readRDS(checkpoint_file)
} else {
  final_tuning_fold_results <- vector("list", 5)
}

final_tuning_fold_results[[fold_id]] <- fold_result

saveRDS(
  final_tuning_fold_results,
  checkpoint_file
)

# 10. Виводимо результат
cat("\nFOLD 1 RESULTS:\n")
cat("Best lambda index:", fold_result$best_lambda_index, "\n")
cat("Best lambda:", enet_fit$lambda[fold_result$best_lambda_index], "\n")
cat("Best validation log-loss:", fold_result$best_logloss, "\n")

cat("\nCHECKPOINT SAVED:", file.exists(checkpoint_file), "\n")
cat("FINAL TUNING FOLD 1: PASS\n")



###############################################################################
# ЕТАП 29. FINAL TRAINING TUNING — FOLD 2
###############################################################################

set.seed(20261106)

fold_id <- 2

# 1. Перевірка checkpoint
checkpoint_file <- "results/models/TCGA_SKCM_final_training_tuning_partial.rds"

stopifnot(file.exists(checkpoint_file))

final_tuning_fold_results <- readRDS(checkpoint_file)

stopifnot(
  !is.null(final_tuning_fold_results[[1]]),
  isTRUE(final_tuning_fold_results[[1]]$qc_pass)
)

# 2. Training та validation індекси
valid_idx <- cv_manifest$cv_folds[[fold_id]]
train_idx <- setdiff(seq_len(ncol(X_train)), valid_idx)

cat("FINAL TUNING — FOLD", fold_id, "\n")
cat("Training samples:", length(train_idx), "\n")
cat("Validation samples:", length(valid_idx), "\n")

print(table(y_train[train_idx]))
print(table(y_train[valid_idx]))

stopifnot(
  length(train_idx) == 294,
  length(valid_idx) == 74,
  length(intersect(train_idx, valid_idx)) == 0
)

# 3. Відбір 5000 CpG тільки на training-частині
X_inner_train <- X_train[, train_idx, drop = FALSE]

cpg_variance <- DelayedMatrixStats::rowVars(
  X_inner_train,
  na.rm = TRUE
)

stopifnot(all(is.finite(cpg_variance)))

top_cpg_idx <- order(
  cpg_variance,
  decreasing = TRUE
)[seq_len(5000)]

selected_cpg_ids <- rownames(X_train)[top_cpg_idx]

stopifnot(
  length(selected_cpg_ids) == 5000,
  !anyDuplicated(selected_cpg_ids)
)

# 4. Підготовка матриць
X_fold_train <- t(
  as.matrix(X_train[top_cpg_idx, train_idx, drop = FALSE])
)

X_fold_valid <- t(
  as.matrix(X_train[top_cpg_idx, valid_idx, drop = FALSE])
)

y_fold_train <- as.integer(
  y_train[train_idx] == "Immune-high"
)

y_fold_valid <- as.integer(
  y_train[valid_idx] == "Immune-high"
)

stopifnot(
  identical(dim(X_fold_train), c(294L, 5000L)),
  identical(dim(X_fold_valid), c(74L, 5000L)),
  !anyNA(X_fold_train),
  !anyNA(X_fold_valid),
  all(is.finite(X_fold_train)),
  all(is.finite(X_fold_valid))
)

# 5. Навчання Elastic Net
enet_fit <- glmnet::glmnet(
  x = X_fold_train,
  y = y_fold_train,
  family = "binomial",
  alpha = 0.5,
  nlambda = 50,
  lambda.min.ratio = 0.01,
  standardize = TRUE
)

stopifnot(length(enet_fit$lambda) == 50)

# 6. Прогноз
valid_probabilities <- predict(
  enet_fit,
  newx = X_fold_valid,
  type = "response"
)

valid_probabilities <- as.matrix(valid_probabilities)

stopifnot(
  nrow(valid_probabilities) == 74,
  ncol(valid_probabilities) == 50,
  all(is.finite(valid_probabilities))
)

# 7. Log-loss
epsilon <- 1e-15

valid_probabilities <- pmin(
  pmax(valid_probabilities, epsilon),
  1 - epsilon
)

fold_logloss <- vapply(
  seq_len(50),
  function(j) {
    p <- valid_probabilities[, j]
    
    -mean(
      y_fold_valid * log(p) +
        (1 - y_fold_valid) * log(1 - p)
    )
  },
  numeric(1)
)

stopifnot(
  length(fold_logloss) == 50,
  all(is.finite(fold_logloss))
)

# 8. Результат Fold 2
fold_result <- list(
  fold_id = fold_id,
  training_indices = train_idx,
  validation_indices = valid_idx,
  training_sample_ids = names(y_train)[train_idx],
  validation_sample_ids = names(y_train)[valid_idx],
  selected_cpg_ids = selected_cpg_ids,
  alpha = 0.5,
  n_features = 5000,
  lambda = enet_fit$lambda,
  logloss = fold_logloss,
  best_lambda_index = which.min(fold_logloss),
  best_logloss = min(fold_logloss),
  qc_pass = TRUE
)

# 9. Оновлення checkpoint
final_tuning_fold_results[[fold_id]] <- fold_result

saveRDS(
  final_tuning_fold_results,
  checkpoint_file
)

# 10. Перевірка збереження
checkpoint_check <- readRDS(checkpoint_file)

stopifnot(
  !is.null(checkpoint_check[[1]]),
  !is.null(checkpoint_check[[2]]),
  isTRUE(checkpoint_check[[1]]$qc_pass),
  isTRUE(checkpoint_check[[2]]$qc_pass),
  length(checkpoint_check[[2]]$logloss) == 50
)

cat("\nFOLD 2 RESULTS:\n")
cat("Best lambda index:", fold_result$best_lambda_index, "\n")
cat("Best lambda:", enet_fit$lambda[fold_result$best_lambda_index], "\n")
cat("Best validation log-loss:", fold_result$best_logloss, "\n")

cat("\nCHECKPOINT SAVED:", file.exists(checkpoint_file), "\n")
cat("FOLD 1 PRESERVED:", !is.null(checkpoint_check[[1]]), "\n")
cat("FINAL TUNING FOLD 2: PASS\n")

###############################################################################
# FINAL TRAINING TUNING — FOLD 3
###############################################################################

set.seed(20261106)

fold_id <- 3

# 1. Завантаження checkpoint
checkpoint_file <- "results/models/TCGA_SKCM_final_training_tuning_partial.rds"

stopifnot(file.exists(checkpoint_file))

final_tuning_fold_results <- readRDS(checkpoint_file)

stopifnot(
  isTRUE(final_tuning_fold_results[[1]]$qc_pass),
  isTRUE(final_tuning_fold_results[[2]]$qc_pass)
)

# 2. Індекси
valid_idx <- cv_manifest$cv_folds[[fold_id]]
train_idx <- setdiff(seq_len(ncol(X_train)), valid_idx)

cat("FINAL TUNING — FOLD", fold_id, "\n")
cat("Training:", length(train_idx), "\n")
cat("Validation:", length(valid_idx), "\n")

print(table(y_train[train_idx]))
print(table(y_train[valid_idx]))

stopifnot(
  length(train_idx) == 294,
  length(valid_idx) == 74,
  length(intersect(train_idx, valid_idx)) == 0
)

# 3. Відбір CpG
X_inner_train <- X_train[, train_idx, drop = FALSE]

cpg_variance <- DelayedMatrixStats::rowVars(
  X_inner_train,
  na.rm = TRUE
)

stopifnot(all(is.finite(cpg_variance)))

top_cpg_idx <- order(
  cpg_variance,
  decreasing = TRUE
)[seq_len(5000)]

selected_cpg_ids <- rownames(X_train)[top_cpg_idx]

stopifnot(
  length(selected_cpg_ids) == 5000,
  !anyDuplicated(selected_cpg_ids)
)

# 4. Матриці
X_fold_train <- t(
  as.matrix(X_train[top_cpg_idx, train_idx, drop = FALSE])
)

X_fold_valid <- t(
  as.matrix(X_train[top_cpg_idx, valid_idx, drop = FALSE])
)

y_fold_train <- as.integer(
  y_train[train_idx] == "Immune-high"
)

y_fold_valid <- as.integer(
  y_train[valid_idx] == "Immune-high"
)

stopifnot(
  identical(dim(X_fold_train), c(294L, 5000L)),
  identical(dim(X_fold_valid), c(74L, 5000L)),
  !anyNA(X_fold_train),
  !anyNA(X_fold_valid),
  all(is.finite(X_fold_train)),
  all(is.finite(X_fold_valid))
)

# 5. Elastic Net
enet_fit <- glmnet::glmnet(
  x = X_fold_train,
  y = y_fold_train,
  family = "binomial",
  alpha = 0.5,
  nlambda = 50,
  lambda.min.ratio = 0.01,
  standardize = TRUE
)

stopifnot(length(enet_fit$lambda) == 50)

# 6. Прогнози
valid_probabilities <- predict(
  enet_fit,
  newx = X_fold_valid,
  type = "response"
)

valid_probabilities <- as.matrix(valid_probabilities)

stopifnot(
  nrow(valid_probabilities) == 74,
  ncol(valid_probabilities) == 50,
  all(is.finite(valid_probabilities))
)

# 7. Log-loss
epsilon <- 1e-15

valid_probabilities <- pmin(
  pmax(valid_probabilities, epsilon),
  1 - epsilon
)

fold_logloss <- vapply(
  seq_len(50),
  function(j) {
    p <- valid_probabilities[, j]
    
    -mean(
      y_fold_valid * log(p) +
        (1 - y_fold_valid) * log(1 - p)
    )
  },
  numeric(1)
)

stopifnot(
  length(fold_logloss) == 50,
  all(is.finite(fold_logloss))
)

# 8. Збереження результату
fold_result <- list(
  fold_id = fold_id,
  training_indices = train_idx,
  validation_indices = valid_idx,
  training_sample_ids = names(y_train)[train_idx],
  validation_sample_ids = names(y_train)[valid_idx],
  selected_cpg_ids = selected_cpg_ids,
  alpha = 0.5,
  n_features = 5000,
  lambda = enet_fit$lambda,
  logloss = fold_logloss,
  best_lambda_index = which.min(fold_logloss),
  best_logloss = min(fold_logloss),
  qc_pass = TRUE
)

final_tuning_fold_results[[fold_id]] <- fold_result

saveRDS(
  final_tuning_fold_results,
  checkpoint_file
)

# 9. Перевірка checkpoint
checkpoint_check <- readRDS(checkpoint_file)

stopifnot(
  all(vapply(
    checkpoint_check[1:3],
    function(x) !is.null(x) && isTRUE(x$qc_pass),
    logical(1)
  )),
  all(vapply(
    checkpoint_check[1:3],
    function(x) length(x$logloss) == 50,
    logical(1)
  ))
)

# 10. Результати
cat("\nFOLD 3 RESULTS:\n")
cat("Best lambda index:", fold_result$best_lambda_index, "\n")
cat("Best lambda:", enet_fit$lambda[fold_result$best_lambda_index], "\n")
cat("Best validation log-loss:", fold_result$best_logloss, "\n")

cat("\nCHECKPOINT SAVED:", file.exists(checkpoint_file), "\n")
cat("FOLDS 1-3 PRESERVED: TRUE\n")
cat("FINAL TUNING FOLD 3: PASS\n")

###############################################################################
# FINAL TRAINING TUNING — FOLD 4
###############################################################################

set.seed(20261106)

fold_id <- 4

# 1. Завантаження та перевірка checkpoint
checkpoint_file <- "results/models/TCGA_SKCM_final_training_tuning_partial.rds"

stopifnot(file.exists(checkpoint_file))

final_tuning_fold_results <- readRDS(checkpoint_file)

stopifnot(
  all(vapply(
    final_tuning_fold_results[1:3],
    function(x) !is.null(x) && isTRUE(x$qc_pass),
    logical(1)
  ))
)

# 2. Індекси
valid_idx <- cv_manifest$cv_folds[[fold_id]]
train_idx <- setdiff(seq_len(ncol(X_train)), valid_idx)

cat("FINAL TUNING — FOLD", fold_id, "\n")
cat("Training:", length(train_idx), "\n")
cat("Validation:", length(valid_idx), "\n")

print(table(y_train[train_idx]))
print(table(y_train[valid_idx]))

stopifnot(
  length(train_idx) == 294,
  length(valid_idx) == 74,
  length(intersect(train_idx, valid_idx)) == 0
)

# 3. Відбір 5000 CpG на training-частині
X_inner_train <- X_train[, train_idx, drop = FALSE]

cpg_variance <- DelayedMatrixStats::rowVars(
  X_inner_train,
  na.rm = TRUE
)

stopifnot(all(is.finite(cpg_variance)))

top_cpg_idx <- order(
  cpg_variance,
  decreasing = TRUE
)[seq_len(5000)]

selected_cpg_ids <- rownames(X_train)[top_cpg_idx]

stopifnot(
  length(selected_cpg_ids) == 5000,
  !anyDuplicated(selected_cpg_ids)
)

# 4. Матриці
X_fold_train <- t(
  as.matrix(X_train[top_cpg_idx, train_idx, drop = FALSE])
)

X_fold_valid <- t(
  as.matrix(X_train[top_cpg_idx, valid_idx, drop = FALSE])
)

y_fold_train <- as.integer(
  y_train[train_idx] == "Immune-high"
)

y_fold_valid <- as.integer(
  y_train[valid_idx] == "Immune-high"
)

stopifnot(
  identical(dim(X_fold_train), c(294L, 5000L)),
  identical(dim(X_fold_valid), c(74L, 5000L)),
  !anyNA(X_fold_train),
  !anyNA(X_fold_valid),
  all(is.finite(X_fold_train)),
  all(is.finite(X_fold_valid))
)

# 5. Elastic Net
enet_fit <- glmnet::glmnet(
  x = X_fold_train,
  y = y_fold_train,
  family = "binomial",
  alpha = 0.5,
  nlambda = 50,
  lambda.min.ratio = 0.01,
  standardize = TRUE
)

stopifnot(length(enet_fit$lambda) == 50)

# 6. Прогнози
valid_probabilities <- predict(
  enet_fit,
  newx = X_fold_valid,
  type = "response"
)

valid_probabilities <- as.matrix(valid_probabilities)

stopifnot(
  nrow(valid_probabilities) == 74,
  ncol(valid_probabilities) == 50,
  all(is.finite(valid_probabilities))
)

# 7. Log-loss
epsilon <- 1e-15

valid_probabilities <- pmin(
  pmax(valid_probabilities, epsilon),
  1 - epsilon
)

fold_logloss <- vapply(
  seq_len(50),
  function(j) {
    p <- valid_probabilities[, j]
    
    -mean(
      y_fold_valid * log(p) +
        (1 - y_fold_valid) * log(1 - p)
    )
  },
  numeric(1)
)

stopifnot(
  length(fold_logloss) == 50,
  all(is.finite(fold_logloss))
)

# 8. Результат Fold 4
fold_result <- list(
  fold_id = fold_id,
  training_indices = train_idx,
  validation_indices = valid_idx,
  training_sample_ids = names(y_train)[train_idx],
  validation_sample_ids = names(y_train)[valid_idx],
  selected_cpg_ids = selected_cpg_ids,
  alpha = 0.5,
  n_features = 5000,
  lambda = enet_fit$lambda,
  logloss = fold_logloss,
  best_lambda_index = which.min(fold_logloss),
  best_logloss = min(fold_logloss),
  qc_pass = TRUE
)

# 9. Збереження checkpoint
final_tuning_fold_results[[fold_id]] <- fold_result

saveRDS(
  final_tuning_fold_results,
  checkpoint_file
)

# 10. Перевірка збереження
checkpoint_check <- readRDS(checkpoint_file)

stopifnot(
  all(vapply(
    checkpoint_check[1:4],
    function(x) !is.null(x) && isTRUE(x$qc_pass),
    logical(1)
  )),
  all(vapply(
    checkpoint_check[1:4],
    function(x) length(x$logloss) == 50,
    logical(1)
  ))
)

cat("\nFOLD 4 RESULTS:\n")
cat("Best lambda index:", fold_result$best_lambda_index, "\n")
cat("Best lambda:", enet_fit$lambda[fold_result$best_lambda_index], "\n")
cat("Best validation log-loss:", fold_result$best_logloss, "\n")

cat("\nCHECKPOINT SAVED:", file.exists(checkpoint_file), "\n")
cat("FOLDS 1-4 PRESERVED: TRUE\n")
cat("FINAL TUNING FOLD 4: PASS\n")


###############################################################################
# FINAL TRAINING TUNING — FOLD 5
###############################################################################

set.seed(20261106)

fold_id <- 5

# 1. Завантаження та перевірка checkpoint
checkpoint_file <- "results/models/TCGA_SKCM_final_training_tuning_partial.rds"

stopifnot(file.exists(checkpoint_file))

final_tuning_fold_results <- readRDS(checkpoint_file)

stopifnot(
  all(vapply(
    final_tuning_fold_results[1:4],
    function(x) !is.null(x) && isTRUE(x$qc_pass),
    logical(1)
  ))
)

# 2. Індекси
valid_idx <- cv_manifest$cv_folds[[fold_id]]
train_idx <- setdiff(seq_len(ncol(X_train)), valid_idx)

cat("FINAL TUNING — FOLD", fold_id, "\n")
cat("Training:", length(train_idx), "\n")
cat("Validation:", length(valid_idx), "\n")

print(table(y_train[train_idx]))
print(table(y_train[valid_idx]))

stopifnot(
  length(train_idx) == 296,
  length(valid_idx) == 72,
  length(intersect(train_idx, valid_idx)) == 0
)

# 3. Відбір 5000 CpG тільки на training-частині
X_inner_train <- X_train[, train_idx, drop = FALSE]

cpg_variance <- DelayedMatrixStats::rowVars(
  X_inner_train,
  na.rm = TRUE
)

stopifnot(all(is.finite(cpg_variance)))

top_cpg_idx <- order(
  cpg_variance,
  decreasing = TRUE
)[seq_len(5000)]

selected_cpg_ids <- rownames(X_train)[top_cpg_idx]

stopifnot(
  length(selected_cpg_ids) == 5000,
  !anyDuplicated(selected_cpg_ids)
)

# 4. Матриці
X_fold_train <- t(
  as.matrix(X_train[top_cpg_idx, train_idx, drop = FALSE])
)

X_fold_valid <- t(
  as.matrix(X_train[top_cpg_idx, valid_idx, drop = FALSE])
)

y_fold_train <- as.integer(
  y_train[train_idx] == "Immune-high"
)

y_fold_valid <- as.integer(
  y_train[valid_idx] == "Immune-high"
)

stopifnot(
  identical(dim(X_fold_train), c(296L, 5000L)),
  identical(dim(X_fold_valid), c(72L, 5000L)),
  !anyNA(X_fold_train),
  !anyNA(X_fold_valid),
  all(is.finite(X_fold_train)),
  all(is.finite(X_fold_valid))
)

# 5. Elastic Net
enet_fit <- glmnet::glmnet(
  x = X_fold_train,
  y = y_fold_train,
  family = "binomial",
  alpha = 0.5,
  nlambda = 50,
  lambda.min.ratio = 0.01,
  standardize = TRUE
)

stopifnot(length(enet_fit$lambda) == 50)

# 6. Прогнози
valid_probabilities <- predict(
  enet_fit,
  newx = X_fold_valid,
  type = "response"
)

valid_probabilities <- as.matrix(valid_probabilities)

stopifnot(
  nrow(valid_probabilities) == 72,
  ncol(valid_probabilities) == 50,
  all(is.finite(valid_probabilities))
)

# 7. Log-loss
epsilon <- 1e-15

valid_probabilities <- pmin(
  pmax(valid_probabilities, epsilon),
  1 - epsilon
)

fold_logloss <- vapply(
  seq_len(50),
  function(j) {
    p <- valid_probabilities[, j]
    
    -mean(
      y_fold_valid * log(p) +
        (1 - y_fold_valid) * log(1 - p)
    )
  },
  numeric(1)
)

stopifnot(
  length(fold_logloss) == 50,
  all(is.finite(fold_logloss))
)

# 8. Результат Fold 5
fold_result <- list(
  fold_id = fold_id,
  training_indices = train_idx,
  validation_indices = valid_idx,
  training_sample_ids = names(y_train)[train_idx],
  validation_sample_ids = names(y_train)[valid_idx],
  selected_cpg_ids = selected_cpg_ids,
  alpha = 0.5,
  n_features = 5000,
  lambda = enet_fit$lambda,
  logloss = fold_logloss,
  best_lambda_index = which.min(fold_logloss),
  best_logloss = min(fold_logloss),
  qc_pass = TRUE
)

# 9. Збереження checkpoint
final_tuning_fold_results[[fold_id]] <- fold_result

saveRDS(
  final_tuning_fold_results,
  checkpoint_file
)

# 10. Перевірка всіх п'яти фолдів
checkpoint_check <- readRDS(checkpoint_file)

stopifnot(
  length(checkpoint_check) == 5,
  all(vapply(
    checkpoint_check,
    function(x) !is.null(x) && isTRUE(x$qc_pass),
    logical(1)
  )),
  all(vapply(
    checkpoint_check,
    function(x) length(x$logloss) == 50,
    logical(1)
  ))
)

cat("\nFOLD 5 RESULTS:\n")
cat("Best lambda index:", fold_result$best_lambda_index, "\n")
cat("Best lambda:", enet_fit$lambda[fold_result$best_lambda_index], "\n")
cat("Best validation log-loss:", fold_result$best_logloss, "\n")

cat("\nCHECKPOINT SAVED:", file.exists(checkpoint_file), "\n")
cat("ALL FIVE FOLDS PRESERVED: TRUE\n")
cat("FINAL TRAINING TUNING: PASS\n")


###############################################################################
# ЕТАП 33. AGGREGATE FINAL TRAINING CV RESULTS
###############################################################################

# 1. Завантаження checkpoint
checkpoint_file <- "results/models/TCGA_SKCM_final_training_tuning_partial.rds"

stopifnot(file.exists(checkpoint_file))

final_tuning_fold_results <- readRDS(checkpoint_file)

# 2. Перевірка п'яти фолдів
stopifnot(
  length(final_tuning_fold_results) == 5,
  all(vapply(
    final_tuning_fold_results,
    function(x) !is.null(x) && isTRUE(x$qc_pass),
    logical(1)
  ))
)

# 3. Формуємо матрицю log-loss: 5 фолдів × 50 lambda
logloss_matrix <- do.call(
  rbind,
  lapply(
    final_tuning_fold_results,
    function(x) x$logloss
  )
)

stopifnot(
  identical(dim(logloss_matrix), c(5L, 50L)),
  all(is.finite(logloss_matrix))
)

rownames(logloss_matrix) <- paste0("Fold_", 1:5)
colnames(logloss_matrix) <- paste0("Lambda_", 1:50)

# 4. Середній log-loss та SD для кожної позиції lambda
mean_logloss <- colMeans(logloss_matrix)

sd_logloss <- apply(
  logloss_matrix,
  2,
  sd
)

# 5. Вибір позиції з мінімальним середнім log-loss
best_lambda_index <- which.min(mean_logloss)

best_mean_logloss <- mean_logloss[best_lambda_index]
best_sd_logloss <- sd_logloss[best_lambda_index]

# 6. Таблиця результатів
final_tuning_summary <- data.frame(
  Lambda_Index = seq_len(50),
  Mean_LogLoss = mean_logloss,
  SD_LogLoss = sd_logloss
)

final_tuning_summary <- final_tuning_summary[
  order(final_tuning_summary$Mean_LogLoss),
]

# 7. Перевірки
stopifnot(
  length(mean_logloss) == 50,
  length(sd_logloss) == 50,
  best_lambda_index >= 1,
  best_lambda_index <= 50,
  is.finite(best_mean_logloss),
  is.finite(best_sd_logloss),
  best_mean_logloss == min(mean_logloss)
)

# 8. Збереження результатів
saveRDS(
  list(
    logloss_matrix = logloss_matrix,
    mean_logloss = mean_logloss,
    sd_logloss = sd_logloss,
    best_lambda_index = best_lambda_index,
    best_mean_logloss = best_mean_logloss,
    best_sd_logloss = best_sd_logloss,
    fold_results = final_tuning_fold_results,
    alpha = 0.5,
    n_features = 5000,
    n_folds = 5
  ),
  "results/models/TCGA_SKCM_final_training_tuning.rds"
)

write.csv(
  final_tuning_summary,
  "results/tables/TCGA_SKCM_final_training_tuning_summary.csv",
  row.names = FALSE
)

# 9. Виведення результатів
cat("\nFINAL TRAINING TUNING RESULTS:\n")

cat(
  "Selected lambda index:",
  best_lambda_index,
  "of 50\n"
)

cat(
  "Mean validation log-loss:",
  best_mean_logloss,
  "\n"
)

cat(
  "SD validation log-loss:",
  best_sd_logloss,
  "\n"
)

cat("\nBEST 10 LAMBDA POSITIONS:\n")
print(head(final_tuning_summary, 10), row.names = FALSE)

# 10. Фінальний QC
final_tuning_qc <- c(
  five_folds_present = length(final_tuning_fold_results) == 5,
  fifty_lambda_positions = ncol(logloss_matrix) == 50,
  all_logloss_finite = all(is.finite(logloss_matrix)),
  best_index_valid = best_lambda_index %in% 1:50,
  summary_saved = file.exists(
    "results/tables/TCGA_SKCM_final_training_tuning_summary.csv"
  ),
  model_results_saved = file.exists(
    "results/models/TCGA_SKCM_final_training_tuning.rds"
  )
)

cat("\nFINAL TUNING QC:\n")
print(final_tuning_qc)

stopifnot(all(final_tuning_qc))

cat("\nFINAL TRAINING TUNING AGGREGATION: PASS\n")


###############################################################################
# ЕТАП 34. FINAL ELASTIC NET MODEL — TRAINING ON 368 SAMPLES
###############################################################################

# 1. Перевірка вхідних даних та результатів tuning

stopifnot(
  ncol(X_train) == 368,
  length(y_train) == 368,
  length(best_lambda_index) == 1,
  best_lambda_index == 34,
  file.exists("results/models/TCGA_SKCM_final_training_tuning.rds")
)

cat("Training samples:", ncol(X_train), "\n")
cat("Training labels:\n")
print(table(y_train))

# 2. Відбір 5000 найбільш варіабельних CpG
# Використовуємо виключно training cohort

cat("\nCalculating CpG variance across all training samples...\n")

final_cpg_variance <- DelayedMatrixStats::rowVars(
  X_train,
  na.rm = TRUE
)

stopifnot(
  length(final_cpg_variance) == nrow(X_train),
  all(is.finite(final_cpg_variance))
)

final_top_cpg_idx <- order(
  final_cpg_variance,
  decreasing = TRUE
)[seq_len(5000)]

final_selected_cpg_ids <- rownames(X_train)[final_top_cpg_idx]

stopifnot(
  length(final_selected_cpg_ids) == 5000,
  !anyDuplicated(final_selected_cpg_ids),
  !anyNA(final_selected_cpg_ids)
)

cat("Selected CpG features:", length(final_selected_cpg_ids), "\n")

# 3. Формування training matrix

X_final_train <- t(
  as.matrix(
    X_train[final_top_cpg_idx, , drop = FALSE]
  )
)

y_final_train <- as.integer(
  y_train == "Immune-high"
)

stopifnot(
  identical(dim(X_final_train), c(368L, 5000L)),
  length(y_final_train) == 368,
  !anyNA(X_final_train),
  all(is.finite(X_final_train)),
  !anyNA(y_final_train)
)

cat(
  "Final training matrix:",
  paste(dim(X_final_train), collapse = " x "),
  "\n"
)

# 4. Навчання Elastic Net

set.seed(20261106)

final_enet_fit <- glmnet::glmnet(
  x = X_final_train,
  y = y_final_train,
  family = "binomial",
  alpha = 0.5,
  nlambda = 50,
  lambda.min.ratio = 0.01,
  standardize = TRUE
)

stopifnot(
  length(final_enet_fit$lambda) == 50,
  best_lambda_index <= length(final_enet_fit$lambda)
)

# 5. Вибір lambda за результатами 5-fold tuning

final_selected_lambda <- final_enet_fit$lambda[
  best_lambda_index
]

cat("\nFINAL MODEL PARAMETERS:\n")
cat("Alpha:", 0.5, "\n")
cat("Number of CpG:", length(final_selected_cpg_ids), "\n")
cat("Selected lambda index:", best_lambda_index, "\n")
cat("Selected lambda value:", final_selected_lambda, "\n")

# 6. Збереження фінальної моделі
# Holdout test поки не використовується

final_model <- list(
  model = final_enet_fit,
  selected_cpg_ids = final_selected_cpg_ids,
  selected_cpg_indices = final_top_cpg_idx,
  lambda_index = best_lambda_index,
  lambda = final_selected_lambda,
  alpha = 0.5,
  n_features = 5000,
  n_training_samples = 368,
  training_sample_ids = colnames(X_train),
  training_labels = y_train,
  tuning_mean_logloss = best_mean_logloss,
  tuning_sd_logloss = best_sd_logloss,
  tuning_results_file =
    "results/models/TCGA_SKCM_final_training_tuning.rds",
  seed = 20261106
)

final_model_file <- "results/models/TCGA_SKCM_final_ElasticNet_model.rds"

saveRDS(
  final_model,
  final_model_file
)

# 7. Перевірка збереженої моделі

final_model_check <- readRDS(final_model_file)

final_model_qc <- c(
  model_file_exists = file.exists(final_model_file),
  correct_training_samples =
    final_model_check$n_training_samples == 368,
  correct_feature_count =
    length(final_model_check$selected_cpg_ids) == 5000,
  correct_lambda_index =
    final_model_check$lambda_index == 34,
  lambda_is_finite =
    is.finite(final_model_check$lambda),
  model_has_coefficients =
    nrow(as.matrix(coef(
      final_model_check$model,
      s = final_model_check$lambda
    ))) == 5001,
  training_IDs_unique =
    !anyDuplicated(final_model_check$training_sample_ids)
)

cat("\nFINAL MODEL QC:\n")
print(final_model_qc)

stopifnot(all(final_model_qc))

cat("\nFINAL ELASTIC NET MODEL: PASS\n")
cat("HOLDOUT TEST USED: NO\n")


###############################################################################
# ЕТАП 35. HOLDOUT TEST PREPARATION AND QC
###############################################################################

# 1. Завантаження фінальної моделі та split manifest

final_model_file <- "results/models/TCGA_SKCM_final_ElasticNet_model.rds"
split_file <- "results/models/TCGA_SKCM_ML_train_test_split.rds"

stopifnot(
  file.exists(final_model_file),
  file.exists(split_file)
)

final_model_check <- readRDS(final_model_file)
split_manifest <- readRDS(split_file)

# 2. Перевірка розмірів holdout

stopifnot(
  ncol(X_test) == 92,
  length(y_test) == 92,
  length(split_manifest$test_sample_ids) == 92
)

cat("Holdout samples:", ncol(X_test), "\n")
cat("\nHoldout class distribution:\n")
print(table(y_test))

# 3. Перевірка sample IDs

test_ids_match <- identical(
  as.character(colnames(X_test)),
  as.character(split_manifest$test_sample_ids)
)

test_labels_match <- identical(
  as.character(y_test),
  as.character(split_manifest$test_labels)
)

# 4. Перевірка відсутності перетину з training

train_test_overlap <- intersect(
  final_model_check$training_sample_ids,
  colnames(X_test)
)

# 5. Перевірка відповідності CpG

test_cpg_match <- identical(
  as.character(rownames(X_test)),
  as.character(rownames(X_train))
)

selected_cpg_present <- all(
  final_model_check$selected_cpg_ids %in% rownames(X_test)
)

# 6. Формуємо матрицю holdout
# Використовуємо виключно CpG, вибрані фінальною моделлю

test_cpg_idx <- match(
  final_model_check$selected_cpg_ids,
  rownames(X_test)
)

stopifnot(!anyNA(test_cpg_idx))

X_holdout <- t(
  as.matrix(
    X_test[test_cpg_idx, , drop = FALSE]
  )
)

# 7. Перевірка матриці

holdout_matrix_qc <- c(
  correct_dimensions =
    identical(dim(X_holdout), c(92L, 5000L)),
  no_missing_values =
    !anyNA(X_holdout),
  all_finite =
    all(is.finite(X_holdout)),
  beta_in_valid_range =
    all(X_holdout >= 0 & X_holdout <= 1),
  sample_order_correct =
    identical(rownames(X_holdout), colnames(X_test))
)

# 8. Загальний QC

holdout_qc <- c(
  correct_sample_count =
    ncol(X_test) == 92,
  balanced_classes =
    all(table(y_test) == c(46, 46)),
  sample_IDs_match =
    test_ids_match,
  labels_match =
    test_labels_match,
  no_train_test_overlap =
    length(train_test_overlap) == 0,
  CpG_order_match =
    test_cpg_match,
  all_selected_CpG_present =
    selected_cpg_present,
  holdout_matrix = all(holdout_matrix_qc)
)

cat("\nHOLDOUT MATRIX QC:\n")
print(holdout_matrix_qc)

cat("\nHOLDOUT QC:\n")
print(holdout_qc)

stopifnot(all(holdout_matrix_qc))
stopifnot(all(holdout_qc))

cat("\nHOLDOUT TEST PREPARATION: PASS\n")
cat("MODEL PREDICTIONS GENERATED: NO\n")


###############################################################################
# ЕТАП 36. FINAL HOLDOUT TEST EVALUATION
###############################################################################

# 1. Завантаження фінальної моделі
final_model_check <- readRDS(
  "results/models/TCGA_SKCM_final_ElasticNet_model.rds"
)

stopifnot(
  nrow(X_holdout) == 92,
  ncol(X_holdout) == 5000,
  final_model_check$lambda_index == 34
)

# 2. Прогнозування ймовірностей Immune-high
holdout_prob <- predict(
  final_model_check$model,
  newx = X_holdout,
  s = final_model_check$lambda,
  type = "response"
)

holdout_prob <- as.numeric(holdout_prob)

stopifnot(
  length(holdout_prob) == 92,
  all(is.finite(holdout_prob)),
  all(holdout_prob >= 0 & holdout_prob <= 1)
)

# 3. Класифікація з порогом 0.5
classification_threshold <- 0.5

holdout_pred <- factor(
  ifelse(
    holdout_prob >= classification_threshold,
    "Immune-high",
    "Immune-low"
  ),
  levels = c("Immune-low", "Immune-high")
)

holdout_actual <- factor(
  as.character(y_test),
  levels = c("Immune-low", "Immune-high")
)

# 4. Confusion matrix
holdout_confusion <- table(
  Actual = holdout_actual,
  Predicted = holdout_pred
)

cat("\nHOLDOUT CONFUSION MATRIX:\n")
print(holdout_confusion)

# 5. Обчислення основних метрик
TN <- as.numeric(holdout_confusion[1, 1])
FP <- as.numeric(holdout_confusion[1, 2])
FN <- as.numeric(holdout_confusion[2, 1])
TP <- as.numeric(holdout_confusion[2, 2])

accuracy <- (TP + TN) / (TP + TN + FP + FN)

sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)

precision <- if ((TP + FP) > 0) {
  TP / (TP + FP)
} else {
  NA_real_
}

f1_score <- if (is.finite(precision) && (precision + sensitivity) > 0) {
  2 * precision * sensitivity / (precision + sensitivity)
} else {
  NA_real_
}

balanced_accuracy <- (sensitivity + specificity) / 2

# 6. ROC-AUC
holdout_roc <- pROC::roc(
  response = holdout_actual,
  predictor = holdout_prob,
  levels = c("Immune-low", "Immune-high"),
  direction = "<",
  quiet = TRUE
)

roc_auc <- as.numeric(
  pROC::auc(holdout_roc)
)

auc_ci <- as.numeric(
  pROC::ci.auc(holdout_roc)
)

# 7. Додаткові метрики
epsilon <- 1e-15

holdout_prob_clipped <- pmin(
  pmax(holdout_prob, epsilon),
  1 - epsilon
)

actual_numeric <- as.integer(
  holdout_actual == "Immune-high"
)

brier_score <- mean(
  (holdout_prob - actual_numeric)^2
)

holdout_logloss <- -mean(
  actual_numeric * log(holdout_prob_clipped) +
    (1 - actual_numeric) * log(1 - holdout_prob_clipped)
)

# 8. Формування таблиці метрик
holdout_metrics <- data.frame(
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
    roc_auc,
    auc_ci[1],
    auc_ci[3],
    brier_score,
    holdout_logloss
  )
)

# 9. Таблиця прогнозів
holdout_predictions <- data.frame(
  Sample_ID = names(y_test),
  Actual = as.character(holdout_actual),
  Predicted = as.character(holdout_pred),
  Probability_Immune_high = holdout_prob,
  Classification_threshold = classification_threshold
)

stopifnot(
  nrow(holdout_predictions) == 92,
  !anyDuplicated(holdout_predictions$Sample_ID),
  all(is.finite(holdout_metrics$Value))
)

# 10. Збереження результатів
write.csv(
  holdout_predictions,
  "results/tables/TCGA_SKCM_final_holdout_predictions.csv",
  row.names = FALSE
)

write.csv(
  holdout_metrics,
  "results/tables/TCGA_SKCM_final_holdout_metrics.csv",
  row.names = FALSE
)

saveRDS(
  list(
    predictions = holdout_predictions,
    confusion_matrix = holdout_confusion,
    metrics = holdout_metrics,
    roc = holdout_roc,
    roc_auc_ci = auc_ci,
    threshold = classification_threshold,
    model_file = "results/models/TCGA_SKCM_final_ElasticNet_model.rds",
    test_sample_count = 92
  ),
  "results/models/TCGA_SKCM_final_holdout_evaluation.rds"
)

# 11. Фінальний QC
holdout_evaluation_qc <- c(
  predictions_count = nrow(holdout_predictions) == 92,
  unique_sample_IDs = !anyDuplicated(holdout_predictions$Sample_ID),
  confusion_matrix_count = sum(holdout_confusion) == 92,
  metrics_finite = all(is.finite(holdout_metrics$Value)),
  auc_in_valid_range = roc_auc >= 0 && roc_auc <= 1,
  predictions_saved = file.exists(
    "results/tables/TCGA_SKCM_final_holdout_predictions.csv"
  ),
  metrics_saved = file.exists(
    "results/tables/TCGA_SKCM_final_holdout_metrics.csv"
  ),
  evaluation_saved = file.exists(
    "results/models/TCGA_SKCM_final_holdout_evaluation.rds"
  )
)

cat("\nHOLDOUT METRICS:\n")
print(holdout_metrics, row.names = FALSE)

cat("\nHOLDOUT EVALUATION QC:\n")
print(holdout_evaluation_qc)

stopifnot(all(holdout_evaluation_qc))

cat("\nFINAL HOLDOUT EVALUATION: PASS\n")





