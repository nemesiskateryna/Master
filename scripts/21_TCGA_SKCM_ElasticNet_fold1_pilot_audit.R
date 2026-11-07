

###############################################################################
#  ML INPUT DATA QC
###############################################################################

# Завантажуємо фінальний об'єкт
ml_source <- readRDS(
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune.rds"
)

# Перевіряємо основні компоненти
cat("--- COMPONENTS ---\n")
print(names(ml_source))

cat("\n--- TARGET DISTRIBUTION ---\n")
print(table(ml_source$immune_scores$immune_group))

# Перевіряємо матрицю метилювання
cat("\n--- METHYLATION MATRIX ---\n")
print(dim(ml_source$beta_matrix))

# Перевіряємо кількість доступних CpG та зразків
cat(
  "\nUnique CpG:",
  length(unique(rownames(ml_source$beta_matrix))),
  "\n"
)

cat(
  "Unique samples:",
  length(unique(colnames(ml_source$beta_matrix))),
  "\n"
)

# Перевіряємо відповідність цільової змінної
cat(
  "Target/sample alignment:",
  identical(
    ml_source$immune_scores$Sample_ID,
    colnames(ml_source$beta_matrix)
  ),
  "\n"
)

# Перевіряємо доступність HDF5
test_block <- as.matrix(
  ml_source$beta_matrix[1:10, 1:10]
)

cat(
  "HDF5 accessible:",
  all(is.finite(test_block)),
  "\n"
)

cat(
  "Beta values in [0,1]:",
  all(test_block >= 0 & test_block <= 1),
  "\n"
)


###############################################################################
# STEP 3. DEFINE ML TARGET AND PREDICTOR MATRIX
###############################################################################

# 1. Цільова змінна
y <- factor(
  ml_source$immune_scores$immune_group,
  levels = c("Immune-low", "Immune-high")
)

names(y) <- ml_source$immune_scores$Sample_ID

# 2. Матриця предикторів — DNA methylation
# Залишаємо DelayedMatrix/HDF5, не перетворюємо всю матрицю на matrix
X <- ml_source$beta_matrix

# 3. Перевірка відповідності
stopifnot(
  identical(names(y), colnames(X)),
  nlevels(y) == 2,
  !anyNA(y)
)

# 4. Перевірка розмірів
cat("--- ML DATASET ---\n")

cat("Predictors (CpG):", nrow(X), "\n")
cat("Samples:", ncol(X), "\n")

cat("\nTarget distribution:\n")
print(table(y))

cat("\nTarget proportions:\n")
print(round(prop.table(table(y)), 3))

cat("\nClass levels:\n")
print(levels(y))

cat("\nSample alignment:",
    identical(names(y), colnames(X)), "\n")

cat("\nPredictor matrix class:\n")
print(class(X))

cat("\nML INPUT QC:",
    if (
      nrow(X) == 409556 &&
      ncol(X) == 460 &&
      identical(names(y), colnames(X)) &&
      !anyNA(y)
    ) "PASS" else "CHECK",
    "\n"
)

###############################################################################
# STEP 4. PATIENT-LEVEL DUPLICATE AUDIT
###############################################################################

# 1. Отримуємо patient ID з TCGA barcode
sample_ids <- colnames(X)

patient_ids <- substr(sample_ids, 1, 12)

# 2. Формуємо таблицю зразків
sample_patient_table <- data.frame(
  Sample_ID = sample_ids,
  Patient_ID = patient_ids,
  Immune_group = as.character(y),
  stringsAsFactors = FALSE
)

# 3. Групуємо зразки за пацієнтами
patient_groups <- split(
  seq_along(patient_ids),
  patient_ids
)

# 4. Підраховуємо кількість зразків і класів для кожного пацієнта
patient_audit <- data.frame(
  Patient_ID = names(patient_groups),
  
  Number_of_samples = vapply(
    patient_groups,
    length,
    integer(1)
  ),
  
  Number_of_immune_labels = vapply(
    patient_groups,
    function(i) length(unique(y[i])),
    integer(1)
  ),
  
  stringsAsFactors = FALSE
)

# 5. Загальні результати
cat("--- PATIENT-LEVEL AUDIT ---\n")

cat("Total samples:",
    nrow(sample_patient_table), "\n")

cat("Unique patients:",
    nrow(patient_audit), "\n")

cat("Patients with multiple samples:",
    sum(patient_audit$Number_of_samples > 1), "\n")

cat("Maximum samples per patient:",
    max(patient_audit$Number_of_samples), "\n")

cat("\nSamples per patient:\n")
print(table(patient_audit$Number_of_samples))

cat("\nPatients with inconsistent immune labels:\n")
print(
  sum(patient_audit$Number_of_immune_labels > 1)
)

# 6. Переглядаємо пацієнтів із повторними зразками
cat("\n--- PATIENTS WITH MULTIPLE SAMPLES ---\n")

print(
  patient_audit[
    patient_audit$Number_of_samples > 1,
  ],
  row.names = FALSE
)

# 7. Перевірка формату patient ID
cat("\nAll patient IDs have 12 characters:",
    all(nchar(patient_ids) == 12), "\n")

cat("Patient-level audit completed.\n")


###############################################################################
# STEP 5. STRATIFIED TRAIN / TEST SPLIT
###############################################################################

# 1. Фіксуємо випадковий seed
set.seed(2026)

# 2. Створюємо стратифікований training set:
# 80% зразків окремо в кожному класі
train_idx <- unlist(
  lapply(levels(y), function(class_name) {
    
    class_idx <- which(y == class_name)
    
    sample(
      class_idx,
      size = floor(0.8 * length(class_idx))
    )
  }),
  use.names = FALSE
)

# 3. Сортуємо індекси для збереження початкового порядку
train_idx <- sort(train_idx)

# 4. Решта зразків — test set
test_idx <- setdiff(seq_len(ncol(X)), train_idx)

# 5. Формуємо цільові змінні
y_train <- droplevels(y[train_idx])
y_test  <- droplevels(y[test_idx])

# 6. Формуємо відповідні матриці
# Залишаємо DelayedMatrix — без повного завантаження в RAM
X_train <- X[, train_idx, drop = FALSE]
X_test  <- X[, test_idx, drop = FALSE]

# 7. Перевіряємо розміри
cat("--- TRAIN / TEST DIMENSIONS ---\n")

cat("Training samples:", ncol(X_train), "\n")
cat("Test samples:", ncol(X_test), "\n")

cat("Training CpG:", nrow(X_train), "\n")
cat("Test CpG:", nrow(X_test), "\n")

# 8. Перевіряємо баланс класів
cat("\n--- TRAINING CLASS DISTRIBUTION ---\n")
print(table(y_train))
print(round(prop.table(table(y_train)), 3))

cat("\n--- TEST CLASS DISTRIBUTION ---\n")
print(table(y_test))
print(round(prop.table(table(y_test)), 3))

# 9. Перевіряємо відсутність перетину пацієнтів
train_patients <- patient_ids[train_idx]
test_patients  <- patient_ids[test_idx]

patient_overlap <- intersect(
  train_patients,
  test_patients
)

cat("\nPatient overlap:", length(patient_overlap), "\n")

# 10. Перевіряємо відповідність sample IDs
cat(
  "Training alignment:",
  identical(names(y_train), colnames(X_train)),
  "\n"
)

cat(
  "Test alignment:",
  identical(names(y_test), colnames(X_test)),
  "\n"
)

# 11. Підсумковий QC
split_qc <- (
  ncol(X_train) == 368 &&
    ncol(X_test) == 92 &&
    nrow(X_train) == 409556 &&
    nrow(X_test) == 409556 &&
    length(patient_overlap) == 0 &&
    identical(names(y_train), colnames(X_train)) &&
    identical(names(y_test), colnames(X_test)) &&
    all(table(y_train) == c(184, 184)) &&
    all(table(y_test) == c(46, 46))
)

cat(
  "\nTRAIN / TEST SPLIT QC:",
  if (split_qc) "PASS" else "CHECK",
  "\n"
)
###############################################################################
# STEP 6. SAVE TRAIN / TEST SPLIT AND CHECK ML PACKAGES
###############################################################################

# 1. Створюємо таблицю розподілу зразків
split_table <- data.frame(
  Sample_ID = colnames(X),
  Patient_ID = patient_ids,
  Immune_group = as.character(y),
  Set = ifelse(
    seq_len(ncol(X)) %in% train_idx,
    "Training",
    "Test"
  ),
  stringsAsFactors = FALSE
)

# 2. Створюємо об'єкт із параметрами поділу
split_manifest <- list(
  seed = 2026,
  split_ratio = "80/20",
  train_indices = train_idx,
  test_indices = test_idx,
  train_sample_ids = colnames(X_train),
  test_sample_ids = colnames(X_test),
  train_labels = y_train,
  test_labels = y_test,
  patient_overlap = patient_overlap,
  split_table = split_table
)

# 3. Перевіряємо каталоги
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

# 4. Зберігаємо маніфест
saveRDS(
  split_manifest,
  "results/models/TCGA_SKCM_ML_train_test_split.rds"
)

# 5. Зберігаємо таблицю для аудиту
write.csv(
  split_table,
  "results/tables/TCGA_SKCM_ML_train_test_split.csv",
  row.names = FALSE
)

# 6. Перевіряємо наявність потрібних пакетів
required_packages <- c(
  "DelayedArray",
  "DelayedMatrixStats",
  "matrixStats",
  "glmnet",
  "caret",
  "pROC"
)

package_status <- data.frame(
  Package = required_packages,
  Installed = vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
)

cat("--- SAVED FILES ---\n")

cat(
  "Manifest exists:",
  file.exists(
    "results/models/TCGA_SKCM_ML_train_test_split.rds"
  ),
  "\n"
)

cat(
  "CSV exists:",
  file.exists(
    "results/tables/TCGA_SKCM_ML_train_test_split.csv"
  ),
  "\n"
)

cat("\n--- PACKAGE AVAILABILITY ---\n")
print(package_status, row.names = FALSE)

cat(
  "\nSPLIT MANIFEST QC:",
  if (
    file.exists("results/models/TCGA_SKCM_ML_train_test_split.rds") &&
    file.exists("results/tables/TCGA_SKCM_ML_train_test_split.csv")
  ) "PASS" else "CHECK",
  "\n"
)
###############################################################################
# STEP 7. CREATE STRATIFIED 5-FOLD CROSS-VALIDATION
###############################################################################

# 1. Фіксуємо seed
set.seed(2026)

# 2. Кількість folds
k_folds <- 5

# 3. Створюємо список індексів validation для кожного fold
cv_folds <- vector("list", k_folds)

# 4. Окремо розподіляємо кожен клас
for (class_name in levels(y_train)) {
  
  class_indices <- which(y_train == class_name)
  
  assigned_folds <- sample(
    rep(seq_len(k_folds), length.out = length(class_indices))
  )
  
  for (fold_number in seq_len(k_folds)) {
    
    cv_folds[[fold_number]] <- c(
      cv_folds[[fold_number]],
      class_indices[assigned_folds == fold_number]
    )
  }
}

# 5. Сортуємо індекси
cv_folds <- lapply(cv_folds, sort)

# 6. Перевіряємо розподіл класів у кожному fold
cat("--- CROSS-VALIDATION FOLD DISTRIBUTION ---\n")

cv_distribution <- do.call(
  rbind,
  lapply(seq_len(k_folds), function(i) {
    
    fold_labels <- y_train[cv_folds[[i]]]
    
    data.frame(
      Fold = i,
      Samples = length(cv_folds[[i]]),
      Immune_low = sum(fold_labels == "Immune-low"),
      Immune_high = sum(fold_labels == "Immune-high")
    )
  })
)

print(cv_distribution, row.names = FALSE)

# 7. Перевіряємо, що кожен training sample потрапив
# рівно в один validation fold
all_validation_indices <- unlist(cv_folds)

cat(
  "\nTotal validation indices:",
  length(all_validation_indices),
  "\n"
)

cat(
  "Unique validation indices:",
  length(unique(all_validation_indices)),
  "\n"
)

cat(
  "All training samples assigned:",
  setequal(
    all_validation_indices,
    seq_len(nrow(cv_distribution) * 0 + length(y_train))
  ),
  "\n"
)

# 8. Перевіряємо відсутність повторів
cat(
  "No duplicated indices:",
  !anyDuplicated(all_validation_indices),
  "\n"
)

# 9. Зберігаємо CV folds
cv_manifest <- list(
  seed = 2026,
  k_folds = k_folds,
  cv_folds = cv_folds,
  cv_distribution = cv_distribution,
  training_sample_ids = names(y_train)
)

saveRDS(
  cv_manifest,
  "results/models/TCGA_SKCM_ML_5fold_CV_manifest.rds"
)

write.csv(
  cv_distribution,
  "results/tables/TCGA_SKCM_ML_5fold_CV_distribution.csv",
  row.names = FALSE
)

# 10. Фінальний QC
cv_qc <- (
  length(all_validation_indices) == length(y_train) &&
    length(unique(all_validation_indices)) == length(y_train) &&
    !anyDuplicated(all_validation_indices) &&
    all(cv_distribution$Samples > 0) &&
    all(cv_distribution$Immune_low > 0) &&
    all(cv_distribution$Immune_high > 0)
)

cat(
  "\n5-FOLD CV QC:",
  if (cv_qc) "PASS" else "CHECK",
  "\n"
)
###############################################################################
# STEP 8. PILOT FEATURE FILTERING — FOLD 1
###############################################################################

# 1. Завантажуємо збережений CV manifest
cv_manifest <- readRDS(
  "results/models/TCGA_SKCM_ML_5fold_CV_manifest.rds"
)

# 2. Визначаємо validation та inner-training індекси
fold_number <- 1

fold_valid_idx <- cv_manifest$cv_folds[[fold_number]]

fold_train_idx <- setdiff(
  seq_along(y_train),
  fold_valid_idx
)

# 3. Формуємо матрицю тільки для inner-training
X_fold_train <- X_train[
  ,
  fold_train_idx,
  drop = FALSE
]

y_fold_train <- y_train[fold_train_idx]

# 4. Обчислюємо дисперсію CpG
# DelayedMatrix дозволяє обробляти дані без повного
# завантаження вихідної матриці в оперативну пам'ять
cpg_variance <- DelayedMatrixStats::rowVars(
  X_fold_train,
  na.rm = TRUE
)

# 5. Перевіряємо результати обчислення
cat("--- FOLD 1 FEATURE FILTERING ---\n")

cat("Inner-training samples:", ncol(X_fold_train), "\n")
cat("Validation samples:", length(fold_valid_idx), "\n")

cat("Total CpG:", length(cpg_variance), "\n")

cat("NA variances:", sum(is.na(cpg_variance)), "\n")
cat("Non-finite variances:",
    sum(!is.finite(cpg_variance)), "\n")

cat("Zero-variance CpG:",
    sum(cpg_variance == 0, na.rm = TRUE), "\n")

cat("\nVariance summary:\n")
print(summary(cpg_variance))

# 6. Відбираємо 5000 найбільш варіабельних CpG
valid_variance_idx <- which(
  is.finite(cpg_variance) & cpg_variance > 0
)

if (length(valid_variance_idx) < 5000) {
  stop("Fewer than 5000 CpG have valid non-zero variance.")
}

top_cpg_idx <- valid_variance_idx[
  order(
    cpg_variance[valid_variance_idx],
    decreasing = TRUE
  )[seq_len(5000)]
]

top_cpg_ids <- rownames(X_fold_train)[top_cpg_idx]

# 7. Перевіряємо відбір
cat("\nSelected CpG:", length(top_cpg_idx), "\n")

cat("Unique selected CpG:",
    length(unique(top_cpg_ids)), "\n")

cat("Selected CpG IDs available:",
    !anyNA(top_cpg_ids), "\n")

# 8. QC
feature_filter_qc <- (
  length(top_cpg_idx) == 5000 &&
    length(unique(top_cpg_ids)) == 5000 &&
    !anyNA(top_cpg_ids) &&
    all(is.finite(cpg_variance[top_cpg_idx]))
)

cat(
  "\nFOLD 1 VARIANCE FILTER QC:",
  if (feature_filter_qc) "PASS" else "CHECK",
  "\n"
)
###############################################################################
# STEP 9. PREPARE FOLD 1 MATRICES FOR GLMNET
###############################################################################

# 1. Отримуємо матриці 5000 CpG
# Inner-training: 294 samples
# Validation: 74 samples

X_fold_train_selected <- as.matrix(
  X_train[
    top_cpg_idx,
    fold_train_idx,
    drop = FALSE
  ]
)

X_fold_valid_selected <- as.matrix(
  X_train[
    top_cpg_idx,
    fold_valid_idx,
    drop = FALSE
  ]
)

# 2. Транспонуємо:
# rows = samples
# columns = CpG features

X_fold_train_glmnet <- t(X_fold_train_selected)
X_fold_valid_glmnet <- t(X_fold_valid_selected)

# 3. Формуємо числову цільову змінну
# Immune-low = 0
# Immune-high = 1

y_fold_train_numeric <- as.integer(
  y_fold_train == "Immune-high"
)

y_fold_valid_numeric <- as.integer(
  y_train[fold_valid_idx] == "Immune-high"
)

# 4. Перевіряємо розміри
cat("--- GLMNET INPUT DIMENSIONS ---\n")

cat(
  "Training matrix:",
  paste(dim(X_fold_train_glmnet), collapse = " x "),
  "\n"
)

cat(
  "Validation matrix:",
  paste(dim(X_fold_valid_glmnet), collapse = " x "),
  "\n"
)

cat(
  "Training labels:",
  length(y_fold_train_numeric),
  "\n"
)

cat(
  "Validation labels:",
  length(y_fold_valid_numeric),
  "\n"
)

# 5. Перевіряємо відсутність пропусків і нескінченних значень
cat("\n--- NUMERICAL QC ---\n")

cat(
  "Training NA:",
  sum(is.na(X_fold_train_glmnet)),
  "\n"
)

cat(
  "Validation NA:",
  sum(is.na(X_fold_valid_glmnet)),
  "\n"
)

cat(
  "Training non-finite:",
  sum(!is.finite(X_fold_train_glmnet)),
  "\n"
)

cat(
  "Validation non-finite:",
  sum(!is.finite(X_fold_valid_glmnet)),
  "\n"
)

cat(
  "Training beta range:",
  paste(range(X_fold_train_glmnet), collapse = " - "),
  "\n"
)

cat(
  "Validation beta range:",
  paste(range(X_fold_valid_glmnet), collapse = " - "),
  "\n"
)

# 6. Перевіряємо відповідність CpG і зразків
cat("\n--- ALIGNMENT QC ---\n")

cat(
  "Training CpG alignment:",
  identical(
    colnames(X_fold_train_glmnet),
    top_cpg_ids
  ),
  "\n"
)

cat(
  "Validation CpG alignment:",
  identical(
    colnames(X_fold_valid_glmnet),
    top_cpg_ids
  ),
  "\n"
)

cat(
  "Training sample alignment:",
  identical(
    rownames(X_fold_train_glmnet),
    names(y_fold_train)
  ),
  "\n"
)

cat(
  "Validation sample alignment:",
  identical(
    rownames(X_fold_valid_glmnet),
    names(y_train[fold_valid_idx])
  ),
  "\n"
)

# 7. Фінальна перевірка
glmnet_input_qc <- (
  identical(dim(X_fold_train_glmnet), c(294L, 5000L)) &&
    identical(dim(X_fold_valid_glmnet), c(74L, 5000L)) &&
    length(y_fold_train_numeric) == 294 &&
    length(y_fold_valid_numeric) == 74 &&
    all(is.finite(X_fold_train_glmnet)) &&
    all(is.finite(X_fold_valid_glmnet)) &&
    all(X_fold_train_glmnet >= 0 & X_fold_train_glmnet <= 1) &&
    all(X_fold_valid_glmnet >= 0 & X_fold_valid_glmnet <= 1) &&
    identical(colnames(X_fold_train_glmnet), top_cpg_ids) &&
    identical(colnames(X_fold_valid_glmnet), top_cpg_ids) &&
    identical(rownames(X_fold_train_glmnet), names(y_fold_train)) &&
    identical(rownames(X_fold_valid_glmnet), names(y_train[fold_valid_idx]))
)

cat(
  "\nGLMNET INPUT QC:",
  if (glmnet_input_qc) "PASS" else "CHECK",
  "\n"
)
###############################################################################
# STEP 10. ELASTIC NET — FOLD 1 PILOT MODEL
###############################################################################

# 1. Перевірка пакета
stopifnot(requireNamespace("glmnet", quietly = TRUE))

# 2. Формуємо матрицю для glmnet
# glmnet самостійно стандартизує ознаки за training data
set.seed(2026)

# 3. Внутрішня крос-валідація для підбору lambda
cv_enet <- glmnet::cv.glmnet(
  x = X_fold_train_glmnet,
  y = y_fold_train_numeric,
  family = "binomial",
  alpha = 0.5,
  nfolds = 5,
  type.measure = "deviance",
  standardize = TRUE,
  keep = TRUE
)

# 4. Оптимальні параметри
cat("--- ELASTIC NET PARAMETERS ---\n")

cat("Alpha:", 0.5, "\n")
cat("Lambda.min:", cv_enet$lambda.min, "\n")
cat("Lambda.1se:", cv_enet$lambda.1se, "\n")

# 5. Кількість відібраних ознак
coef_min <- as.matrix(
  coef(cv_enet, s = "lambda.min")
)

coef_1se <- as.matrix(
  coef(cv_enet, s = "lambda.1se")
)

selected_min <- rownames(coef_min)[
  coef_min[, 1] != 0
]

selected_1se <- rownames(coef_1se)[
  coef_1se[, 1] != 0
]

selected_min <- setdiff(selected_min, "(Intercept)")
selected_1se <- setdiff(selected_1se, "(Intercept)")

cat("\nSelected CpG at lambda.min:",
    length(selected_min), "\n")

cat("Selected CpG at lambda.1se:",
    length(selected_1se), "\n")

# 6. Прогноз для validation set
pred_prob <- as.numeric(
  predict(
    cv_enet,
    newx = X_fold_valid_glmnet,
    s = "lambda.1se",
    type = "response"
  )
)

pred_class <- ifelse(
  pred_prob >= 0.5,
  1,
  0
)

# 7. Оцінка класифікації
confusion_matrix <- table(
  Actual = y_fold_valid_numeric,
  Predicted = pred_class
)

cat("\n--- VALIDATION CONFUSION MATRIX ---\n")
print(confusion_matrix)

accuracy <- mean(
  pred_class == y_fold_valid_numeric
)

cat("\nValidation accuracy:",
    round(accuracy, 4), "\n")

# 8. Перевірка прогнозів
cat("\n--- PREDICTION QC ---\n")

cat("Predictions:", length(pred_prob), "\n")
cat("NA predictions:", sum(is.na(pred_prob)), "\n")
cat("Probability range:",
    paste(range(pred_prob), collapse = " - "), "\n")

cat(
  "\nELASTIC NET PILOT QC:",
  if (
    length(pred_prob) == length(y_fold_valid_numeric) &&
    all(is.finite(pred_prob)) &&
    all(pred_prob >= 0 & pred_prob <= 1)
  ) "PASS" else "CHECK",
  "\n"
)

###############################################################################
# STEP 11. FOLD 1 METRICS AND PILOT RESULTS
###############################################################################

# 1. Елементи confusion matrix
tn <- sum(y_fold_valid_numeric == 0 & pred_class == 0)
fp <- sum(y_fold_valid_numeric == 0 & pred_class == 1)
fn <- sum(y_fold_valid_numeric == 1 & pred_class == 0)
tp <- sum(y_fold_valid_numeric == 1 & pred_class == 1)

# 2. Основні метрики
sensitivity <- tp / (tp + fn)
specificity <- tn / (tn + fp)
precision <- tp / (tp + fp)

accuracy <- (tp + tn) / (tp + tn + fp + fn)
balanced_accuracy <- (sensitivity + specificity) / 2

f1_score <- 2 * (
  precision * sensitivity
) / (
  precision + sensitivity
)

# 3. AUC за ранговою статистикою
positive_scores <- pred_prob[y_fold_valid_numeric == 1]
negative_scores <- pred_prob[y_fold_valid_numeric == 0]

auc_value <- (
  sum(rank(c(positive_scores, negative_scores))[seq_along(positive_scores)]) -
    length(positive_scores) * (length(positive_scores) + 1) / 2
) / (
  length(positive_scores) * length(negative_scores)
)

# 4. Формуємо таблицю метрик
pilot_metrics <- data.frame(
  Fold = 1,
  Model = "Elastic Net",
  Alpha = 0.5,
  Lambda_1se = cv_enet$lambda.1se,
  Selected_CpG = length(selected_1se),
  Accuracy = accuracy,
  Balanced_accuracy = balanced_accuracy,
  Sensitivity = sensitivity,
  Specificity = specificity,
  Precision = precision,
  F1 = f1_score,
  AUC = auc_value
)

cat("--- FOLD 1 PILOT METRICS ---\n")
print(pilot_metrics, row.names = FALSE)

# 5. Таблиця прогнозів
pilot_predictions <- data.frame(
  Sample_ID = rownames(X_fold_valid_glmnet),
  Actual = y_fold_valid_numeric,
  Actual_group = as.character(y_train[fold_valid_idx]),
  Predicted = pred_class,
  Probability_Immune_high = pred_prob,
  stringsAsFactors = FALSE
)

# 6. Зберігаємо пілотні результати
saveRDS(
  list(
    model = cv_enet,
    metrics = pilot_metrics,
    predictions = pilot_predictions,
    selected_CpG_lambda_min = selected_min,
    selected_CpG_lambda_1se = selected_1se,
    top_5000_CpG = top_cpg_ids,
    confusion_matrix = confusion_matrix,
    fold_number = 1,
    seed = 2026
  ),
  "results/models/TCGA_SKCM_ElasticNet_fold1_pilot.rds"
)

write.csv(
  pilot_metrics,
  "results/tables/TCGA_SKCM_ElasticNet_fold1_pilot_metrics.csv",
  row.names = FALSE
)

write.csv(
  pilot_predictions,
  "results/tables/TCGA_SKCM_ElasticNet_fold1_pilot_predictions.csv",
  row.names = FALSE
)

# 7. QC
cat("\n--- SAVING QC ---\n")

cat(
  "Pilot model saved:",
  file.exists("results/models/TCGA_SKCM_ElasticNet_fold1_pilot.rds"),
  "\n"
)

cat(
  "Metrics saved:",
  file.exists("results/tables/TCGA_SKCM_ElasticNet_fold1_pilot_metrics.csv"),
  "\n"
)

cat(
  "Predictions saved:",
  file.exists("results/tables/TCGA_SKCM_ElasticNet_fold1_pilot_predictions.csv"),
  "\n"
)















