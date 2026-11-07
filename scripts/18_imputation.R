
###############################################################################
# SCRIPT 15
# TCGA-SKCM: KNN IMPUTATION AND FINAL QUALITY CONTROL
#
# Мета:
# 1. Зафіксувати результати порівняння Median та KNN.
# 2. Виконати blockwise KNN-імпутацію за потреби.
# 3. Перевірити фінальну матрицю.
# 4. Зберегти QC-звіт та метадані.
#
# Вхідна матриця:
# TCGA_SKCM_beta_matrix_460_StandardCpG_20pct.h5
#
# Вихідна матриця:
# TCGA_SKCM_beta_matrix_460_StandardCpG_20pct_KNN.h5
###############################################################################

###############################################################################
# 1. ПАКЕТИ ТА ШЛЯХИ
###############################################################################

library(rhdf5)
library(impute)

input_h5 <- file.path(
  "data", "TCGA", "processed",
  "TCGA_SKCM_beta_matrix_460_StandardCpG_20pct.h5"
)

output_h5 <- file.path(
  "data", "TCGA", "processed",
  "TCGA_SKCM_beta_matrix_460_StandardCpG_20pct_KNN.h5"
)

metadata_path <- file.path(
  "data", "TCGA", "processed",
  "TCGA_SKCM_KNN_row_metadata.rds"
)

results_dir <- file.path("results", "tables")

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(file.exists(input_h5))

###############################################################################
# 2. РЕЗУЛЬТАТИ ПОРІВНЯННЯ MEDIAN ТА KNN
###############################################################################

# Результати п'яти повторів на стратифікованому тестовому наборі.
# Значення відповідають раніше проведеному експерименту.

benchmark_results <- data.frame(
  Method = c("Median", "KNN"),
  
  RMSE_mean = c(
    0.182542306,
    0.144997131
  ),
  
  RMSE_SD = c(
    0.003766858,
    0.002920413
  ),
  
  MAE_mean = c(
    0.111382463,
    0.089038457
  ),
  
  MAE_SD = c(
    0.002143364,
    0.001118388
  ),
  
  Correlation_mean = c(
    0.803447324,
    0.876415205
  ),
  
  Correlation_SD = c(
    0.008834040,
    0.006206067
  )
)

cat("\n--- ПОРІВНЯННЯ МЕТОДІВ ІМПУТАЦІЇ ---\n")
print(benchmark_results)

###############################################################################
# 3. ПАРАМЕТРИ KNN
###############################################################################

k_neighbors <- 10
block_size <- 1000

rowmax_value <- 0.5
colmax_value <- 0.8

###############################################################################
# 4. BLOCKWISE KNN-ІМПУТАЦІЯ
###############################################################################

if (!file.exists(output_h5)) {
  
  cat("\n--- ПОЧАТОК KNN-ІМПУТАЦІЇ ---\n")
  
  # Зчитування ідентифікаторів
  cpg_ids <- h5read(input_h5, "CpG_ID")
  sample_ids <- h5read(input_h5, "Sample_Barcode")
  
  n_cpg <- length(cpg_ids)
  n_samples <- length(sample_ids)
  
  # Кількість блоків
  block_starts <- seq(
    1,
    n_cpg,
    by = block_size
  )
  
  n_blocks <- length(block_starts)
  
  # Створення тимчасового HDF5-файлу
  partial_h5 <- paste0(output_h5, "_partial.h5")
  
  if (file.exists(partial_h5)) {
    stop(
      "Знайдено попередній partial-файл. ",
      "Перевірте його вручну перед продовженням."
    )
  }
  
  h5createFile(partial_h5)
  
  # Створення beta dataset
  h5createDataset(
    file = partial_h5,
    dataset = "beta",
    dims = c(n_cpg, n_samples),
    storage.mode = "double",
    chunk = c(block_size, n_samples),
    level = 6
  )
  
  # Перенесення ідентифікаторів
  h5write(cpg_ids, partial_h5, "CpG_ID")
  h5write(sample_ids, partial_h5, "Sample_Barcode")
  
  total_imputed <- 0
  
  for (i in seq_along(block_starts)) {
    
    start_row <- block_starts[i]
    
    end_row <- min(
      start_row + block_size - 1,
      n_cpg
    )
    
    rows <- start_row:end_row
    
    # Зчитування одного блоку
    beta_block <- h5read(
      input_h5,
      "beta",
      index = list(rows, seq_len(n_samples))
    )
    
    beta_block <- as.matrix(beta_block)
    
    # Кількість пропусків до імпутації
    missing_before <- sum(is.na(beta_block))
    
    # KNN-імпутація
    if (missing_before > 0) {
      
      knn_result <- impute::impute.knn(
        data = beta_block,
        k = k_neighbors,
        rowmax = rowmax_value,
        colmax = colmax_value,
        maxp = block_size
      )
      
      beta_block <- knn_result$data
      
    }
    
    # Контроль результату блоку
    if (anyNA(beta_block)) {
      stop(
        "Залишилися NA у блоці ",
        i
      )
    }
    
    if (any(!is.finite(beta_block))) {
      stop(
        "Виявлено не скінченні значення у блоці ",
        i
      )
    }
    
    if (any(beta_block < 0 | beta_block > 1)) {
      stop(
        "Beta поза діапазоном [0,1] у блоці ",
        i
      )
    }
    
    # Запис імпутованого блоку
    h5write(
      beta_block,
      partial_h5,
      "beta",
      index = list(rows, seq_len(n_samples))
    )
    
    total_imputed <- total_imputed + missing_before
    
    cat(
      "Блок", i, "/", n_blocks,
      "| CpG:", end_row,
      "| Заповнено:", missing_before,
      "\n"
    )
    
    rm(beta_block)
    
    gc()
  }
  
  # Закриття HDF5 та перенесення до фінального файлу
  H5close()
  
  if (file.exists(output_h5)) {
    stop("Фінальний файл уже існує. Перезапис заборонено.")
  }
  
  file.rename(
    partial_h5,
    output_h5
  )
  
  if (!file.exists(output_h5)) {
    stop("Не вдалося створити фінальний HDF5-файл.")
  }
  
  # Метадані імпутації
  knn_metadata <- list(
    Dataset = "TCGA-SKCM",
    Method = "Blockwise KNN",
    K = k_neighbors,
    Block_size = block_size,
    Rowmax = rowmax_value,
    Colmax = colmax_value,
    CpG_count = n_cpg,
    Sample_count = n_samples,
    Number_of_blocks = n_blocks,
    Imputed_values = total_imputed,
    Date = Sys.time()
  )
  
  if (!file.exists(metadata_path)) {
    saveRDS(knn_metadata, metadata_path)
  }
  
  cat("\nKNN-ІМПУТАЦІЮ ЗАВЕРШЕНО\n")
  
} else {
  
  cat(
    "\nФінальна KNN-матриця вже існує. ",
    "Повторна імпутація не виконується.\n"
  )
}

###############################################################################
# 5. ФІНАЛЬНИЙ QC
###############################################################################

original_cpg <- h5read(input_h5, "CpG_ID")
imputed_cpg <- h5read(output_h5, "CpG_ID")

original_samples <- h5read(input_h5, "Sample_Barcode")
imputed_samples <- h5read(output_h5, "Sample_Barcode")

qc_cpg_order <- identical(
  original_cpg,
  imputed_cpg
)

qc_sample_order <- identical(
  original_samples,
  imputed_samples
)

n_cpg <- length(original_cpg)
n_samples <- length(original_samples)

qc_na <- 0
qc_nonfinite <- 0
qc_out_of_range <- 0
qc_changed_observed <- 0

qc_min <- Inf
qc_max <- -Inf

for (start_row in seq(1, n_cpg, by = 5000)) {
  
  end_row <- min(
    start_row + 4999,
    n_cpg
  )
  
  rows <- start_row:end_row
  
  original_block <- h5read(
    input_h5,
    "beta",
    index = list(rows, seq_len(n_samples))
  )
  
  imputed_block <- h5read(
    output_h5,
    "beta",
    index = list(rows, seq_len(n_samples))
  )
  
  qc_na <- qc_na + sum(is.na(imputed_block))
  
  qc_nonfinite <- qc_nonfinite +
    sum(!is.finite(imputed_block))
  
  qc_out_of_range <- qc_out_of_range +
    sum(
      imputed_block < 0 | imputed_block > 1,
      na.rm = TRUE
    )
  
  observed <- !is.na(original_block)
  
  qc_changed_observed <- qc_changed_observed +
    sum(
      imputed_block[observed] !=
        original_block[observed],
      na.rm = TRUE
    )
  
  qc_min <- min(
    qc_min,
    imputed_block,
    na.rm = TRUE
  )
  
  qc_max <- max(
    qc_max,
    imputed_block,
    na.rm = TRUE
  )
  
  rm(original_block, imputed_block, observed)
  
  cat(
    "QC:", end_row, "/", n_cpg, "CpG\n"
  )
}

###############################################################################
# 6. ПІДСУМКОВИЙ QC-СТАТУС
###############################################################################

qc_pass <- (
  qc_cpg_order &&
    qc_sample_order &&
    qc_na == 0 &&
    qc_nonfinite == 0 &&
    qc_out_of_range == 0 &&
    qc_changed_observed == 0
)

qc_report <- data.frame(
  Parameter = c(
    "Dataset",
    "Method",
    "K",
    "Block_size",
    "CpG_count",
    "Sample_count",
    "CpG_order_identical",
    "Sample_order_identical",
    "Residual_NA",
    "Non_finite_values",
    "Beta_outside_0_1",
    "Changed_observed_values",
    "Minimum_beta",
    "Maximum_beta",
    "Final_QC_status"
  ),
  
  Value = as.character(c(
    "TCGA-SKCM",
    "Blockwise KNN",
    k_neighbors,
    block_size,
    n_cpg,
    n_samples,
    qc_cpg_order,
    qc_sample_order,
    qc_na,
    qc_nonfinite,
    qc_out_of_range,
    qc_changed_observed,
    qc_min,
    qc_max,
    ifelse(qc_pass, "PASS", "FAIL")
  ))
)

cat("\n========== FINAL QC REPORT ==========\n")
print(qc_report)

###############################################################################
# 7. ЗБЕРЕЖЕННЯ РЕЗУЛЬТАТІВ
###############################################################################

qc_rds <- file.path(
  results_dir,
  "TCGA_SKCM_KNN_final_QC_report.rds"
)

qc_csv <- file.path(
  results_dir,
  "TCGA_SKCM_KNN_final_QC_report.csv"
)

metadata_copy <- file.path(
  results_dir,
  "TCGA_SKCM_KNN_imputation_metadata.rds"
)

# Не перезаписуємо попередні результати
if (!any(file.exists(c(qc_rds, qc_csv, metadata_copy)))) {
  
  saveRDS(qc_report, qc_rds)
  
  write.csv(
    qc_report,
    qc_csv,
    row.names = FALSE
  )
  
  if (file.exists(metadata_path)) {
    file.copy(
      metadata_path,
      metadata_copy,
      overwrite = FALSE
    )
  }
  
} else {
  
  cat(
    "\nЗвіти вже існують. ",
    "Автоматичний перезапис не виконується.\n"
  )
}

###############################################################################
# 8. ФІНАЛЬНИЙ РЕЗУЛЬТАТ
###############################################################################

cat("\n====================================\n")
cat("TCGA-SKCM KNN PIPELINE COMPLETED\n")
cat("Final QC:", ifelse(qc_pass, "PASS", "FAIL"), "\n")
cat("====================================\n")

if (!qc_pass) {
  stop("Фінальний QC не пройдено.")
}
