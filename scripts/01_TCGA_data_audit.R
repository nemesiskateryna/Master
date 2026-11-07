# Title: 01_TCGA_data_audit.R
# Project: Master


# Опис: Скрипт ініціалізує середовище, перевіряє наявність файлів даних, за потреби
#       копіює їх із кореневої папки у проєкт, завантажує матриці експресії,
#       результати DESeq2, метилювання та клінічні дані, робить їх аудит і зберігає звіти.

# 1. Ініціалізація проєкту
source("scripts/00_setup.R")

# 2. Налаштування шляхів до даних
tcga_dir <- "data/TCGA/raw_saved"
processed_dir <- "data/TCGA/processed"
results_tables_dir <- "results/tables"

# Створюємо необхідні папки, якщо вони ще не існують
dir.create(tcga_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_tables_dir, recursive = TRUE, showWarnings = FALSE)

# Визначення імена файлів
expr_file    <- file.path(tcga_dir, "RNAseq_TCGA_SKCM_VST.rds")
counts_file  <- file.path(tcga_dir, "RNAseq_TCGA_SKCM_normalized_counts.rds")
dds_file     <- file.path(tcga_dir, "DESeq2_TCGA_SKCM.rds")
beta_file    <- file.path(tcga_dir, "beta_matrix_TCGA_SKCM_imputed.rds")
clinical_rds <- file.path(tcga_dir, "clinical_data_TCGA_SKCM_clean.rds")
clinical_csv <- file.path(tcga_dir, "clinical_data_TCGA_SKCM_clean.csv")

# 3. Перевірка наявності та копіювання файлів
source_dir <- "/Users/katie"
files_to_copy <- c(
  "RNAseq_TCGA_SKCM_VST.rds",
  "RNAseq_TCGA_SKCM_normalized_counts.rds",
  "DESeq2_TCGA_SKCM.rds",
  "beta_matrix_TCGA_SKCM_imputed.rds",
  "clinical_data_TCGA_SKCM_clean.rds"
)

for (f in files_to_copy) {
  source_file <- file.path(source_dir, f)
  destination_file <- file.path(tcga_dir, f)
  
  if (file.exists(source_file)) {
    if (!file.exists(destination_file)) {
      file.copy(source_file, destination_file)
      cat("COPIED:", f, "\n")
    } else {
      cat("ALREADY EXISTS:", f, "\n")
    }
  } else {
    cat("SOURCE NOT FOUND:", f, "\n")
  }
}

# 4. Завантаження даних у середовище R
expr.mat      <- readRDS(expr_file)
norm_counts   <- readRDS(counts_file)
dds           <- readRDS(dds_file)
beta_matrix   <- readRDS(beta_file)
clinical_data <- readRDS(clinical_rds)

# Виведення  інформації
cat("\n--- RNA VST matrix ---\n")
print(class(expr.mat))
print(dim(expr.mat))

cat("\n--- Normalized counts ---\n")
print(class(norm_counts))
print(dim(norm_counts))

cat("\n--- DESeq2 object ---\n")
print(class(dds))
print(dim(DESeq2::counts(dds)))

cat("\n--- Beta matrix (Methylation) ---\n")
print(class(beta_matrix))
print(dim(beta_matrix))

cat("\n--- Clinical data ---\n")
print(class(clinical_data))
print(dim(clinical_data))

# 5. Функція для загального аудиту матриць
check_matrix <- function(x, object_name) {
  cat("\nChecking:", object_name, "\n")
  if (!is.matrix(x) && !is.data.frame(x)) {
    warning(object_name, " is not a matrix or data.frame.")
    return(invisible(NULL))
  }
  x_mat <- as.matrix(x)
  cat("Dimensions:", paste(dim(x_mat), collapse = " x "), "\n")
  cat("Missing values:", sum(is.na(x_mat)), "\n")
  cat("Duplicated row names:", sum(duplicated(rownames(x_mat))), "\n")
  cat("Duplicated column names:", sum(duplicated(colnames(x_mat))), "\n")
  
  if (is.numeric(x_mat)) {
    cat("Infinite values:", sum(is.infinite(x_mat)), "\n")
    cat("Finite values:", sum(is.finite(x_mat)), "\n")
    cat("Range:", range(x_mat, na.rm = TRUE), "\n")
  }
  invisible(NULL)
}

# Запуск перевірки для числових матриць
check_matrix(expr.mat, "expr.mat")
check_matrix(norm_counts, "norm_counts")
check_matrix(beta_matrix, "beta_matrix")

# 6. Аналіз зразків та їх перетин (Overlap) між різними оміксними наборами
sample_ids <- list(
  RNA_VST     = colnames(expr.mat),
  RNA_counts  = colnames(norm_counts),
  DESeq2      = colnames(DESeq2::counts(dds)),
  Methylation = colnames(beta_matrix)
)

sample_summary <- data.frame(
  data_type = names(sample_ids),
  n_samples = vapply(sample_ids, length, integer(1))
)
print(sample_summary)

write.csv(sample_summary, file.path(results_tables_dir, "TCGA_sample_summary.csv"), row.names = FALSE)

# Створення матриці перетину зразків (хто є де)
overlap_matrix <- matrix(
  0,
  nrow = length(sample_ids),
  ncol = length(sample_ids),
  dimnames = list(names(sample_ids), names(sample_ids))
)

for (i in seq_along(sample_ids)) {
  for (j in seq_along(sample_ids)) {
    overlap_matrix[i, j] <- length(intersect(sample_ids[[i]], sample_ids[[j]]))
  }
}

print(overlap_matrix)
write.csv(overlap_matrix, file.path(results_tables_dir, "TCGA_sample_overlap.csv"))

# 7. Перевірка клінічних даних
clinical_summary <- data.frame(
  variable = names(clinical_data),
  class = vapply(clinical_data, function(x) paste(class(x), collapse = "/"), character(1)),
  missing = vapply(clinical_data, function(x) sum(is.na(x)), integer(1))
)
print(clinical_summary)

write.csv(clinical_summary, file.path(results_tables_dir, "TCGA_clinical_variable_audit.csv"), row.names = FALSE)

# 8. Збереження результатів перевірки
saveRDS(
  list(
    sample_ids = sample_ids,
    overlap_matrix = overlap_matrix,
    clinical_columns = names(clinical_data)
  ),
  file.path(processed_dir, "TCGA_audit_objects.rds")
)