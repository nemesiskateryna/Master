# Title: 05_GEO_expression.R
# Project: Master


project_dir <- getwd()

expr_dir <- file.path(
  project_dir,
  "data", "GEO", "GSE65904"
)

processed_dir <- file.path(
  project_dir,
  "data", "GEO", "processed"
)

dir.create(expr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# Перевіряємо та завантажуємо пакет GEOquery, якщо він відсутній
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("GEOquery")
}

library(GEOquery)

expression_matrix_file <- file.path(processed_dir, "GSE65904_expression_matrix.rds")

if (!file.exists(expression_matrix_file)) {
  message("Завантаження платформи та матриці експресії GSE65904 через GEOquery...")
  
  # getGEO завантажує Series Matrix напряму з NCBI
  gse <- GEOquery::getGEO("GSE65904", GSEMatrix = TRUE, destdir = expr_dir)
  
  # Отримуємо перший елемент списку (оскільки це один Series)
  gse_expr <- gse[[1]]
  
  # Витягуємо числову матрицю експресії (досліджені значення)
  expr_matrix <- Biobase::exprs(gse_expr)
  
  # Зберігаємо готову матрицю
  saveRDS(expr_matrix, expression_matrix_file)
  
  message("Матрицю експресії успішно завантажено та збережено!")
} else {
  message("Матриця експресії вже була завантажена раніше. Завантажуємо з локального файлу...")
  expr_matrix <- readRDS(expression_matrix_file)
}

cat("\n--- GSE65904 EXPRESSION MATRIX INFO ---\n")
cat("Rows (Genes/Probes):", nrow(expr_matrix), "\n")
cat("Columns (Samples):", ncol(expr_matrix), "\n")


# 2. Load high-confidence matching ---------------------------

matching_file <- file.path(
  project_dir,
  "results", "tables",
  "GEO_sample_matching_high_confidence.csv"
)

if (!file.exists(matching_file)) {
  stop("High-confidence matching file not found.")
}

matching <- read.csv(
  matching_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 3. Inspect expression matrix -------------------------------

cat("\n--- EXPRESSION MATRIX QC ---\n")

cat("Rows (probes):", nrow(expr_matrix), "\n")
cat("Columns (samples):", ncol(expr_matrix), "\n")

cat("\nFirst 10 matrix column names:\n")
print(head(colnames(expr_matrix), 10))

cat("\nFirst 10 expression GSMs from matching:\n")
print(head(matching$GSM_expr, 10))

# 4. Check sample matching ------------------------------------

cat("\n--- SAMPLE MATCHING QC ---\n")

cat(
  "High-confidence pairs:",
  nrow(matching), "\n"
)

cat(
  "Unique expression GSMs:",
  length(unique(matching$GSM_expr)), "\n"
)

cat(
  "Expression GSMs found in matrix:",
  sum(matching$GSM_expr %in% colnames(expr_matrix)), "\n"
)

cat(
  "Expression GSMs absent from matrix:",
  sum(!matching$GSM_expr %in% colnames(expr_matrix)), "\n"
)

cat("\nMissing GSMs:\n")
print(setdiff(matching$GSM_expr, colnames(expr_matrix)))

cat(
  "\nDuplicated expression GSMs in matching:",
  sum(duplicated(matching$GSM_expr)), "\n"
)

# 5. Check expression matrix identifiers ---------------------

cat("\n--- MATRIX IDENTIFIER CHECK ---\n")

cat(
  "Duplicated matrix sample names:",
  sum(duplicated(colnames(expr_matrix))), "\n"
)

cat(
  "Duplicated probe IDs:",
  sum(duplicated(rownames(expr_matrix))), "\n"
)

cat("\nFirst 10 probe IDs:\n")
print(head(rownames(expr_matrix), 10))


# 6. Select high-confidence matched samples -------------------

selected_expr_gsm <- matching$GSM_expr

expr_matrix_high_confidence <- expr_matrix[
  ,
  selected_expr_gsm,
  drop = FALSE
]

# Use unified sample identifiers
colnames(expr_matrix_high_confidence) <- matching$Expression_ID

# Verify selection
cat("\n--- HIGH-CONFIDENCE EXPRESSION MATRIX ---\n")

cat("Dimensions:",
    nrow(expr_matrix_high_confidence), "x",
    ncol(expr_matrix_high_confidence), "\n")

cat("Duplicated sample IDs:",
    sum(duplicated(colnames(expr_matrix_high_confidence))), "\n")

cat("Sample order matches matching table:",
    identical(colnames(expr_matrix_high_confidence),
              matching$Expression_ID), "\n")


# 8. Inspect GEO processing metadata ---------------------------

cat("\n--- GEO SAMPLE METADATA ---\n")

# Inspect available metadata fields
metadata_columns <- colnames(Biobase::pData(gse_expr))

print(metadata_columns)

# Search for processing-related metadata
processing_columns <- grep(
  "process|normal|characteristics|description",
  metadata_columns,
  ignore.case = TRUE,
  value = TRUE
)

cat("\nProcessing-related columns:\n")
print(processing_columns)

# Print processing metadata, if available
if (length(processing_columns) > 0) {
  print(
    head(
      Biobase::pData(gse_expr)[, processing_columns, drop = FALSE],
      5
    )
  )
}

# 9. Check expression distributions by sample -----------------

cat("\n--- SAMPLE-LEVEL EXPRESSION DISTRIBUTION ---\n")

sample_distribution <- data.frame(
  GSM = colnames(expr_matrix_high_confidence),
  min = apply(expr_matrix_high_confidence, 2, min),
  median = apply(expr_matrix_high_confidence, 2, median),
  mean = colMeans(expr_matrix_high_confidence),
  max = apply(expr_matrix_high_confidence, 2, max)
)

print(summary(sample_distribution[, -1]))

# 10. Check quantile consistency -------------------------------

sample_quantiles <- apply(
  expr_matrix_high_confidence,
  2,
  quantile,
  probs = c(0.01, 0.25, 0.50, 0.75, 0.99),
  na.rm = TRUE
)

cat("\nVariation of sample quantiles:\n")

print(
  apply(sample_quantiles, 1, function(x) {
    c(min = min(x), max = max(x), range = diff(range(x)))
  })
)

# 11. Log2 transformation --------------------------------------

# GEO VALUE contains positive normalized signal intensities.
# Apply log2 transformation without additional normalization.

stopifnot(
  all(is.finite(expr_matrix_high_confidence)),
  all(expr_matrix_high_confidence > 0)
)

expr_matrix_log2 <- log2(expr_matrix_high_confidence)

cat("\n--- LOG2 EXPRESSION MATRIX QC ---\n")

cat("Dimensions:",
    nrow(expr_matrix_log2), "x",
    ncol(expr_matrix_log2), "\n")

cat("Minimum:", min(expr_matrix_log2), "\n")
cat("Maximum:", max(expr_matrix_log2), "\n")
cat("Median:", median(expr_matrix_log2), "\n")
cat("Mean:", mean(expr_matrix_log2), "\n")

cat("Missing values:", sum(is.na(expr_matrix_log2)), "\n")
cat("Infinite values:", sum(is.infinite(expr_matrix_log2)), "\n")

cat("\nQuantiles:\n")

print(
  quantile(
    as.vector(expr_matrix_log2),
    probs = c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1)
  )
)

# 12. Save expression matrices --------------------------------

  saveRDS(
    +     expr_matrix_high_confidence,
    +     file.path(
      +         processed_dir,
      +         "GSE65904_expression_matrix_high_confidence_linear.rds"
           )
     )

  saveRDS(
    +     expr_matrix_log2,
    +     file.path(
      +         processed_dir,
      +         "GSE65904_expression_matrix_high_confidence_log2.rds"
           )
    )
 
   cat("\nExpression matrices successfully saved.\n")


















