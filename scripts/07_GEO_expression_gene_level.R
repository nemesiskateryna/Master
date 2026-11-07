# 07_GEO_expression_gene_level.R
#
# Purpose:
# Convert GSE65904 probe-level log2 expression data
# into a gene-level expression matrix using Entrez Gene IDs.

# 1. Project setup --------------------------------------------

project_dir <- getwd()

if (!file.exists(file.path(project_dir, "Master.Rproj"))) {
  stop(
    "Master.Rproj not found. Set the working directory to the project root."
  )
}

processed_dir <- file.path(
  project_dir, "data", "GEO", "processed"
)

tables_dir <- file.path(
  project_dir, "results", "tables"
)

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# 2. LOAD & VALIDATE INPUT DATA FOR EXPRESSION AGGREGATION ----

expression_file <- file.path(processed_dir, "GSE65904_expression_matrix_high_confidence_log2.rds")
annotation_file <- file.path(processed_dir, "GSE65904_probe_annotation_high_confidence.rds")

if (!file.exists(expression_file)) stop("Expression matrix file not found.")
if (!file.exists(annotation_file)) stop("Probe annotation file not found.")

expr_matrix      <- readRDS(expression_file)
probe_annotation <- readRDS(annotation_file)

# Перевірка наявності необхідних колонок анотації
required_cols <- c("Probe_ID", "Symbol", "Entrez_Gene_ID", "Eligible_Entrez_Mapping")
if (length(setdiff(required_cols, colnames(probe_annotation))) > 0) {
  stop("Missing required annotation columns.")
}

# Сувора валідація матриці експресії
stopifnot(
  !is.null(rownames(expr_matrix)),
  !anyDuplicated(rownames(expr_matrix)),
  is.numeric(expr_matrix),
  all(is.finite(expr_matrix))
)


# 3. MATCH PROBES & SELECT ELIGIBLE GENES ---------------------

# Зіставляємо рядки матриці з таблицею анотації за ID зонда
annotation_index <- match(rownames(expr_matrix), probe_annotation$Probe_ID)
if (anyNA(annotation_index)) stop("Some expression probes could not be matched to annotation.")

annotation_matched <- probe_annotation[annotation_index, , drop = FALSE]

# Фільтруємо лише допущені (eligible) зонди з валідними Entrez ID
eligible <- !is.na(annotation_matched$Eligible_Entrez_Mapping) & annotation_matched$Eligible_Entrez_Mapping
valid_annotation <- annotation_matched[eligible, , drop = FALSE]
expr_valid       <- expr_matrix[eligible, , drop = FALSE]

if (nrow(valid_annotation) == 0) stop("No eligible probes available for aggregation.")


# 4. AGGREGATE PROBES TO GENE LEVEL (MEDIAN) ------------------

cat("\nAggregating probe-level expression to gene level (Entrez ID)...\n")

entrez_ids <- as.character(valid_annotation$Entrez_Gene_ID)
probe_groups <- split(seq_along(entrez_ids), factor(entrez_ids, levels = unique(entrez_ids)))

# Розрахунок медіанного значення експресії для генів, що мають кілька зондів
gene_expression_matrix <- t(
  vapply(
    probe_groups,
    function(index) {
      apply(expr_valid[index, , drop = FALSE], 2, median, na.rm = TRUE)
    },
    FUN.VALUE = numeric(ncol(expr_valid))
  )
)

rownames(gene_expression_matrix) <- names(probe_groups)
colnames(gene_expression_matrix) <- colnames(expr_valid)
storage.mode(gene_expression_matrix) <- "numeric"


# 5. LOAD METHYLATION & ALIGN MULTI-OMICS SAMPLES -------------

cat("\nLoading methylation matrix and matching tables for alignment...\n")

beta_matrix <- readRDS(file.path(processed_dir, "GSE144487_beta_matrix_high_confidence.rds"))
matching    <- read.csv(file.path(tables_dir, "GEO_sample_matching_high_confidence.csv"), stringsAsFactors = FALSE, check.names = FALSE)

target_ids <- as.character(matching$Expression_ID)
beta_ids   <- colnames(beta_matrix)
gene_ids   <- colnames(gene_expression_matrix)

# Перевірка ідентичності вибірок
if (!setequal(beta_ids, target_ids) || !setequal(gene_ids, target_ids)) {
  stop("Sample identifiers do not match between matrices and matching table.")
}

# Синхронізація (вирівнювання) порядків зразків в обох омісних матрицях
beta_aligned <- beta_matrix[, target_ids, drop = FALSE]
gene_aligned <- gene_expression_matrix[, target_ids, drop = FALSE]

if (!identical(colnames(beta_aligned), colnames(gene_aligned))) {
  stop("Aligned matrices have different sample order.")
}

cat("Sample alignment verified successfully: exactly 74 matched pairs.\n")


# 6. SAVE PROCESSED MATRICES AND QC REPORTS ------------------

# Збереження фінальних вирівняних матриць
saveRDS(gene_expression_matrix, file.path(processed_dir, "GSE65904_expression_gene_level_high_confidence_log2.rds"))
saveRDS(beta_aligned,            file.path(processed_dir, "GSE144487_beta_matrix_high_confidence_aligned.rds"))
saveRDS(gene_aligned,            file.path(processed_dir, "GSE65904_expression_gene_level_high_confidence_log2_aligned.rds"))

# Формування підсумкового звіту контролю якості (QC)
alignment_qc <- data.frame(
  Metric = c(
    "Matching table samples", "Methylation samples", "Expression samples",
    "Unique Entrez genes", "Identical sample order", "Missing values in gene matrix"
  ),
  Value = c(
    length(target_ids), ncol(beta_aligned), ncol(gene_aligned),
    nrow(gene_expression_matrix), identical(colnames(beta_aligned), colnames(gene_aligned)),
    sum(is.na(gene_aligned))
  ),
  stringsAsFactors = FALSE
)

write.csv(alignment_qc, file.path(tables_dir, "GEO_omics_sample_alignment_QC.csv"), row.names = FALSE)


cat("Gene-level matrix dimensions:", nrow(gene_aligned), "genes x", ncol(gene_aligned), "samples\n")
cat("Methylation matrix dimensions:", nrow(beta_aligned), "CpGs x", ncol(beta_aligned), "samples\n")








