# ============================================================
# GEO STEP: ENTREZ ID -> GENE SYMBOL
# Створення матриці експресії з унікальними Gene Symbols
# ============================================================

cat("\n========================================\n")
cat("GEO GENE SYMBOL MAPPING\n")
cat("========================================\n")

# 1. Завантаження центрального GEO-об'єкта

geo_file <- "data/GEO/processed/GEO_integrated_multiomics_high_confidence.rds"

if (!file.exists(geo_file)) {
  stop("Не знайдено центральний GEO RDS: ", geo_file)
}

geo_integrated <- readRDS(geo_file)

required_components <- c(
  "expression_matrix",
  "gene_annotation"
)

missing_components <- setdiff(
  required_components,
  names(geo_integrated)
)

if (length(missing_components) > 0) {
  stop(
    "Відсутні компоненти: ",
    paste(missing_components, collapse = ", ")
  )
}

expr <- geo_integrated$expression_matrix
gene_annotation <- geo_integrated$gene_annotation

# 2. Перевірка вихідних даних

cat("\n--- INPUT DIMENSIONS ---\n")

cat("Expression matrix:", dim(expr), "\n")
cat("Gene annotation:", dim(gene_annotation), "\n")

if (is.null(rownames(expr))) {
  stop("У матриці експресії відсутні rownames.")
}

if (!all(c("Entrez_Gene_ID", "Symbol") %in%
         names(gene_annotation))) {
  stop("В анотації відсутні Entrez_Gene_ID або Symbol.")
}

expr_ids <- as.character(rownames(expr))
annotation_ids <- as.character(
  gene_annotation$Entrez_Gene_ID
)

if (!identical(expr_ids, annotation_ids)) {
  stop(
    "Порядок Entrez ID в матриці та анотації не збігається."
  )
}

cat("Expression / annotation ID order: PASS\n")

# 3. Підготовка таблиці відповідності

gene_map <- data.frame(
  Entrez_Gene_ID = annotation_ids,
  Symbol = trimws(
    as.character(gene_annotation$Symbol)
  ),
  stringsAsFactors = FALSE
)

# Видаляємо порожні та відсутні символи

gene_map <- gene_map[
  !is.na(gene_map$Symbol) &
    gene_map$Symbol != "" &
    gene_map$Symbol != "NA",
  ,
  drop = FALSE
]

# 4. Перевірка дублікатів Entrez ID та Gene Symbols

cat("\n--- MAPPING DUPLICATE AUDIT ---\n")

duplicate_entrez <- unique(
  gene_map$Entrez_Gene_ID[
    duplicated(gene_map$Entrez_Gene_ID)
  ]
)

symbol_counts <- table(gene_map$Symbol)

ambiguous_symbols <- names(
  symbol_counts[symbol_counts > 1]
)

cat("Rows with valid symbols:", nrow(gene_map), "\n")
cat("Unique Entrez IDs:", length(unique(
  gene_map$Entrez_Gene_ID
)), "\n")

cat(
  "Duplicated Entrez IDs:",
  length(duplicate_entrez),
  "\n"
)

cat(
  "Ambiguous Gene Symbols:",
  length(ambiguous_symbols),
  "\n"
)

# 5. Формування списку однозначно зіставлених генів

unique_symbols <- names(
  symbol_counts[symbol_counts == 1]
)

gene_map_unique <- gene_map[
  gene_map$Symbol %in% unique_symbols &
    !gene_map$Entrez_Gene_ID %in% duplicate_entrez,
  ,
  drop = FALSE
]

# Додаткова перевірка унікальності

stopifnot(
  !anyDuplicated(gene_map_unique$Entrez_Gene_ID)
)

stopifnot(
  !anyDuplicated(gene_map_unique$Symbol)
)

cat("\n--- UNIQUE MAPPING SUMMARY ---\n")

cat(
  "Unique Entrez-to-Symbol mappings:",
  nrow(gene_map_unique),
  "\n"
)

cat(
  "Excluded ambiguous symbols:",
  length(ambiguous_symbols),
  "\n"
)

# 6. Побудова окремої матриці експресії за Gene Symbols

row_index <- match(
  gene_map_unique$Entrez_Gene_ID,
  expr_ids
)

if (anyNA(row_index)) {
  stop("Не всі унікальні Entrez ID знайдені в матриці.")
}

expr_symbol <- expr[row_index, , drop = FALSE]

rownames(expr_symbol) <- gene_map_unique$Symbol

# Перевірки матриці

stopifnot(
  !anyDuplicated(rownames(expr_symbol))
)

stopifnot(
  identical(
    colnames(expr_symbol),
    colnames(expr)
  )
)

stopifnot(
  all(is.finite(expr_symbol))
)

cat("\n--- SYMBOL MATRIX ---\n")

cat("Dimensions:", dim(expr_symbol), "\n")
cat("Unique rownames:", !anyDuplicated(
  rownames(expr_symbol)
), "\n")
cat("Sample order preserved:", identical(
  colnames(expr_symbol),
  colnames(expr)
), "\n")

# 7. Перевірка покриття імунних сигнатур

immune_signatures <- list(
  
  CYT = c(
    "GZMA",
    "PRF1"
  ),
  
  CD8 = c(
    "CD8A",
    "CD8B"
  ),
  
  IFNG = c(
    "IFNG",
    "CXCL9",
    "CXCL10",
    "STAT1",
    "IDO1"
  ),
  
  TIS = c(
    "CCL5",
    "CD27",
    "CD274",
    "CD276",
    "CD8A",
    "CMKLR1",
    "CXCL9",
    "CXCR6",
    "HLA-DQA1",
    "HLA-DRB1",
    "HLA-E",
    "IDO1",
    "LAG3",
    "NKG7",
    "PDCD1LG2",
    "PSMB10",
    "STAT1",
    "TIGIT"
  )
)

signature_coverage <- do.call(
  rbind,
  lapply(names(immune_signatures), function(sig) {
    
    genes <- immune_signatures[[sig]]
    
    present <- intersect(
      genes,
      rownames(expr_symbol)
    )
    
    missing <- setdiff(
      genes,
      rownames(expr_symbol)
    )
    
    data.frame(
      Signature = sig,
      Expected_genes = length(genes),
      Present_genes = length(present),
      Missing_genes = length(missing),
      Coverage_percent = round(
        100 * length(present) / length(genes),
        2
      ),
      Present = paste(present, collapse = ";"),
      Missing = paste(missing, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)

cat("\n--- IMMUNE SIGNATURE COVERAGE ---\n")
print(signature_coverage, row.names = FALSE)

# 8. Фінальний QC

mapping_qc <- c(
  annotation_order_matches_expression =
    identical(expr_ids, annotation_ids),
  
  unique_entrez_ids =
    !anyDuplicated(gene_map_unique$Entrez_Gene_ID),
  
  unique_gene_symbols =
    !anyDuplicated(gene_map_unique$Symbol),
  
  no_missing_expression_values =
    all(is.finite(expr_symbol)),
  
  same_74_sample_columns =
    ncol(expr_symbol) == 74,
  
  sample_order_preserved =
    identical(colnames(expr_symbol), colnames(expr)),
  
  all_TIS_genes_present =
    signature_coverage$Present_genes[
      signature_coverage$Signature == "TIS"
    ] == 18
)

cat("\n--- MAPPING QC ---\n")
print(mapping_qc)

cat(
  "\nAll mapping QC checks passed:",
  all(mapping_qc),
  "\n"
)

# 9. Збереження результатів

dir.create(
  "data/GEO/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

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

saveRDS(
  expr_symbol,
  "data/GEO/processed/GEO_expression_matrix_unique_symbols.rds"
)

saveRDS(
  gene_map_unique,
  "data/GEO/processed/GEO_unique_entrez_symbol_mapping.rds"
)

saveRDS(
  list(
    signatures = immune_signatures,
    coverage = signature_coverage,
    mapping_qc = mapping_qc,
    ambiguous_symbols = ambiguous_symbols,
    duplicate_entrez = duplicate_entrez
  ),
  "results/models/GEO_immune_signature_mapping_audit.rds"
)

write.csv(
  signature_coverage,
  "results/tables/GEO_immune_signature_gene_coverage.csv",
  row.names = FALSE
)

write.csv(
  gene_map_unique,
  "results/tables/GEO_unique_entrez_symbol_mapping.csv",
  row.names = FALSE
)

cat("\nAll mapping outputs saved.\n")
cat("GEO gene symbol mapping completed.\n")


# ============================================================
# GEO: CALCULATE 17-GENE TIS-LIKE IMMUNE SCORE
# ============================================================

# 1. Визначаємо доступні гени TIS
tis_genes_full <- c(
  "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1",
  "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E",
  "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10",
  "STAT1", "TIGIT"
)

# IDO1 відсутній у вихідній анотації
tis_genes_available <- setdiff(tis_genes_full, "IDO1")

# 2. Перевіряємо наявність усіх 17 генів
missing_tis_genes <- setdiff(
  tis_genes_available,
  rownames(expr_symbol)
)

if (length(missing_tis_genes) > 0) {
  stop(
    "Missing TIS genes: ",
    paste(missing_tis_genes, collapse = ", ")
  )
}

cat("Available TIS genes:", length(tis_genes_available), "\n")
cat("Missing TIS genes:", length(missing_tis_genes), "\n")

# 3. Вибираємо експресію 17 генів
tis_expression <- expr_symbol[
  tis_genes_available,
  ,
  drop = FALSE
]

# 4. Z-score стандартизація кожного гена між зразками
tis_zscore <- t(scale(t(tis_expression)))

if (any(!is.finite(tis_zscore))) {
  stop("Non-finite values detected after TIS standardization.")
}

# 5. Розраховуємо середній 17-генний показник
tis_17gene_score <- colMeans(tis_zscore)

# 6. Формуємо таблицю результатів
geo_immune_scores <- data.frame(
  Sample_ID = colnames(expr_symbol),
  TIS_17gene = as.numeric(tis_17gene_score),
  stringsAsFactors = FALSE
)

# 7. Визначаємо Immune-high / Immune-low за медіаною
tis_median <- median(
  geo_immune_scores$TIS_17gene,
  na.rm = TRUE
)

geo_immune_scores$Immune_phenotype <- ifelse(
  geo_immune_scores$TIS_17gene >= tis_median,
  "Immune-high",
  "Immune-low"
)

# 8. Перевірка результатів
cat("\n--- GEO TIS-17 SUMMARY ---\n")

cat("Number of samples:", nrow(geo_immune_scores), "\n")
cat("Number of genes:", length(tis_genes_available), "\n")
cat("TIS-17 median:", tis_median, "\n")

cat("\nPhenotype distribution:\n")
print(table(geo_immune_scores$Immune_phenotype))

cat("\nScore summary:\n")
print(summary(geo_immune_scores$TIS_17gene))

cat("\nSample order preserved:",
    identical(
      geo_immune_scores$Sample_ID,
      colnames(expr_symbol)
    ),
    "\n"
)

# 9. Зберігаємо результати
write.csv(
  geo_immune_scores,
  "results/tables/GEO_TIS_17gene_immune_scores.csv",
  row.names = FALSE
)

saveRDS(
  list(
    scores = geo_immune_scores,
    genes_used = tis_genes_available,
    gene_missing = "IDO1",
    original_signature_size = length(tis_genes_full),
    calculated_signature_size = length(tis_genes_available),
    standardization = "Gene-wise Z-score across GEO samples",
    classification = "Median split within GEO cohort"
  ),
  "results/models/GEO_TIS_17gene_immune_scores.rds"
)

cat("\nSaved:\n")
cat("- results/tables/GEO_TIS_17gene_immune_scores.csv\n")
cat("- results/models/GEO_TIS_17gene_immune_scores.rds\n")


# ============================================================
# GEO: AUXILIARY IMMUNE SIGNATURES AND CONVERGENT CHECK
# ============================================================

# 1. Визначаємо гени допоміжних сигнатур
immune_signatures <- list(
  CYT = c("GZMA", "PRF1"),
  CD8 = c("CD8A", "CD8B"),
  IFNG = c("IFNG", "CXCL9", "CXCL10", "STAT1", "IDO1")
)

# 2. Визначаємо доступні гени
signature_coverage <- do.call(
  rbind,
  lapply(names(immune_signatures), function(signature_name) {
    
    genes <- immune_signatures[[signature_name]]
    available <- intersect(genes, rownames(expr_symbol))
    
    data.frame(
      Signature = signature_name,
      Total_genes = length(genes),
      Available_genes = length(available),
      Missing_genes = paste(
        setdiff(genes, available),
        collapse = ";"
      ),
      Genes_used = paste(available, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)

cat("\n--- AUXILIARY SIGNATURE COVERAGE ---\n")
print(signature_coverage, row.names = FALSE)

# 3. Функція розрахунку середнього gene-wise Z-score
calculate_signature_score <- function(expression_matrix, genes) {
  
  available_genes <- intersect(
    genes,
    rownames(expression_matrix)
  )
  
  if (length(available_genes) == 0) {
    stop("No genes available for signature.")
  }
  
  selected_expression <- expression_matrix[
    available_genes,
    ,
    drop = FALSE
  ]
  
  standardized_expression <- t(
    scale(t(selected_expression))
  )
  
  if (any(!is.finite(standardized_expression))) {
    stop("Non-finite values in standardized signature.")
  }
  
  list(
    score = colMeans(standardized_expression),
    genes_used = available_genes
  )
}

# 4. Розраховуємо допоміжні сигнатури
cyt_result <- calculate_signature_score(
  expr_symbol,
  immune_signatures$CYT
)

cd8_result <- calculate_signature_score(
  expr_symbol,
  immune_signatures$CD8
)

ifng_genes_available <- intersect(
  immune_signatures$IFNG,
  rownames(expr_symbol)
)

ifng_result <- calculate_signature_score(
  expr_symbol,
  ifng_genes_available
)

# 5. Об'єднуємо оцінки
auxiliary_scores <- data.frame(
  Sample_ID = colnames(expr_symbol),
  CYT_score = as.numeric(cyt_result$score),
  CD8_score = as.numeric(cd8_result$score),
  IFNG_score = as.numeric(ifng_result$score),
  stringsAsFactors = FALSE
)

geo_immune_scores_complete <- merge(
  geo_immune_scores,
  auxiliary_scores,
  by = "Sample_ID",
  all = TRUE,
  sort = FALSE
)

# Відновлюємо початковий порядок зразків
geo_immune_scores_complete <- geo_immune_scores_complete[
  match(colnames(expr_symbol),
        geo_immune_scores_complete$Sample_ID),
]

# 6. Кореляції Spearman з TIS-17
auxiliary_names <- c(
  "CYT_score",
  "CD8_score",
  "IFNG_score"
)

correlation_results <- do.call(
  rbind,
  lapply(auxiliary_names, function(score_name) {
    
    test <- cor.test(
      geo_immune_scores_complete$TIS_17gene,
      geo_immune_scores_complete[[score_name]],
      method = "spearman",
      exact = FALSE
    )
    
    data.frame(
      Signature = score_name,
      Spearman_rho = unname(test$estimate),
      P_value = test$p.value,
      N = sum(
        complete.cases(
          geo_immune_scores_complete$TIS_17gene,
          geo_immune_scores_complete[[score_name]]
        )
      )
    )
  })
)

correlation_results$FDR_BH <- p.adjust(
  correlation_results$P_value,
  method = "BH"
)

correlation_results$Significant_FDR_005 <-
  correlation_results$FDR_BH < 0.05

# 7. Виводимо результати
cat("\n--- TIS-17 VS AUXILIARY IMMUNE SIGNATURES ---\n")
print(correlation_results, row.names = FALSE)

cat("\n--- SCORE SUMMARY ---\n")
print(
  summary(
    geo_immune_scores_complete[
      ,
      c(
        "TIS_17gene",
        "CYT_score",
        "CD8_score",
        "IFNG_score"
      )
    ]
  )
)

cat(
  "\nSample order preserved:",
  identical(
    geo_immune_scores_complete$Sample_ID,
    colnames(expr_symbol)
  ),
  "\n"
)

cat(
  "All scores finite:",
  all(is.finite(as.matrix(
    geo_immune_scores_complete[
      ,
      c(
        "TIS_17gene",
        "CYT_score",
        "CD8_score",
        "IFNG_score"
      )
    ]
  ))),
  "\n"
)

# 8. Зберігаємо результати
write.csv(
  signature_coverage,
  "results/tables/GEO_auxiliary_signature_coverage.csv",
  row.names = FALSE
)

write.csv(
  correlation_results,
  "results/tables/GEO_TIS17_auxiliary_correlations.csv",
  row.names = FALSE
)

write.csv(
  geo_immune_scores_complete,
  "results/tables/GEO_complete_immune_scores.csv",
  row.names = FALSE
)

saveRDS(
  list(
    scores = geo_immune_scores_complete,
    signature_coverage = signature_coverage,
    correlations = correlation_results,
    methods = list(
      scoring = "Mean of gene-wise Z-scores across 74 GEO samples",
      TIS = "17 genes; IDO1 unavailable",
      IFNG = "4 genes; IDO1 unavailable",
      correlation = "Spearman",
      multiple_testing = "Benjamini-Hochberg"
    )
  ),
  "results/models/GEO_complete_immune_scores.rds"
)

cat("\nSaved GEO auxiliary immune signature results.\n")



