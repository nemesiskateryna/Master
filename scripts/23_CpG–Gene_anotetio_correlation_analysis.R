 
# ============================================================
# КРОК 1. ВИТЯГУВАННЯ CpG-КОЕФІЦІЄНТІВ ФІНАЛЬНОЇ ELASTIC NET
# ============================================================

cat("\n--- LOADING FINAL MODEL ---\n")

model_file <- "results/models/TCGA_SKCM_final_ElasticNet_model.rds"

stopifnot(file.exists(model_file))

final_model_check <- readRDS(model_file)

# Перевірка структури моделі
cat("Model fields:\n")
print(names(final_model_check))

cat("\nSelected CpG:", length(final_model_check$selected_cpg_ids), "\n")
cat("Lambda:", final_model_check$lambda, "\n")
cat("Lambda index:", final_model_check$lambda_index, "\n")


# ------------------------------------------------------------
# 1. Отримання коефіцієнтів за фінальним lambda
# ------------------------------------------------------------

coef_matrix <- as.matrix(
  stats::coef(
    final_model_check$model,
    s = final_model_check$lambda
  )
)

cat("\nCoefficient matrix dimensions:\n")
print(dim(coef_matrix))

cat("\nFirst coefficients:\n")
print(head(coef_matrix))

# ------------------------------------------------------------
# 2. Формування таблиці коефіцієнтів
# ------------------------------------------------------------

coefficient_table <- data.frame(
  Feature = rownames(coef_matrix),
  Coefficient = as.numeric(coef_matrix[, 1]),
  stringsAsFactors = FALSE
)

# Видалення intercept
coefficient_table <- coefficient_table[
  coefficient_table$Feature != "(Intercept)",
]

# Перевірка відповідності CpG
stopifnot(
  setequal(
    coefficient_table$Feature,
    final_model_check$selected_cpg_ids
  )
)

# ------------------------------------------------------------
# 3. Визначення ненульових коефіцієнтів
# ------------------------------------------------------------

nonzero_cpg <- coefficient_table[
  coefficient_table$Coefficient != 0,
]

nonzero_cpg$Absolute_Coefficient <- abs(
  nonzero_cpg$Coefficient
)

# Напрямок асоціації
nonzero_cpg$Direction <- ifelse(
  nonzero_cpg$Coefficient > 0,
  "Immune-high",
  "Immune-low"
)

# Сортування за абсолютною величиною коефіцієнта
nonzero_cpg <- nonzero_cpg[
  order(
    nonzero_cpg$Absolute_Coefficient,
    decreasing = TRUE
  ),
]

rownames(nonzero_cpg) <- NULL

# ------------------------------------------------------------
# 4. Підсумкова статистика
# ------------------------------------------------------------

n_selected <- nrow(coefficient_table)
n_nonzero <- nrow(nonzero_cpg)
n_zero <- n_selected - n_nonzero

n_positive <- sum(nonzero_cpg$Coefficient > 0)
n_negative <- sum(nonzero_cpg$Coefficient < 0)

cat("\n--- CpG SIGNATURE SUMMARY ---\n")

cat("Total selected CpG:", n_selected, "\n")
cat("Non-zero coefficients:", n_nonzero, "\n")
cat("Zero coefficients:", n_zero, "\n")
cat("Positive coefficients:", n_positive, "\n")
cat("Negative coefficients:", n_negative, "\n")

cat("\n--- TOP 20 CpG FEATURES ---\n")
print(head(nonzero_cpg, 20))

# ------------------------------------------------------------
# 5. QC
# ------------------------------------------------------------

cpg_interpretation_qc <- c(
  correct_selected_feature_count = n_selected == 5000,
  no_intercept_in_feature_table =
    !("(Intercept)" %in% coefficient_table$Feature),
  all_coefficients_finite =
    all(is.finite(coefficient_table$Coefficient)),
  nonzero_coefficients_identified =
    n_nonzero > 0,
  directions_correct =
    all(nonzero_cpg$Direction == ifelse(
      nonzero_cpg$Coefficient > 0,
      "Immune-high",
      "Immune-low"
    )),
  no_duplicate_CpG =
    !anyDuplicated(nonzero_cpg$Feature)
)

cat("\n--- CpG INTERPRETATION QC ---\n")
print(cpg_interpretation_qc)

if (!all(cpg_interpretation_qc)) {
  stop("CpG interpretation QC: FAIL")
}

cat("\nCpG COEFFICIENT EXTRACTION: PASS\n")


cat("\n--- FINAL CpG SIGNATURE SUMMARY ---\n")

cat("Total selected CpG:", n_selected, "\n")
cat("Non-zero CpG:", n_nonzero, "\n")
cat("Zero CpG:", n_zero, "\n")
cat("Positive coefficients (Immune-high):", n_positive, "\n")
cat("Negative coefficients (Immune-low):", n_negative, "\n")

cat("\n--- TOP 20 CpG BY ABSOLUTE COEFFICIENT ---\n")
print(head(nonzero_cpg, 20))



# ============================================================
# КРОК 2. ЗБЕРЕЖЕННЯ ФІНАЛЬНОЇ CpG-СИГНАТУРИ
# ============================================================

# 1. Підготовка повної таблиці коефіцієнтів

coefficient_table_export <- coefficient_table

coefficient_table_export$Absolute_Coefficient <- abs(
  coefficient_table_export$Coefficient
)

coefficient_table_export$Direction <- ifelse(
  coefficient_table_export$Coefficient > 0,
  "Immune-high",
  ifelse(
    coefficient_table_export$Coefficient < 0,
    "Immune-low",
    "Not_selected_by_ElasticNet"
  )
)

coefficient_table_export <- coefficient_table_export[
  order(
    coefficient_table_export$Absolute_Coefficient,
    decreasing = TRUE
  ),
]

rownames(coefficient_table_export) <- NULL

# 2. Підготовка резюме сигнатури

signature_summary <- data.frame(
  Metric = c(
    "Total_variance_selected_CpG",
    "Nonzero_ElasticNet_CpG",
    "Zero_coefficient_CpG",
    "Positive_Immune_high_CpG",
    "Negative_Immune_low_CpG",
    "Nonzero_CpG_percentage",
    "ElasticNet_alpha",
    "Selected_lambda",
    "Lambda_path_index",
    "Training_sample_count"
  ),
  Value = c(
    n_selected,
    n_nonzero,
    n_zero,
    n_positive,
    n_negative,
    round(100 * n_nonzero / n_selected, 2),
    final_model_check$alpha,
    final_model_check$lambda,
    final_model_check$lambda_index,
    final_model_check$n_training_samples
  )
)

# 3. Шляхи збереження

signature_files <- c(
  "results/tables/TCGA_SKCM_ElasticNet_all_5000_CpG_coefficients.csv",
  "results/tables/TCGA_SKCM_ElasticNet_final_140_CpG_signature.csv",
  "results/tables/TCGA_SKCM_ElasticNet_CpG_signature_summary.csv",
  "results/models/TCGA_SKCM_ElasticNet_CpG_signature.rds"
)

# 4. Перевірка каталогів

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/models", recursive = TRUE, showWarnings = FALSE)

# 5. Збереження CSV

write.csv(
  coefficient_table_export,
  signature_files[1],
  row.names = FALSE
)

write.csv(
  nonzero_cpg,
  signature_files[2],
  row.names = FALSE
)

write.csv(
  signature_summary,
  signature_files[3],
  row.names = FALSE
)

# 6. Збереження RDS з усіма результатами

cpg_signature_results <- list(
  all_coefficients = coefficient_table_export,
  nonzero_signature = nonzero_cpg,
  summary = signature_summary,
  selected_cpg_ids = final_model_check$selected_cpg_ids,
  model_lambda = final_model_check$lambda,
  model_alpha = final_model_check$alpha,
  lambda_index = final_model_check$lambda_index,
  training_sample_count = final_model_check$n_training_samples
)

saveRDS(
  cpg_signature_results,
  signature_files[4]
)

# 7. QC збережених файлів

signature_save_qc <- file.exists(signature_files)

names(signature_save_qc) <- basename(signature_files)

cat("\n--- SAVED CpG SIGNATURE FILES ---\n")
print(data.frame(
  File = signature_files,
  Exists = signature_save_qc
))

# Перевірка кількості рядків
stopifnot(
  nrow(coefficient_table_export) == 5000,
  nrow(nonzero_cpg) == 140,
  sum(coefficient_table_export$Coefficient != 0) == 140,
  all(signature_save_qc)
)

cat("\nFINAL CpG SIGNATURE SAVING: PASS\n")




# ============================================================
# КРОК 3.1. ПЕРЕВІРКА ПАКЕТІВ ДЛЯ CpG-АНОТАЦІЇ
# ============================================================

required_annotation_packages <- c(
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "minfi",
  "AnnotationDbi"
)

annotation_package_status <- data.frame(
  Package = required_annotation_packages,
  Installed = vapply(
    required_annotation_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
)

cat("\n--- ANNOTATION PACKAGE STATUS ---\n")
print(annotation_package_status)

cat("\n--- GENOME BUILD ---\n")
cat("Expected annotation genome: hg19 / GRCh37\n")

# Перевірка наявності основної таблиці сигнатури
signature_file <- "results/tables/TCGA_SKCM_ElasticNet_final_140_CpG_signature.csv"

stopifnot(file.exists(signature_file))

signature_140 <- read.csv(
  signature_file,
  stringsAsFactors = FALSE
)

cat("\nSignature dimensions:\n")
print(dim(signature_140))

cat("\nFirst CpG IDs:\n")
print(head(signature_140$Feature, 10))

stopifnot(
  nrow(signature_140) == 140,
  !anyDuplicated(signature_140$Feature),
  all(grepl("^cg[0-9]+$", signature_140$Feature))
)

cat("\nANNOTATION PREPARATION QC: PASS\n")



# ============================================================
# КРОК 3.2. ВИПРАВЛЕНЕ ОТРИМАННЯ ILLUMINA 450K ANNOTATION
# ============================================================

# 1. Підключення пакета анотації
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)

# 2. Отримання об'єкта анотації
annotation_object <- IlluminaHumanMethylation450kanno.ilmn12.hg19

# 3. Отримання таблиці анотації
annotation_450k <- minfi::getAnnotation(annotation_object)

# 4. Перетворення на data.frame
annotation_450k_df <- as.data.frame(annotation_450k)

# 5. Перевірка структури
cat("\n--- ANNOTATION TABLE ---\n")

cat("Class:\n")
print(class(annotation_450k))

cat("\nDimensions:\n")
print(dim(annotation_450k_df))

cat("\nColumn names:\n")
print(colnames(annotation_450k_df))

cat("\nFirst rows:\n")
print(head(annotation_450k_df, 3))

cat("\nAnnotation row names example:\n")
print(head(rownames(annotation_450k_df)))

cat("\nCpG IDs in signature:", length(signature_140$Feature), "\n")
cat("CpG IDs in annotation:", nrow(annotation_450k_df), "\n")

cat("\nANNOTATION EXTRACTION: COMPLETED\n")


# ============================================================
# КРОК 3.3. ЗІСТАВЛЕННЯ CpG-СИГНАТУРИ З АНОТАЦІЄЮ hg19
# ============================================================

# 1. Поля анотації, які потрібні для інтерпретації

annotation_columns <- c(
  "chr",
  "pos",
  "strand",
  "UCSC_RefGene_Name",
  "UCSC_RefGene_Accession",
  "UCSC_RefGene_Group",
  "Relation_to_Island",
  "Islands_Name",
  "Enhancer",
  "Regulatory_Feature_Group",
  "DHS"
)

# Перевірка наявності всіх потрібних полів
stopifnot(
  all(annotation_columns %in% colnames(annotation_450k_df))
)

# 2. Зіставлення за CpG ID (rownames анотації)

annotation_match_index <- match(
  signature_140$Feature,
  rownames(annotation_450k_df)
)

# 3. Підрахунок відповідностей

n_matched <- sum(!is.na(annotation_match_index))
n_unmatched <- sum(is.na(annotation_match_index))

cat("\n--- CpG ANNOTATION MATCHING ---\n")
cat("Total signature CpG:", nrow(signature_140), "\n")
cat("Matched CpG:", n_matched, "\n")
cat("Unmatched CpG:", n_unmatched, "\n")

# 4. Створення анотованої таблиці

annotation_for_signature <- annotation_450k_df[
  annotation_match_index,
  annotation_columns,
  drop = FALSE
]

rownames(annotation_for_signature) <- NULL

cpg_annotated <- cbind(
  signature_140,
  annotation_for_signature
)

# 5. Перевірка відсутніх анотацій

unmatched_cpg <- signature_140$Feature[
  is.na(annotation_match_index)
]

cat("\nUnmatched CpG IDs:\n")
print(unmatched_cpg)

# 6. Попередня статистика анотації

cat("\n--- GENOMIC REGION DISTRIBUTION ---\n")
print(
  sort(
    table(cpg_annotated$UCSC_RefGene_Group, useNA = "ifany"),
    decreasing = TRUE
  )
)

cat("\n--- CpG ISLAND DISTRIBUTION ---\n")
print(
  sort(
    table(cpg_annotated$Relation_to_Island, useNA = "ifany"),
    decreasing = TRUE
  )
)

# 7. QC

annotation_matching_qc <- c(
  correct_signature_size = nrow(cpg_annotated) == 140,
  all_CpG_matched = n_matched == 140,
  no_unmatched_CpG = n_unmatched == 0,
  CpG_ID_order_preserved = identical(
    cpg_annotated$Feature,
    signature_140$Feature
  ),
  all_required_columns_present = all(
    annotation_columns %in% colnames(cpg_annotated)
  ),
  coefficient_values_preserved = isTRUE(all.equal(
    cpg_annotated$Coefficient,
    signature_140$Coefficient
  ))
)

cat("\n--- ANNOTATION MATCHING QC ---\n")
print(annotation_matching_qc)

if (!all(annotation_matching_qc)) {
  stop("Annotation matching QC: FAIL")
}

cat("\nCpG ANNOTATION MATCHING: PASS\n")





# ============================================================
# КРОК 3.4. ФОРМУВАННЯ CpG–GENE ТАБЛИЦІ
# ============================================================

# 1. Витягування зв'язків CpG–Gene

cpg_gene_list <- strsplit(
  ifelse(
    is.na(cpg_annotated$UCSC_RefGene_Name),
    "",
    cpg_annotated$UCSC_RefGene_Name
  ),
  split = ";",
  fixed = TRUE
)

cpg_gene_mapping <- data.frame(
  CpG_ID = rep(cpg_annotated$Feature, lengths(cpg_gene_list)),
  Gene = trimws(unlist(cpg_gene_list)),
  Coefficient = rep(
    cpg_annotated$Coefficient,
    lengths(cpg_gene_list)
  ),
  Direction = rep(
    cpg_annotated$Direction,
    lengths(cpg_gene_list)
  ),
  Genomic_Region = rep(
    cpg_annotated$UCSC_RefGene_Group,
    lengths(cpg_gene_list)
  ),
  stringsAsFactors = FALSE
)

# Видалення порожніх назв генів
cpg_gene_mapping <- cpg_gene_mapping[
  !is.na(cpg_gene_mapping$Gene) &
    cpg_gene_mapping$Gene != "",
]

# Видалення можливих дублікатів відповідностей
cpg_gene_mapping <- unique(cpg_gene_mapping)

rownames(cpg_gene_mapping) <- NULL

# 2. Зведена таблиця генів

gene_summary <- aggregate(
  CpG_ID ~ Gene,
  data = cpg_gene_mapping,
  FUN = function(x) length(unique(x))
)

colnames(gene_summary)[2] <- "Associated_CpG_Count"

# Кількість позитивних і негативних CpG для кожного гена
gene_directions <- aggregate(
  CpG_ID ~ Gene + Direction,
  data = cpg_gene_mapping,
  FUN = function(x) length(unique(x))
)

gene_direction_wide <- reshape(
  gene_directions,
  idvar = "Gene",
  timevar = "Direction",
  direction = "wide"
)

colnames(gene_direction_wide) <- sub(
  "CpG_ID\\.",
  "",
  colnames(gene_direction_wide)
)

gene_summary <- merge(
  gene_summary,
  gene_direction_wide,
  by = "Gene",
  all.x = TRUE
)

# Перейменування колонок
if ("Immune-high" %in% colnames(gene_summary)) {
  colnames(gene_summary)[
    colnames(gene_summary) == "Immune-high"
  ] <- "Immune_high_CpG"
}

if ("Immune-low" %in% colnames(gene_summary)) {
  colnames(gene_summary)[
    colnames(gene_summary) == "Immune-low"
  ] <- "Immune_low_CpG"
}

gene_summary[is.na(gene_summary)] <- 0

gene_summary <- gene_summary[
  order(
    gene_summary$Associated_CpG_Count,
    decreasing = TRUE
  ),
]

rownames(gene_summary) <- NULL

# 3. Підсумкова статистика

cat("\n--- CpG–GENE MAPPING SUMMARY ---\n")

cat("Annotated CpG:", nrow(cpg_annotated), "\n")
cat("CpG–Gene associations:", nrow(cpg_gene_mapping), "\n")
cat("Unique annotated genes:", length(unique(cpg_gene_mapping$Gene)), "\n")
cat("CpG without annotated gene:",
    sum(is.na(cpg_annotated$UCSC_RefGene_Name) |
          cpg_annotated$UCSC_RefGene_Name == ""),
    "\n")

cat("\n--- TOP 20 GENES BY NUMBER OF ASSOCIATED CpG ---\n")
print(head(gene_summary, 20))

# 4. Збереження результатів

annotation_output_files <- c(
  "results/tables/TCGA_SKCM_ElasticNet_140_CpG_annotated_hg19.csv",
  "results/tables/TCGA_SKCM_ElasticNet_CpG_gene_mapping.csv",
  "results/tables/TCGA_SKCM_ElasticNet_gene_summary.csv",
  "results/models/TCGA_SKCM_ElasticNet_annotated_signature.rds"
)

write.csv(
  cpg_annotated,
  annotation_output_files[1],
  row.names = FALSE
)

write.csv(
  cpg_gene_mapping,
  annotation_output_files[2],
  row.names = FALSE
)

write.csv(
  gene_summary,
  annotation_output_files[3],
  row.names = FALSE
)

saveRDS(
  list(
    annotated_CpG = cpg_annotated,
    CpG_gene_mapping = cpg_gene_mapping,
    gene_summary = gene_summary,
    genome_build = "hg19 / GRCh37",
    annotation_package = annotation_package
  ),
  annotation_output_files[4]
)

# 5. QC

annotation_output_qc <- c(
  annotated_CpG_count = nrow(cpg_annotated) == 140,
  all_CpG_have_unique_IDs = !anyDuplicated(cpg_annotated$Feature),
  gene_mapping_has_valid_CpG = all(
    cpg_gene_mapping$CpG_ID %in% cpg_annotated$Feature
  ),
  gene_mapping_has_no_empty_genes = all(
    !is.na(cpg_gene_mapping$Gene) &
      cpg_gene_mapping$Gene != ""
  ),
  gene_summary_has_unique_genes = !anyDuplicated(gene_summary$Gene),
  all_output_files_saved = all(file.exists(annotation_output_files))
)

cat("\n--- ANNOTATION OUTPUT QC ---\n")
print(annotation_output_qc)

if (!all(annotation_output_qc)) {
  stop("Annotation output QC: FAIL")
}

cat("\nANNOTATED CpG–GENE RESULTS: PASS\n")



cat("Checking gene annotation package...\n")

cat(
  "org.Hs.eg.db installed:",
  requireNamespace("org.Hs.eg.db", quietly = TRUE),
  "\n"
)

cat(
  "AnnotationDbi installed:",
  requireNamespace("AnnotationDbi", quietly = TRUE),
  "\n"
)

if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  
  cat("\nAvailable key types:\n")
  
  print(
    AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)
  )
  
  cat("\nAvailable columns:\n")
  
  print(
    AnnotationDbi::columns(org.Hs.eg.db::org.Hs.eg.db)
  )
}


# 1. Отримуємо Ensembl ID з матриці експресії
expression_ensembl <- sub(
  "\\..*$",
  "",
  rownames(expression_matrix)
)

cat("Expression Ensembl IDs:", length(expression_ensembl), "\n")
cat("Unique Ensembl IDs after version removal:",
    length(unique(expression_ensembl)), "\n")

# 2. Отримуємо відповідності SYMBOL -> ENSEMBL
gene_annotation_mapping <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = unique(annotated_genes),
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENSEMBL", "GENENAME")
)

# 3. Залишаємо лише гени, присутні в матриці експресії
gene_annotation_mapping <- gene_annotation_mapping[
  !is.na(gene_annotation_mapping$ENSEMBL) &
    gene_annotation_mapping$ENSEMBL %in% expression_ensembl,
]

gene_annotation_mapping <- unique(gene_annotation_mapping)

# 4. Перевіряємо кількість відповідностей
cat("\n--- GENE MAPPING SUMMARY ---\n")

cat(
  "Annotated genes:",
  length(unique(annotated_genes)),
  "\n"
)

cat(
  "Genes matched to expression:",
  length(unique(gene_annotation_mapping$SYMBOL)),
  "\n"
)

cat(
  "Genes without expression match:",
  length(setdiff(
    annotated_genes,
    unique(gene_annotation_mapping$SYMBOL)
  )),
  "\n"
)

# 5. Перевіряємо, чи є гени з декількома Ensembl ID
mapping_counts <- table(gene_annotation_mapping$SYMBOL)

multiple_ensembl <- mapping_counts[mapping_counts > 1]

cat("\nGenes with multiple Ensembl IDs:\n")
print(multiple_ensembl)

# 6. Переглядаємо отриману таблицю
cat("\n--- MAPPING PREVIEW ---\n")
print(head(gene_annotation_mapping, 20))

# 7. Список генів, які не вдалося зіставити
unmatched_genes <- setdiff(
  annotated_genes,
  unique(gene_annotation_mapping$SYMBOL)
)

cat("\n--- UNMATCHED GENES ---\n")
print(unmatched_genes)



# Перевіряємо незнайдені символи як альтернативні назви
alias_mapping <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = unmatched_genes,
  keytype = "ALIAS",
  columns = c("ALIAS", "SYMBOL", "ENSEMBL", "GENENAME")
)

# Залишаємо лише відповідності з Ensembl ID,
# які присутні в матриці експресії
alias_mapping <- alias_mapping[
  !is.na(alias_mapping$SYMBOL) &
    !is.na(alias_mapping$ENSEMBL) &
    alias_mapping$ENSEMBL %in% expression_ensembl,
]

alias_mapping <- unique(alias_mapping)

cat("\n--- ALIAS MAPPING RESULTS ---\n")
print(alias_mapping)

cat("\nNumber of unmatched aliases recovered:\n")
print(length(unique(alias_mapping$ALIAS)))

cat("\nAliases with multiple possible gene symbols:\n")
alias_counts <- table(alias_mapping$ALIAS)
print(alias_counts[alias_counts > 1])

cat("\nUnmatched genes still unresolved:\n")
unresolved_genes <- setdiff(
  unmatched_genes,
  unique(alias_mapping$ALIAS)
)
print(unresolved_genes)


# 1. Створюємо таблицю потенційних відповідностей
alias_review <- data.frame(
  Original_Gene = alias_mapping$ALIAS,
  Proposed_Symbol = alias_mapping$SYMBOL,
  Ensembl_ID = alias_mapping$ENSEMBL,
  Gene_Name = alias_mapping$GENENAME,
  stringsAsFactors = FALSE
)

# 2. Перевіряємо, чи є запропонований символ
# безпосередньо серед анотованих генів
alias_review$Proposed_Symbol_is_direct_annotation <-
  alias_review$Proposed_Symbol %in% annotated_genes

# 3. Перевіряємо, чи є відповідність однозначною
alias_review$Unique_Ensembl_for_Alias <-
  !duplicated(alias_review$Original_Gene) &
  !duplicated(alias_review$Original_Gene, fromLast = TRUE)

cat("\n--- ALIAS REVIEW TABLE ---\n")
print(alias_review, row.names = FALSE)

# 4. Окремо перевіряємо ADAM18
cat("\n--- ADAM18 CHECK ---\n")

print(
  AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = "ADAM18",
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENSEMBL", "GENENAME")
  )
)

# 5. Підсумок
cat("\n--- MAPPING COUNTS ---\n")
cat("Direct SYMBOL matches:",
    length(unique(gene_annotation_mapping$SYMBOL)), "\n")
cat("Alias candidates:",
    nrow(alias_review), "\n")
cat("Unresolved original symbols:",
    length(unresolved_genes), "\n")



# Перевіряємо прямі та альтернативні Ensembl ID
genes_to_check <- unique(c(
  unmatched_genes,
  alias_mapping$SYMBOL
))

direct_check <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = genes_to_check,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENSEMBL", "GENENAME")
)

direct_check <- unique(direct_check)

# Видаляємо NA та перевіряємо присутність у матриці
direct_check$Present_in_expression <-
  !is.na(direct_check$ENSEMBL) &
  direct_check$ENSEMBL %in% expression_ensembl

# Додаємо інформацію про кількість Ensembl ID на ген
direct_check$Ensembl_ID_count <- ave(
  as.character(direct_check$ENSEMBL),
  direct_check$SYMBOL,
  FUN = function(x) length(unique(x[!is.na(x)]))
)

cat("\n--- DIRECT GENE CHECK ---\n")

print(
  direct_check[
    order(direct_check$SYMBOL),
  ],
  row.names = FALSE
)

cat("\n--- GENES PRESENT IN EXPRESSION ---\n")

print(
  unique(
    direct_check$SYMBOL[
      direct_check$Present_in_expression
    ]
  )
)

cat("\n--- ADAM18 EXPRESSION CHECK ---\n")

adam18_ids <- direct_check$ENSEMBL[
  direct_check$SYMBOL == "ADAM18"
]

print(adam18_ids)

print(
  data.frame(
    Ensembl_ID = adam18_ids,
    Present_in_expression =
      adam18_ids %in% expression_ensembl
  )
)


# 1. Прямі відповідності
direct_mapping_final <- unique(data.frame(
  Original_Gene = gene_annotation_mapping$SYMBOL,
  Expression_Gene = gene_annotation_mapping$SYMBOL,
  ENSEMBL = gene_annotation_mapping$ENSEMBL,
  Mapping_Type = "Direct",
  stringsAsFactors = FALSE
))

# 2. Альтернативні відповідності
# Виключаємо ZNF833 і ADAM18
excluded_aliases <- c("ZNF833", "ADAM18")

alias_mapping_final <- alias_review[
  !alias_review$Original_Gene %in% excluded_aliases,
  c("Original_Gene", "Proposed_Symbol", "Ensembl_ID")
]

names(alias_mapping_final) <- c(
  "Original_Gene",
  "Expression_Gene",
  "ENSEMBL"
)

alias_mapping_final$Mapping_Type <- "Alias-derived"

# 3. Об'єднуємо таблиці
gene_mapping_final <- rbind(
  direct_mapping_final,
  alias_mapping_final
)

gene_mapping_final <- unique(gene_mapping_final)

# 4. Перевіряємо неоднозначні відповідності
mapping_duplicates <- gene_mapping_final[
  duplicated(gene_mapping_final$Original_Gene) |
    duplicated(gene_mapping_final$Original_Gene, fromLast = TRUE),
]

cat("\n--- FINAL GENE MAPPING ---\n")
print(gene_mapping_final, row.names = FALSE)

cat("\n--- MAPPING SUMMARY ---\n")
cat("Total mapping rows:", nrow(gene_mapping_final), "\n")
cat("Unique original genes:",
    length(unique(gene_mapping_final$Original_Gene)), "\n")
cat("Direct mappings:",
    sum(gene_mapping_final$Mapping_Type == "Direct"), "\n")
cat("Alias-derived mappings:",
    sum(gene_mapping_final$Mapping_Type == "Alias-derived"), "\n")

cat("\n--- DUPLICATE ORIGINAL GENES ---\n")
print(mapping_duplicates)

# 5. Перевіряємо наявність усіх Ensembl ID в експресії
gene_mapping_final$Present_in_expression <-
  gene_mapping_final$ENSEMBL %in% expression_ensembl

mapping_qc <- c(
  no_missing_ensembl = !anyNA(gene_mapping_final$ENSEMBL),
  all_ensembl_present = all(gene_mapping_final$Present_in_expression),
  no_duplicate_original_gene =
    !anyDuplicated(gene_mapping_final$Original_Gene),
  no_duplicate_ensembl =
    !anyDuplicated(gene_mapping_final$ENSEMBL)
)

cat("\n--- MAPPING QC ---\n")
print(mapping_qc)



# 1. Перевіряємо вихідну таблицю CpG–gene
cat("Original CpG-gene associations:", nrow(cpg_gene_mapping), "\n")
print(names(cpg_gene_mapping))

# 2. Додаємо відповідності генів
cpg_gene_expression_mapping <- merge(
  cpg_gene_mapping,
  gene_mapping_final,
  by.x = "Gene",
  by.y = "Original_Gene",
  all.x = TRUE,
  sort = FALSE
)

# 3. Визначаємо, які пари можна аналізувати
cpg_gene_expression_mapping$Mapping_Status <- ifelse(
  !is.na(cpg_gene_expression_mapping$ENSEMBL),
  "Matched",
  "Unmatched"
)

# 4. Розділяємо зіставлені та незiставлені пари
matched_cpg_gene_pairs <- cpg_gene_expression_mapping[
  cpg_gene_expression_mapping$Mapping_Status == "Matched",
]

unmatched_cpg_gene_pairs <- cpg_gene_expression_mapping[
  cpg_gene_expression_mapping$Mapping_Status == "Unmatched",
]

# 5. Підсумок
cat("\n--- CpG-GENE MAPPING SUMMARY ---\n")

cat("Original associations:",
    nrow(cpg_gene_mapping), "\n")

cat("Matched associations:",
    nrow(matched_cpg_gene_pairs), "\n")

cat("Unmatched associations:",
    nrow(unmatched_cpg_gene_pairs), "\n")

cat("Unique matched CpG:",
    length(unique(matched_cpg_gene_pairs$CpG_ID)), "\n")

cat("Unique matched genes:",
    length(unique(matched_cpg_gene_pairs$Gene)), "\n")

cat("\n--- MATCHED PAIRS PREVIEW ---\n")
print(head(matched_cpg_gene_pairs, 15))

cat("\n--- UNMATCHED PAIRS ---\n")
print(
  unmatched_cpg_gene_pairs[
    ,
    intersect(
      c("CpG_ID", "Gene", "Coefficient", "Direction"),
      names(unmatched_cpg_gene_pairs)
    )
  ]
)

# 6. Контроль якості
pair_mapping_qc <- c(
  original_pair_count_preserved =
    nrow(matched_cpg_gene_pairs) +
    nrow(unmatched_cpg_gene_pairs) ==
    nrow(cpg_gene_mapping),
  
  matched_pairs_have_ensembl =
    all(!is.na(matched_cpg_gene_pairs$ENSEMBL)),
  
  matched_ensembl_present_in_expression =
    all(matched_cpg_gene_pairs$ENSEMBL %in% expression_ensembl),
  
  no_duplicate_cpg_gene_pairs =
    !anyDuplicated(
      matched_cpg_gene_pairs[, c("CpG_ID", "Gene")]
    )
)

cat("\n--- PAIR MAPPING QC ---\n")
print(pair_mapping_qc)

if (!all(pair_mapping_qc)) {
  stop("CpG-gene mapping QC: FAIL")
}

cat("\nCpG-GENE MAPPING QC: PASS\n")



# 1. Перевірка Ensembl ID після видалення версій
cat("Expression rows:", nrow(expression_matrix), "\n")
cat("Unique Ensembl IDs:", length(unique(expression_ensembl)), "\n")
cat(
  "Duplicated Ensembl IDs:",
  sum(duplicated(expression_ensembl)),
  "\n"
)

# 2. Перевірка CpG та генів
beta_cpg_ids <- rownames(beta_matrix_check)

matched_cpg_gene_pairs$CpG_present <-
  matched_cpg_gene_pairs$CpG_ID %in% beta_cpg_ids

matched_cpg_gene_pairs$Gene_present <-
  matched_cpg_gene_pairs$ENSEMBL %in% expression_ensembl

cat("\n--- MATRIX PRESENCE CHECK ---\n")

cat(
  "CpG IDs present:",
  sum(matched_cpg_gene_pairs$CpG_present),
  "of", nrow(matched_cpg_gene_pairs), "\n"
)

cat(
  "Gene IDs present:",
  sum(matched_cpg_gene_pairs$Gene_present),
  "of", nrow(matched_cpg_gene_pairs), "\n"
)

# 3. Перевірка кількості пар
cat("\nMatched CpG-gene pairs:", nrow(matched_cpg_gene_pairs), "\n")
cat(
  "Unique CpG IDs:",
  length(unique(matched_cpg_gene_pairs$CpG_ID)),
  "\n"
)
cat(
  "Unique genes:",
  length(unique(matched_cpg_gene_pairs$ENSEMBL)),
  "\n"
)

# 4. Підсумкова перевірка
correlation_preparation_qc <- c(
  expression_ensembl_unique =
    !anyDuplicated(expression_ensembl),
  
  all_cpg_present =
    all(matched_cpg_gene_pairs$CpG_present),
  
  all_genes_present =
    all(matched_cpg_gene_pairs$Gene_present),
  
  correct_sample_count =
    ncol(expression_matrix) == 460 &&
    ncol(beta_matrix_check) == 460,
  
  identical_sample_order =
    identical(
      colnames(expression_matrix),
      colnames(beta_matrix_check)
    ),
  
  expected_matched_pair_count =
    nrow(matched_cpg_gene_pairs) == 94
)

cat("\n--- CORRELATION PREPARATION QC ---\n")
print(correlation_preparation_qc)


# ============================================================
# STEP 4.10. SPEARMAN CORRELATION:
# CpG METHYLATION vs GENE EXPRESSION
# ============================================================

# 1. Підготовка індексів
cpg_indices <- match(
  matched_cpg_gene_pairs$CpG_ID,
  rownames(beta_matrix_check)
)

gene_indices <- match(
  matched_cpg_gene_pairs$ENSEMBL,
  expression_ensembl
)

# 2. Функція обчислення кореляції
calculate_spearman <- function(i) {
  
  # Рівень метилювання CpG
  beta_values <- as.numeric(
    beta_matrix_check[cpg_indices[i], ]
  )
  
  # Рівень експресії відповідного гена
  expression_values <- as.numeric(
    expression_matrix[gene_indices[i], ]
  )
  
  # Враховуємо лише повні спостереження
  valid <- is.finite(beta_values) &
    is.finite(expression_values)
  
  beta_values <- beta_values[valid]
  expression_values <- expression_values[valid]
  
  n_samples <- length(beta_values)
  
  # Перевірка кількості та варіабельності
  if (
    n_samples < 10 ||
    length(unique(beta_values)) < 2 ||
    length(unique(expression_values)) < 2
  ) {
    return(data.frame(
      N = n_samples,
      Spearman_rho = NA_real_,
      P_value = NA_real_
    ))
  }
  
  # Spearman correlation test
  test_result <- tryCatch(
    suppressWarnings(
      cor.test(
        beta_values,
        expression_values,
        method = "spearman",
        exact = FALSE
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(test_result)) {
    return(data.frame(
      N = n_samples,
      Spearman_rho = NA_real_,
      P_value = NA_real_
    ))
  }
  
  data.frame(
    N = n_samples,
    Spearman_rho = unname(test_result$estimate),
    P_value = test_result$p.value
  )
}

# 3. Обчислюємо 94 кореляції
correlation_results_list <- lapply(
  seq_len(nrow(matched_cpg_gene_pairs)),
  calculate_spearman
)

correlation_statistics <- do.call(
  rbind,
  correlation_results_list
)

# 4. Формуємо результуючу таблицю
cpg_gene_correlations <- cbind(
  matched_cpg_gene_pairs,
  correlation_statistics
)

# 5. Поправка Benjamini-Hochberg
cpg_gene_correlations$P_value_BH <- p.adjust(
  cpg_gene_correlations$P_value,
  method = "BH"
)

# 6. Напрямок кореляції
cpg_gene_correlations$Correlation_Direction <- ifelse(
  is.na(cpg_gene_correlations$Spearman_rho),
  NA_character_,
  ifelse(
    cpg_gene_correlations$Spearman_rho > 0,
    "Positive",
    ifelse(
      cpg_gene_correlations$Spearman_rho < 0,
      "Negative",
      "Zero"
    )
  )
)

# 7. Статистична значущість
cpg_gene_correlations$Significant_BH_05 <-
  !is.na(cpg_gene_correlations$P_value_BH) &
  cpg_gene_correlations$P_value_BH < 0.05

# 8. Сортуємо за BH-adjusted p-value
cpg_gene_correlations <- cpg_gene_correlations[
  order(
    cpg_gene_correlations$P_value_BH,
    na.last = TRUE
  ),
]

rownames(cpg_gene_correlations) <- NULL

# 9. Перегляд результатів
cat("\n--- SPEARMAN CORRELATION RESULTS ---\n")

print(
  cpg_gene_correlations[
    ,
    c(
      "CpG_ID",
      "Gene",
      "Mapping_Type",
      "N",
      "Spearman_rho",
      "P_value",
      "P_value_BH",
      "Correlation_Direction",
      "Significant_BH_05"
    )
  ],
  row.names = FALSE
)

# 10. Загальна статистика
cat("\n--- CORRELATION SUMMARY ---\n")

cat(
  "Total tested pairs:",
  nrow(cpg_gene_correlations),
  "\n"
)

cat(
  "Successful tests:",
  sum(!is.na(cpg_gene_correlations$Spearman_rho)),
  "\n"
)

cat(
  "Significant before BH (p < 0.05):",
  sum(cpg_gene_correlations$P_value < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Significant after BH (FDR < 0.05):",
  sum(cpg_gene_correlations$Significant_BH_05),
  "\n"
)

cat(
  "Positive correlations:",
  sum(cpg_gene_correlations$Correlation_Direction == "Positive",
      na.rm = TRUE),
  "\n"
)

cat(
  "Negative correlations:",
  sum(cpg_gene_correlations$Correlation_Direction == "Negative",
      na.rm = TRUE),
  "\n"
)





# 1. Топ-20 кореляцій за абсолютною силою
top_correlations <- cpg_gene_correlations[
  order(
    abs(cpg_gene_correlations$Spearman_rho),
    decreasing = TRUE
  ),
]

cat("\n--- TOP 20 CORRELATIONS ---\n")

print(
  top_correlations[
    seq_len(min(20, nrow(top_correlations))),
    c(
      "CpG_ID",
      "Gene",
      "Mapping_Type",
      "N",
      "Spearman_rho",
      "P_value",
      "P_value_BH",
      "Correlation_Direction",
      "Significant_BH_05"
    )
  ],
  row.names = FALSE
)

# 2. Розподіл сили кореляцій
cat("\n--- CORRELATION DISTRIBUTION ---\n")

print(
  summary(cpg_gene_correlations$Spearman_rho)
)

cat(
  "\nMedian absolute correlation:",
  median(abs(cpg_gene_correlations$Spearman_rho)),
  "\n"
)

# 3. Перевірка значущості
cat("\n--- P-VALUE DISTRIBUTION ---\n")

print(
  summary(cpg_gene_correlations$P_value_BH)
)

# 4. Перевірка кількості зразків у тестах
cat("\n--- SAMPLE COUNT CHECK ---\n")

print(
  table(cpg_gene_correlations$N)
)

# 5. Порівняння прямих та alias-derived відповідностей
cat("\n--- RESULTS BY MAPPING TYPE ---\n")

print(
  aggregate(
    cbind(
      Spearman_rho,
      P_value_BH
    ) ~ Mapping_Type,
    data = cpg_gene_correlations,
    FUN = median
  )
)


# ============================================================
# STEP 4.12. STRATIFIED CpG-GENE CORRELATION ANALYSIS
# PRIMARY vs METASTATIC
# ============================================================

# 1. Перевіряємо типи зразків
sample_types <- integrated_knn_immune$immune_scores$sample_type

names(sample_types) <- integrated_knn_immune$immune_scores$Sample_ID

sample_types <- sample_types[colnames(expression_matrix)]

cat("\n--- SAMPLE TYPE DISTRIBUTION ---\n")
print(table(sample_types, useNA = "ifany"))

stopifnot(
  length(sample_types) == 460,
  !anyNA(sample_types),
  identical(
    names(sample_types),
    colnames(expression_matrix)
  )
)

# 2. Функція кореляції для окремої групи
calculate_stratified_spearman <- function(i, group_name) {
  
  group_indices <- which(sample_types == group_name)
  
  beta_values <- as.numeric(
    beta_matrix_check[cpg_indices[i], group_indices]
  )
  
  expression_values <- as.numeric(
    expression_matrix[gene_indices[i], group_indices]
  )
  
  valid <- is.finite(beta_values) &
    is.finite(expression_values)
  
  beta_values <- beta_values[valid]
  expression_values <- expression_values[valid]
  
  n_samples <- length(beta_values)
  
  if (
    n_samples < 10 ||
    length(unique(beta_values)) < 2 ||
    length(unique(expression_values)) < 2
  ) {
    return(data.frame(
      N = n_samples,
      Spearman_rho = NA_real_,
      P_value = NA_real_
    ))
  }
  
  test_result <- tryCatch(
    suppressWarnings(
      cor.test(
        beta_values,
        expression_values,
        method = "spearman",
        exact = FALSE
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(test_result)) {
    return(data.frame(
      N = n_samples,
      Spearman_rho = NA_real_,
      P_value = NA_real_
    ))
  }
  
  data.frame(
    N = n_samples,
    Spearman_rho = unname(test_result$estimate),
    P_value = test_result$p.value
  )
}

# 3. Розрахунок для Primary та Metastatic
stratified_results <- list()

for (group_name in c("Primary", "Metastatic")) {
  
  group_results <- do.call(
    rbind,
    lapply(
      seq_len(nrow(matched_cpg_gene_pairs)),
      function(i) {
        calculate_stratified_spearman(i, group_name)
      }
    )
  )
  
  group_results$P_value_BH <- p.adjust(
    group_results$P_value,
    method = "BH"
  )
  
  group_results$Sample_Type <- group_name
  
  stratified_results[[group_name]] <- cbind(
    matched_cpg_gene_pairs,
    group_results
  )
}

# 4. Об'єднуємо результати
stratified_correlations <- do.call(
  rbind,
  stratified_results
)

rownames(stratified_correlations) <- NULL

stratified_correlations$Significant_BH_05 <-
  !is.na(stratified_correlations$P_value_BH) &
  stratified_correlations$P_value_BH < 0.05

# 5. Підсумок за групами
cat("\n--- STRATIFIED CORRELATION SUMMARY ---\n")

for (group_name in c("Primary", "Metastatic")) {
  
  current <- stratified_correlations[
    stratified_correlations$Sample_Type == group_name,
  ]
  
  cat("\n", group_name, "\n")
  cat("Tested pairs:", nrow(current), "\n")
  cat(
    "Successful tests:",
    sum(!is.na(current$Spearman_rho)),
    "\n"
  )
  cat(
    "Significant after BH:",
    sum(current$Significant_BH_05),
    "\n"
  )
  cat(
    "Median absolute rho:",
    median(abs(current$Spearman_rho), na.rm = TRUE),
    "\n"
  )
}

# 6. Порівняння напрямку та сили між групами
primary_results <- stratified_correlations[
  stratified_correlations$Sample_Type == "Primary",
]

metastatic_results <- stratified_correlations[
  stratified_correlations$Sample_Type == "Metastatic",
]

comparison_correlations <- data.frame(
  CpG_ID = matched_cpg_gene_pairs$CpG_ID,
  Gene = matched_cpg_gene_pairs$Gene,
  Mapping_Type = matched_cpg_gene_pairs$Mapping_Type,
  Pooled_rho = cpg_gene_correlations$Spearman_rho[
    match(
      matched_cpg_gene_pairs$CpG_ID,
      cpg_gene_correlations$CpG_ID
    )
  ],
  Primary_rho = primary_results$Spearman_rho,
  Primary_FDR = primary_results$P_value_BH,
  Metastatic_rho = metastatic_results$Spearman_rho,
  Metastatic_FDR = metastatic_results$P_value_BH
)

comparison_correlations$Same_direction <-
  sign(comparison_correlations$Primary_rho) ==
  sign(comparison_correlations$Metastatic_rho)

comparison_correlations$Significant_in_both <-
  comparison_correlations$Primary_FDR < 0.05 &
  comparison_correlations$Metastatic_FDR < 0.05

cat("\n--- CROSS-GROUP COMPARISON ---\n")

cat(
  "Same correlation direction:",
  sum(comparison_correlations$Same_direction, na.rm = TRUE),
  "of 94\n"
)

cat(
  "Significant in both groups:",
  sum(comparison_correlations$Significant_in_both, na.rm = TRUE),
  "of 94\n"
)

# 7. Топ-20 за абсолютною кореляцією в обох групах
comparison_correlations <- comparison_correlations[
  order(
    -pmax(
      abs(comparison_correlations$Primary_rho),
      abs(comparison_correlations$Metastatic_rho),
      na.rm = TRUE
    )
  ),
]

cat("\n--- TOP 20 STRATIFIED COMPARISON ---\n")

print(
  head(
    comparison_correlations[
      ,
      c(
        "CpG_ID",
        "Gene",
        "Mapping_Type",
        "Pooled_rho",
        "Primary_rho",
        "Primary_FDR",
        "Metastatic_rho",
        "Metastatic_FDR",
        "Same_direction",
        "Significant_in_both"
      )
    ],
    20
  ),
  row.names = FALSE
)




# ============================================================
# STEP 4.13. CORRECT POOLED CORRELATION MATCHING
# ============================================================

# 1. Створюємо унікальний ключ CpG–Gene
make_pair_key <- function(cpg, gene) {
  paste(cpg, gene, sep = "|||")
}

pooled_keys <- make_pair_key(
  cpg_gene_correlations$CpG_ID,
  cpg_gene_correlations$Gene
)

comparison_keys <- make_pair_key(
  comparison_correlations$CpG_ID,
  comparison_correlations$Gene
)

# 2. Правильно зіставляємо pooled rho за CpG + Gene
pooled_match <- match(comparison_keys, pooled_keys)

comparison_correlations$Pooled_rho <- 
  cpg_gene_correlations$Spearman_rho[pooled_match]

comparison_correlations$Pooled_FDR <-
  cpg_gene_correlations$P_value_BH[pooled_match]

# 3. Перевіряємо коректність зіставлення
pooled_matching_qc <- c(
  all_pairs_matched = !anyNA(pooled_match),
  unique_pair_keys = !anyDuplicated(pooled_keys),
  correct_pair_count =
    length(pooled_match) == nrow(comparison_correlations),
  pooled_values_finite =
    all(is.finite(comparison_correlations$Pooled_rho))
)

cat("\n--- POOLED MATCHING QC ---\n")
print(pooled_matching_qc)

stopifnot(all(pooled_matching_qc))

# 4. Перераховуємо однаковий напрямок
comparison_correlations$Same_direction <-
  sign(comparison_correlations$Primary_rho) ==
  sign(comparison_correlations$Metastatic_rho)

# 5. Визначаємо стійкі асоціації
comparison_correlations$Significant_in_both <-
  !is.na(comparison_correlations$Primary_FDR) &
  !is.na(comparison_correlations$Metastatic_FDR) &
  comparison_correlations$Primary_FDR < 0.05 &
  comparison_correlations$Metastatic_FDR < 0.05

comparison_correlations$Robust_association <-
  comparison_correlations$Same_direction &
  comparison_correlations$Significant_in_both

# 6. Підсумок
cat("\n--- CORRECTED ROBUSTNESS SUMMARY ---\n")

cat(
  "Same direction:",
  sum(comparison_correlations$Same_direction, na.rm = TRUE),
  "of 94\n"
)

cat(
  "Significant in both groups:",
  sum(comparison_correlations$Significant_in_both, na.rm = TRUE),
  "of 94\n"
)

cat(
  "Robust associations:",
  sum(comparison_correlations$Robust_association, na.rm = TRUE),
  "of 94\n"
)

# 7. Пари зі зміною напрямку
opposite_direction_pairs <- comparison_correlations[
  !comparison_correlations$Same_direction |
    is.na(comparison_correlations$Same_direction),
  c(
    "CpG_ID", "Gene",
    "Pooled_rho",
    "Primary_rho", "Primary_FDR",
    "Metastatic_rho", "Metastatic_FDR"
  )
]

cat("\n--- OPPOSITE DIRECTION PAIRS ---\n")
print(opposite_direction_pairs, row.names = FALSE)

# 8. Стійкі асоціації: однаковий напрямок + FDR < 0.05 в обох групах
robust_pairs <- comparison_correlations[
  comparison_correlations$Robust_association,
]

robust_pairs <- robust_pairs[
  order(
    pmax(
      robust_pairs$Primary_FDR,
      robust_pairs$Metastatic_FDR
    )
  ),
]

cat("\n--- ROBUST CpG-GENE ASSOCIATIONS ---\n")

print(
  robust_pairs[
    ,
    c(
      "CpG_ID", "Gene", "Mapping_Type",
      "Pooled_rho",
      "Primary_rho", "Primary_FDR",
      "Metastatic_rho", "Metastatic_FDR"
    )
  ],
  row.names = FALSE
)






# ============================================================
# STEP 4.14. BOOTSTRAP COMPARISON OF SPEARMAN CORRELATIONS
# ============================================================

set.seed(20261106)

B <- 1000

primary_indices <- which(sample_types == "Primary")
metastatic_indices <- which(sample_types == "Metastatic")

bootstrap_correlation_difference <- function(i, B = 1000) {
  
  # Дані для відповідної CpG–gene пари
  beta_values <- as.numeric(
    beta_matrix_check[cpg_indices[i], ]
  )
  
  expression_values <- as.numeric(
    expression_matrix[gene_indices[i], ]
  )
  
  # Первинні пухлини
  x_primary <- beta_values[primary_indices]
  y_primary <- expression_values[primary_indices]
  
  # Метастатичні пухлини
  x_metastatic <- beta_values[metastatic_indices]
  y_metastatic <- expression_values[metastatic_indices]
  
  # Видалення невалідних спостережень
  valid_primary <- is.finite(x_primary) &
    is.finite(y_primary)
  
  valid_metastatic <- is.finite(x_metastatic) &
    is.finite(y_metastatic)
  
  x_primary <- x_primary[valid_primary]
  y_primary <- y_primary[valid_primary]
  
  x_metastatic <- x_metastatic[valid_metastatic]
  y_metastatic <- y_metastatic[valid_metastatic]
  
  n_primary <- length(x_primary)
  n_metastatic <- length(x_metastatic)
  
  # Початкові кореляції
  rho_primary <- cor(
    x_primary, y_primary,
    method = "spearman"
  )
  
  rho_metastatic <- cor(
    x_metastatic, y_metastatic,
    method = "spearman"
  )
  
  observed_difference <- rho_metastatic - rho_primary
  
  # Bootstrap
  bootstrap_differences <- replicate(B, {
    
    idx_primary <- sample.int(
      n_primary,
      size = n_primary,
      replace = TRUE
    )
    
    idx_metastatic <- sample.int(
      n_metastatic,
      size = n_metastatic,
      replace = TRUE
    )
    
    rho_p <- suppressWarnings(
      cor(
        x_primary[idx_primary],
        y_primary[idx_primary],
        method = "spearman"
      )
    )
    
    rho_m <- suppressWarnings(
      cor(
        x_metastatic[idx_metastatic],
        y_metastatic[idx_metastatic],
        method = "spearman"
      )
    )
    
    rho_m - rho_p
  })
  
  bootstrap_differences <- bootstrap_differences[
    is.finite(bootstrap_differences)
  ]
  
  if (length(bootstrap_differences) < B * 0.9) {
    return(data.frame(
      N_Primary = n_primary,
      N_Metastatic = n_metastatic,
      Delta_rho = observed_difference,
      CI_95_Lower = NA_real_,
      CI_95_Upper = NA_real_,
      Bootstrap_valid = length(bootstrap_differences)
    ))
  }
  
  ci <- quantile(
    bootstrap_differences,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  data.frame(
    N_Primary = n_primary,
    N_Metastatic = n_metastatic,
    Delta_rho = observed_difference,
    CI_95_Lower = unname(ci[1]),
    CI_95_Upper = unname(ci[2]),
    Bootstrap_valid = length(bootstrap_differences)
  )
}

# Розрахунок для всіх 94 пар
bootstrap_results <- do.call(
  rbind,
  lapply(
    seq_len(nrow(matched_cpg_gene_pairs)),
    function(i) {
      bootstrap_correlation_difference(i, B = B)
    }
  )
)

# Формуємо підсумкову таблицю
bootstrap_comparison <- data.frame(
  CpG_ID = matched_cpg_gene_pairs$CpG_ID,
  Gene = matched_cpg_gene_pairs$Gene,
  Mapping_Type = matched_cpg_gene_pairs$Mapping_Type,
  bootstrap_results
)

# Інтервал не містить нуля
bootstrap_comparison$CI_excludes_zero <-
  !is.na(bootstrap_comparison$CI_95_Lower) &
  (
    bootstrap_comparison$CI_95_Lower > 0 |
      bootstrap_comparison$CI_95_Upper < 0
  )

# Сортування за абсолютною різницею кореляцій
bootstrap_comparison <- bootstrap_comparison[
  order(
    abs(bootstrap_comparison$Delta_rho),
    decreasing = TRUE
  ),
]

rownames(bootstrap_comparison) <- NULL

# Підсумок
cat("\n--- BOOTSTRAP COMPARISON SUMMARY ---\n")

cat("Pairs tested:", nrow(bootstrap_comparison), "\n")

cat(
  "Valid bootstrap intervals:",
  sum(!is.na(bootstrap_comparison$CI_95_Lower)),
  "\n"
)

cat(
  "95% CI excluding zero:",
  sum(bootstrap_comparison$CI_excludes_zero),
  "\n"
)

# Топ-20 найбільших відмінностей
cat("\n--- TOP 20 DIFFERENCES ---\n")

print(
  head(
    bootstrap_comparison[
      ,
      c(
        "CpG_ID",
        "Gene",
        "Mapping_Type",
        "N_Primary",
        "N_Metastatic",
        "Delta_rho",
        "CI_95_Lower",
        "CI_95_Upper",
        "CI_excludes_zero"
      )
    ],
    20
  ),
  row.names = FALSE
)




# ============================================================
# STEP 4.15. FIX SAMPLE SIZES AND CALCULATE FISHER Z
# ============================================================

# 1. Додаємо правильну кількість зразків
comparison_correlations$N_Primary <- 102
comparison_correlations$N_Metastatic <- 358

# 2. Перевірка
stopifnot(
  nrow(comparison_correlations) == 94,
  all(comparison_correlations$N_Primary == 102),
  all(comparison_correlations$N_Metastatic == 358)
)

# 3. Fisher Z transformation
z_primary <- atanh(
  pmax(-0.999999, pmin(0.999999,
                       comparison_correlations$Primary_rho
  ))
)

z_metastatic <- atanh(
  pmax(-0.999999, pmin(0.999999,
                       comparison_correlations$Metastatic_rho
  ))
)

# 4. Standard error
se_difference <- sqrt(
  1 / (comparison_correlations$N_Primary - 3) +
    1 / (comparison_correlations$N_Metastatic - 3)
)

# 5. Z statistic
z_difference <- (
  z_metastatic - z_primary
) / se_difference

# 6. Двостороннє p-value
comparison_correlations$Fisher_Z_p_value <-
  2 * pnorm(
    abs(z_difference),
    lower.tail = FALSE
  )

# 7. BH correction
comparison_correlations$Fisher_Z_FDR <- p.adjust(
  comparison_correlations$Fisher_Z_p_value,
  method = "BH"
)

comparison_correlations$Fisher_Z_significant <-
  comparison_correlations$Fisher_Z_FDR < 0.05

# 8. QC
fisher_qc <- c(
  correct_pair_count =
    length(z_difference) == 94,
  finite_z_statistics =
    all(is.finite(z_difference)),
  finite_p_values =
    all(is.finite(comparison_correlations$Fisher_Z_p_value)),
  valid_p_value_range =
    all(comparison_correlations$Fisher_Z_p_value >= 0 &
          comparison_correlations$Fisher_Z_p_value <= 1),
  finite_adjusted_p_values =
    all(is.finite(comparison_correlations$Fisher_Z_FDR))
)

cat("\n--- FISHER Z QC ---\n")
print(fisher_qc)

stopifnot(all(fisher_qc))

# 9. Підсумок
cat("\n--- FISHER Z COMPARISON SUMMARY ---\n")

cat(
  "Nominal p < 0.05:",
  sum(comparison_correlations$Fisher_Z_p_value < 0.05),
  "\n"
)

cat(
  "BH FDR < 0.05:",
  sum(comparison_correlations$Fisher_Z_significant),
  "\n"
)

cat(
  "Bootstrap CI excludes zero:",
  sum(comparison_correlations$CI_excludes_zero),
  "\n"
)

# 10. Топ-20
fisher_top <- comparison_correlations[
  order(comparison_correlations$Fisher_Z_FDR),
]

cat("\n--- TOP 20 FISHER Z RESULTS ---\n")

print(
  head(
    fisher_top[
      ,
      c(
        "CpG_ID",
        "Gene",
        "Mapping_Type",
        "Primary_rho",
        "Metastatic_rho",
        "Delta_rho",
        "Fisher_Z_p_value",
        "Fisher_Z_FDR",
        "Fisher_Z_significant",
        "CI_excludes_zero"
      )
    ],
    20
  ),
  row.names = FALSE
)

cat("\n--- BOOTSTRAP COMPARISON ---\n")
print(dim(bootstrap_comparison))
print(names(bootstrap_comparison))
print(head(bootstrap_comparison, 3))

cat("\n--- BOOTSTRAP RESULTS ---\n")
print(class(bootstrap_results))
print(dim(bootstrap_results))
print(names(bootstrap_results))

if (is.data.frame(bootstrap_results)) {
  print(head(bootstrap_results, 3))
}

cat("\n--- BOOTSTRAP CORRELATION DIFFERENCE ---\n")
print(class(bootstrap_correlation_difference))
print(dim(bootstrap_correlation_difference))
print(names(bootstrap_correlation_difference))


# ============================================================
# STEP 4.16. MERGE FISHER Z AND BOOTSTRAP RESULTS
# ============================================================

# 1. Перевірка структури bootstrap
stopifnot(
  nrow(bootstrap_results) == 94,
  all(bootstrap_results$N_Primary == 102),
  all(bootstrap_results$N_Metastatic == 358),
  all(bootstrap_results$Bootstrap_valid == 1000)
)

# 2. Перевірка порядку CpG–Gene пар
stopifnot(
  nrow(comparison_correlations) == 94,
  nrow(matched_cpg_gene_pairs) == 94
)

same_order <- (
  comparison_correlations$CpG_ID ==
    matched_cpg_gene_pairs$CpG_ID &
    comparison_correlations$Gene ==
    matched_cpg_gene_pairs$Gene
)

cat("\n--- PAIR ORDER CHECK ---\n")
print(table(same_order))

# 3. Додаємо bootstrap-результати
# За умови однакового порядку 94 пар

comparison_correlations$Delta_rho <-
  bootstrap_results$Delta_rho

comparison_correlations$CI_95_Lower <-
  bootstrap_results$CI_95_Lower

comparison_correlations$CI_95_Upper <-
  bootstrap_results$CI_95_Upper

comparison_correlations$Bootstrap_valid <-
  bootstrap_results$Bootstrap_valid

# 4. Визначаємо, чи виключає CI нуль
comparison_correlations$CI_excludes_zero <-
  comparison_correlations$CI_95_Lower > 0 |
  comparison_correlations$CI_95_Upper < 0

# 5. Перевірка різниці кореляцій
comparison_correlations$Delta_rho_check <-
  comparison_correlations$Metastatic_rho -
  comparison_correlations$Primary_rho

delta_difference <- max(
  abs(
    comparison_correlations$Delta_rho -
      comparison_correlations$Delta_rho_check
  )
)

cat("\nMaximum Delta rho difference:", delta_difference, "\n")

# 6. Підсумкова таблиця
cat("\n--- COMBINED ANALYSIS SUMMARY ---\n")

cat(
  "Total CpG-Gene pairs:",
  nrow(comparison_correlations), "\n"
)

cat(
  "Fisher Z nominal p < 0.05:",
  sum(comparison_correlations$Fisher_Z_p_value < 0.05), "\n"
)

cat(
  "Fisher Z BH FDR < 0.05:",
  sum(comparison_correlations$Fisher_Z_significant), "\n"
)

cat(
  "Bootstrap CI excludes zero:",
  sum(comparison_correlations$CI_excludes_zero), "\n"
)

# 7. QC
combined_qc <- c(
  correct_pair_count =
    nrow(comparison_correlations) == 94,
  
  correct_pair_order =
    all(same_order),
  
  bootstrap_iterations_valid =
    all(comparison_correlations$Bootstrap_valid == 1000),
  
  finite_delta_rho =
    all(is.finite(comparison_correlations$Delta_rho)),
  
  finite_confidence_intervals =
    all(is.finite(comparison_correlations$CI_95_Lower)) &&
    all(is.finite(comparison_correlations$CI_95_Upper)),
  
  confidence_intervals_ordered =
    all(comparison_correlations$CI_95_Lower <=
          comparison_correlations$CI_95_Upper),
  
  delta_rho_consistent =
    delta_difference < 1e-10,
  
  ci_classification_correct =
    all(
      comparison_correlations$CI_excludes_zero ==
        (
          comparison_correlations$CI_95_Lower > 0 |
            comparison_correlations$CI_95_Upper < 0
        )
    )
)

cat("\n--- COMBINED ANALYSIS QC ---\n")
print(combined_qc)

stopifnot(all(combined_qc))

cat("\nCOMBINED FISHER Z + BOOTSTRAP: PASS\n")

# 8. Топ-20 Fisher Z
fisher_top <- comparison_correlations[
  order(comparison_correlations$Fisher_Z_FDR),
]

cat("\n--- TOP 20 FISHER Z RESULTS ---\n")

print(
  head(
    fisher_top[
      ,
      c(
        "CpG_ID",
        "Gene",
        "Mapping_Type",
        "Primary_rho",
        "Metastatic_rho",
        "Delta_rho",
        "Fisher_Z_p_value",
        "Fisher_Z_FDR",
        "Fisher_Z_significant",
        "CI_excludes_zero"
      ),
      drop = FALSE
    ],
    20
  ),
  row.names = FALSE
)




# ============================================================
# STEP 4.18. CORRECT MERGE BY CpG_ID + Gene
# ============================================================

# 1. Створюємо складені ключі
comparison_correlations$Pair_Key <- paste(
  comparison_correlations$CpG_ID,
  comparison_correlations$Gene,
  sep = "__"
)

bootstrap_comparison$Pair_Key <- paste(
  bootstrap_comparison$CpG_ID,
  bootstrap_comparison$Gene,
  sep = "__"
)

# 2. Перевіряємо унікальність ключів
stopifnot(
  !anyDuplicated(comparison_correlations$Pair_Key),
  !anyDuplicated(bootstrap_comparison$Pair_Key)
)

# 3. Знаходимо відповідність
bootstrap_match <- match(
  comparison_correlations$Pair_Key,
  bootstrap_comparison$Pair_Key
)

# 4. Перевіряємо, що всі 94 пари знайдені
stopifnot(
  length(bootstrap_match) == 94,
  !anyNA(bootstrap_match)
)

# 5. Переносимо bootstrap-результати у правильному порядку
bootstrap_columns <- c(
  "N_Primary",
  "N_Metastatic",
  "Delta_rho",
  "CI_95_Lower",
  "CI_95_Upper",
  "Bootstrap_valid",
  "CI_excludes_zero"
)

for (column in bootstrap_columns) {
  comparison_correlations[[column]] <-
    bootstrap_comparison[[column]][bootstrap_match]
}

# 6. Перевірка різниці кореляцій
comparison_correlations$Delta_rho_check <-
  comparison_correlations$Metastatic_rho -
  comparison_correlations$Primary_rho

delta_difference <- max(
  abs(
    comparison_correlations$Delta_rho -
      comparison_correlations$Delta_rho_check
  )
)

# 7. QC
combined_qc <- c(
  correct_pair_count =
    nrow(comparison_correlations) == 94,
  
  unique_pair_keys =
    !anyDuplicated(comparison_correlations$Pair_Key),
  
  all_bootstrap_pairs_matched =
    all(!is.na(bootstrap_match)),
  
  bootstrap_iterations_valid =
    all(comparison_correlations$Bootstrap_valid == 1000),
  
  finite_delta_rho =
    all(is.finite(comparison_correlations$Delta_rho)),
  
  finite_confidence_intervals =
    all(is.finite(comparison_correlations$CI_95_Lower)) &&
    all(is.finite(comparison_correlations$CI_95_Upper)),
  
  confidence_intervals_ordered =
    all(comparison_correlations$CI_95_Lower <=
          comparison_correlations$CI_95_Upper),
  
  delta_rho_consistent =
    delta_difference < 1e-10,
  
  ci_classification_correct =
    all(
      comparison_correlations$CI_excludes_zero ==
        (
          comparison_correlations$CI_95_Lower > 0 |
            comparison_correlations$CI_95_Upper < 0
        )
    )
)

cat("\n--- CORRECTED COMBINED QC ---\n")
print(combined_qc)

cat("\nMaximum Delta rho difference:",
    delta_difference, "\n")

stopifnot(all(combined_qc))

cat("\nCORRECTED FISHER Z + BOOTSTRAP: PASS\n")

# 8. Підсумок
cat("\n--- FINAL COMPARISON SUMMARY ---\n")

cat(
  "Total pairs:",
  nrow(comparison_correlations), "\n"
)

cat(
  "Fisher Z nominal p < 0.05:",
  sum(comparison_correlations$Fisher_Z_p_value < 0.05), "\n"
)

cat(
  "Fisher Z BH FDR < 0.05:",
  sum(comparison_correlations$Fisher_Z_significant), "\n"
)

cat(
  "Bootstrap CI excludes zero:",
  sum(comparison_correlations$CI_excludes_zero), "\n"
)

# 9. Топ-20 за Fisher Z FDR
fisher_top <- comparison_correlations[
  order(comparison_correlations$Fisher_Z_FDR),
]

cat("\n--- TOP 20 FISHER Z RESULTS ---\n")

print(
  head(
    fisher_top[
      ,
      c(
        "CpG_ID",
        "Gene",
        "Mapping_Type",
        "Primary_rho",
        "Metastatic_rho",
        "Delta_rho",
        "Fisher_Z_p_value",
        "Fisher_Z_FDR",
        "Fisher_Z_significant",
        "CI_excludes_zero"
      )
    ],
    20
  ),
  row.names = FALSE
)


# ============================================================
# STEP 4.19.1. CORRECT OVERLAP QC
# ============================================================

overlap_qc <- c(
  total_pairs_94 =
    length(fisher_significant) == 94 &&
    length(bootstrap_significant) == 94,
  
  fisher_significant_7 =
    sum(fisher_significant) == 7,
  
  bootstrap_significant_17 =
    sum(bootstrap_significant) == 17,
  
  both_methods_7 =
    sum(fisher_significant & bootstrap_significant) == 7,
  
  bootstrap_only_10 =
    sum(!fisher_significant & bootstrap_significant) == 10,
  
  fisher_only_0 =
    sum(fisher_significant & !bootstrap_significant) == 0,
  
  neither_77 =
    sum(!fisher_significant & !bootstrap_significant) == 77,
  
  partition_sums_to_94 =
    sum(
      fisher_significant & bootstrap_significant,
      !fisher_significant & bootstrap_significant,
      fisher_significant & !bootstrap_significant,
      !fisher_significant & !bootstrap_significant
    ) == 94
)

cat("\n--- CORRECTED OVERLAP QC ---\n")
print(overlap_qc)

stopifnot(all(overlap_qc))

cat("\nSTEP 4.19 OVERLAP QC: PASS\n")

# ============================================================
# STEP 4.20. SAVE FINAL CpG-EXPRESSION CORRELATION RESULTS
# ============================================================

# 1. Перевірка наявності об'єктів
required_objects <- c(
  "comparison_correlations",
  "supported_pairs",
  "overlap_summary",
  "overlap_qc",
  "cpg_gene_expression_mapping",
  "unmatched_cpg_gene_pairs",
  "stratified_correlations",
  "bootstrap_comparison"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = .GlobalEnv,
    inherits = FALSE
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Missing objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

# 2. Формування таблиці стійких асоціацій
robust_pairs_to_save <- comparison_correlations[
  comparison_correlations$Robust_association %in% TRUE,
]

# 3. Створення директорій
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/models", recursive = TRUE, showWarnings = FALSE)

# 4. Шляхи до файлів
output_files <- c(
  full_correlations =
    "results/tables/TCGA_SKCM_CpG_gene_expression_correlations.csv",
  
  expression_mapping =
    "results/tables/TCGA_SKCM_CpG_gene_expression_mapping.csv",
  
  unmatched_pairs =
    "results/tables/TCGA_SKCM_unmatched_CpG_gene_pairs.csv",
  
  stratified_correlations =
    "results/tables/TCGA_SKCM_stratified_CpG_gene_correlations.csv",
  
  group_comparison =
    "results/tables/TCGA_SKCM_CpG_gene_correlation_group_comparison.csv",
  
  bootstrap =
    "results/tables/TCGA_SKCM_CpG_gene_correlation_bootstrap.csv",
  
  robust_pairs =
    "results/tables/TCGA_SKCM_CpG_gene_correlation_robust_pairs.csv",
  
  supported_pairs =
    "results/tables/TCGA_SKCM_CpG_gene_correlation_supported_by_both_methods.csv",
  
  overlap_summary =
    "results/tables/TCGA_SKCM_CpG_gene_correlation_method_overlap.csv",
  
  overlap_qc =
    "results/tables/TCGA_SKCM_CpG_gene_correlation_overlap_QC.csv",
  
  combined_rds =
    "results/models/TCGA_SKCM_CpG_gene_expression_correlation_analysis.rds"
)

# 5. Збереження CSV
write.csv(
  comparison_correlations,
  output_files["full_correlations"],
  row.names = FALSE
)

write.csv(
  cpg_gene_expression_mapping,
  output_files["expression_mapping"],
  row.names = FALSE
)

write.csv(
  unmatched_cpg_gene_pairs,
  output_files["unmatched_pairs"],
  row.names = FALSE
)

write.csv(
  stratified_correlations,
  output_files["stratified_correlations"],
  row.names = FALSE
)

write.csv(
  comparison_correlations,
  output_files["group_comparison"],
  row.names = FALSE
)

write.csv(
  bootstrap_comparison,
  output_files["bootstrap"],
  row.names = FALSE
)

write.csv(
  robust_pairs_to_save,
  output_files["robust_pairs"],
  row.names = FALSE
)

write.csv(
  supported_pairs,
  output_files["supported_pairs"],
  row.names = FALSE
)

write.csv(
  overlap_summary,
  output_files["overlap_summary"],
  row.names = FALSE
)

write.csv(
  data.frame(
    QC_Check = names(overlap_qc),
    Passed = as.logical(overlap_qc)
  ),
  output_files["overlap_qc"],
  row.names = FALSE
)

# 6. Повний RDS-об'єкт аналізу
correlation_analysis_results <- list(
  pooled_and_stratified_comparison = comparison_correlations,
  robust_pairs = robust_pairs_to_save,
  supported_by_both_methods = supported_pairs,
  overlap_summary = overlap_summary,
  overlap_qc = overlap_qc,
  cpg_gene_expression_mapping = cpg_gene_expression_mapping,
  unmatched_cpg_gene_pairs = unmatched_cpg_gene_pairs,
  stratified_correlations = stratified_correlations,
  bootstrap_comparison = bootstrap_comparison,
  metadata = list(
    project = "TCGA-SKCM",
    total_samples = 460,
    primary_samples = 102,
    metastatic_samples = 358,
    tested_cpg_gene_pairs = 94,
    bootstrap_iterations = 1000,
    correlation_method = "Spearman",
    multiple_testing_correction = "Benjamini-Hochberg",
    fisher_z_interpretation = "Approximate supplementary comparison",
    bootstrap_intervals = "Individual 95% percentile confidence intervals"
  )
)

saveRDS(
  correlation_analysis_results,
  output_files["combined_rds"]
)

# 7. Перевірка фактичного збереження
save_qc <- file.exists(output_files)

cat("\n--- SAVED FILES ---\n")

print(
  data.frame(
    File = unname(output_files),
    Exists = as.logical(save_qc)
  ),
  row.names = FALSE
)

cat("\n--- FINAL OBJECT DIMENSIONS ---\n")
cat("All tested pairs:", nrow(comparison_correlations), "\n")
cat("Robust pairs:", nrow(robust_pairs_to_save), "\n")
cat("Supported by both methods:", nrow(supported_pairs), "\n")

# 8. Фінальний QC
final_save_qc <- c(
  all_files_exist = all(save_qc),
  full_table_has_94_pairs =
    nrow(comparison_correlations) == 94,
  robust_pairs_have_47 =
    nrow(robust_pairs_to_save) == 47,
  supported_pairs_have_7 =
    nrow(supported_pairs) == 7,
  overlap_qc_pass =
    all(overlap_qc),
  combined_rds_exists =
    file.exists(output_files["combined_rds"])
)

cat("\n--- FINAL SAVE QC ---\n")
print(final_save_qc)

stopifnot(all(final_save_qc))

cat("\nSTEP 4.20: RESULTS SAVED SUCCESSFULLY\n")



