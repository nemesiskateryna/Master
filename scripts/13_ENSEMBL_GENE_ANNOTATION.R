# ============================================================
# TCGA-SKCM: ENSEMBL GENE ANNOTATION
# Ensembl ID -> Gene Symbol
# Мета:
# 1. Перевірити RNA-seq матрицю експресії.
# 2. Анотувати Ensembl gene IDs символами генів.
# 3. Оцінити повноту та однозначність анотації.
# 4. Перевірити дублікати Gene Symbols.
# ============================================================


# ============================================================
# 0. ПІДГОТОВКА СЕРЕДОВИЩА
# ============================================================

# Перевірка робочої директорії
cat("\n--- PROJECT DIRECTORY ---\n")
cat(getwd(), "\n")

# Створення каталогів для результатів
dir.create(
  "data/TCGA/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

# Перевірка необхідних пакетів
required_packages <- c(
  "dplyr",
  "AnnotationDbi",
  "org.Hs.eg.db"
)

package_status <- data.frame(
  Package = required_packages,
  Installed = vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
)

cat("\n--- PACKAGE STATUS ---\n")
print(package_status, row.names = FALSE)

missing_packages <- package_status$Package[
  !package_status$Installed
]

if (length(missing_packages) > 0) {
  stop(
    paste(
      "Не встановлені необхідні пакети:",
      paste(missing_packages, collapse = ", ")
    )
  )
}

library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)


# ============================================================
# 1. ЗАВАНТАЖЕННЯ ІНТЕГРОВАНОЇ КОГОРТИ
# ============================================================

integrated <- readRDS(
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics.rds"
)

# Отримання матриці RNA-seq експресії
rna_matrix <- integrated$expression_matrix

if (is.null(rna_matrix)) {
  stop("У інтегрованому об'єкті відсутня expression_matrix.")
}

# Перевірка структури матриці
if (is.null(rownames(rna_matrix)) ||
    is.null(colnames(rna_matrix))) {
  stop("Матриця RNA-seq не має назв рядків або стовпців.")
}

cat("\n--- RNA-seq MATRIX DIMENSIONS ---\n")

cat("Кількість генів:", nrow(rna_matrix), "\n")
cat("Кількість зразків:", ncol(rna_matrix), "\n")


# ============================================================
# 2. КОНТРОЛЬ ЯКОСТІ ІДЕНТИФІКАТОРІВ ГЕНІВ
# ============================================================

gene_ids <- rownames(rna_matrix)

cat("\n--- GENE ID QUALITY CONTROL ---\n")

missing_gene_ids <- sum(
  is.na(gene_ids) | gene_ids == ""
)

duplicated_gene_ids <- sum(
  duplicated(gene_ids)
)

ensembl_gene_count <- sum(
  grepl("^ENSG[0-9]+", gene_ids)
)

other_gene_id_count <- sum(
  !grepl("^ENSG[0-9]+", gene_ids)
)

cat("Відсутні ідентифікатори:", missing_gene_ids, "\n")
cat("Дубльовані ідентифікатори:", duplicated_gene_ids, "\n")
cat("Ensembl ID:", ensembl_gene_count, "\n")
cat("Інші ідентифікатори:", other_gene_id_count, "\n")

if (missing_gene_ids > 0) {
  stop("Виявлено відсутні ідентифікатори генів.")
}

if (duplicated_gene_ids > 0) {
  stop("Виявлено дубльовані Ensembl ID у RNA-seq матриці.")
}

if (other_gene_id_count > 0) {
  warning(
    "Не всі ідентифікатори відповідають формату Ensembl."
  )
}

# Приклади ідентифікаторів
cat("\n--- FIRST 15 GENE IDS ---\n")
print(head(gene_ids, 15))


# ============================================================
# 3. КОНТРОЛЬ ЯКОСТІ ЗНАЧЕНЬ ЕКСПРЕСІЇ
# ============================================================

cat("\n--- EXPRESSION VALUE QUALITY CONTROL ---\n")

missing_values <- sum(is.na(rna_matrix))

infinite_values <- sum(is.infinite(rna_matrix))

negative_values <- sum(
  rna_matrix < 0,
  na.rm = TRUE
)

minimum_expression <- min(
  rna_matrix,
  na.rm = TRUE
)

maximum_expression <- max(
  rna_matrix,
  na.rm = TRUE
)

cat("Відсутні значення:", missing_values, "\n")
cat("Нескінченні значення:", infinite_values, "\n")
cat("Від'ємні значення:", negative_values, "\n")
cat("Мінімальне значення:", minimum_expression, "\n")
cat("Максимальне значення:", maximum_expression, "\n")

if (missing_values > 0 ||
    infinite_values > 0 ||
    negative_values > 0) {
  stop(
    "Виявлено некоректні значення в RNA-seq матриці."
  )
}


# ============================================================
# 4. ОЦІНКА ВАРІАБЕЛЬНОСТІ ЕКСПРЕСІЇ ГЕНІВ
# ============================================================

gene_variance <- apply(
  rna_matrix,
  1,
  var
)

cat("\n--- GENE EXPRESSION VARIABILITY ---\n")

zero_variance_genes <- sum(
  gene_variance == 0,
  na.rm = TRUE
)

nonfinite_variance_genes <- sum(
  !is.finite(gene_variance)
)

cat(
  "Гени з нульовою дисперсією:",
  zero_variance_genes,
  "\n"
)

cat(
  "Гени з некоректною дисперсією:",
  nonfinite_variance_genes,
  "\n"
)

if (nonfinite_variance_genes > 0) {
  stop("Виявлено некоректні значення дисперсії.")
}


# ============================================================
# 5. ПІДГОТОВКА ENSEMBL ID ДО АНОТАЦІЇ
# ============================================================

# Вихідні ідентифікатори містять версійні суфікси.
# Наприклад:
# ENSG00000000003.15
#
# Для пошуку в org.Hs.eg.db суфікс версії видаляється.
# Оригінальні ідентифікатори матриці не змінюються.

ensembl_ids_clean <- sub(
  "\\.[0-9]+$",
  "",
  gene_ids
)

ensembl_unique <- unique(
  ensembl_ids_clean
)

cat("\n--- ENSEMBL ID PREPARATION ---\n")

cat(
  "Оригінальні ідентифікатори:",
  length(gene_ids),
  "\n"
)

cat(
  "Унікальні Ensembl ID без версії:",
  length(ensembl_unique),
  "\n"
)

if (anyDuplicated(ensembl_ids_clean) > 0) {
  warning(
    "Після видалення версійних суфіксів з'явилися повтори Ensembl ID."
  )
}


# ============================================================
# 6. ОТРИМАННЯ АНOТАЦІЇ З org.Hs.eg.db
# ============================================================

gene_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_unique,
  keytype = "ENSEMBL",
  columns = c(
    "SYMBOL",
    "GENENAME",
    "ENTREZID"
  )
)

cat("\n--- RAW GENE ANNOTATION ---\n")

cat("Рядків у таблиці анотації:", nrow(gene_annotation), "\n")

cat(
  "Ensembl ID з непорожнім Gene Symbol:",
  sum(
    !is.na(gene_annotation$SYMBOL) &
      gene_annotation$SYMBOL != ""
  ),
  "\n"
)

cat(
  "Відсутні Gene Symbols у рядках анотації:",
  sum(
    is.na(gene_annotation$SYMBOL) |
      gene_annotation$SYMBOL == ""
  ),
  "\n"
)

cat(
  "Відсутні Entrez ID у рядках анотації:",
  sum(is.na(gene_annotation$ENTREZID)),
  "\n"
)

# Збереження початкового результату анотації
saveRDS(
  gene_annotation,
  "data/TCGA/processed/TCGA_SKCM_gene_annotation_raw.rds"
)

write.csv(
  gene_annotation,
  "results/tables/TCGA_SKCM_gene_annotation_raw.csv",
  row.names = FALSE
)


# ============================================================
# 7. АУДИТ ВІДПОВІДНОСТІ ENSEMBL ID
# ============================================================

# Видалення рядків без валідного Gene Symbol
annotation_valid <- gene_annotation[
  !is.na(gene_annotation$ENSEMBL) &
    !is.na(gene_annotation$SYMBOL) &
    gene_annotation$SYMBOL != "",
  ,
  drop = FALSE
]

# Видалення рядків без валідного Entrez ID
entrez_valid <- gene_annotation[
  !is.na(gene_annotation$ENSEMBL) &
    !is.na(gene_annotation$ENTREZID) &
    gene_annotation$ENTREZID != "",
  ,
  drop = FALSE
]

# Кількість унікальних Gene Symbols для кожного Ensembl ID
symbol_counts <- tapply(
  annotation_valid$SYMBOL,
  annotation_valid$ENSEMBL,
  function(x) length(unique(x))
)

# Кількість унікальних Entrez ID для кожного Ensembl ID
entrez_counts <- tapply(
  entrez_valid$ENTREZID,
  entrez_valid$ENSEMBL,
  function(x) length(unique(x))
)

# Створення таблиці аудиту:
# один рядок на кожен вихідний Ensembl ID

annotation_audit <- data.frame(
  ENSEMBL = ensembl_unique,
  stringsAsFactors = FALSE
)

annotation_audit$Symbol_count <- as.integer(
  symbol_counts[annotation_audit$ENSEMBL]
)

annotation_audit$Entrez_count <- as.integer(
  entrez_counts[annotation_audit$ENSEMBL]
)

annotation_audit$Symbol_count[
  is.na(annotation_audit$Symbol_count)
] <- 0

annotation_audit$Entrez_count[
  is.na(annotation_audit$Entrez_count)
] <- 0

# Класифікація статусу анотації
annotation_audit$Mapping_status <- dplyr::case_when(
  annotation_audit$Symbol_count == 0 ~ "No gene symbol",
  annotation_audit$Symbol_count == 1 ~ "Unique gene symbol",
  annotation_audit$Symbol_count > 1 ~ "Multiple gene symbols"
)

cat("\n--- ENSEMBL MAPPING AUDIT ---\n")

print(
  table(annotation_audit$Mapping_status)
)

cat(
  "Унікальних Ensembl ID:",
  nrow(annotation_audit),
  "\n"
)

cat(
  "ID з одним Gene Symbol:",
  sum(annotation_audit$Symbol_count == 1),
  "\n"
)

cat(
  "ID з декількома Gene Symbols:",
  sum(annotation_audit$Symbol_count > 1),
  "\n"
)

cat(
  "ID без Gene Symbol:",
  sum(annotation_audit$Symbol_count == 0),
  "\n"
)


# ============================================================
# 8. ВИЗНАЧЕННЯ НЕАНОТОВАНИХ ТА НЕОДНОЗНАЧНИХ ID
# ============================================================

unmapped_ids <- annotation_audit$ENSEMBL[
  annotation_audit$Symbol_count == 0
]

ambiguous_ids <- annotation_audit$ENSEMBL[
  annotation_audit$Symbol_count > 1
]

ambiguous_annotation <- annotation_valid[
  annotation_valid$ENSEMBL %in% ambiguous_ids,
  ,
  drop = FALSE
]

cat("\n--- UNMAPPED ENSEMBL ID EXAMPLES ---\n")
print(head(unmapped_ids, 15))

cat("\n--- AMBIGUOUS ANNOTATION EXAMPLES ---\n")
print(head(ambiguous_annotation, 20))


# ============================================================
# 9. ФОРМУВАННЯ ФІНАЛЬНОЇ ТАБЛИЦІ АНОТАЦІЇ
# ============================================================

# Функція об'єднання унікальних значень.
# Всі альтернативні відповідності зберігаються через ";".

collapse_unique <- function(x) {
  
  x <- unique(as.character(x))
  
  x <- x[
    !is.na(x) &
      x != ""
  ]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(
    sort(x),
    collapse = ";"
  )
}

# Розподіл анотації за Ensembl ID
annotation_split <- split(
  gene_annotation,
  gene_annotation$ENSEMBL
)

# Створення таблиці з одним рядком на вихідний ген
gene_annotation_final <- data.frame(
  ENSEMBL_original = gene_ids,
  ENSEMBL = ensembl_ids_clean,
  stringsAsFactors = FALSE
)

# Збереження всіх кандидатів Gene Symbol
gene_annotation_final$SYMBOL_candidates <- vapply(
  ensembl_ids_clean,
  function(id) {
    
    x <- annotation_split[[id]]
    
    if (is.null(x)) {
      return(NA_character_)
    }
    
    collapse_unique(x$SYMBOL)
  },
  character(1)
)

# Збереження всіх кандидатів Gene Name
gene_annotation_final$GENENAME_candidates <- vapply(
  ensembl_ids_clean,
  function(id) {
    
    x <- annotation_split[[id]]
    
    if (is.null(x)) {
      return(NA_character_)
    }
    
    collapse_unique(x$GENENAME)
  },
  character(1)
)

# Збереження всіх кандидатів Entrez ID
gene_annotation_final$ENTREZID_candidates <- vapply(
  ensembl_ids_clean,
  function(id) {
    
    x <- annotation_split[[id]]
    
    if (is.null(x)) {
      return(NA_character_)
    }
    
    collapse_unique(x$ENTREZID)
  },
  character(1)
)

# Додавання результатів аудиту
audit_index <- match(
  gene_annotation_final$ENSEMBL,
  annotation_audit$ENSEMBL
)

gene_annotation_final$Mapping_status <-
  annotation_audit$Mapping_status[audit_index]

gene_annotation_final$Symbol_count <-
  annotation_audit$Symbol_count[audit_index]

gene_annotation_final$Entrez_count <-
  annotation_audit$Entrez_count[audit_index]

# Створення основної колонки Gene_Symbol:
# символ присвоюється лише за однозначної відповідності

gene_annotation_final$Gene_Symbol <- ifelse(
  gene_annotation_final$Mapping_status == "Unique gene symbol",
  gene_annotation_final$SYMBOL_candidates,
  NA_character_
)

rownames(gene_annotation_final) <- NULL


# ============================================================
# 10. КОНТРОЛЬ ЯКОСТІ ФІНАЛЬНОЇ АНОТАЦІЇ
# ============================================================

cat("\n--- FINAL GENE ANNOTATION QC ---\n")

cat(
  "Загальна кількість генів:",
  nrow(gene_annotation_final),
  "\n"
)

cat(
  "Унікальні оригінальні Ensembl ID:",
  length(unique(gene_annotation_final$ENSEMBL_original)),
  "\n"
)

cat(
  "Гени з однозначним Gene Symbol:",
  sum(!is.na(gene_annotation_final$Gene_Symbol)),
  "\n"
)

cat(
  "Неоднозначні відповідності:",
  sum(
    gene_annotation_final$Mapping_status ==
      "Multiple gene symbols"
  ),
  "\n"
)

cat(
  "Неанотовані гени:",
  sum(
    gene_annotation_final$Mapping_status ==
      "No gene symbol"
  ),
  "\n"
)

cat(
  "Початковий порядок генів збережено:",
  identical(
    gene_annotation_final$ENSEMBL_original,
    gene_ids
  ),
  "\n"
)

if (!identical(
  gene_annotation_final$ENSEMBL_original,
  gene_ids
)) {
  stop("Порядок генів у фінальній анотації не збігається з RNA-seq матрицею.")
}


# ============================================================
# 11. ПЕРЕВІРКА ПОВТОРНИХ GENE SYMBOLS
# ============================================================

unique_symbols <- gene_annotation_final$Gene_Symbol[
  !is.na(gene_annotation_final$Gene_Symbol)
]

duplicated_symbols <- unique(
  unique_symbols[duplicated(unique_symbols)]
)

duplicate_symbol_table <- gene_annotation_final[
  !is.na(gene_annotation_final$Gene_Symbol) &
    gene_annotation_final$Gene_Symbol %in% duplicated_symbols,
  c(
    "ENSEMBL_original",
    "Gene_Symbol",
    "GENENAME_candidates",
    "ENTREZID_candidates",
    "Mapping_status"
  ),
  drop = FALSE
]

duplicate_symbol_table <- duplicate_symbol_table[
  order(duplicate_symbol_table$Gene_Symbol),
  ,
  drop = FALSE
]

cat("\n--- DUPLICATED GENE SYMBOLS ---\n")

cat(
  "Кількість повторних Gene Symbols:",
  length(duplicated_symbols),
  "\n"
)

print(
  duplicate_symbol_table,
  row.names = FALSE
)


# ============================================================
# 12. ПОРІВНЯННЯ ПРОФІЛІВ ЕКСПРЕСІЇ ПОВТОРНИХ SYMBOLS
# ============================================================

# Перевіряємо, чи однакові профілі експресії
# Ensembl ID, що отримали однаковий Gene Symbol.
#
# Важливо: однаковий Gene Symbol не означає,
# що профілі експресії двох Ensembl ID ідентичні.

if (length(duplicated_symbols) > 0) {
  
  duplicate_expression_qc <- do.call(
    rbind,
    lapply(
      duplicated_symbols,
      function(symbol) {
        
        idx <- which(
          !is.na(gene_annotation_final$Gene_Symbol) &
            gene_annotation_final$Gene_Symbol == symbol
        )
        
        ids <- gene_annotation_final$ENSEMBL_original[idx]
        
        # Перевірка кількості Ensembl ID
        if (length(ids) != 2) {
          stop(
            paste(
              "Для Gene Symbol",
              symbol,
              "очікувалося 2 Ensembl ID, отримано:",
              length(ids)
            )
          )
        }
        
        # Перевірка наявності ID у матриці
        if (!all(ids %in% rownames(rna_matrix))) {
          stop(
            paste(
              "Ensembl ID відсутній у RNA-seq матриці:",
              symbol
            )
          )
        }
        
        # Отримання профілів експресії
        expression_profiles <- rna_matrix[
          ids,
          ,
          drop = FALSE
        ]
        
        profile_1 <- as.numeric(
          expression_profiles[1, ]
        )
        
        profile_2 <- as.numeric(
          expression_profiles[2, ]
        )
        
        # Кореляція Пірсона
        profile_correlation <- cor(
          profile_1,
          profile_2,
          method = "pearson"
        )
        
        data.frame(
          Gene_Symbol = symbol,
          Ensembl_ID_1 = ids[1],
          Ensembl_ID_2 = ids[2],
          Identical_profiles = identical(
            profile_1,
            profile_2
          ),
          Correlation = profile_correlation,
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  rownames(duplicate_expression_qc) <- NULL
  
} else {
  
  duplicate_expression_qc <- data.frame(
    Gene_Symbol = character(),
    Ensembl_ID_1 = character(),
    Ensembl_ID_2 = character(),
    Identical_profiles = logical(),
    Correlation = numeric()
  )
}

cat("\n--- DUPLICATE EXPRESSION PROFILE QC ---\n")

print(
  duplicate_expression_qc,
  row.names = FALSE
)


# ============================================================
# 13. ВЕРСІЯ БАЗИ АНОТАЦІЇ
# ============================================================

cat("\n--- ANNOTATION DATABASE VERSION ---\n")

annotation_metadata <- AnnotationDbi::metadata(
  org.Hs.eg.db
)

print(annotation_metadata)


# ============================================================
# 14. ФОРМУВАННЯ ПІДСУМКОВОГО QC-ЗВІТУ
# ============================================================

gene_annotation_qc <- list(
  
  # RNA-seq матриця
  total_genes = nrow(rna_matrix),
  total_samples = ncol(rna_matrix),
  
  # Ідентифікатори
  missing_gene_ids = missing_gene_ids,
  duplicated_gene_ids = duplicated_gene_ids,
  unique_ensembl_ids = length(ensembl_unique),
  
  # Експресія
  missing_expression_values = missing_values,
  infinite_expression_values = infinite_values,
  negative_expression_values = negative_values,
  minimum_expression = minimum_expression,
  maximum_expression = maximum_expression,
  
  # Варіабельність
  zero_variance_genes = zero_variance_genes,
  nonfinite_variance_genes = nonfinite_variance_genes,
  
  # Анотація
  uniquely_annotated_genes = sum(
    gene_annotation_final$Mapping_status ==
      "Unique gene symbol"
  ),
  
  ambiguous_mappings = sum(
    gene_annotation_final$Mapping_status ==
      "Multiple gene symbols"
  ),
  
  unmapped_genes = sum(
    gene_annotation_final$Mapping_status ==
      "No gene symbol"
  ),
  
  duplicated_gene_symbols = duplicate_expression_qc,
  
  # База даних
  annotation_database = "org.Hs.eg.db",
  annotation_metadata = annotation_metadata,
  
  # Метод
  annotation_method = "Ensembl ID to Gene Symbol mapping",
  version_suffix_removed_for_annotation = TRUE
)


# ============================================================
# 15. ЗБЕРЕЖЕННЯ РЕЗУЛЬТАТІВ
# ============================================================

# Фінальна таблиця анотації
saveRDS(
  gene_annotation_final,
  "data/TCGA/processed/TCGA_SKCM_gene_annotation_final.rds"
)

write.csv(
  gene_annotation_final,
  "results/tables/TCGA_SKCM_gene_annotation_final.csv",
  row.names = FALSE
)

# Таблиця аудиту Ensembl ID
saveRDS(
  annotation_audit,
  "data/TCGA/processed/TCGA_SKCM_gene_annotation_audit.rds"
)

write.csv(
  annotation_audit,
  "results/tables/TCGA_SKCM_gene_annotation_audit.csv",
  row.names = FALSE
)

# Таблиця повторних Gene Symbols
write.csv(
  duplicate_symbol_table,
  "results/tables/TCGA_SKCM_duplicate_gene_symbols.csv",
  row.names = FALSE
)

# Порівняння профілів експресії повторних символів
saveRDS(
  duplicate_expression_qc,
  "data/TCGA/processed/TCGA_SKCM_duplicate_gene_symbol_QC.rds"
)

write.csv(
  duplicate_expression_qc,
  "results/tables/TCGA_SKCM_duplicate_gene_symbol_QC.csv",
  row.names = FALSE
)

# Загальний QC-звіт
saveRDS(
  gene_annotation_qc,
  "data/TCGA/processed/TCGA_SKCM_gene_annotation_QC_final.rds"
)


# ============================================================
# 16. ПЕРЕВІРКА ЗБЕРЕЖЕННЯ ФАЙЛІВ
# ============================================================

output_files <- c(
  "data/TCGA/processed/TCGA_SKCM_gene_annotation_final.rds",
  "data/TCGA/processed/TCGA_SKCM_gene_annotation_audit.rds",
  "data/TCGA/processed/TCGA_SKCM_duplicate_gene_symbol_QC.rds",
  "data/TCGA/processed/TCGA_SKCM_gene_annotation_QC_final.rds",
  "results/tables/TCGA_SKCM_gene_annotation_final.csv",
  "results/tables/TCGA_SKCM_gene_annotation_audit.csv",
  "results/tables/TCGA_SKCM_duplicate_gene_symbols.csv",
  "results/tables/TCGA_SKCM_duplicate_gene_symbol_QC.csv"
)

output_status <- data.frame(
  File = output_files,
  Exists = file.exists(output_files),
  stringsAsFactors = FALSE
)

cat("\n--- SAVED OUTPUT FILES ---\n")

print(
  output_status,
  row.names = FALSE
)

if (!all(output_status$Exists)) {
  warning("Не всі вихідні файли було збережено.")
}


# ============================================================
# 17. ЗАВЕРШЕННЯ
# ============================================================

cat("\n============================================================\n")
cat("TCGA-SKCM RNA-seq QC AND GENE ANNOTATION COMPLETED\n")
cat("============================================================\n")

cat("Генів у матриці:", nrow(rna_matrix), "\n")
cat(
  "Однозначно анотованих генів:",
  sum(!is.na(gene_annotation_final$Gene_Symbol)),
  "\n"
)
cat(
  "Неоднозначних Ensembl ID:",
  sum(gene_annotation_final$Mapping_status == "Multiple gene symbols"),
  "\n"
)
cat(
  "Неанотованих Ensembl ID:",
  sum(gene_annotation_final$Mapping_status == "No gene symbol"),
  "\n"
)
cat(
  "Повторних Gene Symbols:",
  length(duplicated_symbols),
  "\n"
)

cat("\nУсі основні результати збережено.\n")

# ============================================================

# 1. Load required packages

library(AnnotationDbi)
library(org.Hs.eg.db)

# 2. Extract Ensembl IDs

gene_ids <- rownames(rna_matrix)

# Remove Ensembl version suffix for annotation only

ensembl_ids_clean <- sub(
  "\\.[0-9]+$",
  "",
  gene_ids
)

# 3. Create annotation table

gene_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_ids_clean),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME", "ENTREZID")
)

# 4. Inspect annotation dimensions

cat("\n--- GENE ANNOTATION DIMENSIONS ---\n")

cat("Original RNA-seq genes:", length(gene_ids), "\n")
cat("Unique Ensembl IDs:", length(unique(ensembl_ids_clean)), "\n")
cat("Annotation rows:", nrow(gene_annotation), "\n")

# 5. Check annotation completeness

cat("\n--- ANNOTATION COMPLETENESS ---\n")

cat(
  "Ensembl IDs with gene symbols:",
  sum(!is.na(gene_annotation$SYMBOL) & gene_annotation$SYMBOL != ""),
  "\n"
)

cat(
  "Missing gene symbols:",
  sum(is.na(gene_annotation$SYMBOL) | gene_annotation$SYMBOL == ""),
  "\n"
)

cat(
  "Missing Entrez IDs:",
  sum(is.na(gene_annotation$ENTREZID)),
  "\n"
)

# 6. Check duplicate Ensembl mappings

cat("\n--- MAPPING QC ---\n")

cat(
  "Duplicated Ensembl IDs in annotation:",
  sum(duplicated(gene_annotation$ENSEMBL)),
  "\n"
)

cat(
  "Ensembl IDs with multiple annotation rows:",
  sum(duplicated(gene_annotation$ENSEMBL) |
        duplicated(gene_annotation$ENSEMBL, fromLast = TRUE)),
  "\n"
)

# 7. Inspect first records

cat("\n--- ANNOTATION PREVIEW ---\n")

print(head(gene_annotation, 10))

# 8. Save original annotation result

saveRDS(
  gene_annotation,
  "data/TCGA/processed/TCGA_SKCM_gene_annotation_raw.rds"
)

write.csv(
  gene_annotation,
  "results/tables/TCGA_SKCM_gene_annotation_raw.csv",
  row.names = FALSE
)

cat("\nAnnotation saved:",
    file.exists("data/TCGA/processed/TCGA_SKCM_gene_annotation_raw.rds"),
    "\n"
)

# ============================================================
# TCGA-SKCM: ENSEMBL ANNOTATION MAPPING AUDIT
# ============================================================

# 1. Unique Ensembl IDs from the original RNA-seq matrix

ensembl_unique <- unique(ensembl_ids_clean)

# 2. Remove empty and missing symbols from annotation

annotation_valid <- gene_annotation[
  !is.na(gene_annotation$ENSEMBL) &
    !is.na(gene_annotation$SYMBOL) &
    gene_annotation$SYMBOL != "",
  ,
  drop = FALSE
]

# 3. Count unique symbols and Entrez IDs per Ensembl ID

symbol_counts <- tapply(
  annotation_valid$SYMBOL,
  annotation_valid$ENSEMBL,
  function(x) length(unique(x))
)

entrez_valid <- gene_annotation[
  !is.na(gene_annotation$ENSEMBL) &
    !is.na(gene_annotation$ENTREZID) &
    gene_annotation$ENTREZID != "",
  ,
  drop = FALSE
]

entrez_counts <- tapply(
  entrez_valid$ENTREZID,
  entrez_valid$ENSEMBL,
  function(x) length(unique(x))
)

# 4. Create one-row-per-Ensembl-ID audit table

annotation_audit <- data.frame(
  ENSEMBL = ensembl_unique,
  stringsAsFactors = FALSE
)

annotation_audit$Symbol_count <- as.integer(
  symbol_counts[annotation_audit$ENSEMBL]
)

annotation_audit$Entrez_count <- as.integer(
  entrez_counts[annotation_audit$ENSEMBL]
)

annotation_audit$Symbol_count[
  is.na(annotation_audit$Symbol_count)
] <- 0

annotation_audit$Entrez_count[
  is.na(annotation_audit$Entrez_count)
] <- 0

# 5. Classify mapping status

annotation_audit$Mapping_status <- dplyr::case_when(
  annotation_audit$Symbol_count == 0 ~ "No gene symbol",
  annotation_audit$Symbol_count == 1 ~ "Unique gene symbol",
  annotation_audit$Symbol_count > 1 ~ "Multiple gene symbols"
)

# 6. Summary

cat("\n--- UNIQUE ENSEMBL MAPPING SUMMARY ---\n")

print(
  table(
    annotation_audit$Mapping_status
  )
)

cat("\nTotal unique Ensembl IDs:",
    nrow(annotation_audit), "\n")

cat(
  "IDs with at least one gene symbol:",
  sum(annotation_audit$Symbol_count > 0), "\n"
)

cat(
  "IDs without gene symbol:",
  sum(annotation_audit$Symbol_count == 0), "\n"
)

cat(
  "IDs with multiple gene symbols:",
  sum(annotation_audit$Symbol_count > 1), "\n"
)

cat(
  "IDs with multiple Entrez IDs:",
  sum(annotation_audit$Entrez_count > 1), "\n"
)

# 7. Inspect ambiguous mappings

ambiguous_ids <- annotation_audit$ENSEMBL[
  annotation_audit$Symbol_count > 1
]

ambiguous_annotation <- annotation_valid[
  annotation_valid$ENSEMBL %in% ambiguous_ids,
  ,
  drop = FALSE
]

cat("\n--- AMBIGUOUS ANNOTATION EXAMPLES ---\n")

print(
  head(
    ambiguous_annotation,
    20
  )
)

# 8. Inspect unmapped Ensembl IDs

unmapped_ids <- annotation_audit$ENSEMBL[
  annotation_audit$Symbol_count == 0
]

cat("\n--- UNMAPPED ENSEMBL ID EXAMPLES ---\n")

print(
  head(unmapped_ids, 20)
)

# 9. Annotation database version

cat("\n--- ANNOTATION DATABASE VERSION ---\n")

print(
  AnnotationDbi::metadata(org.Hs.eg.db)
)

# 10. Save audit

saveRDS(
  annotation_audit,
  "data/TCGA/processed/TCGA_SKCM_gene_annotation_audit.rds"
)

write.csv(
  annotation_audit,
  "results/tables/TCGA_SKCM_gene_annotation_audit.csv",
  row.names = FALSE
)

cat("\nAudit saved:",
    file.exists(
      "data/TCGA/processed/TCGA_SKCM_gene_annotation_audit.rds"
    ),
    "\n"
)
# ============================================================
# TCGA-SKCM: FINAL GENE ANNOTATION TABLE
# ============================================================

# 1. Prepare annotation lookup

annotation_split <- split(
  gene_annotation,
  gene_annotation$ENSEMBL
)

# 2. Function to preserve all annotation candidates

collapse_unique <- function(x) {
  
  x <- unique(as.character(x))
  
  x <- x[
    !is.na(x) &
      x != ""
  ]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  paste(sort(x), collapse = ";")
}

# 3. Build annotation for every original RNA-seq gene

gene_annotation_final <- data.frame(
  ENSEMBL_original = gene_ids,
  ENSEMBL = ensembl_ids_clean,
  stringsAsFactors = FALSE
)

gene_annotation_final$SYMBOL_candidates <- vapply(
  ensembl_ids_clean,
  function(id) {
    x <- annotation_split[[id]]
    
    if (is.null(x)) return(NA_character_)
    
    collapse_unique(x$SYMBOL)
  },
  character(1)
)

gene_annotation_final$GENENAME_candidates <- vapply(
  ensembl_ids_clean,
  function(id) {
    x <- annotation_split[[id]]
    
    if (is.null(x)) return(NA_character_)
    
    collapse_unique(x$GENENAME)
  },
  character(1)
)

gene_annotation_final$ENTREZID_candidates <- vapply(
  ensembl_ids_clean,
  function(id) {
    x <- annotation_split[[id]]
    
    if (is.null(x)) return(NA_character_)
    
    collapse_unique(x$ENTREZID)
  },
  character(1)
)

# 4. Add mapping status

gene_annotation_final <- merge(
  gene_annotation_final,
  annotation_audit[
    ,
    c("ENSEMBL", "Mapping_status",
      "Symbol_count", "Entrez_count")
  ],
  by = "ENSEMBL",
  all.x = TRUE,
  sort = FALSE
)

# Restore original RNA-seq gene order

gene_annotation_final <- gene_annotation_final[
  match(gene_ids, gene_annotation_final$ENSEMBL_original),
  ,
  drop = FALSE
]

rownames(gene_annotation_final) <- NULL

# 5. Create unambiguous gene symbol column

gene_annotation_final$Gene_Symbol <- ifelse(
  gene_annotation_final$Mapping_status == "Unique gene symbol",
  gene_annotation_final$SYMBOL_candidates,
  NA_character_
)

# 6. Check final annotation

cat("\n--- FINAL ANNOTATION QC ---\n")

cat("Total rows:", nrow(gene_annotation_final), "\n")

cat(
  "Unique original Ensembl IDs:",
  length(unique(gene_annotation_final$ENSEMBL_original)),
  "\n"
)

cat(
  "Unique gene symbols:",
  sum(!is.na(gene_annotation_final$Gene_Symbol)),
  "\n"
)

cat(
  "Ambiguous mappings:",
  sum(gene_annotation_final$Mapping_status == "Multiple gene symbols"),
  "\n"
)

cat(
  "Unmapped genes:",
  sum(gene_annotation_final$Mapping_status == "No gene symbol"),
  "\n"
)

cat(
  "Original gene order preserved:",
  identical(
    gene_annotation_final$ENSEMBL_original,
    gene_ids
  ),
  "\n"
)

# 7. Check duplicate symbols among uniquely annotated genes

unique_symbols <- gene_annotation_final$Gene_Symbol[
  !is.na(gene_annotation_final$Gene_Symbol)
]

cat(
  "Duplicated unique Gene Symbols:",
  sum(duplicated(unique_symbols)),
  "\n"
)

# 8. Save final annotation

saveRDS(
  gene_annotation_final,
  "data/TCGA/processed/TCGA_SKCM_gene_annotation_final.rds"
)

write.csv(
  gene_annotation_final,
  "results/tables/TCGA_SKCM_gene_annotation_final.csv",
  row.names = FALSE
)

cat(
  "\nFinal annotation saved:",
  file.exists(
    "data/TCGA/processed/TCGA_SKCM_gene_annotation_final.rds"
  ),
  "\n"
)

# Identify duplicated Gene Symbols

duplicated_symbols <- unique(
  unique_symbols[duplicated(unique_symbols)]
)

duplicate_symbol_table <- gene_annotation_final[
  !is.na(gene_annotation_final$Gene_Symbol) &
    gene_annotation_final$Gene_Symbol %in% duplicated_symbols,
  c(
    "ENSEMBL_original",
    "Gene_Symbol",
    "GENENAME_candidates",
    "ENTREZID_candidates",
    "Mapping_status"
  ),
  drop = FALSE
]

duplicate_symbol_table <- duplicate_symbol_table[
  order(duplicate_symbol_table$Gene_Symbol),
]

cat("\n--- DUPLICATED GENE SYMBOLS ---\n")

print(
  duplicate_symbol_table,
  row.names = FALSE
)

cat(
  "\nNumber of duplicated Gene Symbols:",
  length(duplicated_symbols),
  "\n"
)

