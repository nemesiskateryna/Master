# ============================================================
# 06_GEO_expression_annotation.R
# 
#
# Purpose: Annotate GSE65904 Illumina probes using GPL10558 platform 
#          and perform an exhaustive audit of probe-to-gene mapping.
# ============================================================


# 1. PROJECT SETUP & DIRECTORIES ------------------------------

project_dir <- getwd()

processed_dir <- file.path(project_dir, "data", "GEO", "processed")
geo_dir       <- file.path(project_dir, "data", "GEO", "GSE65904")
tables_dir    <- file.path(project_dir, "results", "tables")

# Створюємо необхідні директорії, якщо вони відсутні
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n--- PROJECT SETUP COMPLETED ---\n")
cat("Project directory:", project_dir, "\n")


# 2. LOAD REQUIRED PACKAGES -----------------------------------

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("GEOquery")
}

library(GEOquery)


# 3. LOAD HIGH-CONFIDENCE EXPRESSION MATRIX -------------------

expression_file <- file.path(
  processed_dir,
  "GSE65904_expression_matrix_high_confidence_log2.rds"
)

if (!file.exists(expression_file)) {
  stop(
    "Expression matrix not found: ", expression_file,
    "\nPlease ensure script 05 has successfully saved the log2 matrix."
  )
}

expr_matrix_log2 <- readRDS(expression_file)

cat("\n--- EXPRESSION MATRIX SUMMARY ---\n")
cat("Dimensions (Probes x Samples):", nrow(expr_matrix_log2), "x", ncol(expr_matrix_log2), "\n")
cat("Unique probe IDs:", length(unique(rownames(expr_matrix_log2))), "\n")
cat("Unique sample IDs:", length(unique(colnames(expr_matrix_log2))), "\n")

# Сувора перевірка цілісності даних за допомогою stopifnot
stopifnot(
  !is.null(rownames(expr_matrix_log2)),
  !is.null(colnames(expr_matrix_log2)),
  !anyDuplicated(rownames(expr_matrix_log2)),
  !anyDuplicated(colnames(expr_matrix_log2))
)


# 4. DOWNLOAD OR LOAD GPL10558 PLATFORM ANNOTATION -----------

annotation_file <- file.path(processed_dir, "GPL10558_annotation_full.rds")

if (file.exists(annotation_file)) {
  cat("\nLoading previously saved full GPL10558 annotation...\n")
  gpl_table_full <- readRDS(annotation_file)
} else {
  cat("\nDownloading full GPL10558 platform annotation from NCBI GEO...\n")
  gpl <- GEOquery::getGEO("GPL10558", destdir = geo_dir)
  gpl_table_full <- GEOquery::Table(gpl)
  
  saveRDS(gpl_table_full, annotation_file)
  cat("Full annotation saved to:", annotation_file, "\n")
}

cat("\n--- FULL GPL10558 ANNOTATION INFO ---\n")
cat("Dimensions:", nrow(gpl_table_full), "x", ncol(gpl_table_full), "\n")


# 5. CLEANING & MATCHING PROBES TO ANNOTATION -----------------

# Допоміжна функція для очищення текстових полів (видалення пробілів, перетворення порожніх рядків у NA)
clean_field <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

gpl_table_full$Symbol_clean     <- clean_field(gpl_table_full$Symbol)
gpl_table_full$Entrez_clean     <- clean_field(gpl_table_full$Entrez_Gene_ID)
gpl_table_full$ILMN_Gene_clean  <- clean_field(gpl_table_full$ILMN_Gene)
gpl_table_full$Probe_Type_clean <- clean_field(gpl_table_full$Probe_Type)

# Зіставляємо зонди з матриці експресії з повною таблицею платформи
expression_probe_ids <- rownames(expr_matrix_log2)
annotation_match <- match(expression_probe_ids, gpl_table_full$ID)

matched_annotation <- gpl_table_full[annotation_match, , drop = FALSE]

cat("\n--- ANNOTATION MATCHING AUDIT ---\n")
cat("Expression probes:", length(expression_probe_ids), "\n")
cat("Matched in GPL annotation:", sum(!is.na(annotation_match)), "\n")
cat("Unmatched probes:", sum(is.na(annotation_match)), "\n")


# 6. MAPPING STATISTICS & QUALITY CONTROL ---------------------

mapping_audit <- data.frame(
  Probe_ID   = as.character(matched_annotation$ID),
  Symbol     = matched_annotation$Symbol_clean,
  Entrez_ID  = matched_annotation$Entrez_clean,
  ILMN_Gene  = matched_annotation$ILMN_Gene_clean,
  Source     = clean_field(matched_annotation$Source),
  Probe_Type = matched_annotation$Probe_Type_clean,
  Species    = clean_field(matched_annotation$Species),
  stringsAsFactors = FALSE
)

cat("\n--- PROBE-TO-GENE MAPPING SUMMARY ---\n")
cat("Total probes:", nrow(mapping_audit), "\n")
cat("Probes with Symbol:", sum(!is.na(mapping_audit$Symbol)), "\n")
cat("Probes with Entrez ID:", sum(!is.na(mapping_audit$Entrez_ID)), "\n")
cat("Probes with both Symbol and Entrez ID:", sum(!is.na(mapping_audit$Symbol) & !is.na(mapping_audit$Entrez_ID)), "\n")
cat("Probes with neither:", sum(is.na(mapping_audit$Symbol) & is.na(mapping_audit$Entrez_ID)), "\n")


# 7. IDENTIFYING AMBIGUITIES (1-to-many / many-to-1) ----------

valid_mapping <- mapping_audit[
  !is.na(mapping_audit$Symbol) &
    !is.na(mapping_audit$Entrez_ID) &
    grepl("^[0-9]+$", mapping_audit$Entrez_ID),
  , drop = FALSE
]

# Перевірка: чи відповідає один Entrez ID кільком різним символам
symbols_per_entrez <- aggregate(Symbol ~ Entrez_ID, data = valid_mapping, FUN = function(x) length(unique(x)))
multi_symbol_entrez <- symbols_per_entrez[symbols_per_entrez$Symbol > 1, , drop = FALSE]

# Перевірка: чи відповідає один символ кільком різним Entrez ID
entrez_per_symbol <- aggregate(Entrez_ID ~ Symbol, data = valid_mapping, FUN = function(x) length(unique(x)))
multi_entrez_symbols <- entrez_per_symbol[entrez_per_symbol$Entrez_ID > 1, , drop = FALSE]

cat("\n--- AMBIGUITY CHECK ---\n")
cat("Entrez IDs with multiple symbols:", nrow(multi_symbol_entrez), "\n")
cat("Symbols with multiple Entrez IDs:", nrow(multi_entrez_symbols), "\n")

# Формуємо таблиці амбівалентних мапінгів для збереження
ambiguous_symbol_mapping <- unique(
  valid_mapping[valid_mapping$Symbol %in% multi_entrez_symbols$Symbol, c("Symbol", "Entrez_ID"), drop = FALSE]
)
ambiguous_symbol_mapping <- ambiguous_symbol_mapping[order(ambiguous_symbol_mapping$Symbol, ambiguous_symbol_mapping$Entrez_ID), ]


# 8. SAVE AUDIT & AMBIGUITY TABLES ----------------------------

write.csv(
  ambiguous_symbol_mapping,
  file.path(tables_dir, "GSE65904_ambiguous_symbol_to_entrez.csv"),
  row.names = FALSE
)

write.csv(
  multi_symbol_entrez,
  file.path(tables_dir, "GSE65904_ambiguous_entrez_to_symbol.csv"),
  row.names = FALSE
)

cat("\n============================================\n")
cat("EXPRESSION ANNOTATION AUDIT COMPLETED SUCCESSFULLY\n")
cat("Saved tables to directory:", tables_dir, "\n")
cat("============================================\n")


