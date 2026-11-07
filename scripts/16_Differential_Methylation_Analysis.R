# ============================================================
# STEP 1. TCGA-SKCM METHYLATION DATA QC
# Перевірка об'єктів та відповідності зразків
# ============================================================

# 1. Inspect integrated object

cat("\n--- INTEGRATED OBJECT ---\n")

print(names(integrated))

# 2. Inspect matrix objects inside integrated list

matrix_fields <- names(integrated)[
  vapply(integrated, is.matrix, logical(1))
]

cat("\n--- MATRIX OBJECTS ---\n")

for (nm in matrix_fields) {
  
  cat(
    nm, ": ",
    paste(dim(integrated[[nm]]), collapse = " x "),
    "\n"
  )
}

# 3. Check immune score object

cat("\n--- IMMUNE SCORES ---\n")

cat("Number of samples:", nrow(immune_scores), "\n")

cat(
  "Duplicated sample IDs:",
  sum(duplicated(immune_scores$Sample_ID)),
  "\n"
)

cat(
  "Missing TIS scores:",
  sum(is.na(immune_scores$TIS_mean_zscore)),
  "\n"
)

# 4. Check sample type

cat("\n--- SAMPLE TYPE DISTRIBUTION ---\n")

print(table(immune_scores$sample_type))

# 5. Identify methylation matrix candidate

# Перевіримо назви об'єктів, не припускаючи назву матриці

cat("\n--- POSSIBLE METHYLATION OBJECTS ---\n")

for (nm in matrix_fields) {
  
  x <- integrated[[nm]]
  
  cat(
    nm,
    "| Rows:", nrow(x),
    "| Columns:", ncol(x),
    "| CpG-like rownames:",
    mean(grepl("^cg[0-9]+$", rownames(x))),
    "\n"
  )
}
# ============================================================
# STEP 2. SAMPLE ID MATCHING QC
# Перевірка відповідності зразків метилювання та TIS
# ============================================================

# 1. Use methylation matrix from integrated object

beta_matrix <- integrated$beta_matrix

# 2. Extract sample IDs

beta_ids <- colnames(beta_matrix)

immune_ids <- as.character(
  immune_scores$Sample_ID
)

# 3. Inspect first sample IDs

cat("\n--- METHYLATION SAMPLE IDS ---\n")
print(head(beta_ids, 10))

cat("\n--- IMMUNE SCORE SAMPLE IDS ---\n")
print(head(immune_ids, 10))

# 4. Check duplicates

cat("\n--- DUPLICATE ID CHECK ---\n")

cat(
  "Duplicated methylation IDs:",
  sum(duplicated(beta_ids)),
  "\n"
)

cat(
  "Duplicated immune IDs:",
  sum(duplicated(immune_ids)),
  "\n"
)

# 5. Check exact matching

cat("\n--- EXACT MATCHING ---\n")

cat(
  "Same number of samples:",
  length(beta_ids) == length(immune_ids),
  "\n"
)

cat(
  "Same set of IDs:",
  setequal(beta_ids, immune_ids),
  "\n"
)

cat(
  "Methylation IDs absent in immune scores:",
  length(setdiff(beta_ids, immune_ids)),
  "\n"
)

cat(
  "Immune IDs absent in methylation:",
  length(setdiff(immune_ids, beta_ids)),
  "\n"
)

# 6. Show unmatched IDs, if any

cat("\n--- UNMATCHED METHYLATION IDS ---\n")
print(head(setdiff(beta_ids, immune_ids), 10))

cat("\n--- UNMATCHED IMMUNE IDS ---\n")
print(head(setdiff(immune_ids, beta_ids), 10))

# ============================================================
# STEP 3. ALIGN SAMPLE ORDER
# Вирівнювання порядку зразків
# ============================================================

# 1. Find matching positions

immune_index <- match(
  beta_ids,
  immune_ids
)

# 2. Check matching positions

stopifnot(!anyNA(immune_index))

# 3. Reorder immune scores

immune_scores_aligned <- immune_scores[
  immune_index,
  ,
  drop = FALSE
]

# 4. Verify sample alignment

cat("\n--- SAMPLE ALIGNMENT QC ---\n")

cat(
  "Same sample IDs:",
  identical(
    colnames(beta_matrix),
    immune_scores_aligned$Sample_ID
  ),
  "\n"
)

cat(
  "Number of samples:",
  nrow(immune_scores_aligned),
  "\n"
)

cat(
  "Sample order preserved:",
  identical(
    colnames(beta_matrix),
    immune_scores_aligned$Sample_ID
  ),
  "\n"
)

# 5. Verify phenotype variables

cat("\n--- IMMUNE GROUPS ---\n")

print(
  table(immune_scores_aligned$immune_group)
)

cat("\n--- SAMPLE TYPES ---\n")

print(
  table(immune_scores_aligned$sample_type)
)

# 6. Stop if alignment is incorrect

stopifnot(
  identical(
    colnames(beta_matrix),
    immune_scores_aligned$Sample_ID
  )
)

stopifnot(
  nrow(immune_scores_aligned) == ncol(beta_matrix)
)
# ============================================================
# STEP 4. METHYLATION MATRIX QUALITY CONTROL
# Контроль якості матриці метилювання
# ============================================================

# 1. Matrix dimensions

cat("\n--- MATRIX DIMENSIONS ---\n")

cat("CpG sites:", nrow(beta_matrix), "\n")
cat("Samples:", ncol(beta_matrix), "\n")

# 2. Check data type

cat("\n--- DATA TYPE ---\n")

print(class(beta_matrix))

# 3. Check missing values

cat("\n--- MISSING VALUES ---\n")

cat(
  "NA values:",
  sum(is.na(beta_matrix)),
  "\n"
)

cat(
  "Infinite values:",
  sum(!is.finite(beta_matrix)),
  "\n"
)

# 4. Check beta-value range

cat("\n--- BETA VALUE RANGE ---\n")

print(
  range(beta_matrix, na.rm = TRUE)
)

# 5. Check values outside valid range

cat(
  "Values below 0:",
  sum(beta_matrix < 0, na.rm = TRUE),
  "\n"
)

cat(
  "Values above 1:",
  sum(beta_matrix > 1, na.rm = TRUE),
  "\n"
)

# 6. Check exact boundary values

cat("\n--- BOUNDARY VALUES ---\n")

cat(
  "Values equal to 0:",
  sum(beta_matrix == 0, na.rm = TRUE),
  "\n"
)

cat(
  "Values equal to 1:",
  sum(beta_matrix == 1, na.rm = TRUE),
  "\n"
)

# 7. Check CpG identifiers

cpg_ids <- rownames(beta_matrix)

cat("\n--- CpG IDENTIFIERS ---\n")

cat(
  "Unique CpG IDs:",
  length(unique(cpg_ids)),
  "\n"
)

cat(
  "Duplicated CpG IDs:",
  sum(duplicated(cpg_ids)),
  "\n"
)

cat(
  "Valid cg identifiers:",
  sum(grepl("^cg[0-9]+$", cpg_ids)),
  "\n"
)

cat(
  "Other identifiers:",
  sum(!grepl("^cg[0-9]+$", cpg_ids)),
  "\n"
)

# 8. Check zero-variance CpGs

cpg_variance <- matrixStats::rowVars(
  beta_matrix
)

cat("\n--- VARIANCE QC ---\n")

cat(
  "Zero-variance CpGs:",
  sum(cpg_variance == 0),
  "\n"
)

cat(
  "Non-finite variance:",
  sum(!is.finite(cpg_variance)),
  "\n"
)
# ============================================================
# STEP 5. ПЕРЕВІРКА НЕСТАНДАРТНИХ ІДЕНТИФІКАТОРІВ CpG
# ============================================================

# 1. Виділення ідентифікаторів, які не відповідають формату cg...
non_cpg_ids <- cpg_ids[
  !grepl("^cg[0-9]+$", cpg_ids)
]

# 2. Загальна кількість
cat("\n--- НЕСТАНДАРТНІ ІДЕНТИФІКАТОРИ ---\n")

cat(
  "Кількість:",
  length(non_cpg_ids),
  "\n"
)

# 3. Перші 50 ідентифікаторів
cat("\n--- ПЕРШІ 50 ID ---\n")

print(
  head(non_cpg_ids, 50)
)

# 4. Перевірка типів ідентифікаторів
cat("\n--- СТРУКТУРА ІДЕНТИФІКАТОРІВ ---\n")

id_prefixes <- sub(
  "^([A-Za-z]+).*",
  "\\1",
  non_cpg_ids
)

print(
  sort(
    table(id_prefixes),
    decreasing = TRUE
  )
)

# 5. Перевірка, чи є серед них контрольні зонди
cat("\n--- МОЖЛИВІ КОНТРОЛЬНІ ЗОНДИ ---\n")

control_ids <- non_cpg_ids[
  grepl(
    "control|ctl|ch\\.|rs[0-9]+",
    non_cpg_ids,
    ignore.case = TRUE
  )
]

cat(
  "Кількість потенційних контрольних/інших зондів:",
  length(control_ids),
  "\n"
)

print(
  head(control_ids, 30)
)
# ============================================================
# STEP 6. ДЕТАЛЬНА ПЕРЕВІРКА CH-ЗОНДІВ
# ============================================================

# 1. Перевірка загального формату
cat("\n--- ФОРМАТ CH-ID ---\n")

cat(
  "Відповідають формату ch.1.позиціяF/R:",
  sum(grepl("^ch\\.[0-9]+\\.[0-9]+[FR]$", non_cpg_ids)),
  "\n"
)

# 2. Перевірка суфіксів
cat("\n--- РОЗПОДІЛ СУФІКСІВ ---\n")

print(
  table(
    sub("^.*([FR])$", "\\1", non_cpg_ids)
  )
)

# 3. Перевірка другого елемента ID
cat("\n--- ЕЛЕМЕНТ ПІСЛЯ ch. ---\n")

print(
  sort(
    table(
      sub("^ch\\.([0-9]+)\\..*$", "\\1", non_cpg_ids)
    ),
    decreasing = TRUE
  )
)

# 4. Перевірка ідентифікаторів, які не відповідають
# очікуваному формату
unexpected_ch_ids <- non_cpg_ids[
  !grepl("^ch\\.[0-9]+\\.[0-9]+[FR]$", non_cpg_ids)
]

cat("\n--- НЕОЧІКУВАНІ CH-ID ---\n")

cat(
  "Кількість:",
  length(unexpected_ch_ids),
  "\n"
)

print(
  head(unexpected_ch_ids, 30)
)
