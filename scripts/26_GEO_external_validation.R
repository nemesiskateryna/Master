# ============================================================
# GEO EXTERNAL VALIDATION
# STEP 1: TCGA-DERIVED 140 CpG SIGNATURE PLATFORM OVERLAP
# ============================================================

# 1. Завантажуємо фінальну TCGA-сигнатуру
tcga_signature_file <-
  "results/tables/TCGA_SKCM_ElasticNet_140_CpG_annotated_hg19.csv"

if (!file.exists(tcga_signature_file)) {
  stop("TCGA signature file not found: ", tcga_signature_file)
}

tcga_signature <- read.csv(
  tcga_signature_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"Feature" %in% colnames(tcga_signature)) {
  stop("Column 'Feature' not found in TCGA signature.")
}

tcga_cpg_ids <- unique(
  trimws(as.character(tcga_signature$Feature))
)

tcga_cpg_ids <- tcga_cpg_ids[
  !is.na(tcga_cpg_ids) & tcga_cpg_ids != ""
]

# 2. Завантажуємо інтегрований GEO-об'єкт
geo_file <-
  "data/GEO/processed/GEO_integrated_multiomics_high_confidence.rds"

if (!file.exists(geo_file)) {
  stop("GEO integrated file not found: ", geo_file)
}

geo_integrated <- readRDS(geo_file)

if (!"beta_matrix" %in% names(geo_integrated)) {
  stop("beta_matrix not found in GEO integrated object.")
}

geo_beta <- geo_integrated$beta_matrix

# 3. Перевіряємо базові характеристики
cat("\n--- INPUT OBJECTS ---\n")

cat("TCGA signature CpGs:", length(tcga_cpg_ids), "\n")
cat("GEO beta matrix dimensions:",
    paste(dim(geo_beta), collapse = " x "), "\n")

cat("TCGA CpG IDs unique:",
    !anyDuplicated(tcga_cpg_ids), "\n")

cat("GEO CpG rownames unique:",
    !anyDuplicated(rownames(geo_beta)), "\n")

cat("GEO samples:", ncol(geo_beta), "\n")

# 4. Перевіряємо точне перекриття CpG
geo_cpg_ids <- rownames(geo_beta)

if (is.null(geo_cpg_ids)) {
  stop("GEO beta matrix has no CpG rownames.")
}

shared_cpg_ids <- intersect(
  tcga_cpg_ids,
  geo_cpg_ids
)

tcga_cpg_not_in_geo <- setdiff(
  tcga_cpg_ids,
  geo_cpg_ids
)

# 5. Формуємо таблицю перекриття
cpg_overlap_qc <- data.frame(
  Metric = c(
    "TCGA signature CpGs",
    "GEO EPIC CpGs",
    "Shared CpGs",
    "TCGA CpGs absent from GEO",
    "Overlap percentage"
  ),
  Value = c(
    length(tcga_cpg_ids),
    length(geo_cpg_ids),
    length(shared_cpg_ids),
    length(tcga_cpg_not_in_geo),
    round(
      100 * length(shared_cpg_ids) /
        length(tcga_cpg_ids),
      2
    )
  )
)

cat("\n--- TCGA / GEO CpG OVERLAP ---\n")
print(cpg_overlap_qc, row.names = FALSE)

# 6. Витягуємо доступні TCGA-derived CpG з GEO
if (length(shared_cpg_ids) == 0) {
  stop("No shared CpGs found. Check genome build and CpG identifiers.")
}

geo_signature_beta <- geo_beta[
  shared_cpg_ids,
  ,
  drop = FALSE
]

# 7. Перевіряємо якість доступних CpG у GEO
cpg_quality <- data.frame(
  CpG_ID = rownames(geo_signature_beta),
  N_samples = ncol(geo_signature_beta),
  Missing_values = rowSums(is.na(geo_signature_beta)),
  Non_missing_samples = rowSums(!is.na(geo_signature_beta)),
  Mean_beta = rowMeans(geo_signature_beta, na.rm = TRUE),
  SD_beta = apply(geo_signature_beta, 1, sd, na.rm = TRUE),
  stringsAsFactors = FALSE
)

cpg_quality$Complete_74_samples <-
  cpg_quality$Missing_values == 0

cpg_quality$Variable_beta <-
  is.finite(cpg_quality$SD_beta) &
  cpg_quality$SD_beta > 0

cat("\n--- GEO CpG QUALITY ---\n")

cat("Shared CpGs:", nrow(cpg_quality), "\n")

cat("CpGs complete across 74 samples:",
    sum(cpg_quality$Complete_74_samples), "\n")

cat("CpGs with non-zero variability:",
    sum(cpg_quality$Variable_beta), "\n")

cat("CpGs with zero variability:",
    sum(!cpg_quality$Variable_beta), "\n")

cat("CpGs with missing values:",
    sum(cpg_quality$Missing_values > 0), "\n")

# 8. Перевіряємо порядок та ідентифікатори
cat("\n--- FINAL QC ---\n")

cat("All extracted CpGs belong to TCGA signature:",
    all(rownames(geo_signature_beta) %in% tcga_cpg_ids),
    "\n")

cat("All extracted CpGs belong to GEO matrix:",
    all(rownames(geo_signature_beta) %in% geo_cpg_ids),
    "\n")

cat("GEO sample count preserved:",
    ncol(geo_signature_beta) == 74,
    "\n")

# 9. Зберігаємо результати
dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/models",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  cpg_overlap_qc,
  "results/tables/GEO_TCGA_140_CpG_overlap_QC.csv",
  row.names = FALSE
)

write.csv(
  cpg_quality,
  "results/tables/GEO_TCGA_shared_CpG_quality.csv",
  row.names = FALSE
)

write.csv(
  data.frame(CpG_ID = tcga_cpg_not_in_geo),
  "results/tables/GEO_TCGA_CpGs_not_found_in_GEO.csv",
  row.names = FALSE
)

saveRDS(
  list(
    tcga_signature = tcga_signature,
    shared_cpg_ids = shared_cpg_ids,
    missing_from_geo = tcga_cpg_not_in_geo,
    geo_signature_beta = geo_signature_beta,
    cpg_quality = cpg_quality,
    overlap_qc = cpg_overlap_qc
  ),
  "results/models/GEO_TCGA_140_CpG_overlap.rds"
)

cat("\nSaved files:\n")
cat("- results/tables/GEO_TCGA_140_CpG_overlap_QC.csv\n")
cat("- results/tables/GEO_TCGA_shared_CpG_quality.csv\n")
cat("- results/tables/GEO_TCGA_CpGs_not_found_in_GEO.csv\n")
cat("- results/models/GEO_TCGA_140_CpG_overlap.rds\n")


# ============================================================
# GEO EXTERNAL VALIDATION: TCGA SIGNATURE STRUCTURE AUDIT
# ============================================================

# 1. Завантажуємо збережений результат CpG overlap
overlap_file <- "results/models/GEO_TCGA_140_CpG_overlap.rds"

if (!file.exists(overlap_file)) {
  stop("Overlap RDS file not found: ", overlap_file)
}

overlap_object <- readRDS(overlap_file)

# 2. Перевіряємо компоненти об'єкта
cat("\n--- OVERLAP OBJECT COMPONENTS ---\n")
print(names(overlap_object))

# 3. Витягуємо TCGA-сигнатуру
if (!"tcga_signature" %in% names(overlap_object)) {
  stop("tcga_signature component is missing.")
}

tcga_signature <- overlap_object$tcga_signature

cat("\n--- TCGA SIGNATURE STRUCTURE ---\n")
cat("Class:", class(tcga_signature), "\n")

if (is.data.frame(tcga_signature) || is.matrix(tcga_signature)) {
  
  cat("Dimensions:", dim(tcga_signature), "\n")
  cat("\nColumn names:\n")
  print(colnames(tcga_signature))
  
  cat("\nFirst 10 rows:\n")
  print(
    head(as.data.frame(tcga_signature), 10),
    row.names = FALSE
  )
  
  cat("\nData types:\n")
  print(
    data.frame(
      Column = names(as.data.frame(tcga_signature)),
      Class = vapply(
        as.data.frame(tcga_signature),
        function(x) paste(class(x), collapse = "/"),
        character(1)
      )
    ),
    row.names = FALSE
  )
  
} else if (is.list(tcga_signature)) {
  
  cat("\nList components:\n")
  print(names(tcga_signature))
  
  str(tcga_signature, max.level = 2)
  
} else {
  
  str(tcga_signature)
}

# 4. Перевіряємо спільні CpG
cat("\n--- SHARED CpG CHECK ---\n")

shared_cpg_ids <- overlap_object$shared_cpg_ids

cat("Shared CpGs:", length(shared_cpg_ids), "\n")
cat("Unique shared CpGs:", length(unique(shared_cpg_ids)), "\n")

cat(
  "GEO beta matrix dimensions:",
  dim(overlap_object$geo_signature_beta),
  "\n"
)

cat(
  "All shared CpGs present in GEO beta matrix:",
  all(shared_cpg_ids %in% rownames(
    overlap_object$geo_signature_beta
  )),
  "\n"
)

# 5. Перевіряємо наявність коефіцієнтів або напрямку
if (is.data.frame(tcga_signature)) {
  
  signature_columns <- names(tcga_signature)
  
  possible_columns <- grep(
    "cpg|coefficient|coef|weight|estimate|direction|rho|correlation|gene|symbol",
    signature_columns,
    ignore.case = TRUE,
    value = TRUE
  )
  
  cat("\n--- POTENTIALLY RELEVANT COLUMNS ---\n")
  print(possible_columns)
  
}

cat("\nSignature structure audit completed.\n")


# ============================================================
# GEO EXTERNAL VALIDATION
# TCGA-DERIVED CpG SIGNATURE: DIRECTIONAL VALIDATION
# ============================================================

# 1. Перевірка необхідних об'єктів
required_objects <- c(
  "overlap_object",
  "tcga_signature",
  "geo_immune_scores_complete"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = .GlobalEnv
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Missing objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

# 2. Отримуємо GEO beta matrix та phenotype
geo_beta_signature <- overlap_object$geo_signature_beta

geo_pheno <- geo_immune_scores_complete

if (!all(c("Sample_ID", "Immune_phenotype") %in% names(geo_pheno))) {
  stop("Required phenotype columns are missing.")
}

if (!setequal(colnames(geo_beta_signature), geo_pheno$Sample_ID)) {
  stop("GEO sample IDs do not match between beta matrix and phenotype.")
}

# Відновлюємо однаковий порядок зразків
geo_pheno <- geo_pheno[
  match(colnames(geo_beta_signature), geo_pheno$Sample_ID),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(geo_beta_signature),
    geo_pheno$Sample_ID
  )
)

# 3. Перевіряємо групи
geo_pheno$Immune_phenotype <- factor(
  geo_pheno$Immune_phenotype,
  levels = c("Immune-low", "Immune-high")
)

if (anyNA(geo_pheno$Immune_phenotype)) {
  stop("Unexpected or missing immune phenotype labels.")
}

cat("\n--- GEO GROUP SIZES ---\n")
print(table(geo_pheno$Immune_phenotype))

# 4. Зіставляємо CpG із TCGA-сигнатурою
tcga_shared <- tcga_signature[
  tcga_signature$Feature %in% rownames(geo_beta_signature),
  ,
  drop = FALSE
]

tcga_shared <- tcga_shared[
  match(rownames(geo_beta_signature), tcga_shared$Feature),
  ,
  drop = FALSE
]

if (anyNA(tcga_shared$Feature)) {
  stop("Some shared CpGs could not be matched to TCGA signature.")
}

stopifnot(
  identical(
    tcga_shared$Feature,
    rownames(geo_beta_signature)
  )
)

if (any(!is.finite(tcga_shared$Coefficient))) {
  stop("Non-finite TCGA coefficients detected.")
}

# 5. Перевіряємо відповідність Direction та знаку коефіцієнта
expected_direction <- ifelse(
  tcga_shared$Coefficient > 0,
  "Immune-high",
  ifelse(
    tcga_shared$Coefficient < 0,
    "Immune-low",
    "Zero"
  )
)

direction_consistency <- sum(
  tcga_shared$Direction == expected_direction
)

cat("\n--- TCGA DIRECTION QC ---\n")

cat(
  "Direction consistent with coefficient sign:",
  direction_consistency,
  "of",
  nrow(tcga_shared),
  "\n"
)

if (any(tcga_shared$Coefficient == 0)) {
  stop("Zero TCGA coefficients detected; inspect signature.")
}

if (direction_consistency != nrow(tcga_shared)) {
  stop("TCGA Direction does not consistently match coefficient signs.")
}

# 6. Розраховуємо GEO ефекти для кожного CpG
geo_effects <- lapply(
  seq_len(nrow(geo_beta_signature)),
  function(i) {
    
    cpg <- rownames(geo_beta_signature)[i]
    
    beta_values <- as.numeric(
      geo_beta_signature[i, ]
    )
    
    high_values <- beta_values[
      geo_pheno$Immune_phenotype == "Immune-high"
    ]
    
    low_values <- beta_values[
      geo_pheno$Immune_phenotype == "Immune-low"
    ]
    
    test <- wilcox.test(
      high_values,
      low_values,
      exact = FALSE
    )
    
    mean_high <- mean(high_values)
    mean_low <- mean(low_values)
    
    median_high <- median(high_values)
    median_low <- median(low_values)
    
    delta_beta <- mean_high - mean_low
    
    data.frame(
      CpG_ID = cpg,
      TCGA_Coefficient = tcga_shared$Coefficient[i],
      TCGA_Direction = tcga_shared$Direction[i],
      GEO_Mean_Immune_high = mean_high,
      GEO_Mean_Immune_low = mean_low,
      GEO_Median_Immune_high = median_high,
      GEO_Median_Immune_low = median_low,
      GEO_Delta_beta = delta_beta,
      GEO_Wilcoxon_p = test$p.value,
      stringsAsFactors = FALSE
    )
  }
)

geo_cpg_validation <- do.call(rbind, geo_effects)

# 7. Корекція множинного тестування
geo_cpg_validation$GEO_FDR_BH <- p.adjust(
  geo_cpg_validation$GEO_Wilcoxon_p,
  method = "BH"
)

geo_cpg_validation$GEO_Direction <- ifelse(
  geo_cpg_validation$GEO_Delta_beta > 0,
  "Immune-high",
  ifelse(
    geo_cpg_validation$GEO_Delta_beta < 0,
    "Immune-low",
    "No difference"
  )
)

# 8. Перевірка напрямку між TCGA та GEO
geo_cpg_validation$Direction_concordant <-
  geo_cpg_validation$TCGA_Direction ==
  geo_cpg_validation$GEO_Direction

geo_cpg_validation$GEO_significant_FDR_005 <-
  geo_cpg_validation$GEO_FDR_BH < 0.05

geo_cpg_validation$Significant_and_concordant <-
  geo_cpg_validation$GEO_significant_FDR_005 &
  geo_cpg_validation$Direction_concordant

# 9. Загальна статистика
n_total <- nrow(geo_cpg_validation)

n_concordant <- sum(
  geo_cpg_validation$Direction_concordant,
  na.rm = TRUE
)

n_significant <- sum(
  geo_cpg_validation$GEO_significant_FDR_005,
  na.rm = TRUE
)

n_significant_concordant <- sum(
  geo_cpg_validation$Significant_and_concordant,
  na.rm = TRUE
)

cat("\n--- DIRECTIONAL VALIDATION SUMMARY ---\n")

cat("Shared CpGs tested:", n_total, "\n")

cat(
  "Directionally concordant:",
  n_concordant,
  "/",
  n_total,
  "\n"
)

cat(
  "Concordance percentage:",
  round(100 * n_concordant / n_total, 2),
  "%\n"
)

cat(
  "GEO significant CpGs (FDR < 0.05):",
  n_significant,
  "\n"
)

cat(
  "Significant and directionally concordant:",
  n_significant_concordant,
  "\n"
)

# 10. Кореляція TCGA coefficient vs GEO Delta beta
coefficient_correlation <- cor.test(
  geo_cpg_validation$TCGA_Coefficient,
  geo_cpg_validation$GEO_Delta_beta,
  method = "spearman",
  exact = FALSE
)

cat("\n--- TCGA COEFFICIENT VS GEO DELTA BETA ---\n")

print(
  data.frame(
    Spearman_rho = unname(
      coefficient_correlation$estimate
    ),
    P_value = coefficient_correlation$p.value,
    N = n_total
  ),
  row.names = FALSE
)

# 11. Перевіряємо розподіл напрямків
cat("\n--- DIRECTIONAL CONTINGENCY TABLE ---\n")

print(
  table(
    TCGA = geo_cpg_validation$TCGA_Direction,
    GEO = geo_cpg_validation$GEO_Direction
  )
)

# 12. Топ CpG за GEO FDR
geo_cpg_validation <- geo_cpg_validation[
  order(
    geo_cpg_validation$GEO_FDR_BH,
    geo_cpg_validation$GEO_Wilcoxon_p
  ),
]

cat("\n--- TOP 15 GEO CpGs ---\n")

print(
  head(
    geo_cpg_validation[
      ,
      c(
        "CpG_ID",
        "TCGA_Coefficient",
        "TCGA_Direction",
        "GEO_Delta_beta",
        "GEO_Direction",
        "GEO_Wilcoxon_p",
        "GEO_FDR_BH",
        "Direction_concordant"
      )
    ],
    15
  ),
  row.names = FALSE
)

# 13. Збереження
write.csv(
  geo_cpg_validation,
  "results/tables/GEO_TCGA_CpG_directional_validation.csv",
  row.names = FALSE
)

saveRDS(
  list(
    validation_results = geo_cpg_validation,
    coefficient_correlation = coefficient_correlation,
    summary = data.frame(
      Total_CpGs = n_total,
      Directionally_concordant = n_concordant,
      Concordance_percent = 100 * n_concordant / n_total,
      GEO_significant_FDR05 = n_significant,
      Significant_and_concordant = n_significant_concordant
    ),
    phenotype_definition = "GEO TIS-17 median split",
    test = "Wilcoxon rank-sum",
    multiple_testing = "Benjamini-Hochberg"
  ),
  "results/models/GEO_TCGA_CpG_directional_validation.rds"
)

cat("\nSaved directional validation results.\n")



# ============================================================
# GEO EXTERNAL VALIDATION
# SAMPLE-LEVEL PERMUTATION TEST
# ============================================================

set.seed(20261106)

# 1. Підготовка даних
beta_matrix_perm <- overlap_object$geo_signature_beta

phenotype_perm <- geo_pheno$Immune_phenotype

tcga_coef_perm <- tcga_shared$Coefficient

tcga_sign_perm <- sign(tcga_coef_perm)

n_samples_perm <- ncol(beta_matrix_perm)
n_cpg_perm <- nrow(beta_matrix_perm)

n_permutations <- 9999

# Перевірки
stopifnot(
  n_samples_perm == 74,
  n_cpg_perm == 121,
  length(phenotype_perm) == n_samples_perm,
  length(tcga_sign_perm) == n_cpg_perm
)

# 2. Функція розрахунку GEO Delta beta
calculate_delta_beta <- function(beta_matrix, labels) {
  
  high <- labels == "Immune-high"
  low <- labels == "Immune-low"
  
  rowMeans(beta_matrix[, high, drop = FALSE]) -
    rowMeans(beta_matrix[, low, drop = FALSE])
}

# 3. Спостережувані статистики
observed_delta <- calculate_delta_beta(
  beta_matrix_perm,
  phenotype_perm
)

observed_concordance <- mean(
  sign(observed_delta) == tcga_sign_perm
)

observed_rho <- cor(
  tcga_coef_perm,
  observed_delta,
  method = "spearman"
)

cat("\n--- OBSERVED STATISTICS ---\n")

cat(
  "Directional concordance:",
  round(observed_concordance, 4),
  "\n"
)

cat(
  "Spearman rho:",
  round(observed_rho, 4),
  "\n"
)

# 4. Пермутації групових міток
permutation_concordance <- numeric(n_permutations)
permutation_rho <- numeric(n_permutations)

for (i in seq_len(n_permutations)) {
  
  permuted_labels <- sample(phenotype_perm)
  
  permuted_delta <- calculate_delta_beta(
    beta_matrix_perm,
    permuted_labels
  )
  
  permutation_concordance[i] <- mean(
    sign(permuted_delta) == tcga_sign_perm
  )
  
  permutation_rho[i] <- cor(
    tcga_coef_perm,
    permuted_delta,
    method = "spearman"
  )
}

# 5. Односторонні permutation P-values
p_concordance <- (
  1 + sum(
    permutation_concordance >= observed_concordance
  )
) / (n_permutations + 1)

p_rho <- (
  1 + sum(
    permutation_rho >= observed_rho
  )
) / (n_permutations + 1)

# 6. Підсумкова таблиця
permutation_summary <- data.frame(
  Statistic = c(
    "Directional concordance",
    "Spearman correlation"
  ),
  Observed = c(
    observed_concordance,
    observed_rho
  ),
  Permutation_null_mean = c(
    mean(permutation_concordance),
    mean(permutation_rho)
  ),
  Permutation_null_SD = c(
    sd(permutation_concordance),
    sd(permutation_rho)
  ),
  Permutation_P_value = c(
    p_concordance,
    p_rho
  ),
  N_permutations = n_permutations
)

cat("\n--- PERMUTATION TEST RESULTS ---\n")
print(permutation_summary, row.names = FALSE)

# 7. Квантилі нульових розподілів
cat("\n--- NULL DISTRIBUTION QUANTILES ---\n")

print(
  data.frame(
    Statistic = c(
      "Directional concordance",
      "Spearman correlation"
    ),
    Null_Q025 = c(
      quantile(permutation_concordance, 0.025),
      quantile(permutation_rho, 0.025)
    ),
    Null_Median = c(
      median(permutation_concordance),
      median(permutation_rho)
    ),
    Null_Q975 = c(
      quantile(permutation_concordance, 0.975),
      quantile(permutation_rho, 0.975)
    )
  ),
  row.names = FALSE
)

# 8. Збереження результатів
saveRDS(
  list(
    summary = permutation_summary,
    observed_delta_beta = observed_delta,
    observed_concordance = observed_concordance,
    observed_spearman_rho = observed_rho,
    permutation_concordance = permutation_concordance,
    permutation_rho = permutation_rho,
    n_permutations = n_permutations,
    seed = 20261106,
    method = paste(
      "Sample-level permutation of immune phenotype labels;",
      "group sizes preserved;",
      "CpG correlation structure retained."
    )
  ),
  "results/models/GEO_TCGA_CpG_permutation_validation.rds"
)

write.csv(
  permutation_summary,
  "results/tables/GEO_TCGA_CpG_permutation_summary.csv",
  row.names = FALSE
)

cat("\nPermutation validation results saved.\n")



# ============================================================
# GEO EXTERNAL VALIDATION — FINAL RESULTS INVENTORY
# ============================================================

cat("\n============================================\n")
cat("GEO EXTERNAL VALIDATION: FILE INVENTORY\n")
cat("============================================\n")

required_files <- c(
  
  # Immune phenotype
  "results/tables/GEO_TIS_17gene_immune_scores.csv",
  "results/tables/GEO_auxiliary_signature_coverage.csv",
  "results/tables/GEO_TIS17_auxiliary_correlations.csv",
  "results/tables/GEO_complete_immune_scores.csv",
  
  # TCGA-GEO CpG overlap
  "results/tables/GEO_TCGA_140_CpG_overlap_QC.csv",
  "results/tables/GEO_TCGA_shared_CpG_quality.csv",
  "results/tables/GEO_TCGA_CpGs_not_found_in_GEO.csv",
  
  # Directional validation
  "results/tables/GEO_TCGA_CpG_directional_validation.csv",
  
  # Permutation validation
  "results/tables/GEO_TCGA_CpG_permutation_summary.csv",
  
  # RDS objects
  "results/models/GEO_TIS_17gene_immune_scores.rds",
  "results/models/GEO_complete_immune_scores.rds",
  "results/models/GEO_TCGA_140_CpG_overlap.rds",
  "results/models/GEO_TCGA_CpG_directional_validation.rds",
  "results/models/GEO_TCGA_CpG_permutation_validation.rds"
)

file_inventory <- data.frame(
  File = required_files,
  Exists = file.exists(required_files),
  Size_KB = ifelse(
    file.exists(required_files),
    round(file.info(required_files)$size / 1024, 2),
    NA
  )
)

print(file_inventory, row.names = FALSE)

cat("\nTotal files:", nrow(file_inventory), "\n")
cat("Existing files:", sum(file_inventory$Exists), "\n")
cat("Missing files:", sum(!file_inventory$Exists), "\n")

# ------------------------------------------------------------
# CHECK IMPORTANT CSV STRUCTURES
# ------------------------------------------------------------

csv_files <- required_files[
  grepl("\\.csv$", required_files)
]

cat("\n============================================\n")
cat("CSV STRUCTURE CHECK\n")
cat("============================================\n")

for (f in csv_files) {
  
  if (file.exists(f)) {
    
    x <- tryCatch(
      read.csv(f, check.names = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(x)) {
      
      cat("\nFile:", f, "\n")
      cat("Dimensions:", nrow(x), "x", ncol(x), "\n")
      cat(
        "Columns:",
        paste(names(x), collapse = " | "),
        "\n"
      )
      
    } else {
      
      cat("\nCould not read:", f, "\n")
    }
  }
}

# ------------------------------------------------------------
# CHECK SAVED RDS OBJECT STRUCTURES
# ------------------------------------------------------------

rds_files <- required_files[
  grepl("\\.rds$", required_files)
]

cat("\n============================================\n")
cat("RDS OBJECT STRUCTURE CHECK\n")
cat("============================================\n")

for (f in rds_files) {
  
  if (file.exists(f)) {
    
    obj <- tryCatch(
      readRDS(f),
      error = function(e) NULL
    )
    
    cat("\nFile:", f, "\n")
    
    if (is.null(obj)) {
      
      cat("Could not read RDS object\n")
      
    } else {
      
      cat("Class:", paste(class(obj), collapse = ", "), "\n")
      
      if (is.list(obj)) {
        cat(
          "Components:",
          paste(names(obj), collapse = " | "),
          "\n"
        )
      } else {
        cat("Dimensions:", paste(dim(obj), collapse = " x "), "\n")
      }
    }
  }
}

cat("\n============================================\n")
cat("INVENTORY COMPLETED\n")
cat("============================================\n")





# ============================================================
# GEO EXTERNAL VALIDATION — FINAL INTEGRATED SUMMARY
# ============================================================

cat("\n============================================\n")
cat("GEO EXTERNAL VALIDATION: FINAL SUMMARY\n")
cat("============================================\n")

# ------------------------------------------------------------
# 1. Завантаження результатів
# ------------------------------------------------------------

tis_scores <- read.csv(
  "results/tables/GEO_TIS_17gene_immune_scores.csv",
  stringsAsFactors = FALSE
)

complete_scores <- read.csv(
  "results/tables/GEO_complete_immune_scores.csv",
  stringsAsFactors = FALSE
)

aux_coverage <- read.csv(
  "results/tables/GEO_auxiliary_signature_coverage.csv",
  stringsAsFactors = FALSE
)

aux_correlations <- read.csv(
  "results/tables/GEO_TIS17_auxiliary_correlations.csv",
  stringsAsFactors = FALSE
)

overlap_qc <- read.csv(
  "results/tables/GEO_TCGA_140_CpG_overlap_QC.csv",
  stringsAsFactors = FALSE
)

cpg_quality <- read.csv(
  "results/tables/GEO_TCGA_shared_CpG_quality.csv",
  stringsAsFactors = FALSE
)

missing_cpg <- read.csv(
  "results/tables/GEO_TCGA_CpGs_not_found_in_GEO.csv",
  stringsAsFactors = FALSE
)

directional <- read.csv(
  "results/tables/GEO_TCGA_CpG_directional_validation.csv",
  stringsAsFactors = FALSE
)

permutation <- read.csv(
  "results/tables/GEO_TCGA_CpG_permutation_summary.csv",
  stringsAsFactors = FALSE
)

overlap_object <- readRDS(
  "results/models/GEO_TCGA_140_CpG_overlap.rds"
)

permutation_object <- readRDS(
  "results/models/GEO_TCGA_CpG_permutation_validation.rds"
)

# ------------------------------------------------------------
# 2. Основні характеристики когорти
# ------------------------------------------------------------

n_samples <- nrow(tis_scores)

n_high <- sum(
  tis_scores$Immune_phenotype == "Immune-high",
  na.rm = TRUE
)

n_low <- sum(
  tis_scores$Immune_phenotype == "Immune-low",
  na.rm = TRUE
)

n_tcga_cpg <- nrow(overlap_object$tcga_signature)

n_shared_cpg <- length(overlap_object$shared_cpg_ids)

n_missing_cpg <- length(overlap_object$missing_from_geo)

overlap_percent <- 100 * n_shared_cpg / n_tcga_cpg

n_complete_cpg <- sum(
  cpg_quality$Complete_74_samples,
  na.rm = TRUE
)

n_variable_cpg <- sum(
  cpg_quality$Variable_beta,
  na.rm = TRUE
)

# ------------------------------------------------------------
# 3. Результати напрямкової валідації
# ------------------------------------------------------------

n_tested <- nrow(directional)

n_concordant <- sum(
  directional$Direction_concordant,
  na.rm = TRUE
)

concordance_percent <- 100 * n_concordant / n_tested

n_significant <- sum(
  directional$GEO_significant_FDR_005,
  na.rm = TRUE
)

n_significant_concordant <- sum(
  directional$Significant_and_concordant,
  na.rm = TRUE
)

# ------------------------------------------------------------
# 4. Результати permutation test
# ------------------------------------------------------------

concordance_perm <- permutation[
  permutation$Statistic == "Directional concordance",
]

rho_perm <- permutation[
  permutation$Statistic == "Spearman correlation",
]

stopifnot(
  nrow(concordance_perm) == 1,
  nrow(rho_perm) == 1
)

# ------------------------------------------------------------
# 5. Формування інтегрованої таблиці
# ------------------------------------------------------------

final_summary <- data.frame(
  Section = character(),
  Metric = character(),
  Value = character(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

add_metric <- function(section, metric, value, interpretation) {
  
  final_summary <<- rbind(
    final_summary,
    data.frame(
      Section = section,
      Metric = metric,
      Value = as.character(value),
      Interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  )
}

# A. GEO cohort

add_metric(
  "GEO cohort",
  "Total matched samples",
  n_samples,
  "Final analytical GEO cohort"
)

add_metric(
  "GEO cohort",
  "Immune-high samples",
  n_high,
  "Samples classified using the GEO TIS-17gene median"
)

add_metric(
  "GEO cohort",
  "Immune-low samples",
  n_low,
  "Samples classified using the GEO TIS-17gene median"
)

add_metric(
  "GEO cohort",
  "Immune phenotype signature",
  "TIS-17gene",
  "17 of 18 original genes available; IDO1 unavailable"
)

# B. Auxiliary immune signatures

for (i in seq_len(nrow(aux_coverage))) {
  
  add_metric(
    "Auxiliary immune signatures",
    paste(
      aux_coverage$Signature[i],
      "gene coverage"
    ),
    paste0(
      aux_coverage$Available_genes[i],
      "/",
      aux_coverage$Total_genes[i]
    ),
    paste(
      "Genes used:",
      aux_coverage$Genes_used[i]
    )
  )
}

for (i in seq_len(nrow(aux_correlations))) {
  
  add_metric(
    "Auxiliary immune signatures",
    paste(
      aux_correlations$Signature[i],
      "Spearman correlation with TIS-17gene"
    ),
    paste0(
      "rho = ",
      round(aux_correlations$Spearman_rho[i], 3),
      "; P = ",
      format(
        aux_correlations$P_value[i],
        scientific = TRUE,
        digits = 3
      ),
      "; FDR = ",
      format(
        aux_correlations$FDR_BH[i],
        scientific = TRUE,
        digits = 3
      )
    ),
    "Convergent association; signatures are not independent"
  )
}

# C. TCGA-GEO CpG overlap

add_metric(
  "CpG overlap",
  "Original TCGA CpG signature",
  n_tcga_cpg,
  "Number of CpGs in the TCGA-derived signature"
)

add_metric(
  "CpG overlap",
  "Shared CpGs",
  n_shared_cpg,
  "CpGs present in both TCGA signature and GEO EPIC matrix"
)

add_metric(
  "CpG overlap",
  "Overlap percentage",
  paste0(round(overlap_percent, 2), "%"),
  "Proportion of TCGA signature represented in GEO"
)

add_metric(
  "CpG overlap",
  "CpGs missing from GEO",
  n_missing_cpg,
  "TCGA signature CpGs unavailable in GEO"
)

add_metric(
  "CpG overlap",
  "CpGs complete across 74 samples",
  n_complete_cpg,
  "Shared CpGs without missing beta values"
)

add_metric(
  "CpG overlap",
  "CpGs with non-zero variability",
  n_variable_cpg,
  "Shared CpGs with variable beta values"
)

# D. Individual CpG validation

add_metric(
  "Directional validation",
  "CpGs tested",
  n_tested,
  "Shared CpGs tested between GEO immune groups"
)

add_metric(
  "Directional validation",
  "Directionally concordant CpGs",
  paste0(
    n_concordant,
    "/",
    n_tested,
    " (",
    round(concordance_percent, 2),
    "%)"
  ),
  "GEO delta-beta direction agrees with TCGA coefficient direction"
)

add_metric(
  "Directional validation",
  "CpGs significant after BH correction",
  n_significant,
  "Number with GEO FDR < 0.05"
)

add_metric(
  "Directional validation",
  "Significant and concordant CpGs",
  n_significant_concordant,
  "CpGs both significant after BH correction and directionally concordant"
)

# E. Global permutation validation

add_metric(
  "Permutation validation",
  "Observed directional concordance",
  round(permutation$Observed[
    permutation$Statistic == "Directional concordance"
  ], 4),
  "Observed proportion of directionally concordant CpGs"
)

add_metric(
  "Permutation validation",
  "Directional concordance permutation P",
  round(concordance_perm$Permutation_P_value, 4),
  "Not significant at alpha = 0.05"
)

add_metric(
  "Permutation validation",
  "Observed Spearman correlation",
  round(permutation$Observed[
    permutation$Statistic == "Spearman correlation"
  ], 4),
  "Correlation between TCGA coefficients and GEO delta-beta"
)

add_metric(
  "Permutation validation",
  "Spearman permutation P",
  round(rho_perm$Permutation_P_value, 4),
  "Borderline trend; not significant at alpha = 0.05"
)

add_metric(
  "Permutation validation",
  "Number of permutations",
  permutation$N_permutations[1],
  "Sample-level permutation preserving group sizes"
)

# ------------------------------------------------------------
# 6. Перевірка цілісності підсумкової таблиці
# ------------------------------------------------------------

cat("\n--- FINAL SUMMARY TABLE ---\n")

print(
  final_summary,
  row.names = FALSE
)

cat("\n--- FINAL QC ---\n")

final_qc <- c(
  all_74_samples_present = n_samples == 74,
  immune_groups_sum_to_74 = n_high + n_low == 74,
  tcga_signature_has_140_cpgs = n_tcga_cpg == 140,
  shared_cpg_count_is_121 = n_shared_cpg == 121,
  missing_cpg_count_is_19 = n_missing_cpg == 19,
  all_shared_cpgs_complete = n_complete_cpg == 121,
  all_shared_cpgs_variable = n_variable_cpg == 121,
  all_121_cpgs_tested = n_tested == 121,
  permutation_has_9999_iterations =
    all(permutation$N_permutations == 9999)
)

print(final_qc)

cat(
  "\nAll final QC checks passed:",
  all(final_qc),
  "\n"
)

stopifnot(all(final_qc))

# ------------------------------------------------------------
# 7. Збереження результатів
# ------------------------------------------------------------

write.csv(
  final_summary,
  "results/tables/GEO_external_validation_final_summary.csv",
  row.names = FALSE
)

saveRDS(
  list(
    final_summary = final_summary,
    cohort = list(
      n_samples = n_samples,
      n_immune_high = n_high,
      n_immune_low = n_low
    ),
    cpg_overlap = list(
      original_tcga_cpgs = n_tcga_cpg,
      shared_cpgs = n_shared_cpg,
      missing_cpgs = n_missing_cpg,
      overlap_percent = overlap_percent
    ),
    directional_validation = list(
      n_tested = n_tested,
      n_concordant = n_concordant,
      concordance_percent = concordance_percent,
      n_significant = n_significant
    ),
    permutation_validation = permutation,
    qc = final_qc
  ),
  "results/models/GEO_external_validation_final.rds"
)

cat("\n--- SAVED FILES ---\n")

final_files <- c(
  "results/tables/GEO_external_validation_final_summary.csv",
  "results/models/GEO_external_validation_final.rds"
)

print(data.frame(
  File = final_files,
  Exists = file.exists(final_files)
))

cat("\nGEO EXTERNAL VALIDATION SUMMARY COMPLETED.\n")










