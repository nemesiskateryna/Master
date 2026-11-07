
# ============================================================
# 13_TCGA_SKCM_integrated_multiomics.R
#
# TCGA-SKCM: ФОРМУВАННЯ ІНТЕГРОВАНОЇ MULTI-OMICS КОГОРТИ
#
# Мета:
# 1. Завантажити попередньо оброблені дані метилювання,
#    RNA-seq та клінічну інформацію.
# 2. Перевірити відповідність ідентифікаторів зразків.
# 3. Провести контроль якості клінічних даних.
# 4. Підготувати аналітичні змінні загального виживання.
# 5. Сформувати узгоджені матриці multi-omics.
# 6. Зберегти інтегровану когорту та результати аудиту.
#
# Когорта: TCGA-SKCM
# Призначення: Discovery cohort
# ============================================================


# ============================================================
# 1. ПІДГОТОВКА РОБОЧОГО СЕРЕДОВИЩА
# ============================================================

# Очищення робочого середовища
rm(list = ls())

# Встановлення кореневої директорії проєкту
project_dir <- getwd()

# Перевірка необхідних пакетів
required_packages <- c(
  "SummarizedExperiment"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Відсутні необхідні пакети: ",
    paste(missing_packages, collapse = ", ")
  )
}

# Визначення директорій
raw_dir <- file.path(
  project_dir,
  "data/TCGA/raw_saved"
)

processed_dir <- file.path(
  project_dir,
  "data/TCGA/processed"
)

results_tables_dir <- file.path(
  project_dir,
  "results/tables"
)

# Створення директорій, якщо вони відсутні
dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. ЗАВАНТАЖЕННЯ ПОПЕРЕДНЬО ОБРОБЛЕНИХ ДАНИХ
# ============================================================

# Матриця DNA-метилювання (beta-values)
tcga_beta <- readRDS(
  file.path(
    raw_dir,
    "beta_matrix_TCGA_SKCM_imputed.rds"
  )
)

# Матриця нормалізованої експресії RNA-seq (VST)
tcga_rna_vst <- readRDS(
  file.path(
    raw_dir,
    "RNAseq_TCGA_SKCM_VST.rds"
  )
)

# Об'єкт DESeq2 із метаданими RNA-seq
tcga_dds <- readRDS(
  file.path(
    raw_dir,
    "DESeq2_TCGA_SKCM.rds"
  )
)

# Очищені клінічні дані
clinical <- readRDS(
  file.path(
    raw_dir,
    "clinical_data_TCGA_SKCM_clean.rds"
  )
)

# Перевірка структури завантажених об'єктів
cat("\n--- РОЗМІРИ ВХІДНИХ ДАНИХ ---\n")

cat(
  "Метилювання:",
  nrow(tcga_beta), "CpG x",
  ncol(tcga_beta), "зразків\n"
)

cat(
  "RNA-seq:",
  nrow(tcga_rna_vst), "генів x",
  ncol(tcga_rna_vst), "зразків\n"
)

cat(
  "Клінічні дані:",
  nrow(clinical), "записів x",
  ncol(clinical), "змінних\n"
)


# ============================================================
# 3. ПІДГОТОВКА ІДЕНТИФІКАТОРІВ ЗРАЗКІВ
# ============================================================

# Повні TCGA barcode зразків метилювання
methylation_ids <- colnames(tcga_beta)

# Метадані RNA-seq
rna_coldata <- as.data.frame(
  SummarizedExperiment::colData(tcga_dds)
)

# Повні barcode RNA-seq зразків
rna_full_ids <- as.character(
  rna_coldata$barcode
)

# Ідентифікатори колонок об'єкта DESeq2
rna_dds_ids <- colnames(
  SummarizedExperiment::assay(tcga_dds)
)

# Перевірка відповідності RNA-seq матриці та метаданих
stopifnot(
  length(rna_full_ids) == length(rna_dds_ids),
  identical(
    rna_dds_ids,
    rownames(rna_coldata)
  )
)

# Перевірка унікальності barcode
stopifnot(
  !anyNA(methylation_ids),
  !anyNA(rna_full_ids),
  !anyDuplicated(methylation_ids),
  !anyDuplicated(rna_full_ids)
)

# Виділення 12-символьного ідентифікатора пацієнта
meth_patient_ids <- substr(
  methylation_ids,
  1,
  12
)

rna_patient_ids <- substr(
  rna_full_ids,
  1,
  12
)

# Виділення типу зразка:
# 01 — первинна пухлина
# 06 — метастатичний зразок
meth_sample_type <- substr(
  methylation_ids,
  14,
  15
)

rna_sample_type <- substr(
  rna_full_ids,
  14,
  15
)


# ============================================================
# 4. ФОРМУВАННЯ ТА ЗІСТАВЛЕННЯ МЕТАДАНИХ
# ============================================================

# Метадані метилювання
meth_matching <- data.frame(
  Methylation_barcode = methylation_ids,
  Patient_ID = meth_patient_ids,
  Methylation_sample_type = meth_sample_type,
  stringsAsFactors = FALSE
)

# Метадані RNA-seq
rna_matching <- data.frame(
  RNA_matrix_ID = rna_dds_ids,
  RNA_barcode = rna_full_ids,
  Patient_ID = rna_patient_ids,
  RNA_sample_type = rna_sample_type,
  stringsAsFactors = FALSE
)

# Точне зіставлення за повним TCGA barcode
tcga_matching <- merge(
  meth_matching,
  rna_matching,
  by.x = "Methylation_barcode",
  by.y = "RNA_barcode",
  all = FALSE,
  sort = FALSE
)

# Зберігаємо порядок RNA-seq зразків
tcga_matching <- tcga_matching[
  match(
    rna_full_ids,
    tcga_matching$Methylation_barcode
  ),
  ,
  drop = FALSE
]

# Видалення можливих незнайдених записів
tcga_matching <- tcga_matching[
  !is.na(tcga_matching$Methylation_barcode),
  ,
  drop = FALSE
]

# Повні barcode інтегрованої когорти
matched_ids <- as.character(
  tcga_matching$Methylation_barcode
)

# Контроль повноти зіставлення
cat("\n--- ЗІСТАВЛЕННЯ МЕТИЛЮВАННЯ ТА RNA-seq ---\n")

cat(
  "Зразків метилювання:",
  length(methylation_ids), "\n"
)

cat(
  "Зразків RNA-seq:",
  length(rna_full_ids), "\n"
)

cat(
  "Спільних зразків:",
  length(matched_ids), "\n"
)

cat(
  "Зразків без відповідності:",
  length(rna_full_ids) - length(matched_ids),
  "\n"
)

# Перевірка точного збігу barcode
stopifnot(
  identical(
    matched_ids,
    rna_full_ids[
      rna_full_ids %in% methylation_ids
    ]
  )
)


# ============================================================
# 5. АНАЛІЗ ТИПІВ ЗРАЗКІВ
# ============================================================

matched_sample_types <- substr(
  matched_ids,
  14,
  15
)

sample_type_summary <- as.data.frame(
  table(matched_sample_types)
)

colnames(sample_type_summary) <- c(
  "Тип_зразка",
  "Кількість"
)

cat("\n--- РОЗПОДІЛ ТИПІВ ЗРАЗКІВ ---\n")

print(
  sample_type_summary,
  row.names = FALSE
)

# Перевірка очікуваної структури когорти
stopifnot(
  length(matched_ids) == 460
)


# ============================================================
# 6. КОНТРОЛЬ ЯКОСТІ КЛІНІЧНИХ ДАНИХ
# ============================================================

# Перевірка наявності необхідних змінних
required_clinical_variables <- c(
  "barcode",
  "patient",
  "age_at_index",
  "sex_at_birth",
  "vital_status",
  "follow_ups_disease_response",
  "days_to_death",
  "days_to_last_follow_up",
  "overall_survival_time",
  "overall_survival_event"
)

stopifnot(
  all(
    required_clinical_variables %in%
      names(clinical)
  )
)

# Підрахунок пропущених значень
missing_summary <- data.frame(
  Variable = names(clinical),
  Missing = vapply(
    clinical,
    function(x) sum(is.na(x)),
    integer(1)
  ),
  Missing_percent = vapply(
    clinical,
    function(x) mean(is.na(x)) * 100,
    numeric(1)
  )
)

cat("\n--- ПРОПУЩЕНІ КЛІНІЧНІ ЗНАЧЕННЯ ---\n")

print(
  missing_summary,
  row.names = FALSE
)

# Перевірка статусу виживання
cat("\n--- СТАТУС ВИЖИВАННЯ ---\n")

print(
  table(
    clinical$vital_status,
    useNA = "ifany"
  )
)

# Перевірка подій загального виживання
cat("\n--- ПОДІЇ ЗАГАЛЬНОГО ВИЖИВАННЯ ---\n")

print(
  table(
    clinical$overall_survival_event,
    useNA = "ifany"
  )
)

# Перевірка клінічної відповіді під час спостереження
cat("\n--- ВІДПОВІДЬ ЗА ДАНИМИ FOLLOW-UP ---\n")

print(
  table(
    clinical$follow_ups_disease_response,
    useNA = "ifany"
  )
)


# ============================================================
# 7. ПЕРЕВІРКА ЧАСУ ЗАГАЛЬНОГО ВИЖИВАННЯ
# ============================================================

# Виявлення від'ємних значень
negative_survival <- clinical[
  !is.na(clinical$overall_survival_time) &
    clinical$overall_survival_time < 0,
  ,
  drop = FALSE
]

cat("\n--- ВІД'ЄМНИЙ ЧАС ВИЖИВАННЯ ---\n")

print(
  negative_survival[
    ,
    c(
      "patient",
      "vital_status",
      "overall_survival_time",
      "overall_survival_event"
    ),
    drop = FALSE
  ]
)

# Виявлення нульових значень
zero_survival <- clinical[
  !is.na(clinical$overall_survival_time) &
    clinical$overall_survival_time == 0,
  ,
  drop = FALSE
]

cat("\n--- НУЛЬОВИЙ ЧАС ВИЖИВАННЯ ---\n")

print(
  zero_survival[
    ,
    c(
      "patient",
      "vital_status",
      "overall_survival_time",
      "overall_survival_event"
    ),
    drop = FALSE
  ]
)

# Перевірка нульового часу до смерті
zero_death_time <- clinical[
  !is.na(clinical$days_to_death) &
    clinical$days_to_death == 0,
  ,
  drop = FALSE
]

cat("\n--- НУЛЬОВИЙ ЧАС ДО СМЕРТІ ---\n")

print(
  zero_death_time[
    ,
    c(
      "patient",
      "vital_status",
      "days_to_death"
    ),
    drop = FALSE
  ]
)

# Перевірка відповідності розрахованого часу
# первинним клінічним змінним
expected_survival_time <- ifelse(
  clinical$vital_status == "Dead",
  clinical$days_to_death,
  clinical$days_to_last_follow_up
)

survival_difference <-
  clinical$overall_survival_time -
  expected_survival_time

cat("\n--- ПЕРЕВІРКА УЗГОДЖЕНОСТІ OS ---\n")

cat(
  "Неузгоджених записів:",
  sum(
    !is.na(survival_difference) &
      survival_difference != 0
  ),
  "\n"
)

stopifnot(
  all(
    survival_difference == 0 |
      is.na(survival_difference)
  )
)


# ============================================================
# 8. ФОРМУВАННЯ АНАЛІТИЧНИХ ЗМІННИХ ВИЖИВАННЯ
# ============================================================

# Створення окремої аналітичної копії
clinical_survival <- clinical

# Збереження вихідних значень
# у вихідному наборі клінічних даних
clinical_survival$OS_time_analysis <-
  clinical_survival$overall_survival_time

clinical_survival$OS_event_analysis <-
  clinical_survival$overall_survival_event

# Від'ємні значення часу замінюємо на NA
# лише в аналітичній змінній
clinical_survival$OS_time_analysis[
  clinical_survival$OS_time_analysis < 0
] <- NA

# Окремо позначаємо нульовий час виживання
clinical_survival$OS_zero_time_flag <-
  !is.na(clinical_survival$OS_time_analysis) &
  clinical_survival$OS_time_analysis == 0

# Перевірка кількості подій
cat("\n--- АНАЛІТИЧНІ ЗМІННІ OS ---\n")

cat(
  "Усього записів:",
  nrow(clinical_survival), "\n"
)

cat(
  "Пропущений OS time:",
  sum(is.na(clinical_survival$OS_time_analysis)),
  "\n"
)

cat(
  "Нульовий OS time:",
  sum(clinical_survival$OS_zero_time_flag),
  "\n"
)

cat(
  "Події смерті:",
  sum(clinical_survival$OS_event_analysis == 1),
  "\n"
)

cat(
  "Цензуровані спостереження:",
  sum(clinical_survival$OS_event_analysis == 0),
  "\n"
)

# Перевірка дублікатів пацієнтів
duplicate_patients <- clinical_survival[
  duplicated(clinical_survival$patient) |
    duplicated(
      clinical_survival$patient,
      fromLast = TRUE
    ),
  ,
  drop = FALSE
]

cat(
  "Кількість записів пацієнтів із повтореннями:",
  nrow(duplicate_patients),
  "\n"
)


# ============================================================
# 9. ТОЧНЕ ЗІСТАВЛЕННЯ КЛІНІЧНИХ ДАНИХ
# ============================================================

# Повні barcode клінічних записів
clinical_ids <- as.character(
  clinical_survival$barcode
)

# Перевірка унікальності клінічних barcode
stopifnot(
  !anyNA(clinical_ids),
  !anyDuplicated(clinical_ids)
)

# Точне зіставлення клінічних записів
# з інтегрованими RNA-seq / methylation зразками
clinical_match_index <- match(
  matched_ids,
  clinical_ids
)

# Перевірка повноти відповідності
stopifnot(
  !anyNA(clinical_match_index)
)

# Вибір клінічних записів у відповідному порядку
clinical_matched <- clinical_survival[
  clinical_match_index,
  ,
  drop = FALSE
]

# Перевірка ідентичності barcode
stopifnot(
  identical(
    as.character(clinical_matched$barcode),
    matched_ids
  )
)

# Перевірка ідентичності пацієнтів
stopifnot(
  all(
    as.character(clinical_matched$patient) ==
      substr(matched_ids, 1, 12)
  )
)

cat("\n--- КЛІНІЧНЕ ЗІСТАВЛЕННЯ ---\n")

cat(
  "Узгоджених зразків:",
  nrow(clinical_matched), "\n"
)

cat(
  "Унікальних пацієнтів:",
  length(unique(clinical_matched$patient)),
  "\n"
)


# ============================================================
# 10. ФОРМУВАННЯ ІНДЕКСІВ ДЛЯ ІНТЕГРАЦІЇ
# ============================================================

# Індекс зразків метилювання за повним barcode
beta_index <- match(
  matched_ids,
  colnames(tcga_beta)
)

# Індекс RNA-seq за ідентифікатором пацієнта
# Матриця VST містить ідентифікатори рівня пацієнта
rna_index <- match(
  substr(matched_ids, 1, 12),
  colnames(tcga_rna_vst)
)

# Перевірка повноти відповідності
stopifnot(
  !anyNA(beta_index),
  !anyNA(rna_index)
)

# Перевірка відсутності повторних індексів
stopifnot(
  !anyDuplicated(beta_index),
  !anyDuplicated(rna_index)
)


# ============================================================
# 11. ФОРМУВАННЯ УЗГОДЖЕНИХ MULTI-OMICS МАТРИЦЬ
# ============================================================

# Матриця метилювання
tcga_beta_matched <- tcga_beta[
  ,
  beta_index,
  drop = FALSE
]

# Матриця експресії RNA-seq
tcga_rna_matched <- tcga_rna_vst[
  ,
  rna_index,
  drop = FALSE
]

# Узгоджена клінічна таблиця
tcga_clinical_matched <- clinical_matched

# Присвоєння однакових повних barcode
# усім трьом компонентам інтегрованої когорти
colnames(tcga_beta_matched) <- matched_ids

colnames(tcga_rna_matched) <- matched_ids

rownames(tcga_clinical_matched) <- matched_ids


# ============================================================
# 12. ФІНАЛЬНА ПЕРЕВІРКА УЗГОДЖЕНОСТІ
# ============================================================

# Перевірка кількості зразків
stopifnot(
  ncol(tcga_beta_matched) == 460,
  ncol(tcga_rna_matched) == 460,
  nrow(tcga_clinical_matched) == 460
)

# Перевірка порядку зразків
stopifnot(
  identical(
    colnames(tcga_beta_matched),
    colnames(tcga_rna_matched)
  )
)

stopifnot(
  identical(
    colnames(tcga_beta_matched),
    rownames(tcga_clinical_matched)
  )
)

# Перевірка кількості унікальних пацієнтів
stopifnot(
  length(
    unique(tcga_clinical_matched$patient)
  ) == 460
)

cat("\n--- ФІНАЛЬНА ПЕРЕВІРКА ---\n")

cat("Метилювання / RNA-seq: TRUE\n")
cat("Omics / клінічні дані: TRUE\n")
cat("460 унікальних пацієнтів: TRUE\n")


# ============================================================
# 13. СТВОРЕННЯ ІНТЕГРОВАНОГО ОБ'ЄКТА TCGA-SKCM
# ============================================================

TCGA_SKCM_integrated <- list(
  
  # Матриця beta-values
  beta_matrix = tcga_beta_matched,
  
  # Матриця нормалізованої експресії
  expression_matrix = tcga_rna_matched,
  
  # Клінічна інформація
  clinical_data = tcga_clinical_matched,
  
  # Метадані дослідження
  metadata = list(
    
    project = "TCGA-SKCM",
    
    n_samples = 460,
    
    n_unique_patients = 460,
    
    original_methylation_samples = 463,
    
    original_rnaseq_samples = 460,
    
    methylation_platform =
      "Illumina Human Methylation 450",
    
    methylation_measure =
      "Imputed beta-values",
    
    expression_data =
      "RNA-seq variance stabilizing transformation",
    
    matching_method =
      "Exact full TCGA barcode for methylation and clinical data; patient ID for RNA-seq matrix",
    
    clinical_alignment =
      "Exact full barcode",
    
    sample_types = list(
      primary = 102,
      metastatic = 358
    ),
    
    cohort_role = "Discovery cohort"
  )
)


# ============================================================
# 14. ЗБЕРЕЖЕННЯ ІНТЕГРОВАНОЇ КОГОРТИ
# ============================================================

integrated_file <- file.path(
  processed_dir,
  "TCGA_SKCM_integrated_multiomics.rds"
)

saveRDS(
  TCGA_SKCM_integrated,
  integrated_file
)

# Збереження таблиці відповідності зразків
write.csv(
  tcga_matching,
  file.path(
    results_tables_dir,
    "TCGA_SKCM_sample_matching.csv"
  ),
  row.names = FALSE
)

# Збереження розподілу типів зразків
write.csv(
  sample_type_summary,
  file.path(
    results_tables_dir,
    "TCGA_SKCM_sample_type_summary.csv"
  ),
  row.names = FALSE
)

# Збереження клінічного QC
write.csv(
  missing_summary,
  file.path(
    results_tables_dir,
    "TCGA_SKCM_clinical_missing_summary.csv"
  ),
  row.names = FALSE
)

# Збереження записів із невалідним часом виживання
write.csv(
  negative_survival,
  file.path(
    results_tables_dir,
    "TCGA_SKCM_negative_survival_records.csv"
  ),
  row.names = FALSE
)

# Збереження записів із нульовим часом виживання
write.csv(
  zero_survival,
  file.path(
    results_tables_dir,
    "TCGA_SKCM_zero_survival_records.csv"
  ),
  row.names = FALSE
)

# Збереження аудиту інтеграції
TCGA_SKCM_integration_audit <- list(
  
  sample_matching = tcga_matching,
  
  sample_type_summary = sample_type_summary,
  
  clinical_missing_summary = missing_summary,
  
  negative_survival_records = negative_survival,
  
  zero_survival_records = zero_survival,
  
  duplicate_patient_records = duplicate_patients,
  
  final_sample_ids = matched_ids
)

saveRDS(
  TCGA_SKCM_integration_audit,
  file.path(
    processed_dir,
    "TCGA_SKCM_integration_audit.rds"
  )
)


# ============================================================
# 15. ПІДСУМКОВИЙ ЗВІТ
# ============================================================

cat("\n")
cat("================================================\n")
cat("ІНТЕГРОВАНА КОГОРТА TCGA-SKCM\n")
cat("================================================\n")

cat(
  "Метилювання:",
  nrow(TCGA_SKCM_integrated$beta_matrix),
  "CpG x",
  ncol(TCGA_SKCM_integrated$beta_matrix),
  "зразків\n"
)

cat(
  "RNA-seq:",
  nrow(TCGA_SKCM_integrated$expression_matrix),
  "генів x",
  ncol(TCGA_SKCM_integrated$expression_matrix),
  "зразків\n"
)

cat(
  "Клінічні дані:",
  nrow(TCGA_SKCM_integrated$clinical_data),
  "зразків x",
  ncol(TCGA_SKCM_integrated$clinical_data),
  "змінних\n"
)

cat(
  "Унікальних пацієнтів:",
  length(
    unique(
      TCGA_SKCM_integrated$clinical_data$patient
    )
  ),
  "\n"
)

cat(
  "Метилювання / RNA-seq узгоджені:",
  identical(
    colnames(TCGA_SKCM_integrated$beta_matrix),
    colnames(TCGA_SKCM_integrated$expression_matrix)
  ),
  "\n"
)

cat(
  "Omics / клінічні дані узгоджені:",
  identical(
    colnames(TCGA_SKCM_integrated$beta_matrix),
    rownames(TCGA_SKCM_integrated$clinical_data)
  ),
  "\n"
)

cat(
  "Файл інтегрованої когорти збережено:",
  file.exists(integrated_file),
  "\n"
)

cat("\nАналіз успішно завершено.\n")
cat("================================================\n")
