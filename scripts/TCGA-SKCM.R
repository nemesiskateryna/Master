
###############################################################################
# TCGA_methylation_reprocessing.R
#
# 
#
# ДАНІ:
# TCGA-SKCM — DNA methylation, Illumina Human Methylation 450
#
# МЕТА:
# 1. Отримання DNA methylation data з GDC.
# 2. Перевірка завантажених файлів.
# 3. Зіставлення файлів із метаданими GDC.
# 4. Контроль якості CpG та Beta values.
# 5. Узгодження метилювання з RNA-seq когортою.
# 6. Аналіз пропущених значень.
# 7. Відбір Standard CpG.
# 8. Створення HDF5-матриці.
# 9. Контроль якості фінальної матриці.
#
# ВАЖЛИВО:
# Імпутація на цьому етапі НЕ проводиться.
# Вихідні GDC-файли не змінюються.

###############################################################################


###############################################################################
# 0. ПІДГОТОВКА СЕРЕДОВИЩА
###############################################################################

required_packages <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "dplyr",
  "rhdf5"
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
    "Не встановлені пакети: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(rhdf5)

# Коренева директорія проєкту
project_dir <- getwd()

# Перевірка структури проєкту
required_dirs <- c(
  "data/TCGA/raw_download",
  "data/TCGA/processed",
  "results/tables"
)

for (d in required_dirs) {
  dir.create(
    file.path(project_dir, d),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

raw_dir <- file.path(
  project_dir,
  "data/TCGA/raw_download/methylation_450K_GDC_2026"
)

processed_dir <- file.path(
  project_dir,
  "data/TCGA/processed"
)

tables_dir <- file.path(
  project_dir,
  "results/tables"
)

dir.create(
  raw_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


###############################################################################
# 1. ЗАПИТ ДО GDC
###############################################################################

query.met <- GDCquery(
  project = "TCGA-SKCM",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450",
  sample.type = c(
    "Primary Tumor",
    "Metastatic"
  )
)

met.results <- getResults(query.met)

cat("\n--- РЕЗУЛЬТАТ GDC QUERY ---\n")

cat("Кількість файлів:", nrow(met.results), "\n")

print(
  table(met.results$sample_type)
)

saveRDS(
  met.results,
  file.path(
    tables_dir,
    "TCGA_SKCM_GDC_methylation_metadata.rds"
  )
)


###############################################################################
# 2. ЗАВАНТАЖЕННЯ ФАЙЛІВ
###############################################################################

# Перевіряємо, чи всі очікувані файли вже є на диску.

expected_names <- as.character(
  met.results$file_name
)

existing_files <- list.files(
  raw_dir,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.txt$"
)

existing_names <- basename(existing_files)

missing_names <- setdiff(
  expected_names,
  existing_names
)

if (length(missing_names) > 0) {
  
  cat(
    "\nПотрібно завантажити файлів:",
    length(missing_names),
    "\n"
  )
  
  GDCdownload(
    query.met,
    method = "api",
    files.per.chunk = 100,
    directory = raw_dir
  )
  
} else {
  
  cat(
    "\nУсі файли вже присутні. Повторне завантаження не потрібне.\n"
  )
}

downloaded_files <- list.files(
  raw_dir,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.txt$"
)

cat(
  "\nФактично знайдено TXT-файлів:",
  length(downloaded_files),
  "\n"
)

stopifnot(
  length(downloaded_files) == nrow(met.results)
)


###############################################################################
# 3. ЗІСТАВЛЕННЯ ФАЙЛІВ З МЕТАДАНИМИ GDC
###############################################################################

local_names <- basename(downloaded_files)

metadata_match <- match(
  local_names,
  met.results$file_name
)

if (anyNA(metadata_match)) {
  stop("Не всі локальні файли зіставлені з метаданими GDC.")
}

sample_mapping <- data.frame(
  File_ID = met.results$id[metadata_match],
  File_Name = local_names,
  Patient_ID = met.results$cases.submitter_id[metadata_match],
  Sample_Barcode = met.results$sample.submitter_id[metadata_match],
  Sample_Type = met.results$sample_type[metadata_match],
  File_Path = downloaded_files,
  stringsAsFactors = FALSE
)

stopifnot(
  !anyDuplicated(sample_mapping$Sample_Barcode),
  !anyDuplicated(sample_mapping$File_ID)
)

cat("\n--- ЗІСТАВЛЕННЯ ФАЙЛІВ ---\n")

print(
  table(sample_mapping$Sample_Type)
)

saveRDS(
  sample_mapping,
  file.path(
    tables_dir,
    "TCGA_SKCM_methylation_sample_mapping.rds"
  )
)

write.csv(
  sample_mapping,
  file.path(
    tables_dir,
    "TCGA_SKCM_methylation_sample_mapping.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 4. ФУНКЦІЯ ЗЧИТУВАННЯ BETA-ФАЙЛІВ
###############################################################################

read_beta_file <- function(path) {
  
  x <- read.delim(
    path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    colClasses = c("character", "numeric"),
    na.strings = c("NA", "")
  )
  
  if (ncol(x) != 2) {
    stop("Неправильна структура файлу: ", path)
  }
  
  colnames(x) <- c(
    "CpG_ID",
    "Beta"
  )
  
  x
}


###############################################################################
# 5. ПОВНИЙ АУДИТ 472 ФАЙЛІВ
###############################################################################

reference_data <- read_beta_file(
  sample_mapping$File_Path[1]
)

reference_cpg <- reference_data$CpG_ID

full_file_audit <- vector(
  "list",
  length = nrow(sample_mapping)
)

for (i in seq_len(nrow(sample_mapping))) {
  
  x <- read_beta_file(
    sample_mapping$File_Path[i]
  )
  
  beta <- x$Beta
  
  full_file_audit[[i]] <- data.frame(
    File_Name = sample_mapping$File_Name[i],
    Sample_Barcode = sample_mapping$Sample_Barcode[i],
    Sample_Type = sample_mapping$Sample_Type[i],
    
    CpG_count = length(x$CpG_ID),
    
    Unique_CpG = length(unique(x$CpG_ID)),
    
    Duplicate_CpG = sum(
      duplicated(x$CpG_ID)
    ),
    
    Missing_CpG_ID = sum(
      is.na(x$CpG_ID) | x$CpG_ID == ""
    ),
    
    Missing_beta = sum(
      is.na(beta)
    ),
    
    Beta_below_0 = sum(
      beta < 0,
      na.rm = TRUE
    ),
    
    Beta_above_1 = sum(
      beta > 1,
      na.rm = TRUE
    ),
    
    Non_finite_beta = sum(
      !is.finite(beta) & !is.na(beta)
    ),
    
    Same_CpG_set_and_order = identical(
      x$CpG_ID,
      reference_cpg
    ),
    
    stringsAsFactors = FALSE
  )
  
  rm(x, beta)
  
  if (i %% 50 == 0) {
    cat("Перевірено файлів:", i, "\n")
    gc()
  }
}

full_file_audit <- bind_rows(
  full_file_audit
)

cat("\n--- ПОВНИЙ АУДИТ ФАЙЛІВ ---\n")

print(
  summary(full_file_audit$Missing_beta)
)

cat(
  "Файлів із невідповідністю CpG:",
  sum(!full_file_audit$Same_CpG_set_and_order),
  "\n"
)

saveRDS(
  full_file_audit,
  file.path(
    tables_dir,
    "TCGA_SKCM_methylation_full_file_audit.rds"
  )
)

write.csv(
  full_file_audit,
  file.path(
    tables_dir,
    "TCGA_SKCM_methylation_full_file_audit.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 6. ЗІСТАВЛЕННЯ З НАЯВНОЮ RNA-seq КОГОРОТОЮ
###############################################################################

integrated_existing <- readRDS(
  file.path(
    processed_dir,
    "TCGA_SKCM_integrated_multiomics.rds"
  )
)

rna_barcodes <- colnames(
  integrated_existing$expression_matrix
)

rna_sample_barcodes <- substr(
  rna_barcodes,
  1,
  16
)

methylation_barcodes <- sample_mapping$Sample_Barcode

common_barcodes <- rna_sample_barcodes[
  rna_sample_barcodes %in% methylation_barcodes
]

if (anyDuplicated(common_barcodes)) {
  stop("Виявлені дублікати barcode у спільній когорті.")
}

matched_mapping <- sample_mapping[
  match(
    common_barcodes,
    sample_mapping$Sample_Barcode
  ),
]

matched_mapping <- matched_mapping[
  match(
    rna_sample_barcodes,
    matched_mapping$Sample_Barcode
  ),
]

stopifnot(
  nrow(matched_mapping) == length(rna_sample_barcodes),
  identical(
    matched_mapping$Sample_Barcode,
    rna_sample_barcodes
  )
)

cat("\n--- RNA-seq / METHYLATION MATCHING ---\n")

cat("RNA-seq зразків:", length(rna_sample_barcodes), "\n")
cat("Метилювання зразків:", nrow(sample_mapping), "\n")
cat("Спільних зразків:", nrow(matched_mapping), "\n")

print(
  table(matched_mapping$Sample_Type)
)

saveRDS(
  matched_mapping,
  file.path(
    processed_dir,
    "TCGA_SKCM_matched_methylation_RNA_samples.rds"
  )
)

write.csv(
  matched_mapping,
  file.path(
    tables_dir,
    "TCGA_SKCM_matched_methylation_RNA_samples.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 7. АНАЛІЗ ПРОПУСКІВ У 460 УЗГОДЖЕНИХ ЗРАЗКАХ
###############################################################################

n_samples <- nrow(matched_mapping)

n_cpg_original <- length(reference_cpg)

cpg_missing_count <- numeric(
  n_cpg_original
)

sample_missing_count <- numeric(
  n_samples
)

for (j in seq_len(n_samples)) {
  
  x <- read_beta_file(
    matched_mapping$File_Path[j]
  )
  
  if (!identical(x$CpG_ID, reference_cpg)) {
    stop("Порядок CpG не збігається у файлі: ", j)
  }
  
  missing_vector <- is.na(x$Beta)
  
  sample_missing_count[j] <- sum(
    missing_vector
  )
  
  cpg_missing_count <- cpg_missing_count +
    as.numeric(missing_vector)
  
  rm(x, missing_vector)
  
  if (j %% 50 == 0) {
    cat("Оброблено зразків:", j, "\n")
    gc()
  }
}

cpg_missing_460 <- data.frame(
  CpG_ID = reference_cpg,
  
  Missing_count = cpg_missing_count,
  
  Missing_percent = 100 *
    cpg_missing_count / n_samples,
  
  stringsAsFactors = FALSE
)

cpg_missing_460$Probe_type <- ifelse(
  grepl("^cg[0-9]+$", cpg_missing_460$CpG_ID),
  "Standard CpG",
  ifelse(
    grepl("^ch[0-9]+$", cpg_missing_460$CpG_ID),
    "Non-CpG (ch)",
    "Other"
  )
)

cat("\n--- ПРОПУСКИ CpG У 460 ЗРАЗКАХ ---\n")

print(
  summary(cpg_missing_460$Missing_percent)
)

saveRDS(
  cpg_missing_460,
  file.path(
    tables_dir,
    "TCGA_SKCM_CpG_missingness_audit_460.rds"
  )
)

write.csv(
  cpg_missing_460,
  file.path(
    tables_dir,
    "TCGA_SKCM_CpG_missingness_audit_460.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 8. КЛАСИФІКАЦІЯ CpG ЗА ЧАСТКОЮ ПРОПУСКІВ
###############################################################################

missing_distribution <- data.frame(
  Missing_interval = c(
    "0%",
    ">0–5%",
    ">5–10%",
    ">10–20%",
    ">20–30%",
    ">30–50%",
    ">50–75%",
    ">75–<100%",
    "100%"
  ),
  
  CpG_count = c(
    sum(cpg_missing_460$Missing_percent == 0),
    sum(cpg_missing_460$Missing_percent > 0 &
          cpg_missing_460$Missing_percent <= 5),
    sum(cpg_missing_460$Missing_percent > 5 &
          cpg_missing_460$Missing_percent <= 10),
    sum(cpg_missing_460$Missing_percent > 10 &
          cpg_missing_460$Missing_percent <= 20),
    sum(cpg_missing_460$Missing_percent > 20 &
          cpg_missing_460$Missing_percent <= 30),
    sum(cpg_missing_460$Missing_percent > 30 &
          cpg_missing_460$Missing_percent <= 50),
    sum(cpg_missing_460$Missing_percent > 50 &
          cpg_missing_460$Missing_percent <= 75),
    sum(cpg_missing_460$Missing_percent > 75 &
          cpg_missing_460$Missing_percent < 100),
    sum(cpg_missing_460$Missing_percent == 100)
  )
)

missing_distribution$Percent_all <- 100 *
  missing_distribution$CpG_count /
  nrow(cpg_missing_460)

print(
  missing_distribution,
  row.names = FALSE
)

write.csv(
  missing_distribution,
  file.path(
    tables_dir,
    "TCGA_SKCM_CpG_missingness_distribution.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 9. ВІДБІР STANDARD CpG З ПРОПУСКАМИ НЕ БІЛЬШЕ 20%
###############################################################################

cpg_keep <- cpg_missing_460$CpG_ID[
  cpg_missing_460$Probe_type == "Standard CpG" &
    cpg_missing_460$Missing_percent <= 20
]

stopifnot(
  !anyDuplicated(cpg_keep)
)

cat("\n--- ВІДБІР CpG ---\n")

cat(
  "Початкова кількість Standard CpG:",
  sum(cpg_missing_460$Probe_type == "Standard CpG"),
  "\n"
)

cat(
  "Відібрано CpG:",
  length(cpg_keep),
  "\n"
)

cat(
  "Без пропусків:",
  sum(
    cpg_missing_460$CpG_ID %in% cpg_keep &
      cpg_missing_460$Missing_percent == 0
  ),
  "\n"
)

cat(
  "З пропусками до 20%:",
  sum(
    cpg_missing_460$CpG_ID %in% cpg_keep &
      cpg_missing_460$Missing_percent > 0
  ),
  "\n"
)

saveRDS(
  cpg_keep,
  file.path(
    processed_dir,
    "TCGA_SKCM_selected_CpG_list.rds"
  )
)

write.csv(
  data.frame(CpG_ID = cpg_keep),
  file.path(
    tables_dir,
    "TCGA_SKCM_selected_CpG_list.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 10. СТВОРЕННЯ HDF5-МАТРИЦІ
###############################################################################

h5_file <- file.path(
  processed_dir,
  "TCGA_SKCM_beta_matrix_460_StandardCpG_20pct.h5"
)

if (!file.exists(h5_file)) {
  
  n_selected_cpg <- length(cpg_keep)
  
  n_selected_samples <- nrow(matched_mapping)
  
  h5createFile(h5_file)
  
  h5createDataset(
    file = h5_file,
    dataset = "beta",
    dims = c(
      n_selected_cpg,
      n_selected_samples
    ),
    storage.mode = "double",
    chunk = c(
      min(10000, n_selected_cpg),
      1
    ),
    level = 6
  )
  
  # Ідентифікатори CpG та зразків
  h5write(
    cpg_keep,
    h5_file,
    "CpG_ID"
  )
  
  h5write(
    as.character(matched_mapping$Sample_Barcode),
    h5_file,
    "Sample_Barcode"
  )
  
  # Послідовне записування кожного зразка
  for (j in seq_len(n_selected_samples)) {
    
    x <- read_beta_file(
      matched_mapping$File_Path[j]
    )
    
    selected_index <- match(
      cpg_keep,
      x$CpG_ID
    )
    
    if (anyNA(selected_index)) {
      stop("Не всі CpG знайдені у файлі: ", j)
    }
    
    beta_values <- matrix(
      x$Beta[selected_index],
      ncol = 1
    )
    
    h5write(
      beta_values,
      h5_file,
      "beta",
      index = list(
        seq_len(n_selected_cpg),
        j
      )
    )
    
    rm(x, selected_index, beta_values)
    
    if (j %% 25 == 0) {
      cat(
        "Записано зразків:",
        j,
        "/",
        n_selected_samples,
        "\n"
      )
      
      gc()
    }
  }
  
  cat("\nHDF5-матрицю створено.\n")
  
} else {
  
  cat(
    "\nHDF5-матриця вже існує. Перезаписування пропущено.\n"
  )
}


###############################################################################
# 11. ПЕРЕВІРКА СТРУКТУРИ HDF5-МАТРИЦІ
###############################################################################

h5_structure <- h5ls(h5_file)

print(h5_structure)

h5_cpg <- h5read(
  h5_file,
  "CpG_ID"
)

h5_samples <- h5read(
  h5_file,
  "Sample_Barcode"
)

stopifnot(
  identical(
    as.character(h5_cpg),
    as.character(cpg_keep)
  ),
  
  identical(
    as.character(h5_samples),
    as.character(rna_sample_barcodes)
  )
)

cat("\n--- ПЕРЕВІРКА ІДЕНТИФІКАТОРІВ ---\n")

cat(
  "CpG відповідають списку:",
  identical(as.character(h5_cpg), as.character(cpg_keep)),
  "\n"
)

cat(
  "Зразки відповідають RNA-seq:",
  identical(as.character(h5_samples), as.character(rna_sample_barcodes)),
  "\n"
)


###############################################################################
# 12. QC ФІНАЛЬНОЇ HDF5-МАТРИЦІ
###############################################################################

n_selected_cpg <- length(cpg_keep)

n_selected_samples <- length(h5_samples)

block_size <- 10000

total_na <- 0

beta_min <- Inf
beta_max <- -Inf

sample_missing_count_h5 <- numeric(
  n_selected_samples
)

for (start_row in seq(
  1,
  n_selected_cpg,
  by = block_size
)) {
  
  end_row <- min(
    start_row + block_size - 1,
    n_selected_cpg
  )
  
  beta_block <- h5read(
    h5_file,
    "beta",
    index = list(
      start_row:end_row,
      seq_len(n_selected_samples)
    )
  )
  
  total_na <- total_na +
    sum(is.na(beta_block))
  
  finite_values <- beta_block[
    is.finite(beta_block)
  ]
  
  if (length(finite_values) > 0) {
    beta_min <- min(beta_min, finite_values)
    beta_max <- max(beta_max, finite_values)
  }
  
  sample_missing_count_h5 <-
    sample_missing_count_h5 +
    colSums(is.na(beta_block))
  
  rm(beta_block, finite_values)
  
  gc()
}

total_values <- n_selected_cpg *
  n_selected_samples

overall_missing_percent <- 100 *
  total_na / total_values

cat("\n--- QC HDF5-МАТРИЦІ ---\n")

cat("CpG:", n_selected_cpg, "\n")
cat("Зразків:", n_selected_samples, "\n")
cat("Усього значень:", total_values, "\n")
cat("NA:", total_na, "\n")
cat("Загальні пропуски (%):", overall_missing_percent, "\n")
cat("Мінімальний Beta:", beta_min, "\n")
cat("Максимальний Beta:", beta_max, "\n")


###############################################################################
# 13. QC ПРОПУСКІВ НА РІВНІ ЗРАЗКІВ
###############################################################################

sample_missingness_final <- data.frame(
  Sample_Barcode = as.character(h5_samples),
  
  Sample_Type = matched_mapping$Sample_Type[
    match(
      h5_samples,
      matched_mapping$Sample_Barcode
    )
  ],
  
  Missing_count = sample_missing_count_h5,
  
  Missing_percent = 100 *
    sample_missing_count_h5 /
    n_selected_cpg,
  
  stringsAsFactors = FALSE
)

cat("\n--- ПРОПУСКИ ЗА ЗРАЗКАМИ ---\n")

print(
  summary(sample_missingness_final$Missing_percent)
)

cat(
  "Зразків із пропусками понад 5%:",
  sum(sample_missingness_final$Missing_percent > 5),
  "\n"
)

cat(
  "Зразків із пропусками понад 10%:",
  sum(sample_missingness_final$Missing_percent > 10),
  "\n"
)

saveRDS(
  sample_missingness_final,
  file.path(
    tables_dir,
    "TCGA_SKCM_sample_missingness_final_460.rds"
  )
)

write.csv(
  sample_missingness_final,
  file.path(
    tables_dir,
    "TCGA_SKCM_sample_missingness_final_460.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 14. ПОРІВНЯННЯ ПРОПУСКІВ МІЖ ТИПАМИ ЗРАЗКІВ
###############################################################################

missingness_by_type <- do.call(
  rbind,
  lapply(
    split(
      sample_missingness_final,
      sample_missingness_final$Sample_Type
    ),
    function(x) {
      
      data.frame(
        Sample_Type = x$Sample_Type[1],
        
        N = nrow(x),
        
        Mean_missing_percent = mean(
          x$Missing_percent
        ),
        
        Median_missing_percent = median(
          x$Missing_percent
        ),
        
        Maximum_missing_percent = max(
          x$Missing_percent
        ),
        
        Samples_over_5_percent = sum(
          x$Missing_percent > 5
        ),
        
        Samples_over_10_percent = sum(
          x$Missing_percent > 10
        )
      )
    }
  )
)

rownames(missingness_by_type) <- NULL

cat("\n--- ПРОПУСКИ ЗА ТИПОМ ЗРАЗКІВ ---\n")

print(
  missingness_by_type,
  row.names = FALSE
)

wilcox_result <- wilcox.test(
  Missing_percent ~ Sample_Type,
  data = sample_missingness_final,
  exact = FALSE
)

cat("\n--- ТЕСТ ВІЛКОКСОНА ---\n")

print(wilcox_result)

write.csv(
  missingness_by_type,
  file.path(
    tables_dir,
    "TCGA_SKCM_missingness_by_sample_type.csv"
  ),
  row.names = FALSE
)

saveRDS(
  wilcox_result,
  file.path(
    tables_dir,
    "TCGA_SKCM_missingness_wilcoxon.rds"
  )
)


###############################################################################
# 15. ФІНАЛЬНЕ РІШЕННЯ ЩОДО ЯКОСТІ ЗРАЗКІВ
###############################################################################

sample_qc_decision <- sample_missingness_final

sample_qc_decision$QC_status <- ifelse(
  sample_qc_decision$Missing_percent > 10,
  
  "High_missingness_check_in_sensitivity_analysis",
  
  ifelse(
    sample_qc_decision$Missing_percent > 5,
    
    "Moderate_missingness_retain",
    
    "Acceptable_missingness_retain"
  )
)

cat("\n--- РІШЕННЯ ЩОДО ЯКОСТІ ЗРАЗКІВ ---\n")

print(
  table(sample_qc_decision$QC_status)
)

cat(
  "\nКількість зразків у когорті:",
  nrow(sample_qc_decision),
  "\n"
)

cat(
  "Зразків, які залишаються в основній когорті:",
  sum(sample_qc_decision$QC_status != "Exclude"),
  "\n"
)

saveRDS(
  sample_qc_decision,
  file.path(
    tables_dir,
    "TCGA_SKCM_sample_QC_decision.rds"
  )
)

write.csv(
  sample_qc_decision,
  file.path(
    tables_dir,
    "TCGA_SKCM_sample_QC_decision.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 16. ФІНАЛЬНИЙ ПІДСУМОК
###############################################################################

cat("\n")
cat("====================================================\n")
cat("TCGA-SKCM DNA METHYLATION: ЕТАП QC ЗАВЕРШЕНО\n")
cat("====================================================\n")

cat("Початкових GDC файлів:", nrow(sample_mapping), "\n")

cat("Узгоджених RNA-seq / methylation зразків:",
    nrow(matched_mapping), "\n")

cat("Відібраних Standard CpG:", length(cpg_keep), "\n")

cat("Розмір HDF5: ",
    length(cpg_keep), " x ", nrow(matched_mapping), "\n",
    sep = "")

cat("Загальна частка пропусків:",
    round(overall_missing_percent, 4), "%\n")

cat("Імпутація: НЕ проводилася\n")

cat("Виключення зразків: НЕ проводилося\n")

cat("HDF5-файл:", h5_file, "\n")

cat("====================================================\n")

gc()
