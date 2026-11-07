
###############################################################################
# TCGA-SKCM
# КОНТРОЛЬ ЯКОСТІ DNA METHYLATION ТА АНАЛІЗ MISSINGNESS
#
# Основна матриця:
# 409 556 Standard CpG × 460 зразків
#
# Контрольна матриця без пропусків:
# 320 046 Standard CpG × 460 зразків
#
# Суворий контрольний набір:
# 404 459 CpG
#
# Імпутація на цьому етапі НЕ виконується.
# Основна HDF5-матриця не змінюється.
###############################################################################


###############################################################################
# STEP 12. ПІДГОТОВКА КОНТРОЛЬНОГО НАБОРУ CpG БЕЗ ПРОПУСКІВ
###############################################################################

library(HDF5Array)
library(rhdf5)
library(ggplot2)

# Шляхи до файлів

project_dir <- getwd()

main_h5_file <- file.path(
  project_dir,
  "data", "TCGA", "processed",
  "TCGA_SKCM_beta_matrix_460_StandardCpG_20pct.h5"
)

cpg_audit_file <- file.path(
  project_dir,
  "results", "tables",
  "TCGA_SKCM_CpG_missingness_audit_460.rds"
)

mapping_file <- file.path(
  project_dir,
  "data", "TCGA", "processed",
  "TCGA_SKCM_matched_methylation_RNA_samples.rds"
)

# Завантаження попередньо збережених результатів

cpg_missing_460 <- readRDS(cpg_audit_file)

matched_mapping <- readRDS(mapping_file)

# CpG без жодного пропуску

cpg_complete <- cpg_missing_460$CpG_ID[
  cpg_missing_460$Probe_type == "Standard CpG" &
    cpg_missing_460$Missing_percent == 0
]

cpg_complete <- as.character(cpg_complete)

# Перевірки

stopifnot(
  length(cpg_complete) == 320046,
  !anyNA(cpg_complete),
  !anyDuplicated(cpg_complete)
)

cat("\n--- STEP 12. КОНТРОЛЬНИЙ НАБІР ---\n")

cat("Кількість повністю заповнених CpG:",
    length(cpg_complete), "\n")

# Збереження списку

saveRDS(
  cpg_complete,
  file.path(
    project_dir,
    "data", "TCGA", "processed",
    "TCGA_SKCM_control_CpG_complete.rds"
  )
)

write.csv(
  data.frame(CpG_ID = cpg_complete),
  file.path(
    project_dir,
    "results", "tables",
    "TCGA_SKCM_control_CpG_complete.csv"
  ),
  row.names = FALSE
)


###############################################################################
# STEP 13. ПЕРЕВІРКА ВІДПОВІДНОСТІ CpG ОСНОВНІЙ HDF5-МАТРИЦІ
###############################################################################

beta_h5 <- HDF5Array(
  main_h5_file,
  "beta"
)

main_cpg <- as.character(
  rhdf5::h5read(main_h5_file, "CpG_ID")
)

main_samples <- as.character(
  rhdf5::h5read(main_h5_file, "Sample_Barcode")
)

complete_indices <- match(
  cpg_complete,
  main_cpg
)

stopifnot(!anyNA(complete_indices))

cat("\n--- STEP 13. ПЕРЕВІРКА CpG ---\n")

cat("CpG знайдено в основній матриці:",
    sum(!is.na(complete_indices)), "\n")

cat("Відповідність порядку CpG:",
    identical(main_cpg[complete_indices], cpg_complete), "\n")

cat("Розмір основної матриці:",
    paste(dim(beta_h5), collapse = " × "), "\n")


###############################################################################
# STEP 14. QC КОНТРОЛЬНОГО НАБОРУ БЕЗ ПРОПУСКІВ
###############################################################################

control_na <- 0
control_min <- Inf
control_max <- -Inf

block_size <- 10000

for (start_row in seq(
  1, length(complete_indices), by = block_size
)) {
  
  end_row <- min(
    start_row + block_size - 1,
    length(complete_indices)
  )
  
  current_indices <- complete_indices[
    start_row:end_row
  ]
  
  beta_block <- as.matrix(
    beta_h5[current_indices, , drop = FALSE]
  )
  
  control_na <- control_na + sum(is.na(beta_block))
  
  control_min <- min(
    control_min,
    min(beta_block, na.rm = TRUE)
  )
  
  control_max <- max(
    control_max,
    max(beta_block, na.rm = TRUE)
  )
  
  rm(beta_block)
  gc()
}

control_qc_report <- data.frame(
  Dataset = "TCGA-SKCM Standard CpG complete control set",
  CpG_count = length(complete_indices),
  Sample_count = length(main_samples),
  Total_beta_values = length(complete_indices) * length(main_samples),
  Missing_values = control_na,
  Missingness_percent = control_na /
    (length(complete_indices) * length(main_samples)) * 100,
  Minimum_beta = control_min,
  Maximum_beta = control_max,
  Imputation_performed = FALSE,
  QC_status = ifelse(control_na == 0, "PASS", "CHECK")
)

cat("\n--- STEP 14. QC КОНТРОЛЬНОГО НАБОРУ ---\n")

print(control_qc_report, row.names = FALSE)

saveRDS(
  control_qc_report,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_control_CpG_QC_report.rds"
  )
)

write.csv(
  control_qc_report,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_control_CpG_QC_report.csv"
  ),
  row.names = FALSE
)


###############################################################################
# STEP 15. РОЗПОДІЛ MISSINGNESS У ОСНОВНОМУ НАБОРІ CpG
###############################################################################

main_cpg_audit <- cpg_missing_460[
  cpg_missing_460$CpG_ID %in%
    cpg_missing_460$CpG_ID[
      cpg_missing_460$Probe_type == "Standard CpG" &
        cpg_missing_460$Missing_percent <= 20
    ],
  ,
  drop = FALSE
]

# Сортування відповідно до основної HDF5-матриці

main_cpg_audit <- main_cpg_audit[
  match(main_cpg, main_cpg_audit$CpG_ID),
  ,
  drop = FALSE
]

stopifnot(
  nrow(main_cpg_audit) == length(main_cpg),
  !anyNA(main_cpg_audit$Missing_percent)
)

missing_distribution <- data.frame(
  Missing_interval = c(
    "0%",
    ">0–5%",
    ">5–10%",
    ">10–15%",
    ">15–20%"
  ),
  CpG_count = c(
    sum(main_cpg_audit$Missing_percent == 0),
    sum(main_cpg_audit$Missing_percent > 0 &
          main_cpg_audit$Missing_percent <= 5),
    sum(main_cpg_audit$Missing_percent > 5 &
          main_cpg_audit$Missing_percent <= 10),
    sum(main_cpg_audit$Missing_percent > 10 &
          main_cpg_audit$Missing_percent <= 15),
    sum(main_cpg_audit$Missing_percent > 15 &
          main_cpg_audit$Missing_percent <= 20)
  )
)

missing_distribution$Percentage <- round(
  missing_distribution$CpG_count /
    nrow(main_cpg_audit) * 100,
  2
)

cat("\n--- STEP 15. РОЗПОДІЛ ПРОПУСКІВ CpG ---\n")

print(missing_distribution, row.names = FALSE)

cat("\nСередня missingness:",
    mean(main_cpg_audit$Missing_percent), "%\n")

cat("Медіанна missingness:",
    median(main_cpg_audit$Missing_percent), "%\n")

cat("Максимальна missingness:",
    max(main_cpg_audit$Missing_percent), "%\n")

saveRDS(
  missing_distribution,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_main_CpG_missingness_distribution.rds"
  )
)

write.csv(
  missing_distribution,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_main_CpG_missingness_distribution.csv"
  ),
  row.names = FALSE
)


###############################################################################
# STEP 16. MISSINGNESS НА РІВНІ ОКРЕМИХ ЗРАЗКІВ
###############################################################################

# Перевірка відповідності зразків

stopifnot(
  identical(
    as.character(matched_mapping$Sample_Barcode),
    main_samples
  )
)

sample_types <- as.character(
  matched_mapping$Sample_Type
)

sample_missing_count <- numeric(length(main_samples))

sample_valid_count <- numeric(length(main_samples))

for (start_row in seq(
  1, nrow(beta_h5), by = block_size
)) {
  
  end_row <- min(
    start_row + block_size - 1,
    nrow(beta_h5)
  )
  
  beta_block <- as.matrix(
    beta_h5[start_row:end_row, , drop = FALSE]
  )
  
  sample_missing_count <- sample_missing_count +
    colSums(is.na(beta_block))
  
  sample_valid_count <- sample_valid_count +
    colSums(!is.na(beta_block))
  
  rm(beta_block)
  gc()
}

sample_missingness_final <- data.frame(
  Sample_Barcode = main_samples,
  Sample_Type = sample_types,
  Missing_values = sample_missing_count,
  Missing_percent = sample_missing_count /
    nrow(beta_h5) * 100
)

cat("\n--- STEP 16. MISSINGNESS ЗА ЗРАЗКАМИ ---\n")

print(summary(sample_missingness_final$Missing_percent))

cat("\nЗразків із missingness понад 5%:",
    sum(sample_missingness_final$Missing_percent > 5), "\n")

cat("Зразків із missingness понад 10%:",
    sum(sample_missingness_final$Missing_percent > 10), "\n")

saveRDS(
  sample_missingness_final,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_sample_missingness_final_460.rds"
  )
)

write.csv(
  sample_missingness_final,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_sample_missingness_final_460.csv"
  ),
  row.names = FALSE
)


###############################################################################
# STEP 17. ПОРІВНЯННЯ MISSINGNESS МІЖ ТИПАМИ ЗРАЗКІВ
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
        Sample_count = nrow(x),
        Mean_missingness = mean(x$Missing_percent),
        Median_missingness = median(x$Missing_percent),
        Maximum_missingness = max(x$Missing_percent),
        Samples_over_5_percent = sum(x$Missing_percent > 5),
        Samples_over_10_percent = sum(x$Missing_percent > 10)
      )
    }
  )
)

wilcox_missingness <- wilcox.test(
  Missing_percent ~ Sample_Type,
  data = sample_missingness_final,
  exact = FALSE
)

cat("\n--- STEP 17. MISSINGNESS ЗА ТИПОМ ЗРАЗКА ---\n")

print(missingness_by_type, row.names = FALSE)

cat("\nТест Вілкоксона:\n")

print(wilcox_missingness)

missingness_comparison <- list(
  Summary = missingness_by_type,
  Wilcoxon_test = wilcox_missingness
)

saveRDS(
  missingness_comparison,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_missingness_by_sample_type.rds"
  )
)

write.csv(
  missingness_by_type,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_missingness_by_sample_type.csv"
  ),
  row.names = FALSE
)


###############################################################################
# STEP 18. MISSINGNESS ОКРЕМО ДЛЯ КОЖНОГО CpG У ДВОХ ГРУПАХ
###############################################################################

met_idx <- which(sample_types == "Metastatic")

primary_idx <- which(sample_types == "Primary Tumor")

stopifnot(
  length(met_idx) == 358,
  length(primary_idx) == 102
)

cpg_missing_by_type <- data.frame(
  CpG_ID = main_cpg,
  Missing_Metastatic = integer(length(main_cpg)),
  Missing_Primary = integer(length(main_cpg)),
  Missing_percent_Metastatic = numeric(length(main_cpg)),
  Missing_percent_Primary = numeric(length(main_cpg))
)

for (start_row in seq(
  1, nrow(beta_h5), by = block_size
)) {
  
  end_row <- min(
    start_row + block_size - 1,
    nrow(beta_h5)
  )
  
  rows <- start_row:end_row
  
  beta_block <- as.matrix(
    beta_h5[rows, , drop = FALSE]
  )
  
  cpg_missing_by_type$Missing_Metastatic[rows] <-
    rowSums(is.na(beta_block[, met_idx, drop = FALSE]))
  
  cpg_missing_by_type$Missing_Primary[rows] <-
    rowSums(is.na(beta_block[, primary_idx, drop = FALSE]))
  
  rm(beta_block)
  gc()
}

cpg_missing_by_type$Missing_percent_Metastatic <-
  cpg_missing_by_type$Missing_Metastatic /
  length(met_idx) * 100

cpg_missing_by_type$Missing_percent_Primary <-
  cpg_missing_by_type$Missing_Primary /
  length(primary_idx) * 100

cat("\n--- STEP 18. MISSINGNESS ЗА ГРУПАМИ ---\n")

cat(
  "Пропуски лише в Primary:",
  sum(
    cpg_missing_by_type$Missing_Metastatic == 0 &
      cpg_missing_by_type$Missing_Primary > 0
  ), "\n"
)

cat(
  "Пропуски лише в Metastatic:",
  sum(
    cpg_missing_by_type$Missing_Primary == 0 &
      cpg_missing_by_type$Missing_Metastatic > 0
  ), "\n"
)

cat(
  "Пропуски в обох групах:",
  sum(
    cpg_missing_by_type$Missing_Primary > 0 &
      cpg_missing_by_type$Missing_Metastatic > 0
  ), "\n"
)


###############################################################################
# STEP 19. РОЗПОДІЛ CpG ЗА ІНТЕРВАЛАМИ MISSINGNESS У КОЖНІЙ ГРУПІ
###############################################################################

classify_missingness <- function(x) {
  cut(
    x,
    breaks = c(-Inf, 0, 5, 10, 20, Inf),
    labels = c(
      "0%",
      ">0–5%",
      ">5–10%",
      ">10–20%",
      ">20%"
    ),
    right = TRUE
  )
}

metastatic_bins <- table(
  classify_missingness(
    cpg_missing_by_type$Missing_percent_Metastatic
  )
)

primary_bins <- table(
  classify_missingness(
    cpg_missing_by_type$Missing_percent_Primary
  )
)

group_missingness_bins <- data.frame(
  Missing_interval = names(metastatic_bins),
  Metastatic = as.integer(metastatic_bins),
  Primary_Tumor = as.integer(primary_bins)
)

cat("\n--- STEP 19. ІНТЕРВАЛИ MISSINGNESS ---\n")

print(group_missingness_bins, row.names = FALSE)

write.csv(
  group_missingness_bins,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_CpG_missingness_bins_by_group.csv"
  ),
  row.names = FALSE
)


###############################################################################
# STEP 20. РІЗНИЦЯ MISSINGNESS МІЖ PRIMARY ТА METASTATIC
###############################################################################

cpg_missing_by_type$Missingness_difference <-
  cpg_missing_by_type$Missing_percent_Primary -
  cpg_missing_by_type$Missing_percent_Metastatic

cpg_missing_by_type$Absolute_difference <-
  abs(cpg_missing_by_type$Missingness_difference)

cat("\n--- STEP 20. РІЗНИЦЯ MISSINGNESS ---\n")

print(
  summary(cpg_missing_by_type$Missingness_difference)
)

cat(
  "\nБільше пропусків у Primary:",
  sum(cpg_missing_by_type$Missingness_difference > 0),
  "\n"
)

cat(
  "Більше пропусків у Metastatic:",
  sum(cpg_missing_by_type$Missingness_difference < 0),
  "\n"
)

cat(
  "Однакова частка пропусків:",
  sum(cpg_missing_by_type$Missingness_difference == 0),
  "\n"
)

top_missingness_difference <- cpg_missing_by_type[
  order(-cpg_missing_by_type$Absolute_difference),
  ,
  drop = FALSE
]

print(
  head(
    top_missingness_difference[
      ,
      c(
        "CpG_ID",
        "Missing_percent_Metastatic",
        "Missing_percent_Primary",
        "Missingness_difference",
        "Absolute_difference"
      )
    ],
    20
  )
)


###############################################################################
# STEP 21. ЗБЕРЕЖЕННЯ АУДИТУ MISSINGNESS ЗА ТИПОМ ЗРАЗКА
###############################################################################

saveRDS(
  cpg_missing_by_type,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_CpG_missingness_by_sample_type.rds"
  )
)

write.csv(
  cpg_missing_by_type,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_CpG_missingness_by_sample_type.csv"
  ),
  row.names = FALSE
)


###############################################################################
# STEP 22. КІЛЬКІСНА ОЦІНКА РІЗНИЦІ ПРОПУСКІВ
###############################################################################

difference_summary <- data.frame(
  Difference_threshold = c(
    ">5 percentage points",
    ">10 percentage points",
    ">15 percentage points",
    ">20 percentage points"
  ),
  CpG_count = c(
    sum(cpg_missing_by_type$Absolute_difference > 5),
    sum(cpg_missing_by_type$Absolute_difference > 10),
    sum(cpg_missing_by_type$Absolute_difference > 15),
    sum(cpg_missing_by_type$Absolute_difference > 20)
  )
)

difference_summary$Percentage <- round(
  difference_summary$CpG_count /
    nrow(cpg_missing_by_type) * 100,
  2
)

cat("\n--- STEP 22. РІЗНИЦЯ ПРОПУСКІВ ---\n")

print(difference_summary, row.names = FALSE)

write.csv(
  difference_summary,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_CpG_missingness_difference_summary.csv"
  ),
  row.names = FALSE
)


###############################################################################
# STEP 23. СТВОРЕННЯ СУВОРОГО КОНТРОЛЬНОГО НАБОРУ
###############################################################################

cpg_strict <- cpg_missing_by_type$CpG_ID[
  cpg_missing_by_type$Missing_percent_Primary <= 20 &
    cpg_missing_by_type$Missing_percent_Metastatic <= 20
]

cpg_flagged <- cpg_missing_by_type$CpG_ID[
  cpg_missing_by_type$Missing_percent_Primary > 20 |
    cpg_missing_by_type$Missing_percent_Metastatic > 20
]

cat("\n--- STEP 23. СУВОРА КОНТРОЛЬНА ВИБІРКА ---\n")

cat("Початковий набір:", nrow(cpg_missing_by_type), "\n")
cat("CpG після суворого фільтра:", length(cpg_strict), "\n")
cat("CpG поза суворим набором:", length(cpg_flagged), "\n")

stopifnot(
  length(cpg_strict) == 404459,
  length(cpg_flagged) == 5097,
  length(intersect(cpg_strict, cpg_flagged)) == 0
)

saveRDS(
  cpg_strict,
  file.path(
    project_dir, "data", "TCGA", "processed",
    "TCGA_SKCM_strict_CpG_list_20pct_by_group.rds"
  )
)

saveRDS(
  cpg_flagged,
  file.path(
    project_dir, "data", "TCGA", "processed",
    "TCGA_SKCM_flagged_CpG_list_20pct_by_group.rds"
  )
)

write.csv(
  data.frame(CpG_ID = cpg_strict),
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_strict_CpG_list_20pct_by_group.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(CpG_ID = cpg_flagged),
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_flagged_CpG_list_20pct_by_group.csv"
  ),
  row.names = FALSE
)


###############################################################################
# STEP 24. ПЕРЕВІРКА ЗБЕРЕЖЕННЯ СУВОРОГО НАБОРУ
###############################################################################

files_to_check <- c(
  "data/TCGA/processed/TCGA_SKCM_strict_CpG_list_20pct_by_group.rds",
  "data/TCGA/processed/TCGA_SKCM_flagged_CpG_list_20pct_by_group.rds",
  "results/tables/TCGA_SKCM_strict_CpG_list_20pct_by_group.csv",
  "results/tables/TCGA_SKCM_flagged_CpG_list_20pct_by_group.csv"
)

file_check <- file.exists(
  file.path(project_dir, files_to_check)
)

cat("\n--- STEP 24. ПЕРЕВІРКА ФАЙЛІВ ---\n")

print(
  data.frame(
    File = files_to_check,
    Saved = file_check
  ),
  row.names = FALSE
)


###############################################################################
# STEP 25. ВИБІРКА 1000 CpG ДЛЯ ПОРІВНЯННЯ BETA-ЗНАЧЕНЬ
###############################################################################

set.seed(2026)

selected_indices <- sample(
  complete_indices,
  1000
)

beta_subset <- as.matrix(
  beta_h5[selected_indices, , drop = FALSE]
)

primary_values <- as.numeric(
  beta_subset[, primary_idx, drop = FALSE]
)

metastatic_values <- as.numeric(
  beta_subset[, met_idx, drop = FALSE]
)

beta_distribution <- data.frame(
  Sample_Type = c("Primary Tumor", "Metastatic"),
  N = c(
    length(primary_values),
    length(metastatic_values)
  ),
  Mean = c(
    mean(primary_values),
    mean(metastatic_values)
  ),
  Median = c(
    median(primary_values),
    median(metastatic_values)
  ),
  SD = c(
    sd(primary_values),
    sd(metastatic_values)
  )
)

cat("\n--- STEP 25. РОЗПОДІЛ BETA-ЗНАЧЕНЬ ---\n")

print(beta_distribution, row.names = FALSE)


###############################################################################
# STEP 26. BOXPLOT РОЗПОДІЛУ BETA-ЗНАЧЕНЬ
###############################################################################

beta_plot_data <- data.frame(
  Beta = c(primary_values, metastatic_values),
  Sample_Type = factor(
    c(
      rep("Primary Tumor", length(primary_values)),
      rep("Metastatic", length(metastatic_values))
    ),
    levels = c("Primary Tumor", "Metastatic")
  )
)

p_beta_distribution <- ggplot(
  beta_plot_data,
  aes(x = Sample_Type, y = Beta)
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA
  ) +
  labs(
    title = "Розподіл β-значень метилювання",
    subtitle = "Випадкова вибірка 1000 повністю заповнених CpG",
    x = "Тип зразка",
    y = "β-значення"
  ) +
  theme_classic() +
  coord_cartesian(ylim = c(0, 1))

print(p_beta_distribution)

ggsave(
  file.path(
    project_dir, "results", "figures",
    "TCGA_SKCM_beta_distribution_by_sample_type.png"
  ),
  plot = p_beta_distribution,
  width = 8,
  height = 6,
  dpi = 300
)


###############################################################################
# STEP 27. СЕРЕДНЄ BETA-ЗНАЧЕННЯ ДЛЯ КОЖНОГО ЗРАЗКА
###############################################################################

sample_beta_sum <- numeric(length(main_samples))
sample_beta_count <- numeric(length(main_samples))

for (start_row in seq(
  1, nrow(beta_h5), by = block_size
)) {
  
  end_row <- min(
    start_row + block_size - 1,
    nrow(beta_h5)
  )
  
  beta_block <- as.matrix(
    beta_h5[start_row:end_row, , drop = FALSE]
  )
  
  sample_beta_sum <- sample_beta_sum +
    colSums(beta_block, na.rm = TRUE)
  
  sample_beta_count <- sample_beta_count +
    colSums(!is.na(beta_block))
  
  rm(beta_block)
  gc()
}

sample_mean_beta <- sample_beta_sum / sample_beta_count

sample_missingness_analysis <- data.frame(
  Sample_Barcode = main_samples,
  Sample_Type = sample_types,
  Missing_percent = sample_missingness_final$Missing_percent,
  Mean_beta = sample_mean_beta
)

stopifnot(
  identical(
    sample_missingness_analysis$Sample_Barcode,
    sample_missingness_final$Sample_Barcode
  )
)


###############################################################################
# STEP 28. ЗАГАЛЬНА КОРЕЛЯЦІЯ MISSINGNESS ТА MEAN BETA
###############################################################################

correlation_test <- cor.test(
  sample_missingness_analysis$Missing_percent,
  sample_missingness_analysis$Mean_beta,
  method = "spearman",
  exact = FALSE
)

cat("\n--- STEP 28. ЗАГАЛЬНА КОРЕЛЯЦІЯ ---\n")

print(correlation_test)

cat("\n--- ОПИСОВА СТАТИСТИКА ---\n")

print(
  summary(
    sample_missingness_analysis[
      , c("Missing_percent", "Mean_beta")
    ]
  )
)


###############################################################################
# STEP 29. КОРЕЛЯЦІЯ ОКРЕМО ДЛЯ PRIMARY ТА METASTATIC
###############################################################################

primary_data <- sample_missingness_analysis[
  sample_missingness_analysis$Sample_Type == "Primary Tumor",
  ,
  drop = FALSE
]

metastatic_data <- sample_missingness_analysis[
  sample_missingness_analysis$Sample_Type == "Metastatic",
  ,
  drop = FALSE
]

primary_correlation <- cor.test(
  primary_data$Missing_percent,
  primary_data$Mean_beta,
  method = "spearman",
  exact = FALSE
)

metastatic_correlation <- cor.test(
  metastatic_data$Missing_percent,
  metastatic_data$Mean_beta,
  method = "spearman",
  exact = FALSE
)

cat("\n--- PRIMARY TUMOR ---\n")

print(primary_correlation)

cat("\n--- METASTATIC ---\n")

print(metastatic_correlation)


###############################################################################
# STEP 30. КОРЕКЦІЯ P-VALUE ТА ЗБЕРЕЖЕННЯ КОРЕЛЯЦІЙ
###############################################################################

correlation_results <- data.frame(
  Sample_Type = c("Primary Tumor", "Metastatic"),
  Sample_count = c(
    nrow(primary_data),
    nrow(metastatic_data)
  ),
  Spearman_rho = c(
    unname(primary_correlation$estimate),
    unname(metastatic_correlation$estimate)
  ),
  P_value = c(
    primary_correlation$p.value,
    metastatic_correlation$p.value
  )
)

correlation_results$Adjusted_P_value <- p.adjust(
  correlation_results$P_value,
  method = "BH"
)

cat("\n--- STEP 30. КОРЕЛЯЦІЯ MISSINGNESS ТА MEAN BETA ---\n")

print(correlation_results, row.names = FALSE)

saveRDS(
  correlation_results,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_missingness_mean_beta_correlations.rds"
  )
)

write.csv(
  correlation_results,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_missingness_mean_beta_correlations.csv"
  ),
  row.names = FALSE
)


###############################################################################
# STEP 31. ГРУПОСПЕЦИФІЧНИЙ QC CpG
###############################################################################

primary_over20 <- sum(
  cpg_missing_by_type$Missing_percent_Primary > 20
)

metastatic_over20 <- sum(
  cpg_missing_by_type$Missing_percent_Metastatic > 20
)

over20_either_group <- sum(
  cpg_missing_by_type$Missing_percent_Primary > 20 |
    cpg_missing_by_type$Missing_percent_Metastatic > 20
)

eligible_both_groups <- sum(
  cpg_missing_by_type$Missing_percent_Primary <= 20 &
    cpg_missing_by_type$Missing_percent_Metastatic <= 20
)

group_filter_summary <- data.frame(
  Category = c(
    "Missingness >20% in Primary Tumor",
    "Missingness >20% in Metastatic",
    "Missingness >20% in either group",
    "Missingness <=20% in both groups"
  ),
  CpG_count = c(
    primary_over20,
    metastatic_over20,
    over20_either_group,
    eligible_both_groups
  )
)

group_filter_summary$Percentage <- round(
  group_filter_summary$CpG_count /
    nrow(cpg_missing_by_type) * 100,
  2
)

cat("\n--- STEP 31. ГРУПОСПЕЦИФІЧНИЙ QC CpG ---\n")

print(group_filter_summary, row.names = FALSE)

saveRDS(
  group_filter_summary,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_group_specific_CpG_QC.rds"
  )
)

write.csv(
  group_filter_summary,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_group_specific_CpG_QC.csv"
  ),
  row.names = FALSE
)


###############################################################################
# STEP 32. ПЕРЕВІРКА ЗБЕРЕЖЕННЯ ОСНОВНИХ QC-ФАЙЛІВ
###############################################################################

files_to_check <- c(
  "data/TCGA/processed/TCGA_SKCM_control_CpG_complete.rds",
  "results/tables/TCGA_SKCM_control_CpG_QC_report.rds",
  "results/tables/TCGA_SKCM_main_CpG_missingness_distribution.rds",
  "results/tables/TCGA_SKCM_sample_missingness_final_460.rds",
  "results/tables/TCGA_SKCM_CpG_missingness_by_sample_type.rds",
  "results/tables/TCGA_SKCM_missingness_mean_beta_correlations.rds",
  "data/TCGA/processed/TCGA_SKCM_strict_CpG_list_20pct_by_group.rds",
  "data/TCGA/processed/TCGA_SKCM_flagged_CpG_list_20pct_by_group.rds"
)

final_file_check <- file.exists(
  file.path(project_dir, files_to_check)
)

cat("\n--- STEP 32. ПЕРЕВІРКА QC-ФАЙЛІВ ---\n")

print(
  data.frame(
    File = files_to_check,
    Saved = final_file_check
  ),
  row.names = FALSE
)


###############################################################################
# STEP 33. ФІНАЛЬНА ПЕРЕВІРКА РОЗМІРУ ОСНОВНОЇ HDF5-МАТРИЦІ
###############################################################################

final_dimensions <- dim(beta_h5)

cat("\n--- STEP 33. ОСНОВНА HDF5-МАТРИЦЯ ---\n")

cat("Кількість CpG:", final_dimensions[1], "\n")
cat("Кількість зразків:", final_dimensions[2], "\n")

cat(
  "Основна матриця збережена без структурних змін.\n"
)

cat(
  "Імпутація не виконувалася.\n"
)


###############################################################################
# STEP 34. ФІНАЛЬНИЙ ЗВІТ QC
###############################################################################

final_qc_summary <- data.frame(
  Parameter = c(
    "Main CpG count",
    "Matched sample count",
    "Complete control CpG count",
    "Strict CpG count",
    "Flagged CpG count",
    "Overall missingness (%)",
    "Complete control missing values",
    "Imputation performed"
  ),
  Result = c(
    nrow(beta_h5),
    ncol(beta_h5),
    length(cpg_complete),
    length(cpg_strict),
    length(cpg_flagged),
    round(
      sum(sample_missing_count) /
        (nrow(beta_h5) * ncol(beta_h5)) * 100,
      4
    ),
    control_na,
    "No"
  )
)

cat("\n--- ФІНАЛЬНИЙ ЗВІТ QC ---\n")

print(final_qc_summary, row.names = FALSE)

saveRDS(
  final_qc_summary,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_methylation_QC_final_summary.rds"
  )
)

write.csv(
  final_qc_summary,
  file.path(
    project_dir, "results", "tables",
    "TCGA_SKCM_methylation_QC_final_summary.csv"
  ),
  row.names = FALSE
)

###############################################################################
# КІНЕЦЬ STEP 12–34
#
# Результат:
# - основна HDF5-матриця збережена;
# - контрольний набір CpG без пропусків створений;
# - суворий контрольний набір створений;
# - missingness перевірена на рівні CpG, зразків і груп;
# - кореляційний аналіз виконаний;
# - усі QC-результати збережені;
# - імпутація НЕ виконувалася.
###############################################################################
```

