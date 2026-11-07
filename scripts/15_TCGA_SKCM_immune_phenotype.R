
# ============================================================
# TCGA-SKCM
# ВИЗНАЧЕННЯ ІМУННОГО ФЕНОТИПУ ПУХЛИН
#
# Мета:
# 1. Проаналізувати розподіл TIS (Tumor Inflammation Signature).
# 2. Визначити групи Immune-high та Immune-low.
# 3. Дослідити зв'язок імунного фенотипу з типом пухлини.
# 4. Порівняти TIS між первинними та метастатичними пухлинами.
# 5. Зберегти результати статистичного аналізу та візуалізації.
#
# Вхідні дані:
# immune_scores — таблиця з розрахованими імунними показниками.
# Необхідні колонки:
#   Sample_ID
#   TIS_mean_zscore
#
# Вихідні файли:
# data/TCGA/processed/
# results/tables/
# results/figures/
# ============================================================


# ============================================================
# 0. ПІДГОТОВКА РОБОЧОГО СЕРЕДОВИЩА
# ============================================================

cat("\n============================================================\n")
cat("TCGA-SKCM: АНАЛІЗ ІМУННОГО ФЕНОТИПУ\n")
cat("============================================================\n")

cat("\nРобоча директорія:\n")
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

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. ПЕРЕВІРКА ВХІДНИХ ДАНИХ
# ============================================================

# Об'єкт immune_scores має бути сформований на попередньому
# етапі аналізу імунних сигнатур.

if (!exists("immune_scores")) {
  stop(
    paste(
      "Об'єкт immune_scores не знайдено.",
      "Спочатку завантажте або сформуйте таблицю імунних показників."
    )
  )
}

# Перевірка наявності необхідних колонок
required_columns <- c(
  "Sample_ID",
  "TIS_mean_zscore"
)

missing_columns <- setdiff(
  required_columns,
  colnames(immune_scores)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Відсутні необхідні колонки:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# Перевірка структури таблиці
cat("\n--- СТРУКТУРА ВХІДНИХ ДАНИХ ---\n")

cat("Кількість зразків:", nrow(immune_scores), "\n")
cat("Кількість змінних:", ncol(immune_scores), "\n")

# Перевірка ідентифікаторів
if (anyNA(immune_scores$Sample_ID) ||
    any(immune_scores$Sample_ID == "")) {
  stop("Виявлено відсутні ідентифікатори зразків.")
}

if (anyDuplicated(immune_scores$Sample_ID) > 0) {
  stop("Виявлено дубльовані Sample_ID.")
}

# Перевірка значень TIS
if (!is.numeric(immune_scores$TIS_mean_zscore)) {
  stop("TIS_mean_zscore повинен бути числовою змінною.")
}

if (any(!is.finite(immune_scores$TIS_mean_zscore))) {
  stop("У TIS_mean_zscore виявлено NA або нескінченні значення.")
}

tis_values <- immune_scores$TIS_mean_zscore


# ============================================================
# 2. ОПИСОВА СТАТИСТИКА TIS
# ============================================================

cat("\n============================================================\n")
cat("2. ОПИСОВА СТАТИСТИКА TIS\n")
cat("============================================================\n")

# Основні статистичні характеристики
tis_summary <- summary(tis_values)

tis_sd <- sd(tis_values)

tis_skewness <- mean(
  (tis_values - mean(tis_values))^3
) / sd(tis_values)^3

tis_median <- median(tis_values)

cat("\nРозподіл значень TIS:\n")
print(tis_summary)

cat("\nСтандартне відхилення:", tis_sd, "\n")
cat("Асиметрія розподілу:", tis_skewness, "\n")
cat("Медіана TIS:", tis_median, "\n")

# Формування таблиці описової статистики
tis_descriptive_statistics <- data.frame(
  Statistic = c(
    "Кількість зразків",
    "Мінімум",
    "Перший квартиль",
    "Медіана",
    "Середнє",
    "Третій квартиль",
    "Максимум",
    "Стандартне відхилення",
    "Асиметрія"
  ),
  Value = c(
    length(tis_values),
    min(tis_values),
    unname(quantile(tis_values, 0.25)),
    median(tis_values),
    mean(tis_values),
    unname(quantile(tis_values, 0.75)),
    max(tis_values),
    tis_sd,
    tis_skewness
  )
)

write.csv(
  tis_descriptive_statistics,
  "results/tables/TCGA_SKCM_TIS_descriptive_statistics.csv",
  row.names = FALSE
)


# ============================================================
# 3. ВІЗУАЛІЗАЦІЯ РОЗПОДІЛУ TIS
# ============================================================

# 3.1. Гістограма
png(
  "results/figures/TCGA_SKCM_TIS_histogram.png",
  width = 1800,
  height = 1400,
  res = 200
)

hist(
  tis_values,
  breaks = 30,
  main = "Розподіл TIS у когорті TCGA-SKCM",
  xlab = "Середній z-показник TIS",
  ylab = "Кількість зразків",
  col = "lightblue",
  border = "white"
)

abline(
  v = tis_median,
  col = "red",
  lwd = 2,
  lty = 2
)

dev.off()


# 3.2. Графік щільності розподілу
png(
  "results/figures/TCGA_SKCM_TIS_density.png",
  width = 1800,
  height = 1400,
  res = 200
)

plot(
  density(tis_values),
  main = "Щільність розподілу TIS у TCGA-SKCM",
  xlab = "Середній z-показник TIS",
  ylab = "Щільність імовірності",
  lwd = 2
)

abline(
  v = tis_median,
  col = "red",
  lwd = 2,
  lty = 2
)

dev.off()


# ============================================================
# 4. КЛАСИФІКАЦІЯ ІМУННОГО ФЕНОТИПУ
# ============================================================

# Як порогове значення використовується медіана TIS
# у досліджуваній когорті.

tis_threshold <- median(
  tis_values
)

cat("\n============================================================\n")
cat("4. КЛАСИФІКАЦІЯ ІМУННОГО ФЕНОТИПУ\n")
cat("============================================================\n")

cat("Порогове значення TIS:", tis_threshold, "\n")

# ВАЖЛИВО:
# Вищий TIS відповідає вищій імунній активності.
#
# TIS > медіани  -> Immune-high
# TIS <= медіани -> Immune-low

immune_scores$immune_group <- ifelse(
  immune_scores$TIS_mean_zscore > tis_threshold,
  "Immune-high",
  "Immune-low"
)

# Фіксація порядку категорій
immune_scores$immune_group <- factor(
  immune_scores$immune_group,
  levels = c(
    "Immune-high",
    "Immune-low"
  )
)

# Перевірка розподілу груп
immune_group_table <- table(
  immune_scores$immune_group
)

cat("\n--- РОЗПОДІЛ ІМУННИХ ГРУП ---\n")
print(immune_group_table)

immune_group_percentages <- round(
  prop.table(immune_group_table) * 100,
  2
)

cat("\nВідсотковий розподіл:\n")
print(immune_group_percentages)

# Таблиця розподілу
immune_group_summary <- data.frame(
  Immune_group = names(immune_group_table),
  Number_of_samples = as.integer(immune_group_table),
  Percentage = as.numeric(immune_group_percentages)
)

write.csv(
  immune_group_summary,
  "results/tables/TCGA_SKCM_immune_group_distribution.csv",
  row.names = FALSE
)


# ============================================================
# 5. ВИЗНАЧЕННЯ ТИПУ ПУХЛИННОГО ЗРАЗКА
# ============================================================

# Вилучення коду типу зразка з повного TCGA barcode.
#
# 01 — первинна пухлина
# 06 — метастатична пухлина

immune_scores$sample_type_code <- substr(
  as.character(immune_scores$Sample_ID),
  14,
  15
)

# Перевірка допустимих кодів
if (!all(
  immune_scores$sample_type_code %in% c("01", "06")
)) {
  stop(
    "Виявлено зразки з невідомими кодами типу пухлини."
  )
}

# Перетворення кодів на текстові категорії
immune_scores$sample_type <- factor(
  immune_scores$sample_type_code,
  levels = c("01", "06"),
  labels = c("Primary", "Metastatic")
)

cat("\n--- РОЗПОДІЛ ТИПІВ ПУХЛИН ---\n")

print(
  table(immune_scores$sample_type)
)


# ============================================================
# 6. ЗВ'ЯЗОК ІМУННОГО ФЕНОТИПУ З ТИПОМ ПУХЛИНИ
# ============================================================

cat("\n============================================================\n")
cat("6. ІМУННИЙ ФЕНОТИП ТА ТИП ПУХЛИНИ\n")
cat("============================================================\n")

# Таблиця спряженості
immune_sample_table <- table(
  immune_scores$immune_group,
  immune_scores$sample_type
)

cat("\n--- ТАБЛИЦЯ СПРЯЖЕНОСТІ ---\n")
print(immune_sample_table)

# Частки всередині кожної імунної групи
immune_sample_percentages <- round(
  prop.table(
    immune_sample_table,
    margin = 1
  ) * 100,
  2
)

cat("\n--- ВІДСОТКИ ВСЕРЕДИНІ ІМУННИХ ГРУП ---\n")
print(immune_sample_percentages)

# Збереження таблиць
write.csv(
  as.data.frame.matrix(immune_sample_table),
  "results/tables/TCGA_SKCM_immune_group_by_sample_type_counts.csv"
)

write.csv(
  as.data.frame.matrix(immune_sample_percentages),
  "results/tables/TCGA_SKCM_immune_group_by_sample_type_percentages.csv"
)

# Критерій хі-квадрат Пірсона з поправкою Єйтса
chi_square_test <- chisq.test(
  immune_sample_table,
  correct = TRUE
)

cat("\n--- КРИТЕРІЙ ХІ-КВАДРАТ ---\n")
print(chi_square_test)

# Очікувані частоти
cat("\nОчікувані частоти:\n")
print(chi_square_test$expected)

# Розмір ефекту: коефіцієнт Phi для таблиці 2 x 2
phi_coefficient <- sqrt(
  as.numeric(chi_square_test$statistic) /
    sum(immune_sample_table)
)

cat("\nКоефіцієнт Phi:", phi_coefficient, "\n")

# Збереження результату тесту
chi_square_summary <- data.frame(
  Test = "Chi-squared test with Yates correction",
  Chi_squared = as.numeric(chi_square_test$statistic),
  Degrees_of_freedom = as.numeric(chi_square_test$parameter),
  P_value = chi_square_test$p.value,
  Phi_coefficient = phi_coefficient
)

write.csv(
  chi_square_summary,
  "results/tables/TCGA_SKCM_immune_sample_type_chisquare.csv",
  row.names = FALSE
)


# ============================================================
# 7. ПОРІВНЯННЯ TIS МІЖ ПЕРВИННИМИ ТА МЕТАСТАТИЧНИМИ ПУХЛИНАМИ
# ============================================================

cat("\n============================================================\n")
cat("7. ПОРІВНЯННЯ TIS ЗА ТИПОМ ПУХЛИНИ\n")
cat("============================================================\n")

# Виділення TIS для двох категорій
tis_primary <- immune_scores$TIS_mean_zscore[
  immune_scores$sample_type == "Primary"
]

tis_metastatic <- immune_scores$TIS_mean_zscore[
  immune_scores$sample_type == "Metastatic"
]

# Описова статистика
tis_by_sample_type <- data.frame(
  Sample_type = c("Primary", "Metastatic"),
  N = c(
    length(tis_primary),
    length(tis_metastatic)
  ),
  Mean = c(
    mean(tis_primary),
    mean(tis_metastatic)
  ),
  Median = c(
    median(tis_primary),
    median(tis_metastatic)
  ),
  SD = c(
    sd(tis_primary),
    sd(tis_metastatic)
  ),
  Minimum = c(
    min(tis_primary),
    min(tis_metastatic)
  ),
  Maximum = c(
    max(tis_primary),
    max(tis_metastatic)
  )
)

cat("\n--- ОПИСОВА СТАТИСТИКА TIS ЗА ТИПОМ ПУХЛИНИ ---\n")
print(tis_by_sample_type, row.names = FALSE)

write.csv(
  tis_by_sample_type,
  "results/tables/TCGA_SKCM_TIS_by_sample_type.csv",
  row.names = FALSE
)

# Непараметричне порівняння двох незалежних груп
# Критерій Вілкоксона (Mann–Whitney U)

tis_sample_type_test <- wilcox.test(
  tis_primary,
  tis_metastatic,
  exact = FALSE
)

cat("\n--- КРИТЕРІЙ ВІЛКОКСОНА ---\n")
print(tis_sample_type_test)

# Збереження результату
wilcoxon_summary <- data.frame(
  Test = "Wilcoxon rank-sum test",
  W_statistic = as.numeric(tis_sample_type_test$statistic),
  P_value = tis_sample_type_test$p.value,
  Primary_median = median(tis_primary),
  Metastatic_median = median(tis_metastatic),
  Primary_mean = mean(tis_primary),
  Metastatic_mean = mean(tis_metastatic)
)

write.csv(
  wilcoxon_summary,
  "results/tables/TCGA_SKCM_TIS_wilcoxon_test.csv",
  row.names = FALSE
)


# ============================================================
# 8. ВІЗУАЛІЗАЦІЯ TIS ЗА ТИПОМ ПУХЛИНИ
# ============================================================

png(
  "results/figures/TCGA_SKCM_TIS_by_sample_type.png",
  width = 1800,
  height = 1400,
  res = 200
)

boxplot(
  TIS_mean_zscore ~ sample_type,
  data = immune_scores,
  main = "TIS за типом пухлинного зразка",
  xlab = "Тип пухлини",
  ylab = "Середній z-показник TIS",
  col = c("lightblue", "lightpink")
)

# Горизонтальна лінія медіани всієї когорти
abline(
  h = tis_threshold,
  lty = 2,
  lwd = 2,
  col = "red"
)

dev.off()


# ============================================================
# 9. ФОРМУВАННЯ ПІДСУМКОВОЇ МЕТАДАНИХ КЛАСИФІКАЦІЇ
# ============================================================

immune_classification_metadata <- list(
  
  # Інформація про когорту
  cohort = "TCGA-SKCM",
  number_of_samples = nrow(immune_scores),
  
  # Показник імунної активності
  primary_score = "TIS_mean_zscore",
  
  # Метод класифікації
  threshold = tis_threshold,
  threshold_method = "Cohort median",
  
  classification = c(
    "Immune-high: TIS > cohort median",
    "Immune-low: TIS <= cohort median"
  ),
  
  # Розподіл груп
  group_distribution = immune_group_table,
  
  # Розподіл за типом пухлини
  sample_type_distribution = table(
    immune_scores$sample_type
  ),
  
  immune_group_by_sample_type = immune_sample_table,
  
  # Статистичні результати
  association_test = "Chi-squared test with Yates correction",
  association_p_value = chi_square_test$p.value,
  association_phi = phi_coefficient,
  
  sample_type_comparison = "Wilcoxon rank-sum test",
  sample_type_comparison_p_value = tis_sample_type_test$p.value,
  
  # Описові характеристики
  TIS_summary = tis_summary,
  TIS_standard_deviation = tis_sd,
  TIS_skewness = tis_skewness
)


# ============================================================
# 10. ЗБЕРЕЖЕННЯ РЕЗУЛЬТАТІВ
# ============================================================

# Збереження таблиці з імунними показниками та групами
saveRDS(
  immune_scores,
  "data/TCGA/processed/TCGA_SKCM_immune_scores.rds"
)

write.csv(
  immune_scores,
  "results/tables/TCGA_SKCM_immune_scores.csv",
  row.names = FALSE
)

# Збереження параметрів класифікації
saveRDS(
  immune_classification_metadata,
  "data/TCGA/processed/TCGA_SKCM_immune_classification_metadata.rds"
)


# ============================================================
# 11. ФІНАЛЬНА ПЕРЕВІРКА
# ============================================================

cat("\n============================================================\n")
cat("ФІНАЛЬНА ПЕРЕВІРКА РЕЗУЛЬТАТІВ\n")
cat("============================================================\n")

# Перевірка відсутності пропущених груп
stopifnot(
  !anyNA(immune_scores$immune_group)
)

# Перевірка кількості зразків
stopifnot(
  nrow(immune_scores) ==
    sum(immune_group_table)
)

# Перевірка відповідності класифікації TIS
stopifnot(
  all(
    as.character(immune_scores$immune_group[
      immune_scores$TIS_mean_zscore > tis_threshold
    ]) == "Immune-high"
  )
)

stopifnot(
  all(
    as.character(immune_scores$immune_group[
      immune_scores$TIS_mean_zscore <= tis_threshold
    ]) == "Immune-low"
  )
)

# Перевірка збережених файлів
output_files <- c(
  "data/TCGA/processed/TCGA_SKCM_immune_scores.rds",
  "data/TCGA/processed/TCGA_SKCM_immune_classification_metadata.rds",
  "results/tables/TCGA_SKCM_immune_scores.csv",
  "results/tables/TCGA_SKCM_TIS_descriptive_statistics.csv",
  "results/tables/TCGA_SKCM_immune_group_distribution.csv",
  "results/tables/TCGA_SKCM_immune_group_by_sample_type_counts.csv",
  "results/tables/TCGA_SKCM_immune_group_by_sample_type_percentages.csv",
  "results/tables/TCGA_SKCM_immune_sample_type_chisquare.csv",
  "results/tables/TCGA_SKCM_TIS_by_sample_type.csv",
  "results/tables/TCGA_SKCM_TIS_wilcoxon_test.csv",
  "results/figures/TCGA_SKCM_TIS_histogram.png",
  "results/figures/TCGA_SKCM_TIS_density.png",
  "results/figures/TCGA_SKCM_TIS_by_sample_type.png"
)

output_status <- data.frame(
  File = output_files,
  Exists = file.exists(output_files)
)

cat("\n--- СТАТУС ЗБЕРЕЖЕННЯ ФАЙЛІВ ---\n")
print(output_status, row.names = FALSE)

if (!all(output_status$Exists)) {
  warning("Деякі вихідні файли не були збережені.")
}

cat("\n--- ПІДСУМКОВИЙ РОЗПОДІЛ ІМУННИХ ГРУП ---\n")
print(immune_group_table)

cat("\nПорогове значення TIS:", tis_threshold, "\n")

cat(
  "P-value (зв'язок імунного фенотипу з типом пухлини):",
  chi_square_test$p.value,
  "\n"
)

cat(
  "P-value (порівняння TIS між типами пухлин):",
  tis_sample_type_test$p.value,
  "\n"
)

cat("\n============================================================\n")
cat("АНАЛІЗ ІМУННОГО ФЕНОТИПУ TCGA-SKCM ЗАВЕРШЕНО\n")
cat("============================================================\n")
