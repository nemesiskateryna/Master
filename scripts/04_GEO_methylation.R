# Title: 04_GEO_methylation.R
# Project: Master


# 1. Project directory ----------------------------------------

project_dir <- getwd()

geo_dir <- file.path(
  project_dir,
  "data", "GEO", "GSE144487"
)

dir.create(
  geo_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. Download normalized beta matrix -------------------------

beta_file <- file.path(
  geo_dir,
  "GSE144487_Tumor_beta_normalized.txt.gz"
)

beta_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE144nnn/GSE144487/suppl/",
  "GSE144487_Tumor_beta_normalized.txt.gz"
)

if (!file.exists(beta_file)) {
  download.file(
    url = beta_url,
    destfile = beta_file,
    mode = "wb",
    method = "auto"
  )
}

# 3. Check downloaded file ------------------------------------

if (!file.exists(beta_file)) {
  stop("Beta matrix was not downloaded.")
}

file_size_gb <- file.info(beta_file)$size / 1024^3

cat(
  "Downloaded file:",
  basename(beta_file), "\n",
  "File size:",
  round(file_size_gb, 2), "GB\n"
)

# 4. Inspect file structure -----------------------------------

con <- gzfile(beta_file, open = "rt")

first_lines <- readLines(
  con,
  n = 5,
  warn = FALSE
)

close(con)

cat("\nFirst lines of the beta matrix:\n")
cat(first_lines, sep = "\n")


# 5. Match Illumina array IDs using high-confidence matching

metadata_file <- file.path(
  project_dir,
  "data", "GEO", "processed",
  "GSE144487_metadata_platform_1.csv"
)

matching_file <- file.path(
  project_dir,
  "results", "tables",
  "GEO_sample_matching_high_confidence.csv"
)

if (!file.exists(metadata_file)) {
  stop("GSE144487 metadata CSV not found.")
}

if (!file.exists(matching_file)) {
  stop("High-confidence matching file not found.")
}

metadata_met <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

GEO_sample_matching_high_confidence <- read.csv(
  matching_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Extract Illumina array ID from metadata title
metadata_met$array_id <- sub(
  ":.*$",
  "",
  metadata_met$title
)

# Extract Tumor ID
metadata_met$Tumor_ID <- sub(
  "^.*(Tum[0-9]+).*$",
  "\\1",
  metadata_met$source_name_ch1
)

# Read beta matrix column names
con <- gzfile(beta_file, open = "rt")

beta_header <- readLines(
  con,
  n = 1,
  warn = FALSE
)

close(con)

beta_array_ids <- strsplit(
  beta_header,
  "\t",
  fixed = TRUE
)[[1]][-1]

# Check high-confidence matching structure
cat("\n--- HIGH-CONFIDENCE MATCHING QC ---\n")

cat(
  "High-confidence pairs:",
  nrow(GEO_sample_matching_high_confidence), "\n"
)

cat(
  "Unique methylation GSM:",
  length(unique(GEO_sample_matching_high_confidence$GSM_met)), "\n"
)

cat(
  "Unique Tumor IDs:",
  length(unique(GEO_sample_matching_high_confidence$Tumor_ID)), "\n"
)

# Connect high-confidence GSMs to Illumina array IDs
high_confidence_met <- merge(
  GEO_sample_matching_high_confidence,
  metadata_met[
    ,
    c("geo_accession", "array_id", "Tumor_ID")
  ],
  by.x = "GSM_met",
  by.y = "geo_accession",
  all.x = TRUE,
  sort = FALSE
)

# Check that every matched sample has an array ID
if (any(is.na(high_confidence_met$array_id))) {
  stop("Some high-confidence GSMs have no matching Illumina array ID.")
}

# Check that all array IDs exist in beta matrix
if (any(!high_confidence_met$array_id %in% beta_array_ids)) {
  stop("Some high-confidence array IDs are absent from beta matrix.")
}

# Check duplicates
if (anyDuplicated(high_confidence_met$GSM_met)) {
  stop("Duplicated methylation GSMs in high-confidence matching.")
}

if (anyDuplicated(high_confidence_met$array_id)) {
  stop("Duplicated array IDs in high-confidence matching.")
}

# Reorder samples according to the beta matrix
high_confidence_met <- high_confidence_met[
  match(high_confidence_met$array_id, beta_array_ids),
]

cat(
  "High-confidence samples successfully mapped:",
  nrow(high_confidence_met), "\n"
)

cat("\n--- FIRST 10 MAPPED SAMPLES ---\n")
print(head(high_confidence_met, 10))

# 6. Import normalized beta matrix ---------------------------

cat("\n--- IMPORTING BETA MATRIX ---\n")

# Перевіряємо наявність бібліотеки data.table для швидкого та безпечного читання
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Пакет 'data.table' не встановлено. Будь ласка, встановіть його: install.packages('data.table')")
}

# Використовуємо data.table::fread для швидкого імпорту великих файлів (>1 ГБ)
# `data.table = FALSE` перетворює результат на звичайний DataFrame (якщо це потрібно для сумісності з іншими пакетами)
beta_data <- data.table::fread(
  beta_file,
  header = TRUE,
  sep = "\t",
  data.table = FALSE
)

# Якщо перший стовпчик містить назви CpG-проб, робимо його назвами рядків (row.names)
if ("CpG" %in% colnames(beta_data)[1] || is.character(beta_data[, 1])) {
  rownames(beta_data) <- beta_data[, 1]
  beta_data <- beta_data[, -1, drop = FALSE]
} else {
  # Якщо назви проб знаходяться у першій колонці без імені
  rownames(beta_data) <- beta_data[, 1]
  beta_data <- beta_data[, -1, drop = FALSE]
}

cat(
  "Full beta matrix dimensions:",
  nrow(beta_data), "CpG probes x",
  ncol(beta_data), "samples\n"
)
# Convert to numeric matrix
beta_matrix <- as.matrix(beta_data)
storage.mode(beta_matrix) <- "numeric"

rm(beta_data)
gc()

# 7. Select high-confidence samples --------------------------

selected_array_ids <- high_confidence_met$array_id

beta_matrix_high_confidence <- beta_matrix[
  ,
  selected_array_ids,
  drop = FALSE
]

colnames(beta_matrix_high_confidence) <-
  high_confidence_met$Expression_ID

cat(
  "\nSelected beta matrix dimensions:",
  nrow(beta_matrix_high_confidence), "CpG probes x",
  ncol(beta_matrix_high_confidence), "samples\n"
)

# 8. Beta matrix quality control -----------------------------

cat("\n--- BETA MATRIX QC ---\n")

# Dimensions
cat(
  "Matrix dimensions:",
  nrow(beta_matrix_high_confidence), "x",
  ncol(beta_matrix_high_confidence), "\n"
)

# CpG probe checks
cat(
  "Duplicated CpG probes:",
  sum(duplicated(rownames(beta_matrix_high_confidence))), "\n"
)

cat(
  "Missing values:",
  sum(is.na(beta_matrix_high_confidence)), "\n"
)

# Beta value range
cat(
  "Minimum beta value:",
  min(beta_matrix_high_confidence, na.rm = TRUE), "\n"
)

cat(
  "Maximum beta value:",
  max(beta_matrix_high_confidence, na.rm = TRUE), "\n"
)

# Values outside expected beta range
cat(
  "Values below 0:",
  sum(beta_matrix_high_confidence < 0, na.rm = TRUE), "\n"
)

cat(
  "Values above 1:",
  sum(beta_matrix_high_confidence > 1, na.rm = TRUE), "\n"
)

# Missing values per sample
missing_per_sample <- colSums(
  is.na(beta_matrix_high_confidence)
)

cat(
  "Samples with missing values:",
  sum(missing_per_sample > 0), "\n"
)

# Summary of beta values
beta_summary <- summary(
  as.vector(beta_matrix_high_confidence)
)

print(beta_summary)

# Sample-level QC table
sample_qc <- data.frame(
  Expression_ID = colnames(beta_matrix_high_confidence),
  GSM_met = high_confidence_met$GSM_met,
  Tumor_ID = high_confidence_met$Tumor_ID,
  missing_values = missing_per_sample,
  missing_percent = round(
    missing_per_sample / nrow(beta_matrix_high_confidence) * 100,
    4
  )
)

print(head(sample_qc, 10))


# 9. Save processed methylation data -------------------------

processed_dir <- file.path(
  project_dir,
  "data", "GEO", "processed"
)

tables_dir <- file.path(
  project_dir,
  "results", "tables"
)

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Save beta matrix
saveRDS(
  beta_matrix_high_confidence,
  file = file.path(
    processed_dir,
    "GSE144487_beta_matrix_high_confidence.rds"
  )
)

# Save matching information
write.csv(
  high_confidence_met,
  file = file.path(
    tables_dir,
    "GEO_high_confidence_methylation_mapping.csv"
  ),
  row.names = FALSE
)

# Save sample-level QC
write.csv(
  sample_qc,
  file = file.path(
    tables_dir,
    "GSE144487_methylation_sample_QC.csv"
  ),
  row.names = FALSE
)

# Save overall QC summary
qc_summary <- data.frame(
  metric = c(
    "CpG_probes",
    "Samples",
    "Duplicated_CpG",
    "Missing_values",
    "Samples_with_missing_values",
    "Minimum_beta",
    "Maximum_beta",
    "Mean_beta",
    "Median_beta"
  ),
  value = c(
    nrow(beta_matrix_high_confidence),
    ncol(beta_matrix_high_confidence),
    sum(duplicated(rownames(beta_matrix_high_confidence))),
    sum(is.na(beta_matrix_high_confidence)),
    sum(missing_per_sample > 0),
    min(beta_matrix_high_confidence),
    max(beta_matrix_high_confidence),
    mean(beta_matrix_high_confidence),
    median(beta_matrix_high_confidence)
  )
)

write.csv(
  qc_summary,
  file = file.path(
    tables_dir,
    "GSE144487_methylation_QC_summary.csv"
  ),
  row.names = FALSE
)

cat("\n--- SAVING COMPLETED ---\n")
cat("Beta matrix and QC files saved successfully.\n")
