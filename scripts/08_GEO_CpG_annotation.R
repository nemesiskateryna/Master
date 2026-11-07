# 08_GEO_CpG_annotation.R

#Purpose: Load aligned methylation matrix (GSE144487), match CpG probes 
#          with Illumina EPIC hg19 annotation, audit coverage, and save results.


# 1. PROJECT SETUP & DIRECTORIES ------------------------------

project_dir <- getwd()

if (!file.exists(file.path(project_dir, "Master.Rproj"))) {
  stop("Master.Rproj not found. Set the working directory to the project root.")
}

processed_dir <- file.path(project_dir, "data", "GEO", "processed")
tables_dir    <- file.path(project_dir, "results", "tables")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir,    recursive = TRUE, showWarnings = FALSE)

cat("\n--- PROJECT SETUP COMPLETED ---\n")


# 2. LOAD ALIGNED METHYLATION MATRIX --------------------------

aligned_beta_file <- file.path(processed_dir, "GSE144487_beta_matrix_high_confidence_aligned.rds")

if (!file.exists(aligned_beta_file)) {
  stop("Aligned methylation matrix not found: ", aligned_beta_file)
}

cat("Loading aligned methylation matrix...\n")
beta_matrix <- readRDS(aligned_beta_file)

cpg_ids <- rownames(beta_matrix)

cat("\n--- METHYLATION MATRIX SUMMARY ---\n")
cat("Total CpG probes:", length(cpg_ids), "\n")
cat("Total samples:   ", ncol(beta_matrix), "\n")

# Перевірка формату ідентифікаторів зондів (усі мають починатися з "cg")
invalid_ids <- sum(!grepl("^cg[0-9]+$", cpg_ids))
if (invalid_ids > 0) {
  warning(sprintf("Detected %d non-standard CpG probe IDs.", invalid_ids))
} else {
  cat("[OK] All probe IDs follow the standard 'cgXXXX' format.\n")
}


# 3. LOAD EPIC hg19 ANNOTATION PACKAGE ------------------------

package_name <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"

if (!requireNamespace(package_name, quietly = TRUE)) {
  cat("\nInstalling required EPIC annotation package...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(package_name, ask = FALSE, update = FALSE)
}

library(minfi)
library(package_name, character.only = TRUE)

cat("\nExtracting EPIC annotation (hg19)...\n")
epic_annotation <- as.data.frame(minfi::getAnnotation(get(package_name)))

cat("Annotation dimensions:", nrow(epic_annotation), "probes x", ncol(epic_annotation), "columns\n")


# 4. MATCH METHYLATION PROBES WITH EPIC ANNOTATION ------------

annotation_ids <- rownames(epic_annotation)

matched_ids <- intersect(cpg_ids, annotation_ids)
missing_ids <- setdiff(cpg_ids, annotation_ids)

cat("\n--- PROBE MATCHING AUDIT ---\n")
cat("Matrix CpGs:       ", length(cpg_ids), "\n")
cat("Annotation CpGs:   ", length(annotation_ids), "\n")
cat("Successfully matched:", length(matched_ids), "\n")
cat("Unmatched CpGs:    ", length(missing_ids), "\n")
cat(sprintf("Match rate:        %.3f%%\n", 100 * length(matched_ids) / length(cpg_ids)))


# 5. FILTER & EXTRACT RELEVANT ANNOTATION FIELDS --------------

required_columns <- c(
  "chr", "pos", "strand", 
  "UCSC_RefGene_Name", "UCSC_RefGene_Group", 
  "Islands_Name", "Relation_to_Island"
)

# Перевірка наявності всіх необхідних колонок
missing_cols <- setdiff(required_columns, colnames(epic_annotation))
if (length(missing_cols) > 0) {
  stop("Missing required columns in EPIC annotation: ", paste(missing_cols, collapse = ", "))
}

# Формуємо відфільтровану таблицю анотації в точному порядку зондів матриці
matched_annotation <- epic_annotation[match(cpg_ids, annotation_ids), required_columns, drop = FALSE]

# Перейменовуємо рядок, щоб явно містив Probe_ID
matched_annotation$Probe_ID <- cpg_ids


# 6. SAVE ANNOTATION AND QC METRICS ---------------------------

annotation_output_file <- file.path(processed_dir, "GSE144487_CpG_annotation_matched.rds")
saveRDS(matched_annotation, annotation_output_file)

# Формуємо підсумковий звіт контролю якості (QC) для дипломної роботи
cpg_qc_summary <- data.frame(
  Metric = c(
    "Total input CpG probes",
    "Matched with EPIC annotation",
    "Unmatched probes",
    "Matching percentage",
    "Genomic coordinates available",
    "RefGene annotation available"
  ),
  Value = c(
    length(cpg_ids),
    length(matched_ids),
    length(missing_ids),
    round(100 * length(matched_ids) / length(cpg_ids), 3),
    sum(!is.na(matched_annotation$chr) & !is.na(matched_annotation$pos)),
    sum(!is.na(matched_annotation$UCSC_RefGene_Name) & matched_annotation$UCSC_RefGene_Name != "")
  ),
  stringsAsFactors = FALSE
)

write.csv(cpg_qc_summary, file.path(tables_dir, "GSE144487_CpG_annotation_QC.csv"), row.names = FALSE)

cat("\n============================================\n")
cat("CpG ANNOTATION AND AUDIT COMPLETED SUCCESSFULLY\n")
cat("Matched annotation saved to:", annotation_output_file, "\n")
cat("QC summary report saved to:  results/tables/GSE144487_CpG_annotation_QC.csv\n")
cat("============================================\n")








