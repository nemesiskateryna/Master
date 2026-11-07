# 10_GEO_Final_Integration.R
# 
# Purpose: Final integration and QC of multi-omics data 
#          (DNA Methylation GSE144487 + Gene Expression GSE65904)
#          for the 74 matched melanoma patient samples.
# ============================================================


# 1. PROJECT SETUP & DIRECTORIES ------------------------------

project_dir <- getwd()

if (!file.exists(file.path(project_dir, "Master.Rproj"))) {
  stop("Master.Rproj not found. Set the working directory to the project root.")
}

processed_dir <- file.path(project_dir, "data", "GEO", "processed")
tables_dir    <- file.path(project_dir, "results", "tables")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n--- PROJECT SETUP COMPLETED ---\n")


# 2. LOAD GEO INTEGRATED DATA COMPONENTS ----------------------

cat("Loading processed matrices, annotations, and metadata...\n")

# А. ДНК-метилювання та його анотація
beta_matrix <- readRDS(file.path(processed_dir, "GSE144487_beta_matrix_high_confidence_aligned.rds"))
methylation_annotation <- readRDS(file.path(processed_dir, "GSE144487_CpG_annotation.rds"))

# Б. Експресія генів та її анотація
expr_matrix <- readRDS(file.path(processed_dir, "GSE65904_expression_gene_level_high_confidence_log2_aligned.rds"))
gene_annotation <- readRDS(file.path(processed_dir, "GSE65904_gene_annotation_high_confidence.rds"))

# В. Таблиця вирівнювання зразків
matching_aligned <- read.csv(file.path(tables_dir, "GEO_sample_matching_high_confidence_aligned.csv"), stringsAsFactors = FALSE)

# Г. Клінічні метадані
meta_met  <- readRDS(file.path(processed_dir, "GSE144487_all_metadata.rds"))
meta_expr <- readRDS(file.path(processed_dir, "GSE65904_all_metadata.rds"))


# 3. SUBSET CLINICAL METADATA FOR THE 74 MATCHED SAMPLES ----

matched_met <- meta_met[
  match(matching_aligned$GSM_met, meta_met$geo_accession),
  , drop = FALSE
]

matched_expr <- meta_expr[
  match(matching_aligned$GSM_expr, meta_expr$geo_accession),
  , drop = FALSE
]

cat("Matched methylation clinical records:", nrow(matched_met), "\n")
cat("Matched expression clinical records: ", nrow(matched_expr), "\n")


# 4. COMPREHENSIVE ALIGNMENT & QC AUDIT -----------------------

cat("\n--- RUNNING INTEGRATION QC AUDIT ---\n")

cat("Methylation matrix dimensions:", paste(dim(beta_matrix), collapse = " x "), "\n")
cat("Expression matrix dimensions: ", paste(dim(expr_matrix), collapse = " x "), "\n")

# Перевірка ідентичності колонок (зразків)
samples_identical <- identical(colnames(beta_matrix), colnames(expr_matrix))
matching_identical <- identical(colnames(beta_matrix), matching_aligned$Expression_ID)

cat("Methylation and expression sample IDs identical:", samples_identical, "\n")
cat("Sample IDs identical to matching table IDs:     ", matching_identical, "\n")

if (!samples_identical || !matching_identical) {
  stop("CRITICAL ERROR: Sample IDs or column order do not match between modalities!")
}

cat("Missing GSE144487 metadata entries:", sum(is.na(matched_met$geo_accession)), "\n")
cat("Missing GSE65904 metadata entries: ", sum(is.na(matched_expr$geo_accession)), "\n")


# 5. COMPILE CENTRAL MULTI-OMICS OBJECT -----------------------

geo_integrated <- list(
  beta_matrix              = beta_matrix,
  expression_matrix        = expr_matrix,
  methylation_annotation   = methylation_annotation,
  gene_annotation          = gene_annotation,
  clinical_GSE144487       = matched_met,
  clinical_GSE65904        = matched_expr,
  sample_matching          = matching_aligned,
  metadata = list(
    project              = "Melanoma Multi-omics",
    methylation_dataset  = "GSE144487",
    expression_dataset   = "GSE65904",
    n_samples            = 74,
    methylation_platform = "Illumina Infinium HumanMethylationEPIC",
    expression_platform  = "Illumina HumanHT-12 V4.0",
    genome_build         = "hg19",
    matching_status      = "Pre-matched high-confidence cohort",
    integration_status   = "Sample IDs and matrix order fully verified"
  )
)

# Збереження фінального об'єкта
output_rds_path <- file.path(processed_dir, "GEO_integrated_multiomics_high_confidence.rds")
saveRDS(geo_integrated, file = output_rds_path)

cat("\n[OK] Central GEO multi-omics object successfully saved to:\n", output_rds_path, "\n")


# 6. FINAL INTEGRATION SUMMARY --------------------------------

cat("\n============================================\n")
cat("FINAL INTEGRATION SUMMARY REPORT\n")
cat("Total samples integrated:       ", geo_integrated$metadata$n_samples, "\n")
cat("CpG probes (GSE144487):         ", nrow(geo_integrated$beta_matrix), "\n")
cat("Genes / Probes (GSE65904):      ", nrow(geo_integrated$expression_matrix), "\n")
cat("Saved object file size (MB):    ", round(file.info(output_rds_path)$size / (1024^2), 2), "\n")
cat("============================================\n")