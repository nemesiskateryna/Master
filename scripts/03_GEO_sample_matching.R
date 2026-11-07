# ============================================================
# 03_GEO_sample_matching.R
#
# Project: Master
# Study: Melanoma Multi-omics
#
# Purpose:
# Match DNA methylation samples from GSE144487 with
# gene expression samples from GSE65904 using numeric IDs,
# followed by clinical quality control, stage harmonization,
# and high-confidence sample filtering.
#
# Matching assumption:
# Tum1 <-> LUM1
# Tum2 <-> LUM2
# ...
# Tum196 <-> LUM196
#
# IMPORTANT:
# Matching is based on numeric suffixes and is treated as
# an analytical assumption, validated via clinical stages.
#
# Input:
# data/GEO/processed/GSE144487_metadata_platform_1.csv
# data/GEO/processed/GSE65904_metadata_platform_1.csv
#
# Output:
# results/tables/GEO_sample_matching.csv
# results/tables/GEO_sample_matching.rds
# results/tables/GEO_unmatched_expression_samples.csv
# results/tables/GEO_unmatched_methylation_samples.csv
# results/tables/GEO_sample_matching_QC.csv
# results/tables/GEO_sample_matching_high_confidence.csv
# results/tables/GEO_sample_matching_high_confidence.rds
# ============================================================


# 1. PROJECT DIRECTORIES --------------------------------------

project_dir <- getwd()

if (!dir.exists(file.path(project_dir, "data", "GEO", "processed"))) {
  stop(
    "GEO processed directory not found. ",
    "Please set the working directory to the project root."
  )
}

input_dir <- file.path(
  project_dir,
  "data", "GEO", "processed"
)

output_dir <- file.path(
  project_dir,
  "results", "tables"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. INPUT FILES ----------------------------------------------

methylation_file <- file.path(
  input_dir,
  "GSE144487_metadata_platform_1.csv"
)

expression_file <- file.path(
  input_dir,
  "GSE65904_metadata_platform_1.csv"
)

if (!file.exists(methylation_file)) {
  stop("Methylation metadata file not found: ", methylation_file)
}

if (!file.exists(expression_file)) {
  stop("Expression metadata file not found: ", expression_file)
}


# 3. LOAD GEO METADATA ----------------------------------------

metadata_met <- read.csv(
  methylation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata_expr <- read.csv(
  expression_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

message("GSE144487 metadata loaded: ", nrow(metadata_met), " samples")
message("GSE65904 metadata loaded: ", nrow(metadata_expr), " samples")


# 4. CHECK REQUIRED COLUMNS -----------------------------------

required_met <- c(
  "geo_accession",
  "source_name_ch1",
  "title"
)

required_expr <- c(
  "geo_accession",
  "title"
)

missing_met <- setdiff(
  required_met,
  colnames(metadata_met)
)

missing_expr <- setdiff(
  required_expr,
  colnames(metadata_expr)
)

if (length(missing_met) > 0) {
  stop(
    "Missing methylation metadata columns: ",
    paste(missing_met, collapse = ", ")
  )
}

if (length(missing_expr) > 0) {
  stop(
    "Missing expression metadata columns: ",
    paste(missing_expr, collapse = ", ")
  )
}


# 5. EXTRACT SAMPLE IDENTIFIERS -------------------------------

# DNA methylation:
# Extract numeric suffix from Tum1, Tum2, ..., Tum196.

metadata_met$Tumor_ID <- trimws(
  as.character(metadata_met$source_name_ch1)
)

met_number_match <- regexec(
  "^Tum[[:space:]_-]*([0-9]+)$",
  metadata_met$Tumor_ID,
  ignore.case = TRUE
)

met_number_parts <- regmatches(
  metadata_met$Tumor_ID,
  met_number_match
)

metadata_met$sample_number <- vapply(
  met_number_parts,
  function(x) {
    if (length(x) == 2L) {
      as.integer(x[2])
    } else {
      NA_integer_
    }
  },
  integer(1)
)


# Gene expression:
# Extract patient number from titles such as
# "Melanoma tumor sample from patient 1".
#
# If a patient column with LUM identifiers exists,
# use it as the preferred source.

expr_number_from_title <- regexec(
  "patient[[:space:]_-]*([0-9]+)",
  metadata_expr$title,
  ignore.case = TRUE
)

expr_title_parts <- regmatches(
  metadata_expr$title,
  expr_number_from_title
)

expr_number_title <- vapply(
  expr_title_parts,
  function(x) {
    if (length(x) == 2L) {
      as.integer(x[2])
    } else {
      NA_integer_
    }
  },
  integer(1)
)

expr_number_patient <- rep(
  NA_integer_,
  nrow(metadata_expr)
)

if ("patient" %in% colnames(metadata_expr)) {
  
  patient_ids <- trimws(
    as.character(metadata_expr$patient)
  )
  
  patient_match <- regexec(
    "^LUM[[:space:]_-]*([0-9]+)$",
    patient_ids,
    ignore.case = TRUE
  )
  
  patient_parts <- regmatches(
    patient_ids,
    patient_match
  )
  
  expr_number_patient <- vapply(
    patient_parts,
    function(x) {
      if (length(x) == 2L) {
        as.integer(x[2])
      } else {
        NA_integer_
      }
    },
    integer(1)
  )
}

metadata_expr$sample_number <- ifelse(
  !is.na(expr_number_patient),
  expr_number_patient,
  expr_number_title
)

metadata_expr$sample_number <- as.integer(
  metadata_expr$sample_number
)

metadata_expr$Expression_ID <- ifelse(
  !is.na(metadata_expr$sample_number),
  paste0("LUM", metadata_expr$sample_number),
  NA_character_
)


# 6. VALIDATE SAMPLE IDENTIFIERS ------------------------------

if (anyNA(metadata_met$sample_number)) {
  stop(
    "Some methylation samples do not have valid Tum<number> IDs."
  )
}

if (anyNA(metadata_expr$sample_number)) {
  stop(
    "Some expression samples do not have valid patient/LUM IDs."
  )
}

if (anyDuplicated(metadata_met$sample_number)) {
  stop("Duplicated methylation sample numbers detected.")
}

if (anyDuplicated(metadata_expr$sample_number)) {
  stop("Duplicated expression sample numbers detected.")
}

if (anyDuplicated(metadata_met$geo_accession)) {
  stop("Duplicated methylation GEO accessions detected.")
}

if (anyDuplicated(metadata_expr$geo_accession)) {
  stop("Duplicated expression GEO accessions detected.")
}


# 7. CREATE MAPPING TABLES ------------------------------------

met_mapping <- data.frame(
  GSM_met = as.character(metadata_met$geo_accession),
  Tumor_ID = metadata_met$Tumor_ID,
  sample_number = metadata_met$sample_number,
  stringsAsFactors = FALSE
)

expr_mapping <- data.frame(
  GSM_expr = as.character(metadata_expr$geo_accession),
  Expression_ID = metadata_expr$Expression_ID,
  sample_number = metadata_expr$sample_number,
  stringsAsFactors = FALSE
)


# 8. MATCH SAMPLES BY NUMERIC SUFFIX --------------------------

GEO_sample_matching <- merge(
  met_mapping,
  expr_mapping,
  by = "sample_number",
  all = FALSE,
  sort = TRUE
)

GEO_sample_matching <- GEO_sample_matching[
  order(GEO_sample_matching$sample_number),
]

rownames(GEO_sample_matching) <- NULL

GEO_sample_matching$matching_method <- "Numeric suffix"
GEO_sample_matching$matching_status <- "Assumed match"


# 9. IDENTIFY UNMATCHED SAMPLES -------------------------------

unmatched_expression <- expr_mapping[
  !expr_mapping$sample_number %in%
    GEO_sample_matching$sample_number,
  ,
  drop = FALSE
]

unmatched_methylation <- met_mapping[
  !met_mapping$sample_number %in%
    GEO_sample_matching$sample_number,
  ,
  drop = FALSE
]

rownames(unmatched_expression) <- NULL
rownames(unmatched_methylation) <- NULL


# 10. QUALITY CONTROL -----------------------------------------

n_met <- nrow(met_mapping)
n_expr <- nrow(expr_mapping)
n_matched <- nrow(GEO_sample_matching)
n_unmatched_expr <- nrow(unmatched_expression)
n_unmatched_met <- nrow(unmatched_methylation)

matching_QC <- data.frame(
  dataset = c(
    "GSE144487",
    "GSE65904",
    "Matched pairs",
    "Unmatched expression",
    "Unmatched methylation"
  ),
  number_of_samples = c(
    n_met,
    n_expr,
    n_matched,
    n_unmatched_expr,
    n_unmatched_met
  ),
  stringsAsFactors = FALSE
)

cat("\n============================================\n")
cat("GEO SAMPLE MATCHING — QUALITY CONTROL\n")
cat("============================================\n")

print(matching_QC)

cat("\nFirst 10 matched pairs:\n")
print(head(GEO_sample_matching, 10))

cat("\nLast 10 matched pairs:\n")
print(tail(GEO_sample_matching, 10))

cat("\nUnmatched expression samples:\n")
print(unmatched_expression)


# 10.1. CLINICAL QUALITY CONTROL & HARMONIZATION ------------

clinical_qc_check <- GEO_sample_matching

if ("Stage:ch1" %in% colnames(metadata_met)) {
  clinical_qc_check$met_stage_raw <- metadata_met$`Stage:ch1`[match(clinical_qc_check$GSM_met, metadata_met$geo_accession)]
}

if ("tumor stage:ch1" %in% colnames(metadata_expr)) {
  clinical_qc_check$expr_stage_raw <- metadata_expr$`tumor stage:ch1`[match(clinical_qc_check$GSM_expr, metadata_expr$geo_accession)]
}

# Harmonize textual staging categories for robust validation
if (all(c("met_stage_raw", "expr_stage_raw") %in% colnames(clinical_qc_check))) {
  
  clinical_qc_check$met_stage_clean <- dplyr::case_when(
    grepl("Regional", clinical_qc_check$met_stage_raw, ignore.case = TRUE) ~ "Regional",
    grepl("Primary",  clinical_qc_check$met_stage_raw, ignore.case = TRUE) ~ "Primary",
    grepl("Distant",  clinical_qc_check$met_stage_raw, ignore.case = TRUE) ~ "Distant",
    grepl("Local",    clinical_qc_check$met_stage_raw, ignore.case = TRUE) ~ "Local",
    TRUE ~ as.character(clinical_qc_check$met_stage_raw)
  )
  
  clinical_qc_check$expr_stage_clean <- dplyr::case_when(
    grepl("Regional", clinical_qc_check$expr_stage_raw, ignore.case = TRUE) ~ "Regional",
    grepl("Primary",  clinical_qc_check$expr_stage_raw, ignore.case = TRUE) ~ "Primary",
    grepl("Distant",  clinical_qc_check$expr_stage_raw, ignore.case = TRUE) ~ "Distant",
    grepl("Local",    clinical_qc_check$expr_stage_raw, ignore.case = TRUE) ~ "Local",
    TRUE ~ as.character(clinical_qc_check$expr_stage_raw)
  )
  
  harmonized_matches <- sum(clinical_qc_check$met_stage_clean == clinical_qc_check$expr_stage_clean, na.rm = TRUE)
  total_checked <- sum(!is.na(clinical_qc_check$met_stage_clean) & !is.na(clinical_qc_check$expr_stage_clean))
  
  cat(sprintf("\n[QC] Harmonized Stage matches: %d out of %d checked pairs\n", harmonized_matches, total_checked))
} else {
  cat("\n[QC] Stage columns not available for harmonization.\n")
}

cat("Note: Age column is present in expression metadata ('age:ch1'), but missing in methylation metadata. Age correlation skipped.\n")


# 11. CHECK EXPECTED COHORT STRUCTURE -------------------------

if (n_met != 196) {
  stop("Expected 196 methylation samples, found: ", n_met)
}

if (n_expr != 214) {
  stop("Expected 214 expression samples, found: ", n_expr)
}

if (n_matched != 196) {
  stop("Expected 196 matched pairs, found: ", n_matched)
}

if (n_unmatched_expr != 18) {
  stop(
    "Expected 18 unmatched expression samples, found: ",
    n_unmatched_expr
  )
}

if (n_unmatched_met != 0) {
  stop(
    "Unexpected unmatched methylation samples: ",
    n_unmatched_met
  )
}

if (!identical(
  GEO_sample_matching$sample_number,
  as.integer(1:196)
)) {
  stop("Matched sample numbers are not exactly 1:196.")
}

if (!identical(
  unmatched_expression$sample_number,
  as.integer(197:214)
)) {
  stop("Unmatched expression samples are not exactly 197:214.")
}


# 12. SAVE RESULTS --------------------------------------------

write.csv(
  GEO_sample_matching,
  file.path(output_dir, "GEO_sample_matching.csv"),
  row.names = FALSE
)

saveRDS(
  GEO_sample_matching,
  file.path(output_dir, "GEO_sample_matching.rds")
)

write.csv(
  unmatched_expression,
  file.path(
    output_dir,
    "GEO_unmatched_expression_samples.csv"
  ),
  row.names = FALSE
)

write.csv(
  unmatched_methylation,
  file.path(
    output_dir,
    "GEO_unmatched_methylation_samples.csv"
  ),
  row.names = FALSE
)

write.csv(
  matching_QC,
  file.path(output_dir, "GEO_sample_matching_QC.csv"),
  row.names = FALSE
)


# 12.1. FILTER AND SAVE HIGH-CONFIDENCE MATCHED SAMPLES -------

GEO_sample_matching_filtered <- clinical_qc_check %>%
  dplyr::filter(
    !is.na(met_stage_clean) & 
      !is.na(expr_stage_clean) & 
      (met_stage_clean == expr_stage_clean)
  )

cat(sprintf("\nHigh-confidence matched pairs saved after stage filtering: %d out of %d\n", 
            nrow(GEO_sample_matching_filtered), nrow(GEO_sample_matching)))

write.csv(
  GEO_sample_matching_filtered,
  file.path(output_dir, "GEO_sample_matching_high_confidence.csv"),
  row.names = FALSE
)

saveRDS(
  GEO_sample_matching_filtered,
  file.path(output_dir, "GEO_sample_matching_high_confidence.rds")
)


# 13. FINAL REPORT --------------------------------------------

cat("\n============================================\n")
cat("GEO SAMPLE MATCHING COMPLETED\n")
cat("============================================\n")

cat("Methylation samples:       ", n_met, "\n")
cat("Expression samples:        ", n_expr, "\n")
cat("Assumed matched pairs:     ", n_matched, "\n")
cat("High-confidence pairs:     ", nrow(GEO_sample_matching_filtered), "\n")
cat("Unmatched expression:      ", n_unmatched_expr, "\n")
cat("Unmatched methylation:     ", n_unmatched_met, "\n")

cat("\nMatching rule:\n")
cat("TumN <-> LUMN (numeric suffix assumption, validated by clinical stages)\n")

cat("\nResults saved in:\n")
cat(output_dir, "\n")

cat("\nGEO sample matching script finished successfully.\n")
