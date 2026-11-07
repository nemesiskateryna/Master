# TCGA-SKCM: TUMOR SAMPLE TYPE QC





# 1. Extract sample type codes from full TCGA barcodes

sample_barcodes <- colnames(integrated$beta_matrix)

sample_type <- substr(
  sample_barcodes,
  14,
  15
)

# 2. Create sample type labels

sample_type_label <- dplyr::case_when(
  sample_type == "01" ~ "Primary Tumor",
  sample_type == "06" ~ "Metastatic Tumor",
  TRUE ~ "Other / Unknown"
)

# 3. Distribution

cat("\n--- SAMPLE TYPE DISTRIBUTION ---\n")

print(
  data.frame(
    Sample_type_code = sample_type,
    Sample_type = sample_type_label
  ) |>
    dplyr::count(
      Sample_type_code,
      Sample_type,
      name = "N"
    )
)

# 4. Check whether all samples are primary or metastatic

cat("\n--- SAMPLE TYPE VALIDATION ---\n")

cat(
  "All samples are primary or metastatic:",
  all(sample_type %in% c("01", "06")),
  "\n"
)

cat(
  "Number of primary tumors:",
  sum(sample_type == "01"),
  "\n"
)

cat(
  "Number of metastatic tumors:",
  sum(sample_type == "06"),
  "\n"
)

# 5. Check unique patients per sample type

patient_ids <- substr(sample_barcodes, 1, 12)

sample_type_patient_table <- data.frame(
  Patient = patient_ids,
  Sample_type = sample_type_label
) |>
  dplyr::distinct() |>
  dplyr::count(Sample_type, name = "Unique_patients")

cat("\n--- UNIQUE PATIENTS BY SAMPLE TYPE ---\n")

print(sample_type_patient_table)

# 6. Save sample type annotation

sample_type_annotation <- data.frame(
  barcode = sample_barcodes,
  patient = patient_ids,
  sample_type_code = sample_type,
  sample_type = sample_type_label,
  stringsAsFactors = FALSE
)

saveRDS(
  sample_type_annotation,
  "data/TCGA/processed/TCGA_SKCM_sample_type_annotation.rds"
)

write.csv(
  sample_type_annotation,
  "results/tables/TCGA_SKCM_sample_type_annotation.csv",
  row.names = FALSE
)

cat("\nSample type annotation saved:",
    file.exists(
      "data/TCGA/processed/TCGA_SKCM_sample_type_annotation.rds"
    ),
    "\n"
)
