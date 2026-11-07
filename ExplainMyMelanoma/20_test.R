###############################################################################
# SCRIPT 20 — TEST OBSERVER INITIALIZATION
###############################################################################


###############################################################################
# SCRIPT 24B — CREATE PERMANENT VALID TEST CSV
###############################################################################

cat("\n========== SCRIPT 24B: SAVE TEST CSV ==========\n")

project_root <- "/Users/katie/master/Master"

# 1. Define paths
model_file <- file.path(
  project_root,
  "results/models/TCGA_SKCM_MLP_final_model.rds"
)

integrated_file <- file.path(
  project_root,
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune.rds"
)

split_file <- file.path(
  project_root,
  "results/models/TCGA_SKCM_ML_train_test_split.rds"
)

output_dir <- file.path(
  project_root,
  "app/ExplainMyMelanoma/data/reference"
)

output_file <- file.path(
  output_dir,
  "ExplainMyMelanoma_test_profile.csv"
)

# 2. Safety checks
stopifnot(
  file.exists(model_file),
  file.exists(integrated_file),
  file.exists(split_file),
  dir.exists(output_dir)
)

if (file.exists(output_file)) {
  stop("Test CSV already exists. No file was overwritten.")
}

# 3. Load original objects
mlp <- readRDS(model_file)

tcga_object <- readRDS(integrated_file)

split <- readRDS(split_file)

beta <- tcga_object$beta_matrix

# 4. Select first holdout sample
test_id <- split$test_sample_ids[1]

cat("Selected holdout sample:", test_id, "\n")

# 5. Check feature availability
required_cpgs <- mlp$selected_cpgs

if (!all(required_cpgs %in% rownames(beta))) {
  stop("Some required CpGs are missing from the Beta matrix.")
}

if (!test_id %in% colnames(beta)) {
  stop("Holdout sample is missing from the Beta matrix.")
}

# 6. Reconstruct profile in exact model order
test_profile <- as.numeric(
  beta[required_cpgs, test_id, drop = TRUE]
)

names(test_profile) <- required_cpgs

stopifnot(
  length(test_profile) == 5000,
  !is.null(names(test_profile)),
  identical(names(test_profile), required_cpgs),
  all(is.finite(test_profile)),
  all(test_profile >= 0 & test_profile <= 1)
)

# 7. Create CSV
valid_csv_data <- data.frame(
  CpG = names(test_profile),
  Beta = unname(test_profile),
  stringsAsFactors = FALSE
)

write.csv(
  valid_csv_data,
  file = output_file,
  row.names = FALSE,
  quote = FALSE
)

# 8. Read back and verify
saved_csv <- read.csv(
  output_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n--- CSV VERIFICATION ---\n")

cat("File exists:", file.exists(output_file), "\n")
cat("File size:", file.info(output_file)$size, "bytes\n")
cat("Rows:", nrow(saved_csv), "\n")
cat("Columns:", paste(names(saved_csv), collapse = ", "), "\n")
cat("CpG IDs unique:", !anyDuplicated(saved_csv$CpG), "\n")
cat("CpG order correct:", identical(saved_csv$CpG, required_cpgs), "\n")
cat("Beta values finite:", all(is.finite(saved_csv$Beta)), "\n")
cat("Beta minimum:", min(saved_csv$Beta), "\n")
cat("Beta maximum:", max(saved_csv$Beta), "\n")

# 9. Final checks
csv_ok <- (
  nrow(saved_csv) == 5000 &&
    identical(names(saved_csv), c("CpG", "Beta")) &&
    identical(saved_csv$CpG, required_cpgs) &&
    !anyDuplicated(saved_csv$CpG) &&
    all(is.finite(saved_csv$Beta)) &&
    all(saved_csv$Beta >= 0 & saved_csv$Beta <= 1)
)

cat("\nCSV CREATION PASSED:", csv_ok, "\n")

if (csv_ok) {
  cat("\nPERMANENT TEST FILE:\n")
  cat(output_file, "\n")
}

###############################################################################
