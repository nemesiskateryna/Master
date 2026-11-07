
###############################################################################
# ExplainMyMelanoma AI
# SCRIPT 12 — CSV INPUT VALIDATION MODULE
# Creates a reusable CSV validation function and tests it.
# READ-ONLY for existing thesis data and model files.
###############################################################################

project_root <- "/Users/katie/master/Master"

app_r_dir <- file.path(
  project_root,
  "app/ExplainMyMelanoma/R"
)

output_file <- file.path(
  app_r_dir,
  "02_csv_input_validation.R"
)

if (!dir.exists(app_r_dir)) {
  stop("Application R directory does not exist.")
}

if (file.exists(output_file)) {
  stop("Output file already exists. No file was overwritten.")
}

###############################################################################
# 1. DEFINE CSV VALIDATION FUNCTION
###############################################################################

validation_code <- '
validate_methylation_csv <- function(file, mlp) {

  # Check uploaded file
  if (length(file) != 1 || !file.exists(file)) {
    stop("CSV file does not exist.")
  }

  # Read CSV
  input_data <- tryCatch(
    read.csv(
      file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) {
      stop("Unable to read CSV: ", conditionMessage(e))
    }
  )

  # Check empty file
  if (nrow(input_data) == 0) {
    stop("CSV file contains no data rows.")
  }

  # Check required columns
  required_columns <- c("CpG", "Beta")

  if (!all(required_columns %in% colnames(input_data))) {
    stop(
      "CSV must contain two required columns: CpG and Beta."
    )
  }

  # Keep required columns
  input_data <- input_data[, required_columns, drop = FALSE]

  # Normalize CpG identifiers
  input_data$CpG <- trimws(as.character(input_data$CpG))

  if (any(is.na(input_data$CpG) | input_data$CpG == "")) {
    stop("CpG identifiers contain missing or empty values.")
  }

  # Check duplicated CpGs
  if (anyDuplicated(input_data$CpG)) {
    stop("CSV contains duplicated CpG identifiers.")
  }

  # Convert Beta values to numeric
  beta_values <- suppressWarnings(
    as.numeric(as.character(input_data$Beta))
  )

  if (any(!is.finite(beta_values))) {
    stop("Beta values contain missing or non-numeric values.")
  }

  if (any(beta_values < 0 | beta_values > 1)) {
    stop("Beta values must be within the interval [0, 1].")
  }

  input_data$Beta <- beta_values

  # Check required model features
  required_cpgs <- mlp$selected_cpgs

  missing_cpgs <- setdiff(
    required_cpgs,
    input_data$CpG
  )

  if (length(missing_cpgs) > 0) {
    stop(
      paste0(
        "Missing ",
        length(missing_cpgs),
        " required CpGs. First missing IDs: ",
        paste(head(missing_cpgs, 5), collapse = ", ")
      )
    )
  }

  extra_cpgs <- setdiff(
    input_data$CpG,
    required_cpgs
  )

  # Run the existing MLP inference function
  prediction <- predict_mlp_profile(
    input_data,
    mlp
  )

  # Return structured result
  list(
    status = "success",
    prediction = prediction,
    input_summary = list(
      n_rows_uploaded = nrow(input_data),
      n_required_cpgs = length(required_cpgs),
      n_extra_cpgs = length(extra_cpgs),
      n_cpgs_used = prediction$n_cpgs_used,
      beta_min = min(input_data$Beta),
      beta_max = max(input_data$Beta)
    ),
    input_data = input_data
  )
}
'

writeLines(validation_code, output_file)

cat("Validation module created:\n")
cat(output_file, "\n")

###############################################################################
# 2. LOAD MODULE
###############################################################################

source(output_file)

cat("\n========== FUNCTION AUDIT ==========\n")

cat(
  "Function available:",
  exists("validate_methylation_csv"),
  "\n"
)

###############################################################################
# 3. LOAD MODEL AND TEST DATA
###############################################################################

mlp_file <- file.path(
  project_root,
  "results/models/TCGA_SKCM_MLP_final_model.rds"
)

tcga_file <- file.path(
  project_root,
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune.rds"
)

split_file <- file.path(
  project_root,
  "results/models/TCGA_SKCM_ML_train_test_split.rds"
)

mlp <- readRDS(mlp_file)
tcga_object <- readRDS(tcga_file)
split <- readRDS(split_file)

beta <- tcga_object$beta_matrix

###############################################################################
# 4. CREATE A TEMPORARY VALID CSV
###############################################################################

test_id <- split$test_sample_ids[1]

test_profile <- as.numeric(
  beta[mlp$selected_cpgs, test_id, drop = TRUE]
)

valid_csv_data <- data.frame(
  CpG = mlp$selected_cpgs,
  Beta = test_profile
)

temp_csv <- tempfile(fileext = ".csv")

write.csv(
  valid_csv_data,
  temp_csv,
  row.names = FALSE
)

###############################################################################
# 5. TEST VALID CSV UPLOAD
###############################################################################

cat("\n========== VALID CSV TEST ==========\n")

valid_result <- validate_methylation_csv(
  temp_csv,
  mlp
)

cat("Status:", valid_result$status, "\n")

cat(
  "Uploaded rows:",
  valid_result$input_summary$n_rows_uploaded,
  "\n"
)

cat(
  "Required CpGs:",
  valid_result$input_summary$n_required_cpgs,
  "\n"
)

cat(
  "CpGs used by MLP:",
  valid_result$input_summary$n_cpgs_used,
  "\n"
)

cat(
  "Predicted class:",
  valid_result$prediction$predicted_class,
  "\n"
)

cat(
  "Immune-high probability:",
  valid_result$prediction$probability_immune_high,
  "\n"
)

cat(
  "Immune-low probability:",
  valid_result$prediction$probability_immune_low,
  "\n"
)

cat(
  "Beta range:",
  valid_result$input_summary$beta_min,
  "-",
  valid_result$input_summary$beta_max,
  "\n"
)

###############################################################################
# 6. TEST MISSING CpG IN CSV
###############################################################################

cat("\n========== MISSING CpG CSV TEST ==========\n")

missing_csv <- tempfile(fileext = ".csv")

write.csv(
  valid_csv_data[-1, ],
  missing_csv,
  row.names = FALSE
)

missing_test <- tryCatch(
  validate_methylation_csv(missing_csv, mlp),
  error = function(e) conditionMessage(e)
)

print(missing_test)

###############################################################################
# 7. TEST INVALID BETA IN CSV
###############################################################################

cat("\n========== INVALID BETA CSV TEST ==========\n")

invalid_csv_data <- valid_csv_data
invalid_csv_data$Beta[1] <- 1.5

invalid_csv <- tempfile(fileext = ".csv")

write.csv(
  invalid_csv_data,
  invalid_csv,
  row.names = FALSE
)

invalid_test <- tryCatch(
  validate_methylation_csv(invalid_csv, mlp),
  error = function(e) conditionMessage(e)
)

print(invalid_test)

###############################################################################
# 8. CLEAN TEMPORARY FILES
###############################################################################

unlink(c(temp_csv, missing_csv, invalid_csv))

cat("\nSCRIPT 12 — COMPLETED\n")
###############################################################################

