
###############################################################################
# SCRIPT 27A — CORRECT GEO DATASET REFERENCES IN APP
###############################################################################

cat("\n========== SCRIPT 27A: GEO REFERENCES ==========\n")

project_root <- "/Users/katie/master/Master"

app_dir <- file.path(
  project_root,
  "app/ExplainMyMelanoma"
)

app_file <- file.path(app_dir, "app.R")

backup_file <- file.path(
  app_dir,
  "app.R.backup_before_GEO_correction"
)

temp_file <- file.path(
  app_dir,
  "app.R.GEO_correction_tmp"
)

# 1. Safety checks
if (!file.exists(app_file)) {
  stop("app.R not found.")
}

if (file.exists(backup_file)) {
  stop("Backup already exists. No files were modified.")
}

if (file.exists(temp_file)) {
  stop("Temporary file already exists. No files were modified.")
}

# 2. Read current app
app_lines <- readLines(
  app_file,
  warn = FALSE
)

old_pattern <- 'GEO dataset GSE91061'

old_hits <- grep(
  old_pattern,
  app_lines,
  fixed = TRUE
)

cat("Incorrect GEO references:", length(old_hits), "\n")

if (length(old_hits) != 1) {
  stop("Expected exactly one incorrect GEO reference. No changes made.")
}

# 3. Replace with correct dataset descriptions
new_lines <- app_lines[-old_hits]

new_geo_lines <- c(
  '      shiny::p(shiny::strong("GEO DNA methylation validation: "), "GSE144487 — external validation of the CpG signature."),',
  '      shiny::p(shiny::strong("GEO gene-expression analysis: "), "GSE65904 — immune phenotype and TIS-17 analysis."),'
)

new_lines <- append(
  new_lines,
  new_geo_lines,
  after = old_hits - 1
)

# 4. Write temporary version
writeLines(
  new_lines,
  con = temp_file,
  useBytes = TRUE
)

# 5. Validate syntax before replacing original
syntax_ok <- tryCatch(
  {
    parse(file = temp_file)
    TRUE
  },
  error = function(e) {
    cat("Syntax error:", conditionMessage(e), "\n")
    FALSE
  }
)

if (!syntax_ok) {
  unlink(temp_file)
  stop("Temporary app.R has invalid syntax. Original preserved.")
}

# 6. Create backup
backup_ok <- file.copy(
  from = app_file,
  to = backup_file,
  overwrite = FALSE
)

if (!backup_ok) {
  unlink(temp_file)
  stop("Backup creation failed. Original preserved.")
}

# 7. Replace app.R
replace_ok <- file.copy(
  from = temp_file,
  to = app_file,
  overwrite = TRUE
)

if (!replace_ok) {
  stop("Replacement failed. Original backup is available.")
}

unlink(temp_file)

# 8. Final audit
corrected_lines <- readLines(
  app_file,
  warn = FALSE
)

cat("\n--- FINAL GEO AUDIT ---\n")

cat(
  "GSE91061:",
  sum(grepl("GSE91061", corrected_lines, fixed = TRUE)),
  "\n"
)

cat(
  "GSE144487:",
  sum(grepl("GSE144487", corrected_lines, fixed = TRUE)),
  "\n"
)

cat(
  "GSE65904:",
  sum(grepl("GSE65904", corrected_lines, fixed = TRUE)),
  "\n"
)

cat("Backup exists:", file.exists(backup_file), "\n")

final_syntax <- tryCatch(
  {
    parse(file = app_file)
    TRUE
  },
  error = function(e) FALSE
)

cat("Final syntax valid:", final_syntax, "\n")

cat(
  "SCRIPT 27A PASSED:",
  final_syntax &&
    file.exists(backup_file) &&
    !any(grepl("GSE91061", corrected_lines, fixed = TRUE)) &&
    any(grepl("GSE144487", corrected_lines, fixed = TRUE)) &&
    any(grepl("GSE65904", corrected_lines, fixed = TRUE)),
  "\n"
)

###############################################################################


###############################################################################
# SCRIPT 27B — PREDICTION OBJECT STRUCTURE AUDIT
###############################################################################

cat("\n========== SCRIPT 27B: PREDICTION AUDIT ==========\n")

project_root <- "/Users/katie/master/Master"

app_dir <- file.path(
  project_root,
  "app/ExplainMyMelanoma"
)

model_file <- file.path(
  project_root,
  "results/models/TCGA_SKCM_MLP_final_model.rds"
)

csv_file <- file.path(
  app_dir,
  "data/reference/ExplainMyMelanoma_test_profile.csv"
)

# 1. Check inputs
stopifnot(
  file.exists(model_file),
  file.exists(csv_file)
)

# 2. Load model and validation function
mlp <- readRDS(model_file)

source(
  file.path(app_dir, "R/01_mlp_inference.R")
)

source(
  file.path(app_dir, "R/02_csv_input_validation.R")
)

# 3. Validate test CSV
validated <- validate_methylation_csv(
  csv_file,
  mlp
)

# 4. Inspect prediction structure
cat("\n--- PREDICTION STRUCTURE ---\n")

str(
  validated$prediction,
  max.level = 3
)

# 5. Inspect input summary
cat("\n--- INPUT SUMMARY STRUCTURE ---\n")

str(
  validated$input_summary,
  max.level = 2
)

# 6. Print available prediction values
cat("\n--- PREDICTION VALUES ---\n")

print(validated$prediction)

cat("\nSCRIPT 27B COMPLETED\n")

###############################################################################


