###############################################################################
# SCRIPT 31 — LAUNCH AND TEST PREDICTION OVERVIEW
###############################################################################

cat("\n========== SCRIPT 31: APPLICATION LAUNCH ==========\n")

project_root <- "/Users/katie/master/Master"

app_dir <- file.path(
  project_root,
  "app/ExplainMyMelanoma"
)

required_files <- c(
  file.path(app_dir, "app.R"),
  file.path(app_dir, "R/01_mlp_inference.R"),
  file.path(app_dir, "R/02_csv_input_validation.R"),
  file.path(app_dir, "R/03_upload_module.R"),
  file.path(app_dir, "R/04_prediction_overview_module.R"),
  file.path(
    project_root,
    "results/models/TCGA_SKCM_MLP_final_model.rds"
  ),
  file.path(
    app_dir,
    "data/reference/ExplainMyMelanoma_test_profile.csv"
  )
)

cat("\n--- FILE CHECK ---\n")

file_check <- file.exists(required_files)

print(data.frame(
  File = basename(required_files),
  Exists = file_check
))

if (!all(file_check)) {
  stop("One or more required files are missing.")
}

# Parse all application files
cat("\n--- SYNTAX CHECK ---\n")

invisible(lapply(
  required_files[1:5],
  parse
))

cat("All application R files parsed successfully.\n")

# Launch application
cat("\n--- STARTING SHINY APPLICATION ---\n")

shiny::runApp(
  appDir = app_dir,
  launch.browser = TRUE
)

###############################################################################