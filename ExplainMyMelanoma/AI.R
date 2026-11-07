
###############################################################################
# SCRIPT 92 — LAUNCH SHINY APPLICATION
# ExplainMyMelanoma AI
#
# Purpose:
#   1. Validate app files and saved objects.
#   2. Validate Top 50 CpG importance features.
#   3. Launch Shiny application.
#
# This script does not modify project files.
###############################################################################

cat("\n========== SCRIPT 92: SHINY LAUNCH ==========\n")

###############################################################################
# 1. PROJECT PATHS
###############################################################################

project_root <- "/Users/katie/master/Master"

app_dir <- file.path(
  project_root,
  "app",
  "ExplainMyMelanoma"
)

app_file <- file.path(app_dir, "app.R")

module_file <- file.path(
  app_dir,
  "R",
  "06_dynamic_epigenetic_visualizations.R"
)

model_file <- file.path(
  project_root,
  "results",
  "models",
  "TCGA_SKCM_MLP_final_model.rds"
)

importance_file <- file.path(
  project_root,
  "results",
  "models",
  "TCGA_SKCM_MLP_CpG_permutation_importance_5fold.rds"
)

###############################################################################
# 2. FILE VALIDATION
###############################################################################

cat("\n--- 1. FILE VALIDATION ---\n")

required_files <- c(
  app_file,
  module_file,
  model_file,
  importance_file
)

file_status <- file.exists(required_files)

print(data.frame(
  File = basename(required_files),
  Exists = file_status
))

if (!all(file_status)) {
  stop("Required application or model files are missing.")
}

###############################################################################
# 3. SYNTAX VALIDATION
###############################################################################

cat("\n--- 2. SYNTAX VALIDATION ---\n")

parse(file = app_file)
parse(file = module_file)

cat("app.R syntax: PASS\n")
cat("Module syntax: PASS\n")

###############################################################################
# 4. MODEL AND IMPORTANCE VALIDATION
###############################################################################

cat("\n--- 3. MODEL VALIDATION ---\n")

mlp_launch_check <- readRDS(model_file)

importance_launch_check <- readRDS(importance_file)

if (length(mlp_launch_check$selected_cpgs) != 5000) {
  stop("MLP model does not contain exactly 5000 CpG features.")
}

importance_table_launch <- importance_launch_check$importance_summary

required_columns <- c(
  "CpG",
  "Mean_Delta_LogLoss",
  "Importance_Rank"
)

if (!all(required_columns %in% colnames(importance_table_launch))) {
  stop("Importance summary is missing required columns.")
}

top50_launch <- importance_table_launch[
  order(importance_table_launch$Importance_Rank),
  ,
  drop = FALSE
]

top50_launch <- head(top50_launch, 50)

if (nrow(top50_launch) != 50) {
  stop("Could not identify exactly 50 top-ranked CpGs.")
}

if (!all(as.character(top50_launch$CpG) %in%
         as.character(mlp_launch_check$selected_cpgs))) {
  stop("Top 50 CpGs do not match the MLP feature set.")
}

cat("MLP features:", length(mlp_launch_check$selected_cpgs), "\n")
cat("Top 50 CpGs: PASS\n")

###############################################################################
# 5. SHINY PACKAGE VALIDATION
###############################################################################

cat("\n--- 4. SHINY PACKAGE VALIDATION ---\n")

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Package 'shiny' is not installed.")
}

cat("Shiny version:",
    as.character(utils::packageVersion("shiny")),
    "\n")

###############################################################################
# 6. SET WORKING DIRECTORY
###############################################################################

cat("\n--- 5. APPLICATION DIRECTORY ---\n")

setwd(app_dir)

cat("Working directory:\n")
cat(getwd(), "\n")

###############################################################################
# 7. LAUNCH APPLICATION
###############################################################################

cat("\n--- 6. LAUNCHING SHINY ---\n")

cat(
  "The application will open in your browser.\n",
  "Keep this R session running while testing the application.\n"
)

shiny::runApp(
  appDir = app_dir,
  launch.browser = TRUE,
  display.mode = "normal"
)

###############################################################################
# SCRIPT 92 END
###############################################################################

