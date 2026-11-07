
###############################################################################
# ExplainMyMelanoma AI
# SCRIPT 01 — PROJECT SETUP AND MODEL AUDIT
###############################################################################

cat("\n========================================\n")
cat("EXPLAINMYMELANOMA AI — STEP 1\n")
cat("========================================\n\n")

# 1. Define existing thesis project
project_root <- "/Users/katie/master/Master"

if (!dir.exists(project_root)) {
  stop("Thesis project directory was not found. Check project_root.")
}

setwd(project_root)

cat("Thesis project directory:\n")
cat(getwd(), "\n\n")


# 2. Create separate application directory
app_dir <- file.path(
  project_root,
  "app",
  "ExplainMyMelanoma"
)

app_subdirs <- c(
  "R",
  "www",
  "data/reference",
  "results/reports"
)

if (!dir.exists(app_dir)) {
  dir.create(app_dir, recursive = TRUE)
  cat("Created application directory.\n")
} else {
  cat("Application directory already exists. Preserving its contents.\n")
}

for (folder in app_subdirs) {
  
  folder_path <- file.path(app_dir, folder)
  
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    cat("Created:", folder, "\n")
  } else {
    cat("Already exists:", folder, "\n")
  }
}


# 3. Check essential model files
model_dir <- file.path(project_root, "results", "models")

required_files <- c(
  "TCGA_SKCM_MLP_final_model.rds",
  "TCGA_SKCM_final_ElasticNet_model.rds",
  "TCGA_SKCM_ElasticNet_CpG_signature.rds",
  "TCGA_SKCM_ML_train_test_split.rds"
)

model_paths <- file.path(model_dir, required_files)

model_audit <- data.frame(
  File = required_files,
  Exists = file.exists(model_paths),
  Size_bytes = ifelse(
    file.exists(model_paths),
    file.info(model_paths)$size,
    NA_real_
  )
)

cat("\n========== MODEL FILE AUDIT ==========\n")
print(model_audit, row.names = FALSE)


# 4. Check saved model object structure
cat("\n========== MODEL OBJECT STRUCTURE ==========\n")

for (i in seq_along(model_paths)) {
  
  if (file.exists(model_paths[i])) {
    
    cat("\nFILE:", required_files[i], "\n")
    
    obj <- readRDS(model_paths[i])
    
    cat("Class:", paste(class(obj), collapse = ", "), "\n")
    
    if (is.list(obj)) {
      cat("List elements:\n")
      print(names(obj))
    }
    
    rm(obj)
  }
}


# 5. Check required R packages
required_packages <- c(
  "shiny",
  "bslib",
  "ggplot2",
  "dplyr",
  "glmnet",
  "nnet"
)

package_audit <- data.frame(
  Package = required_packages,
  Installed = vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
)

cat("\n========== R PACKAGE AUDIT ==========\n")
print(package_audit, row.names = FALSE)


# 6. Display application structure
cat("\n========== APPLICATION STRUCTURE ==========\n")

for (folder in c("", app_subdirs)) {
  
  path <- file.path(app_dir, folder)
  
  cat(
    ifelse(folder == "", "ExplainMyMelanoma/", paste0("  ", folder, "/")),
    ifelse(dir.exists(path), " [OK]", " [MISSING]"),
    "\n"
  )
}


# 7. Final checks
cat("\n========== FINAL CHECK ==========\n")

cat(
  "Application directory exists:",
  dir.exists(app_dir), "\n"
)

cat(
  "Required model files found:",
  sum(model_audit$Exists),
  "/",
  nrow(model_audit), "\n"
)

cat(
  "Required packages installed:",
  sum(package_audit$Installed),
  "/",
  nrow(package_audit), "\n"
)

cat("\nSCRIPT 01 — SETUP AND AUDIT COMPLETED\n")
###############################################################################

###############################################################################
# ExplainMyMelanoma AI
# SCRIPT 02 — INSTALL SHINY
###############################################################################

cat("\n========================================\n")
cat("EXPLAINMYMELANOMA AI — STEP 2\n")
cat("========================================\n\n")

# Install Shiny if missing
if (!requireNamespace("shiny", quietly = TRUE)) {
  
  install.packages("shiny")
  
} else {
  
  cat("Shiny is already installed.\n")
  
}

# Verify required packages
required_packages <- c(
  "shiny",
  "bslib",
  "ggplot2",
  "dplyr",
  "glmnet",
  "nnet"
)

package_audit <- data.frame(
  Package = required_packages,
  Installed = vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
)

cat("\n========== PACKAGE AUDIT ==========\n")
print(package_audit, row.names = FALSE)

cat(
  "\nAll required packages installed:",
  all(package_audit$Installed),
  "\n"
)

cat("\nSCRIPT 02 — COMPLETED\n")

###############################################################################


