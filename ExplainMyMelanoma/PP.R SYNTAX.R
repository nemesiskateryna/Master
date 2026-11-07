
###############################################################################
# SCRIPT 22 — APP.R SYNTAX AND DEPENDENCY AUDIT
###############################################################################

cat("\n========== SCRIPT 22: APP AUDIT ==========\n")

project_root <- "/Users/katie/master/Master"

app_dir <- file.path(
  project_root,
  "app/ExplainMyMelanoma"
)

app_file <- file.path(app_dir, "app.R")

model_file <- file.path(
  project_root,
  "results/models/TCGA_SKCM_MLP_final_model.rds"
)

module_files <- file.path(
  app_dir,
  "R",
  c(
    "01_mlp_inference.R",
    "02_csv_input_validation.R",
    "03_upload_module.R"
  )
)

# 1. FILE EXISTENCE
cat("\n--- FILE CHECK ---\n")

cat("app.R exists:", file.exists(app_file), "\n")
cat("MLP model exists:", file.exists(model_file), "\n")

for (f in module_files) {
  cat(
    basename(f),
    ":",
    file.exists(f),
    "\n"
  )
}

# 2. APP.R SYNTAX
cat("\n--- APP.R SYNTAX ---\n")

syntax_ok <- tryCatch(
  {
    parsed_app <- parse(file = app_file)
    TRUE
  },
  error = function(e) {
    cat("SYNTAX ERROR:", conditionMessage(e), "\n")
    FALSE
  }
)

cat("Syntax valid:", syntax_ok, "\n")

# 3. REQUIRED PACKAGES
cat("\n--- PACKAGE CHECK ---\n")

required_packages <- c("shiny", "bslib")

for (pkg in required_packages) {
  cat(
    pkg,
    ":",
    requireNamespace(pkg, quietly = TRUE),
    "\n"
  )
}

# 4. MODULE SOURCE CHECK
cat("\n--- MODULE CHECK ---\n")

module_env <- new.env(parent = globalenv())

modules_ok <- FALSE

if (all(file.exists(module_files))) {
  
  modules_ok <- tryCatch(
    {
      for (f in module_files) {
        sys.source(f, envir = module_env)
      }
      
      expected_functions <- c(
        "melanomaUploadUI",
        "melanomaUploadServer",
        "validate_methylation_csv",
        "predict_mlp_profile"
      )
      
      for (fn in expected_functions) {
        cat(
          fn,
          "available:",
          exists(fn, envir = module_env, inherits = FALSE),
          "\n"
        )
      }
      
      all(vapply(
        expected_functions,
        exists,
        logical(1),
        envir = module_env,
        inherits = FALSE
      ))
    },
    error = function(e) {
      cat("MODULE ERROR:", conditionMessage(e), "\n")
      FALSE
    }
  )
}

cat("All expected module functions available:", modules_ok, "\n")

# 5. MODEL STRUCTURE CHECK
cat("\n--- MLP MODEL CHECK ---\n")

model_ok <- FALSE

if (file.exists(model_file)) {
  
  model_ok <- tryCatch(
    {
      mlp_audit <- readRDS(model_file)
      
      cat(
        "Selected CpGs:",
        length(mlp_audit$selected_cpgs),
        "\n"
      )
      
      cat(
        "PCA components:",
        ncol(mlp_audit$pca$rotation),
        "\n"
      )
      
      length(mlp_audit$selected_cpgs) == 5000 &&
        ncol(mlp_audit$pca$rotation) == 30
    },
    error = function(e) {
      cat("MODEL ERROR:", conditionMessage(e), "\n")
      FALSE
    }
  )
}

cat("MLP model structure valid:", model_ok, "\n")

# 6. FINAL AUDIT
cat("\n========== FINAL AUDIT ==========\n")

audit_passed <- (
  file.exists(app_file) &&
    syntax_ok &&
    all(file.exists(module_files)) &&
    all(vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)) &&
    modules_ok &&
    model_ok
)

cat("APP AUDIT PASSED:", audit_passed, "\n")

if (audit_passed) {
  cat("\nReady for the next step: controlled Shiny launch.\n")
} else {
  cat("\nAudit failed. Do not launch the application yet.\n")
}

###############################################################################

