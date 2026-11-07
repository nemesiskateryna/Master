###############################################################################
# SCRIPT 69 — SHINY APP VALIDATION WITH CORRECT WORKING DIRECTORY
###############################################################################

cat("\n========== SCRIPT 69: APP VALIDATION ==========\n")

# 1. SAVE CURRENT WORKING DIRECTORY --------------------------------------------

original_wd_69 <- getwd()

cat("Original working directory:", original_wd_69, "\n")
cat("Application directory:", app_dir, "\n")

# 2. SWITCH TO APPLICATION DIRECTORY ------------------------------------------

setwd(app_dir)

cat("Active working directory:", getwd(), "\n")

# 3. PARSE APP.R ----------------------------------------------------------------

app_file_69 <- file.path(app_dir, "app.R")

expressions_69 <- parse(app_file_69)

cat("\napp.R parsed successfully:", TRUE, "\n")
cat("Expressions:", length(expressions_69), "\n")

# 4. CREATE ISOLATED ENVIRONMENT ------------------------------------------------

test_env_69 <- new.env(parent = globalenv())

# 5. EVALUATE APP.R -------------------------------------------------------------

app_object_69 <- NULL

evaluation_error_69 <- NULL

tryCatch(
  {
    for (i in seq_along(expressions_69)) {
      
      result_69 <- eval(
        expressions_69[[i]],
        envir = test_env_69
      )
      
      if (i == length(expressions_69)) {
        app_object_69 <- result_69
      }
    }
  },
  error = function(e) {
    evaluation_error_69 <<- conditionMessage(e)
  }
)

# 6. RESTORE ORIGINAL WORKING DIRECTORY ----------------------------------------

setwd(original_wd_69)

cat("\nWorking directory restored:", getwd(), "\n")

# 7. REPORT EVALUATION ----------------------------------------------------------

cat("\n--- APP EVALUATION ---\n")

if (is.null(evaluation_error_69)) {
  cat("app.R evaluated successfully: TRUE\n")
} else {
  cat("app.R evaluation failed:\n")
  cat(evaluation_error_69, "\n")
}

# 8. VALIDATE SHINY OBJECT ------------------------------------------------------

is_shiny_app_69 <- inherits(
  app_object_69,
  "shiny.appobj"
)

cat("\n--- SHINY APP OBJECT ---\n")
cat("Class:", paste(class(app_object_69), collapse = ", "), "\n")
cat("Valid shiny.appobj:", is_shiny_app_69, "\n")

# 9. VALIDATE APP COMPONENTS ----------------------------------------------------

ui_exists_69 <- exists(
  "ui",
  envir = test_env_69,
  inherits = FALSE
)

server_exists_69 <- exists(
  "server",
  envir = test_env_69,
  inherits = FALSE
)

cat("\n--- APP COMPONENTS ---\n")
cat("UI exists:", ui_exists_69, "\n")
cat("Server exists:", server_exists_69, "\n")

# 10. VALIDATE MODULE FUNCTIONS -------------------------------------------------

module_names_69 <- c(
  "melanomaUploadUI",
  "melanomaUploadServer",
  "predictionOverviewUI",
  "predictionOverviewServer",
  "epigeneticEvidenceUI",
  "epigeneticEvidenceServer",
  "dynamicEpigeneticUI",
  "dynamicEpigeneticServer"
)

module_status_69 <- vapply(
  module_names_69,
  function(nm) {
    exists(
      nm,
      envir = globalenv(),
      mode = "function",
      inherits = FALSE
    )
  },
  logical(1)
)

cat("\n--- MODULE FUNCTIONS ---\n")
print(module_status_69)

# 11. VALIDATE FINAL RESULT -----------------------------------------------------

checks_69 <- c(
  evaluation = is.null(evaluation_error_69),
  shiny_app = is_shiny_app_69,
  ui = ui_exists_69,
  server = server_exists_69,
  modules = all(module_status_69)
)

cat("\n========== FINAL VALIDATION ==========\n")
print(checks_69)

if (all(checks_69)) {
  
  cat("\nSCRIPT 69 COMPLETED SUCCESSFULLY\n")
  cat("Shiny app object: VALID\n")
  cat("All modules: AVAILABLE\n")
  
} else {
  
  cat("\nSCRIPT 69 NOT FULLY PASSED\n")
  cat("Please send the complete output for diagnosis.\n")
}

cat("\nNo project files modified.\n")



