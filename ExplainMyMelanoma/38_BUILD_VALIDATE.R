###############################################################################
# SCRIPT 38 — BUILD AND VALIDATE SHINY APP
#
# Does not launch the browser or start the Shiny server.
# Does not modify any files.
###############################################################################

cat("\n========== SCRIPT 38: APP BUILD TEST ==========\n")

# 1. PATHS --------------------------------------------------------------------

app_file <- file.path(app_dir, "app.R")

if (!file.exists(app_file)) {
  stop("app.R not found: ", app_file)
}

if (!file.exists(file.path(
  app_dir, "R", "05_epigenetic_evidence_module.R"
))) {
  stop("Epigenetic evidence module not found.")
}

# 2. CHECK REQUIRED DATA DIRECTORIES ------------------------------------------

cat("\n--- DATA DIRECTORIES ---\n")

cat("tables_dir available:",
    exists("tables_dir", inherits = TRUE), "\n")

cat("figures_dir available:",
    exists("figures_dir", inherits = TRUE), "\n")

if (!exists("tables_dir", inherits = TRUE) ||
    !exists("figures_dir", inherits = TRUE)) {
  stop("tables_dir or figures_dir is unavailable.")
}

cat("Tables directory exists:", dir.exists(tables_dir), "\n")
cat("Figures directory exists:", dir.exists(figures_dir), "\n")

# 3. SOURCE APP IN AN ISOLATED ENVIRONMENT ------------------------------------

cat("\n--- SOURCING APP.R ---\n")

old_wd <- getwd()
app_env <- new.env(parent = globalenv())

source_result <- tryCatch({
  
  setwd(app_dir)
  
  source(
    app_file,
    local = app_env,
    echo = FALSE
  )
  
}, error = function(e) {
  
  cat("\nAPP BUILD ERROR:\n")
  cat(conditionMessage(e), "\n")
  
  NULL
  
}, finally = {
  
  setwd(old_wd)
  
})

if (is.null(source_result)) {
  stop("App build failed. No files were modified.")
}

# 4. CHECK MODULE FUNCTIONS ---------------------------------------------------

cat("\n--- MODULE FUNCTION CHECK ---\n")

required_functions <- c(
  "epigeneticEvidenceUI",
  "epigeneticEvidenceServer",
  "predictionOverviewUI",
  "predictionOverviewServer",
  "melanomaUploadUI",
  "melanomaUploadServer"
)

for (fn in required_functions) {
  
  available <- exists(
    fn,
    envir = app_env,
    inherits = TRUE
  )
  
  cat(
    if (available) "[OK] " else "[MISSING] ",
    fn,
    "\n"
  )
}

functions_ok <- all(vapply(
  required_functions,
  exists,
  logical(1),
  envir = app_env,
  inherits = TRUE
))

# 5. CHECK UI AND SERVER OBJECTS -----------------------------------------------

cat("\n--- APP OBJECTS ---\n")

ui_exists <- exists(
  "ui",
  envir = app_env,
  inherits = FALSE
)

server_exists <- exists(
  "server",
  envir = app_env,
  inherits = FALSE
)

cat("UI object exists:", ui_exists, "\n")
cat("Server function exists:", server_exists, "\n")

if (ui_exists) {
  cat("UI class:",
      paste(class(app_env$ui), collapse = ", "),
      "\n")
}

if (server_exists) {
  cat("Server class:",
      paste(class(app_env$server), collapse = ", "),
      "\n")
}

# 6. CHECK RETURNED SHINY APPLICATION -----------------------------------------

app_object <- source_result$value

app_object_ok <- inherits(
  app_object,
  "shiny.appobj"
)

cat("\n--- SHINY APPLICATION OBJECT ---\n")
cat("Returned object class:",
    paste(class(app_object), collapse = ", "),
    "\n")

cat("Valid shiny.appobj:", app_object_ok, "\n")

# 7. FINAL STATUS -------------------------------------------------------------

cat("\n========== SCRIPT 38 FINAL STATUS ==========\n")

cat("App source completed:", TRUE, "\n")
cat("Required module functions available:", functions_ok, "\n")
cat("UI object exists:", ui_exists, "\n")
cat("Server function exists:", server_exists, "\n")
cat("Shiny application object valid:", app_object_ok, "\n")

if (!all(
  functions_ok,
  ui_exists,
  server_exists,
  app_object_ok
)) {
  stop("SCRIPT 38 FAILED. Do not launch the app yet.")
}

cat("\nSCRIPT 38 PASSED.\n")
cat("App was built but NOT launched.\n")
cat("No files were modified.\n")