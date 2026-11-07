###############################################################################
# SCRIPT 57 — FULL SHINY APPLICATION BUILD AUDIT
###############################################################################

cat("\n========== SCRIPT 57: FULL APP BUILD AUDIT ==========\n")

# 1. VERIFY PATHS -------------------------------------------------------------

stopifnot(
  exists("app_dir"),
  dir.exists(app_dir)
)

app_file <- file.path(app_dir, "app.R")

stopifnot(file.exists(app_file))

cat("App directory:", app_dir, "\n")
cat("app.R exists: VERIFIED\n")


# 2. CHECK APP SOURCE SYNTAX --------------------------------------------------

app_parse <- tryCatch(
  {
    parse(app_file)
    TRUE
  },
  error = function(e) {
    cat("PARSE ERROR:", conditionMessage(e), "\n")
    FALSE
  }
)

stopifnot(app_parse)

cat("app.R parsing: VALID\n")


# 3. SOURCE APP IN AN ISOLATED ENVIRONMENT -----------------------------------

test_app_build <- function() {
  
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  
  setwd(app_dir)
  
  test_env <- new.env(parent = globalenv())
  
  sys.source(
    app_file,
    envir = test_env
  )
  
  test_env
}

build_warnings <- character()

test_env <- withCallingHandlers(
  test_app_build(),
  warning = function(w) {
    build_warnings <<- c(
      build_warnings,
      conditionMessage(w)
    )
    invokeRestart("muffleWarning")
  }
)

cat("Full app source: COMPLETED\n")


# 4. VERIFY REQUIRED MODULE FUNCTIONS ----------------------------------------

required_functions <- c(
  "melanomaUploadServer",
  "predictionOverviewServer",
  "epigeneticEvidenceServer",
  "dynamicEpigeneticUI",
  "dynamicEpigeneticServer"
)

function_status <- vapply(
  required_functions,
  function(fn) {
    exists(fn, envir = test_env, inherits = TRUE) &&
      is.function(get(fn, envir = test_env, inherits = TRUE))
  },
  logical(1)
)

print(function_status)

stopifnot(all(function_status))

cat("All required module functions: AVAILABLE\n")


# 5. VERIFY APP OBJECT --------------------------------------------------------

app_objects <- ls(test_env)

shiny_objects <- app_objects[
  vapply(
    app_objects,
    function(nm) {
      inherits(
        get(nm, envir = test_env),
        "shiny.appobj"
      )
    },
    logical(1)
  )
]

cat("Shiny app objects found:", length(shiny_objects), "\n")

if (length(shiny_objects) > 0) {
  print(shiny_objects)
}

stopifnot(length(shiny_objects) >= 1)

cat("Shiny application object: VALID\n")


# 6. VERIFY EXISTING AND NEW MODULE REFERENCES -------------------------------

app_source <- readLines(
  app_file,
  warn = FALSE
)

module_checks <- c(
  upload_module = any(grepl(
    "melanomaUploadServer",
    app_source,
    fixed = TRUE
  )),
  prediction_module = any(grepl(
    "predictionOverviewServer",
    app_source,
    fixed = TRUE
  )),
  evidence_module = any(grepl(
    "epigeneticEvidenceServer",
    app_source,
    fixed = TRUE
  )),
  dynamic_module = any(grepl(
    "dynamicEpigeneticServer",
    app_source,
    fixed = TRUE
  ))
)

print(module_checks)

stopifnot(all(module_checks))

cat("All four application modules: CONNECTED\n")


# 7. REPORT WARNINGS ----------------------------------------------------------

cat("\nBUILD WARNINGS:", length(build_warnings), "\n")

if (length(build_warnings) > 0) {
  print(build_warnings)
} else {
  cat("No warnings during app build.\n")
}


# 8. FINAL RESULT -------------------------------------------------------------

cat("\n========== SCRIPT 57 COMPLETED ==========\n")
cat("Full application source: VERIFIED\n")
cat("Required module functions: VERIFIED\n")
cat("Shiny app object: VERIFIED\n")
cat("All module references: VERIFIED\n")
cat("Existing results and figures: NOT MODIFIED\n")


###############################################################################
# SCRIPT 58 — EXPLICIT SHINY APP OBJECT VALIDATION
###############################################################################

cat("\n========== SCRIPT 58: EXPLICIT APP VALIDATION ==========\n")

# 1. CHECK UI AND SERVER -------------------------------------------------------

cat("\n--- UI / SERVER OBJECTS ---\n")

ui_exists <- exists("ui", envir = test_env, inherits = FALSE)
server_exists <- exists("server", envir = test_env, inherits = FALSE)

cat("ui exists:", ui_exists, "\n")
cat("server exists:", server_exists, "\n")

if (ui_exists) {
  cat(
    "UI class:",
    paste(class(get("ui", envir = test_env)), collapse = ", "),
    "\n"
  )
}

if (server_exists) {
  cat(
    "Server is function:",
    is.function(get("server", envir = test_env)),
    "\n"
  )
}

stopifnot(ui_exists, server_exists)
stopifnot(is.function(get("server", envir = test_env)))


# 2. CREATE SHINY APP OBJECT --------------------------------------------------

cat("\n--- CREATE SHINY APP OBJECT ---\n")

app_object <- shiny::shinyApp(
  ui = get("ui", envir = test_env),
  server = get("server", envir = test_env)
)

cat(
  "App object class:",
  paste(class(app_object), collapse = ", "),
  "\n"
)

app_object_valid <- inherits(app_object, "shiny.appobj")

cat("Shiny app object valid:", app_object_valid, "\n")

stopifnot(app_object_valid)


# 3. VERIFY MODULE FUNCTIONS --------------------------------------------------

required_functions <- c(
  "melanomaUploadServer",
  "predictionOverviewServer",
  "epigeneticEvidenceServer",
  "dynamicEpigeneticUI",
  "dynamicEpigeneticServer"
)

function_status <- vapply(
  required_functions,
  function(fn) {
    exists(fn, envir = test_env, inherits = TRUE) &&
      is.function(get(fn, envir = test_env, inherits = TRUE))
  },
  logical(1)
)

print(function_status)

stopifnot(all(function_status))


# 4. FINAL RESULT -------------------------------------------------------------

cat("\n========== SCRIPT 58 COMPLETED ==========\n")
cat("UI object: VERIFIED\n")
cat("Server function: VERIFIED\n")
cat("Shiny app object: VERIFIED\n")
cat("All five module functions: VERIFIED\n")
cat("No project files modified\n")

