###############################################################################
# SCRIPT 56 — INTEGRATE DYNAMIC VISUALIZATION MODULE INTO APP.R
###############################################################################

cat("\n========== SCRIPT 56: APP INTEGRATION ==========\n")

# 1. PATHS --------------------------------------------------------------------

app_file <- file.path(app_dir, "app.R")

backup_file <- paste0(
  app_file,
  ".backup_before_dynamic_visualizations"
)

module_file <- file.path(
  app_dir,
  "R",
  "06_dynamic_epigenetic_visualizations.R"
)

stopifnot(
  file.exists(app_file),
  file.exists(module_file),
  file.exists(file.path(app_dir, "R", "03_upload_module.R"))
)

# 2. READ APP -----------------------------------------------------------------

app_text <- paste(
  readLines(app_file, warn = FALSE),
  collapse = "\n"
)

# 3. DEFINE EXACT INSERTION BLOCKS --------------------------------------------

source_anchor <- paste(
  'source(file.path(app_dir, "R", "05_epigenetic_evidence_module.R"))',
  sep = ""
)

source_new <- paste(
  source_anchor,
  'source(file.path(app_dir, "R", "06_dynamic_epigenetic_visualizations.R"))',
  sep = "\n"
)

ui_anchor <- paste(
  '      # EPIGENETIC EVIDENCE',
  '      shiny::div(',
  '        class = "section-card",',
  '        shiny::h3(class = "section-heading", "Epigenetic evidence"),',
  '        shiny::div(class = "accent-line"),',
  '        epigeneticEvidenceUI("epigenetic")',
  '      ),',
  sep = "\n"
)

ui_new <- paste(
  '      # DYNAMIC EPIGENETIC VISUALIZATIONS',
  '      shiny::div(',
  '        class = "section-card",',
  '        shiny::h3(class = "section-heading", "Dynamic methylation visualizations"),',
  '        shiny::div(class = "accent-line"),',
  '        dynamicEpigeneticUI("dynamic")',
  '      ),',
  '',
  ui_anchor,
  sep = "\n"
)

server_anchor <- paste(
  '      })',
  '    })',
  '  )',
  '  ',
  '  # Epigenetic evidence module',
  '  epigeneticEvidenceServer(',
  '    id = "epigenetic",',
  '    tables_dir = tables_dir,',
  '    figures_dir = figures_dir',
  '  )',
  sep = "\n"
)

server_new <- paste(
  '      })',
  '    })',
  '  )',
  '',
  '  # Dynamic methylation visualizations',
  '  dynamicEpigeneticServer(',
  '    id = "dynamic",',
  '    profile = shiny::reactive({',
  '      current_result <- upload_result$result()',
  '      if (is.null(current_result)) {',
  '        return(NULL)',
  '      }',
  '      current_result$input_data',
  '    }),',
  '    mlp = mlp',
  '  )',
  '',
  '  # Epigenetic evidence module',
  '  epigeneticEvidenceServer(',
  '    id = "epigenetic",',
  '    tables_dir = tables_dir,',
  '    figures_dir = figures_dir',
  '  )',
  sep = "\n"
)

# 4. CHECK EXACT ANCHORS -------------------------------------------------------

count_fixed <- function(pattern, x) {
  hit <- gregexpr(pattern, x, fixed = TRUE)[[1]]
  if (identical(hit, -1L)) return(0L)
  length(hit)
}

source_count <- count_fixed(source_anchor, app_text)
ui_count <- count_fixed(ui_anchor, app_text)
server_count <- count_fixed(server_anchor, app_text)

cat("Source anchor:", source_count, "\n")
cat("UI anchor:", ui_count, "\n")
cat("Server anchor:", server_count, "\n")

if (
  source_count != 1 ||
  ui_count != 1 ||
  server_count != 1
) {
  stop(
    "Expected anchors were not found exactly once. ",
    "No app.R changes were written."
  )
}

if (
  grepl("dynamicEpigeneticUI", app_text, fixed = TRUE) ||
  grepl("dynamicEpigeneticServer", app_text, fixed = TRUE)
) {
  stop(
    "Dynamic module integration appears to exist already. ",
    "No changes were written."
  )
}

# 5. CREATE BACKUP -------------------------------------------------------------

if (file.exists(backup_file)) {
  stop("Backup already exists. No changes were written.")
}

backup_ok <- file.copy(
  app_file,
  backup_file,
  overwrite = FALSE
)

if (!backup_ok) {
  stop("Could not create app.R backup. No changes were written.")
}

cat("Backup created:\n", backup_file, "\n")

# 6. APPLY THREE INSERTIONS ----------------------------------------------------

patched_text <- sub(
  source_anchor,
  source_new,
  app_text,
  fixed = TRUE
)

patched_text <- sub(
  ui_anchor,
  ui_new,
  patched_text,
  fixed = TRUE
)

patched_text <- sub(
  server_anchor,
  server_new,
  patched_text,
  fixed = TRUE
)

# 7. VERIFY PATCHED TEXT -------------------------------------------------------

stopifnot(
  count_fixed(
    'source(file.path(app_dir, "R", "06_dynamic_epigenetic_visualizations.R"))',
    patched_text
  ) == 1,
  count_fixed('dynamicEpigeneticUI("dynamic")', patched_text) == 1,
  count_fixed("dynamicEpigeneticServer(", patched_text) == 1,
  count_fixed("current_result$input_data", patched_text) == 1
)

# 8. PARSE BEFORE WRITING ------------------------------------------------------

parse(text = patched_text)

# 9. WRITE APP ----------------------------------------------------------------

writeLines(
  strsplit(patched_text, "\n", fixed = TRUE)[[1]],
  con = app_file,
  useBytes = TRUE
)

# 10. VERIFY ACTUAL FILE -------------------------------------------------------

written_text <- paste(
  readLines(app_file, warn = FALSE),
  collapse = "\n"
)

stopifnot(
  count_fixed(
    'source(file.path(app_dir, "R", "06_dynamic_epigenetic_visualizations.R"))',
    written_text
  ) == 1,
  count_fixed('dynamicEpigeneticUI("dynamic")', written_text) == 1,
  count_fixed("dynamicEpigeneticServer(", written_text) == 1,
  count_fixed("current_result$input_data", written_text) == 1
)

parse(file = app_file)

# 11. SOURCE ALL APP MODULES IN ISOLATED ENVIRONMENT ---------------------------

integration_env <- new.env(parent = globalenv())

sys.source(
  file.path(app_dir, "R", "01_mlp_inference.R"),
  envir = integration_env
)

sys.source(
  file.path(app_dir, "R", "02_csv_input_validation.R"),
  envir = integration_env
)

sys.source(
  file.path(app_dir, "R", "03_upload_module.R"),
  envir = integration_env
)

sys.source(
  file.path(app_dir, "R", "04_prediction_overview_module.R"),
  envir = integration_env
)

sys.source(
  file.path(app_dir, "R", "05_epigenetic_evidence_module.R"),
  envir = integration_env
)

sys.source(
  module_file,
  envir = integration_env
)

stopifnot(
  exists("melanomaUploadServer", envir = integration_env),
  exists("predictionOverviewUI", envir = integration_env),
  exists("epigeneticEvidenceUI", envir = integration_env),
  exists("dynamicEpigeneticUI", envir = integration_env),
  exists("dynamicEpigeneticServer", envir = integration_env)
)

# 12. FINAL SUMMARY ------------------------------------------------------------

cat("\n========== SCRIPT 56 COMPLETED ==========\n")
cat("Module 06 source in app.R: VERIFIED\n")
cat("Dynamic visualization UI: VERIFIED\n")
cat("Dynamic visualization server: VERIFIED\n")
cat("Uploaded profile connection: VERIFIED\n")
cat("Existing Prediction Overview retained: VERIFIED\n")
cat("Existing Epigenetic Evidence retained: VERIFIED\n")
cat("app.R syntax: VALID\n")
cat("All module functions: AVAILABLE\n")
cat("Backup: PRESERVED\n")