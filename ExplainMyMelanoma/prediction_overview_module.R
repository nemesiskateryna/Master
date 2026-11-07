###############################################################################
# SCRIPT 30B — SAFE PREDICTION OVERVIEW INTEGRATION
###############################################################################

cat("\n========== SCRIPT 30B ==========\n")

project_root <- "/Users/katie/master/Master"

app_dir <- file.path(
  project_root,
  "app/ExplainMyMelanoma"
)

app_file <- file.path(app_dir, "app.R")
upload_file <- file.path(app_dir, "R/03_upload_module.R")
prediction_module <- file.path(
  app_dir,
  "R/04_prediction_overview_module.R"
)

# 1. Confirm all files and previous backups
stopifnot(
  file.exists(app_file),
  file.exists(upload_file),
  file.exists(prediction_module),
  file.exists(paste0(app_file, ".backup_before_prediction_overview")),
  file.exists(paste0(upload_file, ".backup_before_prediction_overview"))
)

# 2. Read files
app_lines <- readLines(app_file, warn = FALSE)
upload_lines <- readLines(upload_file, warn = FALSE)

# 3. Add module source
source_anchor <- 'source(file.path(app_dir, "R", "03_upload_module.R"))'

if (sum(app_lines == source_anchor) != 1) {
  stop("Module 03 source anchor is not unique.")
}

if (any(grepl("04_prediction_overview_module.R", app_lines, fixed = TRUE))) {
  stop("Module 04 source already exists. No changes made.")
}

app_lines <- append(
  app_lines,
  'source(file.path(app_dir, "R", "04_prediction_overview_module.R"))',
  after = which(app_lines == source_anchor)
)

# 4. Insert Prediction overview UI after upload section
upload_anchor <- '      melanomaUploadUI("upload")'

if (sum(app_lines == upload_anchor) != 1) {
  stop("Upload UI anchor is not unique.")
}

upload_index <- which(app_lines == upload_anchor)

if (app_lines[upload_index + 1] != "    ),") {
  stop("Unexpected upload section boundary.")
}

prediction_ui <- c(
  "",
  "    # PREDICTION OVERVIEW",
  "    shiny::div(",
  '      class = "section-card",',
  '      shiny::h3(class = "section-heading", "Prediction overview"),',
  '      shiny::div(class = "accent-line"),',
  '      predictionOverviewUI("prediction")',
  "    ),"
)

app_lines <- append(
  app_lines,
  prediction_ui,
  after = upload_index + 1
)

# 5. Replace server block using the Shiny application marker
server_start <- grep(
  "^server <- function\\(input, output, session\\) \\{$",
  app_lines
)

shiny_marker <- grep(
  "^# Shiny application$",
  app_lines
)

if (length(server_start) != 1 || length(shiny_marker) != 1 ||
    shiny_marker <= server_start) {
  stop("Could not safely identify server section.")
}

new_server <- c(
  "server <- function(input, output, session) {",
  "",
  "  upload_result <- melanomaUploadServer(",
  '    id = "upload",',
  "    mlp = mlp",
  "  )",
  "",
  "  predictionOverviewServer(",
  '    id = "prediction",',
  "    prediction = shiny::reactive({",
  "      current_result <- upload_result$result()",
  "      if (is.null(current_result)) {",
  "        return(NULL)",
  "      }",
  "      current_result$prediction",
  "    })",
  "  )",
  "}",
  ""
)

app_lines <- c(
  app_lines[seq_len(server_start - 1)],
  new_server,
  app_lines[shiny_marker:length(app_lines)]
)

# 6. Remove old prediction UI line (indentation-independent)
old_prediction_pattern <- paste0(
  "^[[:space:]]*shiny::uiOutput\\(",
  "ns\\(['\"]prediction['\"]\\)\\),?[[:space:]]*$"
)

old_prediction_matches <- grep(
  old_prediction_pattern,
  upload_lines,
  perl = TRUE
)

if (length(old_prediction_matches) != 1) {
  stop(
    "Expected exactly one old prediction UI line; found ",
    length(old_prediction_matches),
    ". No files written."
  )
}

upload_lines <- upload_lines[-old_prediction_matches]

# 7. Parse in memory before writing
parse(text = app_lines)
parse(text = upload_lines)

# 8. Write updated files
writeLines(app_lines, app_file)
writeLines(upload_lines, upload_file)

# 9. Verify saved files
parse(file = app_file)
parse(file = upload_file)
parse(file = prediction_module)

updated_app <- readLines(app_file, warn = FALSE)
updated_upload <- readLines(upload_file, warn = FALSE)

checks <- c(
  app_backup_preserved = file.exists(
    paste0(app_file, ".backup_before_prediction_overview")
  ),
  upload_backup_preserved = file.exists(
    paste0(upload_file, ".backup_before_prediction_overview")
  ),
  module_04_loaded = any(grepl(
    "04_prediction_overview_module.R",
    updated_app,
    fixed = TRUE
  )),
  prediction_ui_added = any(grepl(
    'predictionOverviewUI("prediction")',
    updated_app,
    fixed = TRUE
  )),
  prediction_server_added = any(grepl(
    "predictionOverviewServer",
    updated_app,
    fixed = TRUE
  )),
  reactive_prediction_connected = any(grepl(
    "current_result$prediction",
    updated_app,
    fixed = TRUE
  )),
  old_prediction_ui_removed = !any(grepl(
    old_prediction_pattern,
    updated_upload,
    perl = TRUE
  )),
  GEO_GSE144487_preserved = any(grepl(
    "GSE144487",
    updated_app,
    fixed = TRUE
  )),
  GEO_GSE65904_preserved = any(grepl(
    "GSE65904",
    updated_app,
    fixed = TRUE
  )),
  GSE91061_absent = !any(grepl(
    "GSE91061",
    updated_app,
    fixed = TRUE
  ))
)

cat("\n--- INTEGRATION CHECKS ---\n")
print(checks)

cat("\nAll checks passed:", all(checks), "\n")

if (!all(checks)) {
  stop("One or more integration checks failed.")
}

cat("\nAPP.R SIZE:", file.info(app_file)$size, "bytes\n")
cat("UPLOAD MODULE SIZE:", file.info(upload_file)$size, "bytes\n")

cat("\nSCRIPT 30B PASSED\n")

###############################################################################