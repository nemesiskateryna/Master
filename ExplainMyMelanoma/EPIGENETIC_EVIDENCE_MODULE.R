###############################################################################
# SCRIPT 36 — INTEGRATE EPIGENETIC EVIDENCE MODULE
#
# Changes app.R only:
#   1. Sources module 05.
#   2. Adds Epigenetic Evidence UI section.
#   3. Adds epigeneticEvidenceServer() call.
#
# Creates a backup before modifying app.R.
# Does NOT modify analysis results, tables, figures, models or other modules.
###############################################################################

cat("\n========== SCRIPT 36: MODULE INTEGRATION ==========\n")

# 1. CHECK ENVIRONMENT --------------------------------------------------------

required_objects <- c(
  "app_dir",
  "tables_dir",
  "figures_dir"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1))
]

if (length(missing_objects) > 0) {
  stop(
    "Missing required R objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

app_file <- file.path(app_dir, "app.R")

module_file <- file.path(
  app_dir,
  "R",
  "05_epigenetic_evidence_module.R"
)

if (!file.exists(app_file)) {
  stop("app.R not found.")
}

if (!file.exists(module_file)) {
  stop("Module 05 not found.")
}

if (!file.exists(tables_dir) || !file.exists(figures_dir)) {
  stop("Tables or figures directory does not exist.")
}

# 2. CHECK MODULE FUNCTIONS ---------------------------------------------------

module_text <- readLines(module_file, warn = FALSE)

if (!any(grepl(
  "^epigeneticEvidenceUI\\s*<-\\s*function",
  module_text
))) {
  stop("epigeneticEvidenceUI() was not found.")
}

if (!any(grepl(
  "^epigeneticEvidenceServer\\s*<-\\s*function",
  module_text
))) {
  stop("epigeneticEvidenceServer() was not found.")
}

parse(module_file)

# 3. READ app.R ---------------------------------------------------------------

app_lines <- readLines(app_file, warn = FALSE)

if (any(grepl(
  'source\\(file.path\\(app_dir, "R", "05_epigenetic_evidence_module.R"\\)\\)',
  app_lines
))) {
  stop("Module 05 is already sourced. No changes made.")
}

if (any(grepl(
  'epigeneticEvidenceUI\\("epigenetic"\\)',
  app_lines,
  fixed = TRUE
))) {
  stop("Epigenetic Evidence UI already exists. No changes made.")
}

if (any(grepl(
  'epigeneticEvidenceServer\\(',
  app_lines,
  fixed = TRUE
))) {
  stop("Epigenetic Evidence server call already exists. No changes made.")
}

# 4. VERIFY UNIQUE INSERTION ANCHORS ------------------------------------------

source_anchor <- 'source(file.path(app_dir, "R", "04_prediction_overview_module.R"))'
ui_anchor <- "# SCIENTIFIC INFORMATION"
app_marker <- "# Shiny application"

source_idx <- which(trimws(app_lines) == source_anchor)
ui_idx <- which(trimws(app_lines) == ui_anchor)
app_marker_idx <- which(trimws(app_lines) == app_marker)

if (length(source_idx) != 1) {
  stop("Expected exactly one source anchor for module 04.")
}

if (length(ui_idx) != 1) {
  stop("Expected exactly one Scientific Information UI anchor.")
}

if (length(app_marker_idx) != 1) {
  stop("Expected exactly one Shiny application marker.")
}

server_start_idx <- which(grepl(
  "^server\\s*<-\\s*function",
  app_lines
))

if (length(server_start_idx) != 1) {
  stop("Could not uniquely identify server function.")
}

# Find the final closing brace between server() and shinyApp().
server_end_candidates <- which(
  seq_along(app_lines) > server_start_idx &
    seq_along(app_lines) < app_marker_idx &
    trimws(app_lines) == "}"
)

if (length(server_end_candidates) == 0) {
  stop("Could not identify the end of server().")
}

server_end_idx <- max(server_end_candidates)

cat("Insertion anchors verified.\n")
cat("Source insertion after line:", source_idx, "\n")
cat("UI insertion before line:", ui_idx, "\n")
cat("Server insertion before line:", server_end_idx, "\n")

# 5. BUILD CANDIDATE app.R IN MEMORY ------------------------------------------

new_lines <- app_lines

# Add module source.
new_lines <- append(
  new_lines,
  'source(file.path(app_dir, "R", "05_epigenetic_evidence_module.R"))',
  after = source_idx
)

# Recalculate UI anchor after source insertion.
ui_idx_new <- which(trimws(new_lines) == ui_anchor)

# Add Epigenetic Evidence UI section.
ui_block <- c(
  "    # EPIGENETIC EVIDENCE",
  "    shiny::div(",
  '      class = "section-card",',
  '      shiny::h3(class = "section-heading", "Epigenetic evidence"),',
  '      shiny::div(class = "accent-line"),',
  '      epigeneticEvidenceUI("epigenetic")',
  "    ),",
  ""
)

new_lines <- append(
  new_lines,
  ui_block,
  after = ui_idx_new - 1
)

# Recalculate server boundary.
app_marker_idx_new <- which(trimws(new_lines) == app_marker)

server_start_idx_new <- which(grepl(
  "^server\\s*<-\\s*function",
  new_lines
))

server_end_candidates_new <- which(
  seq_along(new_lines) > server_start_idx_new &
    seq_along(new_lines) < app_marker_idx_new &
    trimws(new_lines) == "}"
)

server_end_idx_new <- max(server_end_candidates_new)

# Add server module call before server() closing brace.
server_block <- c(
  "",
  "  # Epigenetic evidence module",
  "  epigeneticEvidenceServer(",
  '    id = "epigenetic",',
  "    tables_dir = tables_dir,",
  "    figures_dir = figures_dir",
  "  )",
  ""
)

new_lines <- append(
  new_lines,
  server_block,
  after = server_end_idx_new - 1
)

# 6. VALIDATE CANDIDATE BEFORE WRITING ----------------------------------------

candidate_ok <- tryCatch({
  parse(text = new_lines)
  TRUE
}, error = function(e) {
  cat("Candidate app.R parse error:", conditionMessage(e), "\n")
  FALSE
})

if (!candidate_ok) {
  stop("Candidate app.R is invalid. Original file remains unchanged.")
}

cat("\nCandidate app.R syntax valid: TRUE\n")

# 7. CREATE BACKUP ------------------------------------------------------------

backup_file <- paste0(
  app_file,
  ".backup_before_epigenetic_evidence"
)

if (file.exists(backup_file)) {
  stop(
    "Backup already exists. No file was modified:\n",
    backup_file
  )
}

backup_ok <- file.copy(
  from = app_file,
  to = backup_file,
  overwrite = FALSE
)

if (!backup_ok || !file.exists(backup_file)) {
  stop("Could not create app.R backup. Original file was not modified.")
}

cat("Backup created:\n", backup_file, "\n")

# 8. WRITE UPDATED app.R -------------------------------------------------------

writeLines(new_lines, app_file)

# 9. POST-WRITE VALIDATION -----------------------------------------------------

final_syntax_ok <- tryCatch({
  parse(app_file)
  TRUE
}, error = function(e) {
  cat("Final app.R parse error:", conditionMessage(e), "\n")
  FALSE
})

final_lines <- readLines(app_file, warn = FALSE)

source_ok <- any(grepl(
  'source\\(file.path\\(app_dir, "R", "05_epigenetic_evidence_module.R"\\)\\)',
  final_lines
))

ui_ok <- any(grepl(
  'epigeneticEvidenceUI\\("epigenetic"\\)',
  final_lines,
  fixed = TRUE
))

server_ok <- any(grepl(
  'epigeneticEvidenceServer\\(',
  final_lines,
  fixed = TRUE
))

cat("\n========== SCRIPT 36 FINAL STATUS ==========\n")
cat("Backup exists:", file.exists(backup_file), "\n")
cat("app.R syntax valid:", final_syntax_ok, "\n")
cat("Module 05 source added:", source_ok, "\n")
cat("Epigenetic Evidence UI added:", ui_ok, "\n")
cat("Epigenetic Evidence server added:", server_ok, "\n")

if (!all(
  final_syntax_ok,
  source_ok,
  ui_ok,
  server_ok
)) {
  stop(
    "Post-write validation failed. Do not launch the app. ",
    "The original backup is available."
  )
}

cat("\nSCRIPT 36 PASSED.\n")
cat("Only app.R was modified.\n")
cat("No existing analysis results were modified.\n")