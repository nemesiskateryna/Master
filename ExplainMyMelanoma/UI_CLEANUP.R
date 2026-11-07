###############################################################################
# SCRIPT 60 — EVIDENCE UI CLEANUP + COLLAPSIBLE SECTION
###############################################################################

cat("\n========== SCRIPT 60: EVIDENCE UI CLEANUP ==========\n")

# 1. PATHS AND BACKUP ---------------------------------------------------------

evidence_file <- file.path(
  app_dir,
  "R",
  "05_epigenetic_evidence_module.R"
)

backup_file <- paste0(
  evidence_file,
  ".backup_before_evidence_cleanup"
)

stopifnot(file.exists(evidence_file))

if (file.exists(backup_file)) {
  stop("Backup already exists. No changes made.")
}

evidence_lines <- readLines(
  evidence_file,
  warn = FALSE
)

ui_start <- grep(
  "epigeneticEvidenceUI\\s*<-\\s*function",
  evidence_lines
)

server_start <- grep(
  "epigeneticEvidenceServer\\s*<-\\s*function",
  evidence_lines
)

stopifnot(
  length(ui_start) == 1,
  length(server_start) == 1,
  ui_start < server_start
)

# 2. DEFINE NEW COLLAPSIBLE UI ------------------------------------------------

new_ui <- c(
  "epigeneticEvidenceUI <- function(id) {",
  "",
  "  ns <- shiny::NS(id)",
  "",
  "  shiny::tags$details(",
  "    class = \"epigenetic-evidence-details\",",
  "",
  "    shiny::tags$summary(",
  "      \"Show epigenetic evidence\"",
  "    ),",
  "",
  "    shiny::div(",
  "      class = \"epigenetic-evidence-content\",",
  "",
  "      # GEO EXTERNAL VALIDATION",
  "      shiny::h4(\"GEO external CpG validation\"),",
  "",
  "      shiny::p(",
  "        \"GSE144487: external assessment of the TCGA-derived \",",
  "        \"DNA methylation signature.\"",
  "      ),",
  "",
  "      shiny::tableOutput(ns(\"geo_validation_table\")),",
  "",
  "      shiny::div(",
  "        class = \"scientific-disclaimer\",",
  "        shiny::strong(\"Interpretation\"),",
  "        shiny::p(",
  "          \"The CpG signature shows substantial overlap and \",",
  "          \"partial directional concordance in GEO. However, \",",
  "          \"the permutation results do not reach the conventional \",",
  "          \"0.05 significance threshold. This is evidence of \",",
  "          \"a biological association, not definitive global replication.\"",
  "        )",
  "      ),",
  "",
  "      shiny::hr(),",
  "",
  "      # GEO IMMUNE SIGNATURES",
  "      shiny::h4(\"GSE65904: immune signature analysis\"),",
  "",
  "      shiny::p(",
  "        \"Association of the TIS-17 gene signature with \",",
  "        \"auxiliary immune-related expression signatures.\"",
  "      ),",
  "",
  "      shiny::tableOutput(ns(\"tis_table\")),",
  "",
  "      shiny::div(",
  "        class = \"scientific-disclaimer\",",
  "        shiny::strong(\"Cohort matching limitation\"),",
  "        shiny::p(",
  "          \"The available matching audit labels all 74 records \",",
  "          \"as assumed matches based on numeric sample suffixes. \",",
  "          \"These matches should not be interpreted as independently \",",
  "          \"verified patient-level matches.\"",
  "        )",
  "      )",
  "    )",
  "  )",
  "}"
)

# 3. BUILD UPDATED FILE --------------------------------------------------------

new_lines <- c(
  evidence_lines[seq_len(ui_start - 1)],
  new_ui,
  "",
  evidence_lines[server_start:length(evidence_lines)]
)

# Preserve original server function exactly.
original_server <- evidence_lines[server_start:length(evidence_lines)]
updated_server <- new_lines[
  (length(new_lines) - length(original_server) + 1):length(new_lines)
]

stopifnot(identical(original_server, updated_server))

# 4. CREATE BACKUP ------------------------------------------------------------

backup_created <- file.copy(
  evidence_file,
  backup_file,
  overwrite = FALSE
)

stopifnot(backup_created)

cat("Backup created:", backup_file, "\n")

# 5. WRITE UPDATED UI ---------------------------------------------------------

writeLines(
  new_lines,
  evidence_file
)

cat("Updated evidence module written.\n")

# 6. VERIFY SYNTAX ------------------------------------------------------------

parse_result <- tryCatch(
  {
    parse(evidence_file)
    TRUE
  },
  error = function(e) {
    cat("PARSE ERROR:", conditionMessage(e), "\n")
    FALSE
  }
)

stopifnot(parse_result)

cat("Module syntax: VALID\n")

# 7. SOURCE MODULE IN ISOLATED ENVIRONMENT -----------------------------------

evidence_test_env <- new.env(parent = globalenv())

sys.source(
  evidence_file,
  envir = evidence_test_env
)

stopifnot(
  is.function(evidence_test_env$epigeneticEvidenceUI),
  is.function(evidence_test_env$epigeneticEvidenceServer)
)

cat("UI and server functions: AVAILABLE\n")

# 8. TEST GENERATED UI --------------------------------------------------------

ui_test <- evidence_test_env$epigeneticEvidenceUI("audit")

ui_html <- as.character(ui_test)

required_content <- c(
  "Show epigenetic evidence",
  "GEO external CpG validation",
  "GSE144487",
  "GSE65904",
  "TIS-17",
  "Cohort matching limitation",
  "geo_validation_table",
  "tis_table"
)

excluded_content <- c(
  "CpG–gene expression associations",
  "cg16890093",
  "cg07016276",
  "TAP1",
  "HLA-F",
  "CpG annotation overview",
  "Genomic region distribution",
  "CpG-island distribution"
)

required_status <- vapply(
  required_content,
  function(x) grepl(x, ui_html, fixed = TRUE),
  logical(1)
)

excluded_status <- vapply(
  excluded_content,
  function(x) !grepl(x, ui_html, fixed = TRUE),
  logical(1)
)

cat("\n--- REQUIRED UI CONTENT ---\n")
print(required_status)

cat("\n--- REMOVED UI CONTENT ---\n")
print(excluded_status)

stopifnot(
  all(required_status),
  all(excluded_status),
  grepl("<details", ui_html, fixed = TRUE),
  grepl("<summary", ui_html, fixed = TRUE)
)

cat("\nCollapsible details element: VERIFIED\n")

# 9. VERIFY SERVER FUNCTION WAS NOT CHANGED ----------------------------------

current_lines <- readLines(
  evidence_file,
  warn = FALSE
)

current_server_start <- grep(
  "epigeneticEvidenceServer\\s*<-\\s*function",
  current_lines
)

stopifnot(length(current_server_start) == 1)

server_unchanged <- identical(
  original_server,
  current_lines[
    current_server_start:length(current_lines)
  ]
)

stopifnot(server_unchanged)

cat("Original server function: PRESERVED\n")

# 10. FINAL -------------------------------------------------------------------

cat("\n========== SCRIPT 60 COMPLETED ==========\n")
cat("GEO validation summary: RETAINED\n")
cat("GSE65904 TIS-17 table and caveat: RETAINED\n")
cat("CpG-gene association UI: REMOVED\n")
cat("TAP1 / HLA-F figures UI: REMOVED\n")
cat("CpG annotation overview UI: REMOVED\n")
cat("Epigenetic Evidence collapsible: VERIFIED\n")
cat("Server function: UNCHANGED\n")
cat("Backup: PRESERVED\n")
cat("Existing results and figures: NOT MODIFIED\n")


###############################################################################
# SCRIPT 66 — CONNECT DYNAMIC VISUALIZATION UI AND SERVER
###############################################################################

app_file_66 <- file.path(app_dir, "app.R")

app_lines_66 <- readLines(
  app_file_66,
  warn = FALSE
)

# 1. DEFINE EXACT INSERTION ANCHORS -------------------------------------------

ui_anchor_66 <- paste(
  "    # EPIGENETIC EVIDENCE",
  "    shiny::div(",
  "      class = \"section-card\",",
  "      shiny::h3(class = \"section-heading\", \"Epigenetic evidence\"),",
  "      shiny::div(class = \"accent-line\"),",
  "      epigeneticEvidenceUI(\"epigenetic\")",
  "    ),",
  sep = "\n"
)

server_anchor_66 <- paste(
  "  # Epigenetic evidence module",
  "  epigeneticEvidenceServer(",
  "    id = \"epigenetic\",",
  "    tables_dir = tables_dir,",
  "    figures_dir = figures_dir",
  "  )",
  sep = "\n"
)

app_text_66 <- paste(app_lines_66, collapse = "\n")

# 2. SAFETY CHECKS ------------------------------------------------------------

if (grepl("dynamicEpigeneticUI\\(", app_text_66)) {
  stop("dynamicEpigeneticUI() already exists. No changes made.")
}

if (grepl("dynamicEpigeneticServer\\(", app_text_66)) {
  stop("dynamicEpigeneticServer() already exists. No changes made.")
}

if (length(gregexpr(
  ui_anchor_66,
  app_text_66,
  fixed = TRUE
)[[1]]) != 1) {
  stop("UI insertion anchor is not unique. No changes made.")
}

if (length(gregexpr(
  server_anchor_66,
  app_text_66,
  fixed = TRUE
)[[1]]) != 1) {
  stop("Server insertion anchor is not unique. No changes made.")
}

# 3. CREATE A NEW BACKUP -------------------------------------------------------

backup_file_66 <- paste0(
  app_file_66,
  ".backup_before_dynamic_ui_connection"
)

if (file.exists(backup_file_66)) {
  stop(
    "Backup already exists: ",
    backup_file_66,
    "\nNo files modified. Use a different backup name."
  )
}

backup_ok_66 <- file.copy(
  app_file_66,
  backup_file_66
)

if (!backup_ok_66) {
  stop("Could not create backup. No changes made.")
}

cat("Backup created:", backup_file_66, "\n")

# 4. INSERT DYNAMIC UI --------------------------------------------------------

dynamic_ui_block_66 <- paste(
  "    # DYNAMIC METHYLATION VISUALIZATIONS",
  "    shiny::div(",
  "      class = \"section-card\",",
  "      shiny::h3(class = \"section-heading\", \"Dynamic methylation visualizations\"),",
  "      shiny::div(class = \"accent-line\"),",
  "      dynamicEpigeneticUI(\"dynamic\")",
  "    ),",
  "",
  sep = "\n"
)

new_text_66 <- sub(
  ui_anchor_66,
  paste0(dynamic_ui_block_66, ui_anchor_66),
  app_text_66,
  fixed = TRUE
)

# 5. INSERT DYNAMIC SERVER ----------------------------------------------------

dynamic_server_block_66 <- paste(
  "  # Dynamic methylation visualizations",
  "  dynamicEpigeneticServer(",
  "    id = \"dynamic\",",
  "    profile = shiny::reactive({",
  "      current_result <- upload_result$result()",
  "      if (is.null(current_result)) {",
  "        return(NULL)",
  "      }",
  "      current_result$input_data",
  "    }),",
  "    mlp = mlp",
  "  )",
  "",
  sep = "\n"
)

new_text_66 <- sub(
  server_anchor_66,
  paste0(dynamic_server_block_66, server_anchor_66),
  new_text_66,
  fixed = TRUE
)

# 6. WRITE UPDATED APP.R -------------------------------------------------------

writeLines(
  strsplit(new_text_66, "\n", fixed = TRUE)[[1]],
  app_file_66
)

# 7. VALIDATE SYNTAX -----------------------------------------------------------

syntax_ok_66 <- tryCatch(
  {
    parse(app_file_66)
    TRUE
  },
  error = function(e) {
    cat("Syntax error:", conditionMessage(e), "\n")
    FALSE
  }
)

if (!syntax_ok_66) {
  file.copy(
    backup_file_66,
    app_file_66,
    overwrite = TRUE
  )
  stop("Syntax validation failed. Original app.R restored from backup.")
}

# 8. VERIFY UI AND SERVER REFERENCES ------------------------------------------

updated_lines_66 <- readLines(
  app_file_66,
  warn = FALSE
)

updated_text_66 <- paste(
  updated_lines_66,
  collapse = "\n"
)

checks_66 <- c(
  dynamic_ui = grepl(
    'dynamicEpigeneticUI\\("dynamic"\\)',
    updated_text_66
  ),
  dynamic_server = grepl(
    'dynamicEpigeneticServer\\(',
    updated_text_66
  ),
  profile_connection = grepl(
    'current_result\\$input_data',
    updated_text_66
  ),
  mlp_connection = grepl(
    'mlp = mlp',
    updated_text_66,
    fixed = TRUE
  ),
  prediction_preserved = grepl(
    'predictionOverviewUI\\("prediction"\\)',
    updated_text_66
  ),
  evidence_preserved = grepl(
    'epigeneticEvidenceUI\\("epigenetic"\\)',
    updated_text_66
  )
)

cat("\n========== SCRIPT 66 VERIFICATION ==========\n")
print(checks_66)

cat("\nUI AND SERVER ORDER:\n")

order_refs_66 <- grep(
  "predictionOverviewUI|dynamicEpigeneticUI|epigeneticEvidenceUI|predictionOverviewServer|dynamicEpigeneticServer|epigeneticEvidenceServer",
  updated_lines_66
)

cat(
  paste(
    sprintf(
      "%04d: %s",
      order_refs_66,
      updated_lines_66[order_refs_66]
    ),
    collapse = "\n"
  ),
  "\n"
)

stopifnot(
  all(checks_66),
  syntax_ok_66
)

cat("\n========== SCRIPT 66 COMPLETED ==========\n")
cat("Dynamic UI connected: YES\n")
cat("Dynamic server connected: YES\n")
cat("Existing sections preserved: YES\n")
cat("app.R syntax valid: YES\n")
cat("Backup preserved: YES\n")






