
###############################################################################
# SCRIPT 91I — TOP 50 CpG PERMUTATION IMPORTANCE HEATMAP INTEGRATION
#
# Modifies:
#   1. app/ExplainMyMelanoma/app.R
#   2. app/ExplainMyMelanoma/R/06_dynamic_epigenetic_visualizations.R
#
# Preserves:
#   - selected_profile()
#   - Beta distribution plot
#   - PCA plot
#   - existing MLP model and prediction logic
###############################################################################

cat("\n========== SCRIPT 91I: START ==========\n")

# ---------------------------------------------------------------------------
# 1. PATHS
# ---------------------------------------------------------------------------

project_root <- "/Users/katie/master/Master"

app_file <- file.path(
  project_root,
  "app",
  "ExplainMyMelanoma",
  "app.R"
)

module_file <- file.path(
  project_root,
  "app",
  "ExplainMyMelanoma",
  "R",
  "06_dynamic_epigenetic_visualizations.R"
)

importance_file <- file.path(
  project_root,
  "results",
  "models",
  "TCGA_SKCM_MLP_CpG_permutation_importance_5fold.rds"
)

stopifnot(
  file.exists(app_file),
  file.exists(module_file),
  file.exists(importance_file)
)

cat("Project files found.\n")

# ---------------------------------------------------------------------------
# 2. READ SOURCE FILES
# ---------------------------------------------------------------------------

app_lines <- readLines(app_file, warn = FALSE)
module_lines <- readLines(module_file, warn = FALSE)

app_original <- app_lines
module_original <- module_lines

# ---------------------------------------------------------------------------
# 3. VALIDATE IMPORTANCE OBJECT
# ---------------------------------------------------------------------------

importance_object_91I <- readRDS(importance_file)

if (is.null(importance_object_91I$importance_summary)) {
  stop("importance_summary is missing from the importance object.")
}

importance_columns_91I <- c(
  "CpG",
  "Mean_Delta_LogLoss",
  "Importance_Rank"
)

if (!all(
  importance_columns_91I %in%
  colnames(importance_object_91I$importance_summary)
)) {
  stop("Required importance columns are missing.")
}

importance_table_91I <- importance_object_91I$importance_summary

if (nrow(importance_table_91I) < 50) {
  stop("Fewer than 50 CpGs are available in importance_summary.")
}

if (anyNA(importance_table_91I$Importance_Rank)) {
  stop("Importance_Rank contains missing values.")
}

cat(
  "Importance object validated:",
  nrow(importance_table_91I),
  "CpGs\n"
)

cat("\nTop 5 CpGs:\n")
print(
  head(
    importance_table_91I[
      order(importance_table_91I$Importance_Rank),
      c("CpG", "Mean_Delta_LogLoss", "Importance_Rank")
    ],
    5
  ),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# 4. INSERT IMPORTANCE OBJECT LOADER INTO app.R
# ---------------------------------------------------------------------------

validation_start_91I <- which(
  trimws(app_lines) ==
    'if (length(mlp$selected_cpgs) != 5000) {'
)

if (length(validation_start_91I) != 1) {
  stop("Could not uniquely locate MLP feature validation.")
}

validation_stop_91I <- which(
  seq_along(app_lines) > validation_start_91I &
    trimws(app_lines) ==
    'stop("Unexpected number of selected CpG features.")'
)

if (length(validation_stop_91I) != 1) {
  stop("Could not uniquely locate the MLP validation stop line.")
}

validation_close_91I <- validation_stop_91I + 1

if (trimws(app_lines[validation_close_91I]) != "}") {
  stop("Unexpected structure after MLP feature validation.")
}

if (any(grepl(
  "importance_object <- readRDS",
  app_lines,
  fixed = TRUE
))) {
  stop("Importance object loader already exists. No changes made.")
}

importance_loader_91I <- c(
  "",
  "# Load 5-fold out-of-fold MLP permutation importance",
  "importance_file <- file.path(",
  "  project_root,",
  '  "results",',
  '  "models",',
  '  "TCGA_SKCM_MLP_CpG_permutation_importance_5fold.rds"',
  ")",
  "",
  "if (!file.exists(importance_file)) {",
  '  stop("MLP permutation importance object not found.")',
  "}",
  "",
  "importance_object <- readRDS(importance_file)",
  "",
  "if (is.null(importance_object$importance_summary)) {",
  '  stop("Permutation importance summary is missing.")',
  "}",
  "",
  "required_importance_columns <- c(",
  '  "CpG",',
  '  "Mean_Delta_LogLoss",',
  '  "Importance_Rank"',
  ")",
  "",
  "if (!all(required_importance_columns %in% ",
  "        colnames(importance_object$importance_summary))) {",
  '  stop("Required permutation importance columns are missing.")',
  "}",
  ""
)

app_lines <- append(
  app_lines,
  importance_loader_91I,
  after = validation_close_91I
)

cat("\nImportance loader prepared in app.R.\n")

# ---------------------------------------------------------------------------
# 5. VALIDATE EXISTING DYNAMIC MODULE CALL
# ---------------------------------------------------------------------------

importance_argument_91I <- which(
  trimws(app_lines) == "importance = importance_object"
)

if (length(importance_argument_91I) != 1) {
  stop(
    "Expected exactly one existing importance = importance_object argument."
  )
}

cat("Dynamic module call already contains importance argument.\n")

# ---------------------------------------------------------------------------
# 6. UPDATE MODULE FUNCTION SIGNATURE
# ---------------------------------------------------------------------------

signature_old_91I <- (
  "dynamicEpigeneticServer <- function(id, profile, mlp) {"
)

signature_new_91I <- (
  "dynamicEpigeneticServer <- function(id, profile, mlp, importance) {"
)

signature_hits_91I <- which(
  trimws(module_lines) == signature_old_91I
)

if (length(signature_hits_91I) != 1) {
  stop("Could not uniquely locate the dynamic module function signature.")
}

module_lines[signature_hits_91I] <- signature_new_91I

# ---------------------------------------------------------------------------
# 7. REPLACE ONLY THE HEATMAP UI BLOCK
# ---------------------------------------------------------------------------

ui_start_91I <- which(
  trimws(module_lines) ==
    'shiny::h4("MLP-selected CpG methylation heatmap"),'
)

if (length(ui_start_91I) != 1) {
  stop("Could not uniquely locate the heatmap UI heading.")
}

height_hits_91I <- which(
  seq_along(module_lines) > ui_start_91I &
    trimws(module_lines) == 'height = "430px"'
)

if (length(height_hits_91I) != 1) {
  stop("Could not uniquely locate the original heatmap height.")
}

ui_end_candidates_91I <- which(
  seq_along(module_lines) > height_hits_91I &
    trimws(module_lines) == "),"
)

if (length(ui_end_candidates_91I) < 1) {
  stop("Could not locate the end of the heatmap plotOutput.")
}

ui_end_91I <- ui_end_candidates_91I[1]

heatmap_ui_new_91I <- c(
  '    shiny::h4("Top 50 CpGs ranked by permutation importance"),',
  "",
  "    shiny::p(",
  '      "5-fold out-of-fold MLP importance; colour indicates ',
  'uploaded-sample beta value."',
  "    ),",
  "",
  "    shiny::plotOutput(",
  '      ns("heatmap"),',
  '      height = "850px"',
  "    ),"
)

module_lines <- c(
  module_lines[seq_len(ui_start_91I - 1)],
  heatmap_ui_new_91I,
  module_lines[(ui_end_91I + 1):length(module_lines)]
)

cat("Heatmap UI block prepared.\n")

# ---------------------------------------------------------------------------
# 8. REPLACE ONLY THE HEATMAP SERVER BLOCK
# ---------------------------------------------------------------------------

heat_start_hits_91I <- which(
  grepl(
    "1. METHYLATION HEATMAP",
    module_lines,
    fixed = TRUE
  )
)

heat_end_hits_91I <- which(
  grepl(
    "2. BETA-VALUE DISTRIBUTION",
    module_lines,
    fixed = TRUE
  )
)

if (length(heat_start_hits_91I) != 1 ||
    length(heat_end_hits_91I) != 1) {
  stop("Could not uniquely locate heatmap server block boundaries.")
}

heat_start_91I <- heat_start_hits_91I
heat_end_91I <- heat_end_hits_91I

if (heat_end_91I <= heat_start_91I) {
  stop("Heatmap server block boundaries are in the wrong order.")
}

heatmap_server_new_91I <- c(
  "    # 1. TOP 50 CpGs BY PERMUTATION IMPORTANCE -------------------------------",
  "",
  "    output$heatmap <- shiny::renderPlot({",
  "",
  "      dat <- selected_profile()",
  "      shiny::req(dat)",
  "      shiny::req(importance)",
  "      shiny::req(importance$importance_summary)",
  "",
  "      importance_table <- importance$importance_summary",
  "",
  '      required_columns <- c("CpG", "Importance_Rank")',
  "      shiny::validate(",
  "        shiny::need(",
  "          all(required_columns %in% colnames(importance_table)),",
  '          "Permutation importance data are incomplete."',
  "        )",
  "      )",
  "",
  "      dat$CpG <- as.character(dat$CpG)",
  "      importance_table$CpG <- as.character(importance_table$CpG)",
  "",
  "      importance_table <- importance_table[",
  "        order(importance_table$Importance_Rank),",
  "        ,",
  "        drop = FALSE",
  "      ]",
  "",
  "      ranked_cpgs <- unique(importance_table$CpG)",
  "      ranked_cpgs <- ranked_cpgs[!is.na(ranked_cpgs)]",
  "",
  "      top_cpgs <- ranked_cpgs[ranked_cpgs %in% dat$CpG]",
  "      top_cpgs <- head(top_cpgs, 50)",
  "",
  "      shiny::validate(",
  "        shiny::need(",
  "          length(top_cpgs) == 50,",
  '          "Fewer than 50 ranked CpGs are available in the uploaded profile."',
  "        )",
  "      )",
  "",
  "      plot_data <- dat[dat$CpG %in% top_cpgs, , drop = FALSE]",
  "      plot_data$CpG <- factor(",
  "        plot_data$CpG,",
  "        levels = rev(top_cpgs)",
  "      )",
  "",
  "      ggplot2::ggplot(",
  "        plot_data,",
  "        ggplot2::aes(x = 1, y = CpG, fill = Beta)",
  "      ) +",
  "        ggplot2::geom_tile(",
  "          width = 0.8,",
  "          height = 0.9,",
  "          colour = \"white\",",
  "          linewidth = 0.2",
  "        ) +",
  "        ggplot2::scale_fill_gradientn(",
  '          colours = c("#191970", "#4169E1", "#D4AF37"),',
  '          limits = c(0, 1),',
  '          name = "Beta value"',
  "        ) +",
  "        ggplot2::scale_x_continuous(",
  "          breaks = NULL,",
  "          expand = ggplot2::expansion(mult = c(0.02, 0.02))",
  "        ) +",
  "        ggplot2::scale_y_discrete(",
  "          drop = FALSE",
  "        ) +",
  "        ggplot2::labs(",
  '          title = "Top 50 CpGs ranked by permutation importance",',
  '          subtitle = paste(',
  '            "5-fold out-of-fold MLP importance;",',
  '            "colour indicates uploaded-sample beta value"',
  "          ),",
  '          x = NULL,',
  '          y = "CpG site"',
  "        ) +",
  "        ggplot2::theme_minimal(base_size = 11) +",
  "        ggplot2::theme(",
  "          axis.text.x = ggplot2::element_blank(),",
  "          axis.ticks.x = ggplot2::element_blank(),",
  "          axis.text.y = ggplot2::element_text(",
  "            size = 8,",
  "            colour = \"#191970\"",
  "          ),",
  "          axis.title.y = ggplot2::element_text(",
  "            colour = \"#191970\",",
  "            face = \"bold\"",
  "          ),",
  "          plot.title = ggplot2::element_text(",
  "            colour = \"#191970\",",
  "            face = \"bold\"",
  "          ),",
  "          plot.subtitle = ggplot2::element_text(",
  "            colour = \"#4169E1\"",
  "          ),",
  "          panel.grid = ggplot2::element_blank(),",
  "          legend.position = \"right\"",
  "        )",
  "    })",
  ""
)

module_lines <- c(
  module_lines[seq_len(heat_start_91I - 1)],
  heatmap_server_new_91I,
  module_lines[heat_end_91I:length(module_lines)]
)

cat("Top 50 heatmap server block prepared.\n")

# ---------------------------------------------------------------------------
# 9. CHECK THAT DISTRIBUTION AND PCA BLOCKS WERE PRESERVED
# ---------------------------------------------------------------------------

if (!any(grepl(
  'output\\$distribution <- shiny::renderPlot',
  module_lines
))) {
  stop("Beta distribution output is missing after editing.")
}

if (!any(grepl(
  'output\\$pca <- shiny::renderPlot',
  module_lines
))) {
  stop("PCA output is missing after editing.")
}

if (!any(grepl(
  "selected_profile <- shiny::reactive",
  module_lines,
  fixed = TRUE
))) {
  stop("selected_profile() was not preserved.")
}

cat("selected_profile(), distribution and PCA references preserved.\n")

# ---------------------------------------------------------------------------
# 10. VALIDATE BOTH FILES BEFORE WRITING
# ---------------------------------------------------------------------------

app_text_91I <- paste(app_lines, collapse = "\n")
module_text_91I <- paste(module_lines, collapse = "\n")

parse(text = app_text_91I)
parse(text = module_text_91I)

cat("Both edited R files passed parse().\n")

# ---------------------------------------------------------------------------
# 11. CREATE UNIQUE BACKUPS
# ---------------------------------------------------------------------------

timestamp_91I <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

app_backup_91I <- paste0(
  app_file,
  ".backup_top50_importance_",
  timestamp_91I
)

module_backup_91I <- paste0(
  module_file,
  ".backup_top50_importance_",
  timestamp_91I
)

if (file.exists(app_backup_91I) ||
    file.exists(module_backup_91I)) {
  stop("Backup filename already exists. No files were modified.")
}

# ---------------------------------------------------------------------------
# 12. WRITE TEMPORARY FILES AND VALIDATE AGAIN
# ---------------------------------------------------------------------------

app_temp_91I <- tempfile(
  pattern = "app_top50_",
  tmpdir = dirname(app_file),
  fileext = ".R"
)

module_temp_91I <- tempfile(
  pattern = "module_top50_",
  tmpdir = dirname(module_file),
  fileext = ".R"
)

writeLines(app_lines, app_temp_91I, useBytes = TRUE)
writeLines(module_lines, module_temp_91I, useBytes = TRUE)

parse(file = app_temp_91I)
parse(file = module_temp_91I)

cat("Temporary files passed parse().\n")

# ---------------------------------------------------------------------------
# 13. BACKUP AND REPLACE WITH ROLLBACK
# ---------------------------------------------------------------------------

app_backup_created_91I <- FALSE
module_backup_created_91I <- FALSE

tryCatch({
  
  if (!file.copy(app_file, app_backup_91I, overwrite = FALSE)) {
    stop("Could not create app.R backup.")
  }
  app_backup_created_91I <- TRUE
  
  if (!file.copy(module_file, module_backup_91I, overwrite = FALSE)) {
    stop("Could not create module backup.")
  }
  module_backup_created_91I <- TRUE
  
  if (!file.copy(app_temp_91I, app_file, overwrite = TRUE)) {
    stop("Could not replace app.R.")
  }
  
  if (!file.copy(module_temp_91I, module_file, overwrite = TRUE)) {
    stop("Could not replace dynamic visualization module.")
  }
  
  # Verify saved files parse successfully.
  parse(file = app_file)
  parse(file = module_file)
  
  cat("\nBoth project files updated and parsed successfully.\n")
  
}, error = function(e) {
  
  cat("\nERROR DURING FILE UPDATE:\n")
  cat(conditionMessage(e), "\n")
  
  if (app_backup_created_91I) {
    file.copy(app_backup_91I, app_file, overwrite = TRUE)
    cat("app.R restored from backup.\n")
  }
  
  if (module_backup_created_91I) {
    file.copy(module_backup_91I, module_file, overwrite = TRUE)
    cat("Visualization module restored from backup.\n")
  }
  
  stop("Integration failed. Original files were restored where possible.")
})

# ---------------------------------------------------------------------------
# 14. CLEAN TEMPORARY FILES
# ---------------------------------------------------------------------------

unlink(c(app_temp_91I, module_temp_91I))

# ---------------------------------------------------------------------------
# 15. FINAL READ-ONLY AUDIT
# ---------------------------------------------------------------------------

app_final_91I <- readLines(app_file, warn = FALSE)
module_final_91I <- readLines(module_file, warn = FALSE)

cat("\n========== FINAL AUDIT ==========\n")

cat(
  "Importance loader:",
  any(grepl(
    "importance_object <- readRDS",
    app_final_91I,
    fixed = TRUE
  )),
  "\n"
)

cat(
  "Importance argument:",
  any(trimws(app_final_91I) == "importance = importance_object"),
  "\n"
)

cat(
  "Updated module signature:",
  any(trimws(module_final_91I) ==
        "dynamicEpigeneticServer <- function(id, profile, mlp, importance) {"),
  "\n"
)

cat(
  "Top 50 heatmap title:",
  any(grepl(
    "Top 50 CpGs ranked by permutation importance",
    module_final_91I,
    fixed = TRUE
  )),
  "\n"
)

cat(
  "850px heatmap height:",
  any(trimws(module_final_91I) == 'height = "850px"'),
  "\n"
)

cat(
  "Distribution preserved:",
  any(grepl(
    'output\\$distribution <- shiny::renderPlot',
    module_final_91I
  )),
  "\n"
)

cat(
  "PCA preserved:",
  any(grepl(
    'output\\$pca <- shiny::renderPlot',
    module_final_91I
  )),
  "\n"
)

cat("\nBackup files:\n")
cat(app_backup_91I, "\n")
cat(module_backup_91I, "\n")

cat("\nSCRIPT 91I COMPLETE\n")




