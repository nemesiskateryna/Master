###############################################################################
# SCRIPT 48 — DYNAMIC PLOT TEST
###############################################################################

cat("\n========== SCRIPT 48: DYNAMIC PLOT TEST ==========\n")

# 1. PREPARE PROFILE IN EXACT MODEL ORDER -------------------------------------

profile_aligned <- profile_test[
  match(mlp_test$selected_cpgs, profile_test$CpG),
  c("CpG", "Beta"),
  drop = FALSE
]

stopifnot(
  identical(profile_aligned$CpG, mlp_test$selected_cpgs),
  nrow(profile_aligned) == 5000
)

# 2. INDEPENDENT PCA PROJECTION CHECK -----------------------------------------

pca_ref <- mlp_test$pca

X_test <- matrix(
  as.numeric(profile_aligned$Beta),
  nrow = 1,
  dimnames = list(
    "Uploaded sample",
    mlp_test$selected_cpgs
  )
)

X_scaled_test <- scale(
  X_test,
  center = pca_ref$center,
  scale = pca_ref$scale
)

new_scores_test <- X_scaled_test %*% pca_ref$rotation

cat("\n--- PCA PROJECTION ---\n")

cat(
  "Projected dimensions:",
  paste(dim(new_scores_test), collapse = " x "),
  "\n"
)

cat(
  "All projected coordinates finite:",
  all(is.finite(new_scores_test)),
  "\n"
)

cat(
  "PC1:",
  as.numeric(new_scores_test[1, "PC1"]),
  "\n"
)

cat(
  "PC2:",
  as.numeric(new_scores_test[1, "PC2"]),
  "\n"
)

stopifnot(
  all(is.finite(new_scores_test)),
  all(c("PC1", "PC2") %in% colnames(new_scores_test))
)

# 3. TEST MODULE SERVER AND PLOT RENDERING -------------------------------------

plot_test_results <- new.env(parent = emptyenv())

server_test <- tryCatch(
  {
    shiny::testServer(
      module_env$dynamicEpigeneticServer,
      args = list(
        profile = shiny::reactive(profile_aligned),
        mlp = mlp_test
      ),
      {
        session$flushReact()
        
        plot_test_results$heatmap <- output$heatmap
        plot_test_results$distribution <- output$distribution
        plot_test_results$pca <- output$pca
      }
    )
    
    TRUE
  },
  error = function(e) {
    plot_test_results$error <- conditionMessage(e)
    FALSE
  }
)

cat("\n--- SERVER TEST ---\n")
cat("testServer completed:", server_test, "\n")

if (!server_test) {
  cat("Server error:\n")
  cat(plot_test_results$error, "\n")
}

# 4. INSPECT PLOT OUTPUTS -------------------------------------------------------

if (server_test) {
  
  plot_names <- c(
    "heatmap",
    "distribution",
    "pca"
  )
  
  plot_status <- vapply(
    plot_names,
    function(nm) {
      !is.null(plot_test_results[[nm]])
    },
    logical(1)
  )
  
  cat("\n--- PLOT OUTPUT CHECK ---\n")
  print(plot_status)
  
  for (nm in plot_names) {
    cat(
      nm,
      "class:",
      paste(class(plot_test_results[[nm]]), collapse = ", "),
      "\n"
    )
  }
  
} else {
  
  plot_status <- c(
    heatmap = FALSE,
    distribution = FALSE,
    pca = FALSE
  )
}

# 5. FINAL SUMMARY -------------------------------------------------------------

cat("\n========== SCRIPT 48 SUMMARY ==========\n")

cat(
  "PCA projection valid:",
  all(is.finite(new_scores_test)),
  "\n"
)

cat(
  "Heatmap output available:",
  plot_status["heatmap"],
  "\n"
)

cat(
  "Distribution output available:",
  plot_status["distribution"],
  "\n"
)

cat(
  "PCA output available:",
  plot_status["pca"],
  "\n"
)

cat("\n========== SCRIPT 48 COMPLETED ==========\n")





###############################################################################
# SCRIPT 49 — FIX DYNAMIC PLOT WARNINGS
###############################################################################

cat("\n========== SCRIPT 49: PLOT WARNING FIX ==========\n")

module_file <- file.path(
  app_dir,
  "R",
  "06_dynamic_epigenetic_visualizations.R"
)

backup_file <- paste0(
  module_file,
  ".backup_before_plot_warning_fix"
)

stopifnot(file.exists(module_file))

if (file.exists(backup_file)) {
  stop("Backup already exists. No files changed:\n", backup_file)
}

# 1. READ AND INSPECT SOURCE ---------------------------------------------------

module_lines <- readLines(
  module_file,
  warn = FALSE
)

original_text <- paste(
  module_lines,
  collapse = "\n"
)

# 2. DEFINE EXACT REPLACEMENTS -------------------------------------------------

old_histogram <- paste(
  'ggplot2::scale_x_continuous(limits = c(0, 1),',
  'breaks = seq(0, 1, by = 0.1))',
  sep = "\n"
)

new_histogram <- paste(
  'ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +',
  'ggplot2::coord_cartesian(xlim = c(0, 1))',
  sep = "\n"
)

old_label <- paste(
  'ggplot2::annotate("label", x = new_sample$PC1,',
  'y = new_sample$PC2, label = "Uploaded sample",',
  'nudge_y = 0.06 * diff(range(training_scores$PC2)),',
  'size = 3, colour = "#191970", fill = "#FFF9E8")',
  sep = "\n"
)

new_label <- paste(
  'ggplot2::geom_label(',
  'data = new_sample,',
  'ggplot2::aes(x = PC1, y = PC2),',
  'inherit.aes = FALSE,',
  'label = "Uploaded sample",',
  'nudge_y = 0.06 * diff(range(training_scores$PC2)),',
  'size = 3,',
  'colour = "#191970",',
  'fill = "#FFF9E8"',
  ')',
  sep = "\n"
)

# 3. VERIFY EXACT MATCHES ------------------------------------------------------

histogram_count <- lengths(
  gregexpr(
    old_histogram,
    original_text,
    fixed = TRUE
  )
)

label_count <- lengths(
  gregexpr(
    old_label,
    original_text,
    fixed = TRUE
  )
)

cat("Histogram target occurrences:", histogram_count, "\n")
cat("PCA label target occurrences:", label_count, "\n")

if (histogram_count != 1 || label_count != 1) {
  stop("Expected source blocks were not found exactly once. No files changed.")
}

# 4. CREATE PATCHED SOURCE IN MEMORY ------------------------------------------

patched_text <- sub(
  old_histogram,
  new_histogram,
  original_text,
  fixed = TRUE
)

patched_text <- sub(
  old_label,
  new_label,
  patched_text,
  fixed = TRUE
)

# 5. VALIDATE BEFORE WRITING ---------------------------------------------------

parse(text = patched_text)

cat("Patched source syntax: VALID\n")

# 6. PRESERVE ORIGINAL MODULE --------------------------------------------------

if (!file.copy(
  module_file,
  backup_file,
  overwrite = FALSE
)) {
  stop("Could not create backup. Original module unchanged.")
}

cat("Backup created:\n", backup_file, "\n")

# 7. WRITE PATCHED MODULE ------------------------------------------------------

writeLines(
  strsplit(patched_text, "\n", fixed = TRUE)[[1]],
  con = module_file,
  useBytes = TRUE
)

# 8. FINAL VALIDATION ----------------------------------------------------------

parse(file = module_file)

check_env <- new.env(parent = globalenv())

sys.source(
  module_file,
  envir = check_env
)

stopifnot(
  exists("dynamicEpigeneticUI", envir = check_env),
  exists("dynamicEpigeneticServer", envir = check_env)
)

cat("\n========== SCRIPT 49 COMPLETED ==========\n")
cat("Histogram warning fix: applied\n")
cat("PCA label positioning fix: applied\n")
cat("Module syntax: VALID\n")
cat("Both module functions: available\n")
cat("Original module backup: preserved\n")






###############################################################################
# SCRIPT 50 — CORRECT AND VERIFY PLOT WARNING FIXES
###############################################################################

cat("\n========== SCRIPT 50: PLOT WARNING FIX RETRY ==========\n")

# 1. PATHS --------------------------------------------------------------------

module_file <- file.path(
  app_dir,
  "R",
  "06_dynamic_epigenetic_visualizations.R"
)

backup_file <- paste0(
  module_file,
  ".backup_before_plot_warning_fix_retry"
)

stopifnot(file.exists(module_file))

# 2. READ MODULE --------------------------------------------------------------

module_lines <- readLines(module_file, warn = FALSE)
module_text <- paste(module_lines, collapse = "\n")

# 3. DEFINE EXACT OLD AND NEW CODE --------------------------------------------

old_histogram <- paste(
  '        ggplot2::scale_x_continuous(',
  '          limits = c(0, 1),',
  '          breaks = seq(0, 1, by = 0.1)',
  '        )',
  sep = "\n"
)

new_histogram <- paste(
  '        ggplot2::scale_x_continuous(',
  '          breaks = seq(0, 1, by = 0.1)',
  '        ) +',
  '',
  '        ggplot2::coord_cartesian(',
  '          xlim = c(0, 1)',
  '        )',
  sep = "\n"
)

old_label <- paste(
  '        ggplot2::annotate(',
  '          "label",',
  '          x = new_sample$PC1,',
  '          y = new_sample$PC2,',
  '          label = "Uploaded sample",',
  '          nudge_y = 0.06 * diff(range(training_scores$PC2)),',
  '          size = 3,',
  '          colour = "#191970",',
  '          fill = "#FFF9E8"',
  '        )',
  sep = "\n"
)

new_label <- paste(
  '        ggplot2::geom_label(',
  '          data = new_sample,',
  '          ggplot2::aes(x = PC1, y = PC2),',
  '          inherit.aes = FALSE,',
  '          label = "Uploaded sample",',
  '          nudge_y = 0.06 * diff(range(training_scores$PC2)),',
  '          size = 3,',
  '          colour = "#191970",',
  '          fill = "#FFF9E8"',
  '        )',
  sep = "\n"
)

# 4. CHECK EXACT MATCHES -------------------------------------------------------

hist_count <- lengths(
  gregexpr(old_histogram, module_text, fixed = TRUE)
)

label_count <- lengths(
  gregexpr(old_label, module_text, fixed = TRUE)
)

cat("Old histogram block found:", hist_count == 1, "\n")
cat("Old PCA label block found:", label_count == 1, "\n")

if (hist_count != 1 || label_count != 1) {
  stop(
    "Expected source blocks were not found exactly once. ",
    "No changes were written."
  )
}

# 5. CREATE ADDITIONAL BACKUP --------------------------------------------------

if (file.exists(backup_file)) {
  stop("Retry backup already exists. No changes were written.")
}

file.copy(
  module_file,
  backup_file,
  overwrite = FALSE
)

cat("Retry backup created:\n", backup_file, "\n")

# 6. APPLY EXACT REPLACEMENTS --------------------------------------------------

patched_text <- gsub(
  old_histogram,
  new_histogram,
  module_text,
  fixed = TRUE
)

patched_text <- gsub(
  old_label,
  new_label,
  patched_text,
  fixed = TRUE
)

# 7. VERIFY PATCHED TEXT BEFORE WRITING ---------------------------------------

stopifnot(
  !grepl(old_histogram, patched_text, fixed = TRUE),
  !grepl(old_label, patched_text, fixed = TRUE),
  grepl(new_histogram, patched_text, fixed = TRUE),
  grepl(new_label, patched_text, fixed = TRUE)
)

# 8. PARSE IN MEMORY -----------------------------------------------------------

parse(text = patched_text)

# 9. WRITE MODULE --------------------------------------------------------------

writeLines(
  strsplit(patched_text, "\n", fixed = TRUE)[[1]],
  con = module_file,
  useBytes = TRUE
)

# 10. VERIFY ACTUAL FILE -------------------------------------------------------

written_text <- paste(
  readLines(module_file, warn = FALSE),
  collapse = "\n"
)

stopifnot(
  !grepl(old_histogram, written_text, fixed = TRUE),
  !grepl(old_label, written_text, fixed = TRUE),
  grepl(new_histogram, written_text, fixed = TRUE),
  grepl(new_label, written_text, fixed = TRUE)
)

parse(file = module_file)

check_env <- new.env(parent = globalenv())

sys.source(
  module_file,
  envir = check_env
)

stopifnot(
  exists("dynamicEpigeneticUI", envir = check_env),
  exists("dynamicEpigeneticServer", envir = check_env)
)

cat("\n========== SCRIPT 50 COMPLETED ==========\n")
cat("Histogram uses coord_cartesian(): VERIFIED\n")
cat("PCA label uses geom_label(): VERIFIED\n")
cat("Old warning-producing code removed: VERIFIED\n")
cat("Module syntax: VALID\n")
cat("Both module functions: AVAILABLE\n")
cat("Retry backup: PRESERVED\n")




###############################################################################
# SCRIPT 51 — RE-TEST DYNAMIC VISUALIZATIONS
###############################################################################

cat("\n========== SCRIPT 51: DYNAMIC PLOT RE-TEST ==========\n")

# 1. CHECK REQUIRED OBJECTS ----------------------------------------------------

required_objects <- c(
  "app_dir",
  "mlp_test",
  "profile_test"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = .GlobalEnv,
    inherits = FALSE
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

module_file <- file.path(
  app_dir,
  "R",
  "06_dynamic_epigenetic_visualizations.R"
)

stopifnot(file.exists(module_file))

# 2. LOAD CURRENT MODULE -------------------------------------------------------

test_env <- new.env(parent = globalenv())

sys.source(
  module_file,
  envir = test_env
)

stopifnot(
  exists("dynamicEpigeneticServer", envir = test_env)
)

# 3. CHECK TEST PROFILE --------------------------------------------------------

cat("\n--- TEST PROFILE ---\n")

stopifnot(
  is.data.frame(profile_test),
  all(c("CpG", "Beta") %in% colnames(profile_test))
)

cat("Profile dimensions:", nrow(profile_test), "rows\n")
cat("Required model CpGs:", length(mlp_test$selected_cpgs), "\n")

profile_positions <- match(
  mlp_test$selected_cpgs,
  profile_test$CpG
)

stopifnot(
  !anyNA(profile_positions)
)

profile_for_test <- profile_test[
  profile_positions,
  c("CpG", "Beta"),
  drop = FALSE
]

stopifnot(
  identical(profile_for_test$CpG, mlp_test$selected_cpgs),
  all(is.finite(profile_for_test$Beta)),
  all(profile_for_test$Beta >= 0 & profile_for_test$Beta <= 1)
)

cat("All 5,000 CpGs matched and validated.\n")

# 4. CAPTURE PLOT WARNINGS -----------------------------------------------------

plot_warnings <- character()
plot_test_results <- list()

warning_handler <- function(w) {
  plot_warnings <<- c(
    plot_warnings,
    conditionMessage(w)
  )
  
  invokeRestart("muffleWarning")
}

# 5. RUN SHINY MODULE TEST -----------------------------------------------------

cat("\n--- SHINY MODULE TEST ---\n")

withCallingHandlers(
  
  shiny::testServer(
    test_env$dynamicEpigeneticServer,
    
    args = list(
      profile = shiny::reactive(profile_for_test),
      mlp = mlp_test
    ),
    
    {
      plot_test_results$heatmap <- output$heatmap
      plot_test_results$distribution <- output$distribution
      plot_test_results$pca <- output$pca
      
      cat("Heatmap output generated:", !is.null(output$heatmap), "\n")
      cat("Distribution output generated:", !is.null(output$distribution), "\n")
      cat("PCA output generated:", !is.null(output$pca), "\n")
    }
  ),
  
  warning = warning_handler
)

# 6. REPORT WARNINGS -----------------------------------------------------------

cat("\n--- WARNING AUDIT ---\n")

if (length(plot_warnings) == 0) {
  cat("No warnings captured during module test.\n")
} else {
  cat("Warnings captured:", length(plot_warnings), "\n")
  
  for (i in seq_along(plot_warnings)) {
    cat(
      sprintf("%02d: %s\n", i, plot_warnings[i])
    )
  }
}

# 7. CHECK SPECIFIC PREVIOUS WARNINGS ------------------------------------------

histogram_warning <- any(
  grepl(
    "Removed .* rows containing missing values or values outside the scale range",
    plot_warnings
  )
)

label_warning <- any(
  grepl(
    "Ignoring unknown parameters: nudge_y",
    plot_warnings,
    fixed = TRUE
  )
)

cat("\nHistogram warning remains:", histogram_warning, "\n")
cat("PCA label warning remains:", label_warning, "\n")

# 8. FINAL SUMMARY -------------------------------------------------------------

cat("\n========== SCRIPT 51 SUMMARY ==========\n")

cat(
  "Heatmap output:",
  !is.null(plot_test_results$heatmap),
  "\n"
)

cat(
  "Distribution output:",
  !is.null(plot_test_results$distribution),
  "\n"
)

cat(
  "PCA output:",
  !is.null(plot_test_results$pca),
  "\n"
)

cat(
  "Histogram warning removed:",
  !histogram_warning,
  "\n"
)

cat(
  "PCA label warning removed:",
  !label_warning,
  "\n"
)

cat(
  "Total warnings:",
  length(plot_warnings),
  "\n"
)

cat("\n========== SCRIPT 51 COMPLETED ==========\n")




###############################################################################
# SCRIPT 52 — FORCE REACTIVE PLOT RENDERING
###############################################################################

cat("\n========== SCRIPT 52: REACTIVE PLOT TEST ==========\n")

# 1. PREPARE PROFILE -----------------------------------------------------------

module_file <- file.path(
  app_dir,
  "R",
  "06_dynamic_epigenetic_visualizations.R"
)

test_env <- new.env(parent = globalenv())

sys.source(module_file, envir = test_env)

positions <- match(
  mlp_test$selected_cpgs,
  profile_test$CpG
)

stopifnot(!anyNA(positions))

profile_for_test <- profile_test[
  positions,
  c("CpG", "Beta"),
  drop = FALSE
]

stopifnot(
  nrow(profile_for_test) == 5000,
  identical(profile_for_test$CpG, mlp_test$selected_cpgs)
)

# 2. INITIALIZE TEST STORAGE ---------------------------------------------------

plot_test_results <- list()
plot_warnings <- character()

# 3. RUN TEST SERVER AND FLUSH REACTIVITY -------------------------------------

withCallingHandlers(
  
  shiny::testServer(
    test_env$dynamicEpigeneticServer,
    
    args = list(
      profile = shiny::reactive(profile_for_test),
      mlp = mlp_test
    ),
    
    {
      # Force pending reactive expressions to execute
      session$flushReact()
      
      # Retrieve rendered outputs
      plot_test_results$heatmap <- output$heatmap
      plot_test_results$distribution <- output$distribution
      plot_test_results$pca <- output$pca
      
      cat("\n--- OUTPUT RETRIEVAL ---\n")
      
      for (nm in c("heatmap", "distribution", "pca")) {
        obj <- plot_test_results[[nm]]
        
        cat(
          nm,
          "| NULL:",
          is.null(obj),
          "| class:",
          paste(class(obj), collapse = ", "),
          "\n"
        )
      }
    }
  ),
  
  warning = function(w) {
    plot_warnings <<- c(
      plot_warnings,
      conditionMessage(w)
    )
    
    invokeRestart("muffleWarning")
  }
)

# 4. FINAL CHECK ---------------------------------------------------------------

cat("\n--- FINAL CHECK ---\n")

for (nm in c("heatmap", "distribution", "pca")) {
  cat(
    nm,
    "available:",
    !is.null(plot_test_results[[nm]]),
    "\n"
  )
}

cat("\nWarnings captured:", length(plot_warnings), "\n")

if (length(plot_warnings) > 0) {
  print(plot_warnings)
}

cat("\n========== SCRIPT 52 COMPLETED ==========\n")



###############################################################################
# SCRIPT 53 — PERSIST AND INSPECT PLOT OUTPUTS
###############################################################################

cat("\n========== SCRIPT 53: PLOT OUTPUT PERSISTENCE ==========\n")

# 1. RELOAD MODULE -------------------------------------------------------------

module_file <- file.path(
  app_dir,
  "R",
  "06_dynamic_epigenetic_visualizations.R"
)

test_env <- new.env(parent = globalenv())
sys.source(module_file, envir = test_env)

# 2. PREPARE PROFILE -----------------------------------------------------------

positions <- match(
  mlp_test$selected_cpgs,
  profile_test$CpG
)

stopifnot(!anyNA(positions))

profile_for_test <- profile_test[
  positions,
  c("CpG", "Beta"),
  drop = FALSE
]

stopifnot(
  nrow(profile_for_test) == 5000,
  identical(profile_for_test$CpG, mlp_test$selected_cpgs)
)

# 3. INITIALIZE EXTERNAL STORAGE -----------------------------------------------

plot_test_results <- list()
plot_warnings <- character()

# 4. RUN TEST SERVER -----------------------------------------------------------

withCallingHandlers(
  
  shiny::testServer(
    test_env$dynamicEpigeneticServer,
    
    args = list(
      profile = shiny::reactive(profile_for_test),
      mlp = mlp_test
    ),
    
    {
      session$flushReact()
      
      plot_test_results$heatmap <<- output$heatmap
      plot_test_results$distribution <<- output$distribution
      plot_test_results$pca <<- output$pca
      
      cat("\n--- OUTPUTS INSIDE TEST SERVER ---\n")
      
      for (nm in c("heatmap", "distribution", "pca")) {
        obj <- plot_test_results[[nm]]
        
        cat(
          nm,
          "| available:",
          !is.null(obj),
          "| class:",
          paste(class(obj), collapse = ", "),
          "\n"
        )
      }
    }
  ),
  
  warning = function(w) {
    plot_warnings <<- c(
      plot_warnings,
      conditionMessage(w)
    )
    
    invokeRestart("muffleWarning")
  }
)

# 5. VERIFY PERSISTED OUTPUTS --------------------------------------------------

cat("\n--- PERSISTED OUTPUTS ---\n")

for (nm in c("heatmap", "distribution", "pca")) {
  cat(
    nm,
    "| available:",
    !is.null(plot_test_results[[nm]]),
    "\n"
  )
}

# 6. INSPECT OUTPUT STRUCTURE --------------------------------------------------

cat("\n--- HEATMAP OUTPUT STRUCTURE ---\n")
str(plot_test_results$heatmap, max.level = 2)

cat("\n--- DISTRIBUTION OUTPUT STRUCTURE ---\n")
str(plot_test_results$distribution, max.level = 2)

cat("\n--- PCA OUTPUT STRUCTURE ---\n")
str(plot_test_results$pca, max.level = 2)

# 7. WARNING AUDIT -------------------------------------------------------------

cat("\nWarnings captured:", length(plot_warnings), "\n")

if (length(plot_warnings) > 0) {
  print(plot_warnings)
}

# 8. FINAL SUMMARY -------------------------------------------------------------

all_outputs_available <- all(
  vapply(
    plot_test_results,
    function(x) !is.null(x),
    logical(1)
  )
)

cat("\nAll three outputs persisted:", all_outputs_available, "\n")
cat("Warnings:", length(plot_warnings), "\n")

cat("\n========== SCRIPT 53 COMPLETED ==========\n")



###############################################################################
# SCRIPT 54 — AUDIT UPLOAD MODULE AND APP INTEGRATION POINTS
###############################################################################

cat("\n========== SCRIPT 54: INTEGRATION AUDIT ==========\n")

# 1. DEFINE PATHS --------------------------------------------------------------

upload_file <- file.path(
  app_dir,
  "R",
  "03_upload_module.R"
)

app_file <- file.path(
  app_dir,
  "app.R"
)

stopifnot(
  file.exists(upload_file),
  file.exists(app_file)
)

# 2. READ FILES ----------------------------------------------------------------

upload_lines <- readLines(upload_file, warn = FALSE)
app_lines <- readLines(app_file, warn = FALSE)

# 3. PRINT UPLOAD MODULE -------------------------------------------------------

cat("\n========== R/03_upload_module.R ==========\n")

cat(
  paste(
    sprintf(
      "%03d: %s",
      seq_along(upload_lines),
      upload_lines
    ),
    collapse = "\n"
  ),
  "\n"
)

# 4. IDENTIFY APP SOURCE LOCATIONS --------------------------------------------

cat("\n========== APP.R: MODULE SOURCES ==========\n")

source_idx <- grep(
  'source\\(file.path\\(app_dir, "R"',
  app_lines
)

if (length(source_idx) > 0) {
  for (i in source_idx) {
    cat(
      sprintf("%03d: %s\n", i, app_lines[i])
    )
  }
} else {
  cat("No matching source statements found.\n")
}

# 5. IDENTIFY APP UI STRUCTURE -------------------------------------------------

cat("\n========== APP.R: UI STRUCTURE ==========\n")

ui_patterns <- c(
  "fluidPage",
  "navbarPage",
  "titlePanel",
  "Upload methylation",
  "Prediction Overview",
  "Epigenetic Evidence",
  "Scientific interpretation",
  "dynamicEpigeneticUI",
  "melanomaUploadUI",
  "predictionOverviewUI"
)

for (pattern in ui_patterns) {
  idx <- grep(pattern, app_lines, fixed = TRUE)
  
  cat("\n", pattern, "| occurrences:", length(idx), "\n")
  
  if (length(idx) > 0) {
    for (i in idx) {
      cat(
        sprintf("%03d: %s\n", i, app_lines[i])
      )
    }
  }
}

# 6. IDENTIFY APP SERVER STRUCTURE ---------------------------------------------

cat("\n========== APP.R: SERVER STRUCTURE ==========\n")

server_patterns <- c(
  "server <-",
  "melanomaUploadServer",
  "predictionOverviewServer",
  "epigeneticEvidenceServer",
  "dynamicEpigeneticServer",
  "upload_result",
  "shinyApp"
)

for (pattern in server_patterns) {
  idx <- grep(pattern, app_lines, fixed = TRUE)
  
  cat("\n", pattern, "| occurrences:", length(idx), "\n")
  
  if (length(idx) > 0) {
    for (i in idx) {
      cat(
        sprintf("%03d: %s\n", i, app_lines[i])
      )
    }
  }
}

# 7. PRINT APP SERVER CONTEXT --------------------------------------------------

server_idx <- grep(
  "server <- function",
  app_lines,
  fixed = TRUE
)

if (length(server_idx) == 1) {
  
  start_idx <- server_idx
  
  end_idx <- min(
    length(app_lines),
    start_idx + 100
  )
  
  cat("\n========== APP.R: SERVER CONTEXT ==========\n")
  
  cat(
    paste(
      sprintf(
        "%03d: %s",
        start_idx:end_idx,
        app_lines[start_idx:end_idx]
      ),
      collapse = "\n"
    ),
    "\n"
  )
  
} else {
  
  cat(
    "\nCould not identify exactly one server function. ",
    "Manual inspection required.\n"
  )
}

# 8. PARSE CHECK ---------------------------------------------------------------

cat("\n========== SYNTAX CHECK ==========\n")

parse(file = upload_file)
parse(file = app_file)

cat("\nUpload module syntax: VALID\n")
cat("app.R syntax: VALID\n")

cat("\n========== SCRIPT 54 COMPLETED ==========\n")


###############################################################################
# SCRIPT 55 — PRESERVE UPLOADED METHYLATION PROFILE
###############################################################################

cat("\n========== SCRIPT 55: PRESERVE INPUT PROFILE ==========\n")

# 1. DEFINE PATHS --------------------------------------------------------------

upload_file <- file.path(
  app_dir,
  "R",
  "03_upload_module.R"
)

backup_file <- paste0(
  upload_file,
  ".backup_before_dynamic_visualizations"
)

stopifnot(file.exists(upload_file))

# 2. READ CURRENT MODULE -------------------------------------------------------

upload_text <- paste(
  readLines(upload_file, warn = FALSE),
  collapse = "\n"
)

# 3. DEFINE EXACT REPLACEMENT --------------------------------------------------

old_block <- paste(
  "        result(list(",
  "          prediction = validated$prediction,",
  "          input_summary = validated$input_summary",
  "        ))",
  sep = "\n"
)

new_block <- paste(
  "        result(list(",
  "          prediction = validated$prediction,",
  "          input_summary = validated$input_summary,",
  "          input_data = validated$input_data",
  "        ))",
  sep = "\n"
)

# 4. VERIFY SOURCE -------------------------------------------------------------

if (!grepl(old_block, upload_text, fixed = TRUE)) {
  stop("Expected result block not found. No changes written.")
}

if (grepl(new_block, upload_text, fixed = TRUE)) {
  stop("New input_data field already exists. No changes written.")
}

# 5. CREATE BACKUP -------------------------------------------------------------

if (file.exists(backup_file)) {
  stop("Backup already exists. No changes written.")
}

backup_ok <- file.copy(
  upload_file,
  backup_file,
  overwrite = FALSE
)

if (!backup_ok) {
  stop("Could not create backup. No changes written.")
}

cat("Backup created:\n", backup_file, "\n")

# 6. APPLY PATCH ---------------------------------------------------------------

patched_text <- sub(
  old_block,
  new_block,
  upload_text,
  fixed = TRUE
)

stopifnot(
  grepl(new_block, patched_text, fixed = TRUE),
  !grepl(old_block, patched_text, fixed = TRUE)
)

# 7. PARSE BEFORE WRITING ------------------------------------------------------

parse(text = patched_text)

# 8. WRITE MODULE --------------------------------------------------------------

writeLines(
  strsplit(patched_text, "\n", fixed = TRUE)[[1]],
  con = upload_file,
  useBytes = TRUE
)

# 9. VERIFY WRITTEN FILE -------------------------------------------------------

written_text <- paste(
  readLines(upload_file, warn = FALSE),
  collapse = "\n"
)

stopifnot(
  grepl(new_block, written_text, fixed = TRUE),
  !grepl(old_block, written_text, fixed = TRUE)
)

parse(file = upload_file)

# 10. SOURCE IN ISOLATED ENVIRONMENT -------------------------------------------

check_env <- new.env(parent = globalenv())

sys.source(
  upload_file,
  envir = check_env
)

stopifnot(
  exists("melanomaUploadUI", envir = check_env),
  exists("melanomaUploadServer", envir = check_env)
)

# 11. DISPLAY MODIFIED SECTION -------------------------------------------------

upload_lines <- readLines(upload_file, warn = FALSE)

result_idx <- grep(
  "input_data = validated\\$input_data",
  upload_lines
)

stopifnot(length(result_idx) == 1)

cat("\n--- UPDATED RESULT BLOCK ---\n")

cat(
  paste(
    upload_lines[(result_idx - 3):(result_idx + 2)],
    collapse = "\n"
  ),
  "\n"
)

# 12. FINAL SUMMARY ------------------------------------------------------------

cat("\n========== SCRIPT 55 COMPLETED ==========\n")
cat("Uploaded input_data preserved: VERIFIED\n")
cat("Prediction calculation changed: NO\n")
cat("CSV validation changed: NO\n")
cat("Module syntax: VALID\n")
cat("Module functions: AVAILABLE\n")
cat("Backup: PRESERVED\n")



