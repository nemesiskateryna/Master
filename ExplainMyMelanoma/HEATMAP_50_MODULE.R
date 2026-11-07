
###############################################################################
# SCRIPT 88 — EXACT HEATMAP MODULE INSPECTION
# Read-only. No files modified.
###############################################################################

module_file <- file.path(
  "/Users/katie/master/Master",
  "app",
  "ExplainMyMelanoma",
  "R",
  "06_dynamic_epigenetic_visualizations.R"
)

module_lines <- readLines(module_file, warn = FALSE)

cat("\n========== FUNCTION SIGNATURES ==========\n")

for (pattern in c(
  "^dynamicEpigeneticUI",
  "^dynamicEpigeneticServer"
)) {
  idx <- grep(pattern, module_lines)
  
  for (i in idx) {
    cat(sprintf("%3d | %s\n", i, module_lines[i]))
  }
}

cat("\n========== HEATMAP CODE ==========\n")

heat_idx <- grep(
  "output\\$heatmap\\s*<-\\s*shiny::renderPlot",
  module_lines
)

stopifnot(length(heat_idx) == 1)

heat_start <- heat_idx[1]
heat_end_candidates <- which(
  seq_along(module_lines) > heat_start &
    grepl("^\\s*\\}\\)", module_lines)
)

heat_end <- if (length(heat_end_candidates)) {
  heat_end_candidates[1]
} else {
  min(length(module_lines), heat_start + 100)
}

cat(
  paste(
    sprintf(
      "%3d | %s",
      heat_start:heat_end,
      module_lines[heat_start:heat_end]
    ),
    collapse = "\n"
  ),
  "\n"
)

cat("\n========== UI HEATMAP OUTPUT ==========\n")

ui_idx <- grep(
  "plotOutput.*heatmap",
  module_lines
)

for (i in ui_idx) {
  lo <- max(1, i - 4)
  hi <- min(length(module_lines), i + 5)
  
  cat(
    paste(
      sprintf("%3d | %s", lo:hi, module_lines[lo:hi]),
      collapse = "\n"
    ),
    "\n"
  )
}

cat("\nSCRIPT 88 COMPLETE. No files modified.\n")




###############################################################################
# SCRIPT 89 — HEATMAP UI AND PROFILE AUDIT
# Read-only. No files modified.
###############################################################################

module_file <- file.path(
  "/Users/katie/master/Master",
  "app",
  "ExplainMyMelanoma",
  "R",
  "06_dynamic_epigenetic_visualizations.R"
)

module_lines <- readLines(module_file, warn = FALSE)

# 1. UI function
cat("\n========== UI FUNCTION ==========\n")

ui_start <- grep(
  "^dynamicEpigeneticUI",
  module_lines
)

stopifnot(length(ui_start) == 1)

ui_end <- grep(
  "^dynamicEpigeneticServer",
  module_lines
)

cat(
  paste(
    sprintf(
      "%3d | %s",
      ui_start:(ui_end - 1),
      module_lines[ui_start:(ui_end - 1)]
    ),
    collapse = "\n"
  ),
  "\n"
)

# 2. Selected profile function
cat("\n========== SELECTED PROFILE FUNCTION ==========\n")

profile_idx <- grep(
  "selected_profile\\s*<-\\s*shiny::reactive",
  module_lines
)

cat("Occurrences:", length(profile_idx), "\n")

if (length(profile_idx) == 1) {
  profile_end <- grep(
    "^\\s*output\\$",
    module_lines
  )
  
  profile_end <- profile_end[profile_end > profile_idx]
  
  end_i <- if (length(profile_end)) {
    profile_end[1] - 1
  } else {
    min(length(module_lines), profile_idx + 80)
  }
  
  cat(
    paste(
      sprintf(
        "%3d | %s",
        profile_idx:end_i,
        module_lines[profile_idx:end_i]
      ),
      collapse = "\n"
    ),
    "\n"
  )
}

# 3. Every plotOutput occurrence
cat("\n========== ALL plotOutput REFERENCES ==========\n")

plot_idx <- grep(
  "plotOutput",
  module_lines
)

cat("Occurrences:", length(plot_idx), "\n")

for (i in plot_idx) {
  lo <- max(1, i - 2)
  hi <- min(length(module_lines), i + 5)
  
  cat(
    paste(
      sprintf(
        "%3d | %s",
        lo:hi,
        module_lines[lo:hi]
      ),
      collapse = "\n"
    ),
    "\n\n"
  )
}

# 4. Server function beginning
cat("\n========== SERVER INITIALIZATION ==========\n")

server_start <- grep(
  "^dynamicEpigeneticServer",
  module_lines
)

stopifnot(length(server_start) == 1)

cat(
  paste(
    sprintf(
      "%3d | %s",
      server_start:min(length(module_lines), server_start + 65),
      module_lines[
        server_start:min(length(module_lines), server_start + 65)
      ]
    ),
    collapse = "\n"
  ),
  "\n"
)

cat("\nSCRIPT 89 COMPLETE. No files modified.\n")





