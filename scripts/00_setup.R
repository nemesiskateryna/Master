# Title: 00_setup.R
# Project: Master

#1. Project directory

  project_dir <- getwd()
setwd(project_dir)
#2. Create directories
directories <- c(
  "data/TCGA/raw_saved",
  "data/TCGA/processed",
  "data/GEO/raw",
  "data/GEO/processed",
  "data/metadata",
  "scripts",
  "results/tables",
  "results/figures",
  "results/models",
  "results/logs"
)
for (dir in directories) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
#3. Reproducibility
set.seed(2026)
#4. Required packages
required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "readr",
  "tibble",
  "stringr",
  "survival",
  "survminer",
  "limma",
  "DESeq2",
  "glmnet",
  "matrixStats",
  "AnnotationDbi",
  "Biobase",
  "GEOquery"
)
#Check installed packages
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]
if (length(missing_packages) > 0) {
  message(
    "Потрібно встановити пакети: ",
    paste(missing_packages, collapse = ", ")
  )
}
#5. Load packages
invisible(
  lapply(
    required_packages[
      vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
      )
    ],
    library,
    character.only = TRUE
  )
)
#6. Logging
log_file <- file.path(
  project_dir,
  "results/logs/analysis_log.txt"
)
writeLines(
  paste(
    "Analysis started:",
    Sys.time(),
    "\nR version:",
    R.version.string
  ),
  log_file
)
message("Project initialized: ", project_dir)
