
###############################################################################
# ExplainMyMelanoma AI
# SCRIPT 05 — TRAINING DATA AND CLASS MAPPING AUDIT
###############################################################################

project_root <- "/Users/katie/master/Master"
model_dir <- file.path(project_root, "results", "models")

mlp <- readRDS(
  file.path(model_dir, "TCGA_SKCM_MLP_final_model.rds")
)

elastic_net <- readRDS(
  file.path(model_dir, "TCGA_SKCM_final_ElasticNet_model.rds")
)

split <- readRDS(
  file.path(model_dir, "TCGA_SKCM_ML_train_test_split.rds")
)

# ---------------------------------------------------------------------------
# 1. Inspect training labels
# ---------------------------------------------------------------------------

cat("\n========== MLP TRAINING LABELS ==========\n")

cat("Unique labels:\n")
print(unique(mlp$training_labels))

cat("\nLabel counts:\n")
print(table(mlp$training_labels))

cat("\nFirst 20 labels:\n")
print(head(mlp$training_labels, 20))

cat("\n========== ELASTIC NET TRAINING LABELS ==========\n")

cat("Unique labels:\n")
print(unique(elastic_net$training_labels))

cat("\nLabel counts:\n")
print(table(elastic_net$training_labels))

cat("\nFirst 20 labels:\n")
print(head(elastic_net$training_labels, 20))

# ---------------------------------------------------------------------------
# 2. Compare saved training labels
# ---------------------------------------------------------------------------

cat("\n========== LABEL CONSISTENCY ==========\n")

cat(
  "MLP and Elastic Net labels identical:",
  identical(
    as.character(mlp$training_labels),
    as.character(elastic_net$training_labels)
  ),
  "\n"
)

cat(
  "MLP and split labels identical:",
  identical(
    as.character(mlp$training_labels),
    as.character(split$train_labels)
  ),
  "\n"
)

cat("\nSplit label counts:\n")
print(table(split$train_labels))

# ---------------------------------------------------------------------------
# 3. Inspect train/test sample identifiers
# ---------------------------------------------------------------------------

cat("\n========== SAMPLE IDENTIFIERS ==========\n")

cat("Training samples:", length(split$train_sample_ids), "\n")
cat("Test samples:", length(split$test_sample_ids), "\n")

cat("First 10 training IDs:\n")
print(head(split$train_sample_ids, 10))

cat("\nFirst 10 test IDs:\n")
print(head(split$test_sample_ids, 10))

cat("\nPatient overlap:\n")
print(split$patient_overlap)

# ---------------------------------------------------------------------------
# 4. Locate processed RDS files
# ---------------------------------------------------------------------------

cat("\n========== PROCESSED DATA FILE INVENTORY ==========\n")

search_dirs <- c(
  file.path(project_root, "data", "processed"),
  file.path(project_root, "results", "models")
)

for (directory in search_dirs) {
  
  cat("\nDIRECTORY:", directory, "\n")
  
  if (dir.exists(directory)) {
    
    files <- list.files(
      directory,
      pattern = "\\.(rds|RData|rda)$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
    
    if (length(files) == 0) {
      
      cat("No RDS/RData files found.\n")
      
    } else {
      
      inventory <- data.frame(
        File = basename(files),
        Relative_path = substring(
          files,
          nchar(project_root) + 2
        ),
        Size_MB = round(
          file.info(files)$size / 1024^2,
          3
        )
      )
      
      print(inventory, row.names = FALSE)
      
    }
    
  } else {
    
    cat("Directory does not exist.\n")
    
  }
}

cat("\nSCRIPT 05 — COMPLETED\n")

###############################################################################

