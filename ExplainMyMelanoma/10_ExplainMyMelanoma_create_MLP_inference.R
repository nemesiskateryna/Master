###############################################################################
# ExplainMyMelanoma AI
# SCRIPT 10 — CREATE MLP INFERENCE FUNCTION
# Creates a NEW app component and tests it on one holdout sample
###############################################################################

project_root <- "/Users/katie/master/Master"

app_r_dir <- file.path(
  project_root,
  "app/ExplainMyMelanoma/R"
)

tcga_file <- file.path(
  project_root,
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune.rds"
)

mlp_file <- file.path(
  project_root,
  "results/models/TCGA_SKCM_MLP_final_model.rds"
)

split_file <- file.path(
  project_root,
  "results/models/TCGA_SKCM_ML_train_test_split.rds"
)

output_file <- file.path(
  app_r_dir,
  "01_mlp_inference.R"
)

stopifnot(dir.exists(app_r_dir))
stopifnot(file.exists(tcga_file))
stopifnot(file.exists(mlp_file))
stopifnot(file.exists(split_file))

if (file.exists(output_file)) {
  stop(
    "Output file already exists. ",
    "No file was overwritten. Choose a new filename."
  )
}

###############################################################################
# 1. CREATE INFERENCE FUNCTION
###############################################################################

inference_code <- c(
  
  "###############################################################################",
  "# ExplainMyMelanoma AI — MLP inference",
  "# Uses the previously trained TCGA-SKCM MLP model.",
  "# Does not retrain the model.",
  "###############################################################################",
  "",
  "predict_mlp_profile <- function(beta_profile, mlp, threshold = 0.5) {",
  "",
  "  if (!is.numeric(threshold) || length(threshold) != 1 ||",
  "      is.na(threshold) || threshold <= 0 || threshold >= 1) {",
  "    stop('Threshold must be a single number between 0 and 1.')",
  "  }",
  "",
  "  required_cpgs <- mlp$selected_cpgs",
  "",
  "  if (is.data.frame(beta_profile)) {",
  "    if (!all(c('CpG', 'Beta') %in% colnames(beta_profile))) {",
  "      stop(\"Input data.frame must contain columns 'CpG' and 'Beta'.\")",
  "    }",
  "    if (anyDuplicated(beta_profile$CpG)) {",
  "      stop('Input contains duplicated CpG identifiers.')",
  "    }",
  "    if (!is.numeric(beta_profile$Beta)) {",
  "      stop('Beta column must be numeric.')",
  "    }",
  "    values <- beta_profile$Beta",
  "    names(values) <- as.character(beta_profile$CpG)",
  "  } else if (is.numeric(beta_profile) && !is.null(names(beta_profile))) {",
  "    values <- beta_profile",
  "  } else {",
  "    stop('Input must be a named numeric vector or a CpG/Beta data.frame.')",
  "  }",
  "",
  "  if (anyDuplicated(names(values))) {",
  "    stop('Input contains duplicated CpG identifiers.')",
  "  }",
  "",
  "  missing_cpgs <- setdiff(required_cpgs, names(values))",
  "",
  "  if (length(missing_cpgs) > 0) {",
  "    stop(",
  "      paste0(",
  "        'Missing ', length(missing_cpgs),",
  "        ' required CpGs. First missing IDs: ',",
  "        paste(head(missing_cpgs, 10), collapse = ', ')",
  "      )",
  "    )",
  "  }",
  "",
  "  values <- values[required_cpgs]",
  "",
  "  if (anyNA(values) || any(!is.finite(values))) {",
  "    stop('Beta values contain missing or non-finite values.')",
  "  }",
  "",
  "  if (any(values < 0 | values > 1)) {",
  "    stop('Beta values must be within the interval [0, 1].')",
  "  }",
  "",
  "  X <- matrix(",
  "    as.numeric(values),",
  "    nrow = 1,",
  "    dimnames = list('Input_sample', required_cpgs)",
  "  )",
  "",
  "  X_scaled <- scale(",
  "    X,",
  "    center = mlp$pca$center,",
  "    scale = mlp$pca$scale",
  "  )",
  "",
  "  pca_scores <- X_scaled %*% mlp$pca$rotation",
  "",
  "  pca_scores_scaled <- scale(",
  "    pca_scores,",
  "    center = mlp$score_means,",
  "    scale = mlp$score_sds",
  "  )",
  "",
  "  probabilities <- predict(",
  "    mlp$model,",
  "    newdata = pca_scores_scaled,",
  "    type = 'raw'",
  "  )",
  "",
  "  if (!identical(colnames(probabilities), mlp$class_levels)) {",
  "    stop('Prediction class order does not match the model class levels.')",
  "  }",
  "",
  "  if (anyNA(probabilities) || any(!is.finite(probabilities))) {",
  "    stop('Model returned invalid probabilities.')",
  "  }",
  "",
  "  if (any(abs(rowSums(probabilities) - 1) > 1e-8)) {",
  "    stop('Predicted probabilities do not sum to 1.')",
  "  }",
  "",
  "  p_low <- as.numeric(probabilities[1, 'Immune-low'])",
  "  p_high <- as.numeric(probabilities[1, 'Immune-high'])",
  "",
  "  predicted_class <- if (p_high >= threshold) {",
  "    'Immune-high'",
  "  } else {",
  "    'Immune-low'",
  "  }",
  "",
  "  list(",
  "    predicted_class = predicted_class,",
  "    probability_immune_low = p_low,",
  "    probability_immune_high = p_high,",
  "    threshold = threshold,",
  "    n_cpgs_used = length(required_cpgs),",
  "    model_name = mlp$model_name,",
  "    interpretation = paste(",
  "      'Prediction of transcriptionally defined immune phenotype;',",
  "      'not a prediction of clinical immunotherapy response.'",
  "    )",
  "  )",
  "}",
  "",
  "###############################################################################"
)

writeLines(inference_code, output_file)

cat("\n========== 1. FUNCTION FILE ==========\n")
cat("Created:", output_file, "\n")
cat("File exists:", file.exists(output_file), "\n")

###############################################################################
# 2. LOAD MODEL AND SOURCE FUNCTION
###############################################################################

mlp <- readRDS(mlp_file)
split <- readRDS(split_file)

source(output_file)

cat("\n========== 2. FUNCTION LOADED ==========\n")
cat(
  "Function available:",
  exists("predict_mlp_profile"),
  "\n"
)

###############################################################################
# 3. PREPARE ONE HOLDOUT SAMPLE
###############################################################################

tcga_object <- readRDS(tcga_file)
beta <- tcga_object$beta_matrix

test_id <- split$test_sample_ids[1]

test_profile <- as.numeric(
  beta[mlp$selected_cpgs, test_id, drop = TRUE]
)

names(test_profile) <- mlp$selected_cpgs

###############################################################################
# 4. RUN NEW INFERENCE FUNCTION
###############################################################################

cat("\n========== 3. NEW FUNCTION PREDICTION ==========\n")

new_prediction <- predict_mlp_profile(
  beta_profile = test_profile,
  mlp = mlp
)

print(new_prediction)

###############################################################################
# 5. INDEPENDENT DIRECT PIPELINE COMPARISON
###############################################################################

cat("\n========== 4. DIRECT PIPELINE COMPARISON ==========\n")

X_test <- matrix(
  test_profile,
  nrow = 1,
  dimnames = list(test_id, mlp$selected_cpgs)
)

X_test_scaled <- scale(
  X_test,
  center = mlp$pca$center,
  scale = mlp$pca$scale
)

direct_pca <- X_test_scaled %*% mlp$pca$rotation

direct_pca_scaled <- scale(
  direct_pca,
  center = mlp$score_means,
  scale = mlp$score_sds
)

direct_prediction <- predict(
  mlp$model,
  newdata = direct_pca_scaled,
  type = "raw"
)

direct_low <- as.numeric(direct_prediction[1, "Immune-low"])
direct_high <- as.numeric(direct_prediction[1, "Immune-high"])

cat("Direct Immune-low probability:", direct_low, "\n")
cat("Function Immune-low probability:",
    new_prediction$probability_immune_low, "\n")

cat("Direct Immune-high probability:", direct_high, "\n")
cat("Function Immune-high probability:",
    new_prediction$probability_immune_high, "\n")

cat(
  "Maximum probability difference:",
  max(
    abs(
      c(direct_low, direct_high) -
        c(
          new_prediction$probability_immune_low,
          new_prediction$probability_immune_high
        )
    )
  ),
  "\n"
)

cat(
  "Predicted class identical:",
  identical(
    colnames(direct_prediction)[max.col(direct_prediction)],
    new_prediction$predicted_class
  ),
  "\n"
)

cat("\nSCRIPT 10 — COMPLETED\n")
###############################################################################




