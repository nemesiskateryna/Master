
###############################################################################
# ExplainMyMelanoma AI
# SCRIPT 03 — MODEL STRUCTURE AUDIT
#
# Purpose:
# Inspect saved MLP and Elastic Net objects before app integration.
# No model training, modification, or overwriting.
###############################################################################

cat("\n========================================\n")
cat("EXPLAINMYMELANOMA AI — STEP 3\n")
cat("MODEL STRUCTURE AUDIT\n")
cat("========================================\n\n")

project_root <- "/Users/katie/master/Master"
model_dir <- file.path(project_root, "results", "models")

# ---------------------------------------------------------------------------
# 1. Load saved objects
# ---------------------------------------------------------------------------

mlp <- readRDS(
  file.path(model_dir, "TCGA_SKCM_MLP_final_model.rds")
)

elastic_net <- readRDS(
  file.path(model_dir, "TCGA_SKCM_final_ElasticNet_model.rds")
)

signature <- readRDS(
  file.path(model_dir, "TCGA_SKCM_ElasticNet_CpG_signature.rds")
)

# ---------------------------------------------------------------------------
# 2. MLP structure
# ---------------------------------------------------------------------------

cat("\n========== MLP MODEL ==========\n")

cat("Model name:\n")
print(mlp$model_name)

cat("\nUnderlying model class:\n")
print(class(mlp$model))

cat("\nUnderlying model components:\n")
print(names(mlp$model))

cat("\nNumber of selected CpGs:\n")
print(length(mlp$selected_cpgs))

cat("\nSelected CpG object class:\n")
print(class(mlp$selected_cpgs))

cat("\nFirst 10 selected CpGs:\n")
print(head(mlp$selected_cpgs, 10))

cat("\nNumber of input features:\n")
print(mlp$n_features)

cat("\nNumber of PCA components:\n")
print(mlp$n_pcs)

cat("\nPCA object class:\n")
print(class(mlp$pca))

cat("\nPCA object components:\n")
print(names(mlp$pca))

cat("\nPCA rotation dimensions:\n")
print(dim(mlp$pca$rotation))

cat("\nScore means dimensions:\n")
print(dim(mlp$score_means))

cat("\nScore SDs dimensions:\n")
print(dim(mlp$score_sds))

cat("\nMLP class levels:\n")
print(mlp$class_levels)

cat("\nMLP training dimensions:\n")
print(mlp$training_dimensions)

cat("\nMLP feature selection method:\n")
print(mlp$feature_selection_method)

cat("\nMLP PCA method:\n")
print(mlp$pca_method)

cat("\nMLP score scaling method:\n")
print(mlp$score_scaling_method)

# ---------------------------------------------------------------------------
# 3. Elastic Net structure
# ---------------------------------------------------------------------------

cat("\n========== ELASTIC NET MODEL ==========\n")

cat("Underlying model class:\n")
print(class(elastic_net$model))

cat("\nUnderlying model components:\n")
print(names(elastic_net$model))

cat("\nNumber of selected CpG IDs:\n")
print(length(elastic_net$selected_cpg_ids))

cat("\nSelected CpG IDs class:\n")
print(class(elastic_net$selected_cpg_ids))

cat("\nFirst 10 selected CpG IDs:\n")
print(head(elastic_net$selected_cpg_ids, 10))

cat("\nNumber of selected CpG indices:\n")
print(length(elastic_net$selected_cpg_indices))

cat("\nNumber of features:\n")
print(elastic_net$n_features)

cat("\nNumber of training samples:\n")
print(elastic_net$n_training_samples)

cat("\nAlpha:\n")
print(elastic_net$alpha)

cat("\nLambda:\n")
print(elastic_net$lambda)

cat("\nLambda index:\n")
print(elastic_net$lambda_index)

# ---------------------------------------------------------------------------
# 4. Elastic Net coefficient structure
# ---------------------------------------------------------------------------

cat("\n========== ELASTIC NET COEFFICIENTS ==========\n")

enet_coef <- tryCatch(
  coef(
    elastic_net$model,
    s = elastic_net$lambda
  ),
  error = function(e) e
)

if (inherits(enet_coef, "error")) {
  
  cat("Coefficient extraction error:\n")
  cat(conditionMessage(enet_coef), "\n")
  
} else {
  
  cat("Coefficient object class:\n")
  print(class(enet_coef))
  
  cat("Coefficient dimensions:\n")
  print(dim(enet_coef))
  
  cat("First 10 coefficients:\n")
  print(head(as.matrix(enet_coef), 10))
  
}

# ---------------------------------------------------------------------------
# 5. CpG signature structure
# ---------------------------------------------------------------------------

cat("\n========== CpG SIGNATURE ==========\n")

cat("Signature components:\n")
print(names(signature))

cat("\nNonzero signature dimensions:\n")
print(dim(signature$nonzero_signature))

cat("\nNonzero signature columns:\n")
print(names(signature$nonzero_signature))

cat("\nNumber of nonzero CpGs:\n")
print(nrow(signature$nonzero_signature))

cat("\nFirst 10 signature rows:\n")
print(head(signature$nonzero_signature, 10))

cat("\nNumber of selected CpG IDs in signature:\n")
print(length(signature$selected_cpg_ids))

# ---------------------------------------------------------------------------
# 6. Compare CpG identifiers
# ---------------------------------------------------------------------------

cat("\n========== CpG IDENTIFIER COMPARISON ==========\n")

mlp_cpgs <- as.character(mlp$selected_cpgs)
enet_cpgs <- as.character(elastic_net$selected_cpg_ids)

cat("MLP selected CpGs:", length(mlp_cpgs), "\n")
cat("Elastic Net selected CpGs:", length(enet_cpgs), "\n")

cat(
  "Duplicated MLP CpGs:",
  sum(duplicated(mlp_cpgs)),
  "\n"
)

cat(
  "Duplicated Elastic Net CpGs:",
  sum(duplicated(enet_cpgs)),
  "\n"
)

cat(
  "Shared CpGs:",
  length(intersect(mlp_cpgs, enet_cpgs)),
  "\n"
)

cat(
  "MLP CpGs missing from Elastic Net:",
  length(setdiff(mlp_cpgs, enet_cpgs)),
  "\n"
)

cat(
  "Elastic Net CpGs missing from MLP:",
  length(setdiff(enet_cpgs, mlp_cpgs)),
  "\n"
)

# ---------------------------------------------------------------------------
# 7. Training labels and class definitions
# ---------------------------------------------------------------------------

cat("\n========== CLASS DEFINITIONS ==========\n")

cat("MLP class levels:\n")
print(mlp$class_levels)

cat("\nMLP training label counts:\n")
print(table(mlp$training_labels))

cat("\nElastic Net training label counts:\n")
print(table(elastic_net$training_labels))

# ---------------------------------------------------------------------------
# 8. Final audit summary
# ---------------------------------------------------------------------------

cat("\n========== FINAL AUDIT SUMMARY ==========\n")

cat(
  "MLP model available:",
  !is.null(mlp$model),
  "\n"
)

cat(
  "MLP PCA available:",
  !is.null(mlp$pca),
  "\n"
)

cat(
  "MLP scaling parameters available:",
  !is.null(mlp$score_means) &&
    !is.null(mlp$score_sds),
  "\n"
)

cat(
  "Elastic Net model available:",
  !is.null(elastic_net$model),
  "\n"
)

cat(
  "Elastic Net CpG IDs available:",
  !is.null(elastic_net$selected_cpg_ids),
  "\n"
)

cat(
  "CpG signature available:",
  !is.null(signature$nonzero_signature),
  "\n"
)

cat("\nSCRIPT 03 — MODEL STRUCTURE AUDIT COMPLETED\n")

###############################################################################

