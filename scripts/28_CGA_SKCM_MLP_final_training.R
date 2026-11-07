# TCGA-SKCM FINAL MLP TRAINING
#
# Purpose:
#   1. Tune MLP hyperparameters using 5-fold CV on 368 training samples.
#   2. Perform CpG feature selection independently inside each CV training fold.
#   3. Fit PCA and scaling only on each CV training fold.
#   4. Select the best configuration by mean validation log-loss.
#   5. Refit the final MLP using all 368 training samples.
#
# IMPORTANT:
#   - The existing train/test split is reused.
#   - The 92 holdout samples are NOT used.
#   - No information from holdout is used for tuning or preprocessing.
#   - Existing nested CV and Elastic Net files are not overwritten.
#
# Output:
#   results/models/TCGA_SKCM_MLP_final_model.rds
#   results/models/TCGA_SKCM_MLP_final_tuning.rds
#   results/models/TCGA_SKCM_MLP_final_cv_predictions.rds
#   results/tables/TCGA_SKCM_MLP_final_tuning.csv
#   results/tables/TCGA_SKCM_MLP_final_cv_summary.csv
#
# Random seed: 20261106


# 0. ENVIRONMENT
###############################################################################

set.seed(20261106)

cat("\n============================================\n")
cat("TCGA-SKCM FINAL MLP TRAINING\n")
cat("============================================\n")

required_packages <- c(
  "nnet",
  "DelayedMatrixStats"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(nnet)
library(DelayedMatrixStats)


###############################################################################
# 1. PROJECT PATHS
###############################################################################

integrated_file <- paste0(
  "data/TCGA/processed/",
  "TCGA_SKCM_integrated_multiomics_KNN_immune.rds"
)

split_file <- paste0(
  "results/models/",
  "TCGA_SKCM_ML_train_test_split.rds"
)

output_model <- paste0(
  "results/models/",
  "TCGA_SKCM_MLP_final_model.rds"
)

output_tuning <- paste0(
  "results/models/",
  "TCGA_SKCM_MLP_final_tuning.rds"
)

output_cv_predictions <- paste0(
  "results/models/",
  "TCGA_SKCM_MLP_final_cv_predictions.rds"
)

output_tuning_csv <- paste0(
  "results/tables/",
  "TCGA_SKCM_MLP_final_tuning.csv"
)

output_cv_csv <- paste0(
  "results/tables/",
  "TCGA_SKCM_MLP_final_cv_summary.csv"
)

required_inputs <- c(
  integrated_file,
  split_file
)

if (!all(file.exists(required_inputs))) {
  stop(
    "Required input files are missing:\n",
    paste(
      required_inputs[!file.exists(required_inputs)],
      collapse = "\n"
    )
  )
}

output_files <- c(
  output_model,
  output_tuning,
  output_cv_predictions,
  output_tuning_csv,
  output_cv_csv
)

if (any(file.exists(output_files))) {
  stop(
    "At least one output file already exists. ",
    "To prevent overwriting, the script has stopped:\n",
    paste(
      output_files[file.exists(output_files)],
      collapse = "\n"
    )
  )
}


###############################################################################
# 2. LOAD EXISTING DATA AND SPLIT MANIFEST
###############################################################################

cat("\n========== LOADING INPUT DATA ==========\n")

ml_source <- readRDS(integrated_file)

split_manifest <- readRDS(split_file)

beta_matrix <- ml_source$beta_matrix

if (is.null(beta_matrix)) {
  stop("beta_matrix not found in integrated object.")
}

train_ids <- split_manifest$train_sample_ids
test_ids <- split_manifest$test_sample_ids

y_train <- factor(
  as.character(split_manifest$train_labels),
  levels = c("Immune-low", "Immune-high")
)

names(y_train) <- train_ids

stopifnot(
  length(train_ids) == 368,
  length(test_ids) == 92,
  length(y_train) == 368,
  !anyDuplicated(train_ids),
  !anyDuplicated(test_ids),
  length(intersect(train_ids, test_ids)) == 0,
  !anyNA(y_train)
)

cat("Beta matrix dimensions:\n")
print(dim(beta_matrix))

cat("\nTraining samples:\n")
print(table(y_train))

cat("\nHoldout samples recorded in manifest:\n")
cat(length(test_ids), "\n")


###############################################################################
# 3. ALIGN TRAINING SAMPLES
###############################################################################

matrix_sample_ids <- colnames(beta_matrix)

train_positions <- match(
  train_ids,
  matrix_sample_ids
)

if (anyNA(train_positions)) {
  stop("Some training sample IDs are missing from beta_matrix.")
}

stopifnot(
  identical(
    colnames(beta_matrix)[train_positions],
    train_ids
  )
)

# Only training columns are extracted for model development.
X_train <- beta_matrix[
  ,
  train_positions,
  drop = FALSE
]

stopifnot(
  ncol(X_train) == 368,
  identical(colnames(X_train), train_ids)
)

cat("\nTraining matrix dimensions:\n")
print(dim(X_train))



# 4. DEFINE HYPERPARAMETER GRID
###############################################################################

parameter_grid <- expand.grid(
  n_pcs = c(10, 20, 30),
  size = c(3, 5, 8),
  decay = c(0.01, 0.1),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

parameter_grid$Config_ID <- seq_len(
  nrow(parameter_grid)
)

parameter_grid <- parameter_grid[
  ,
  c("Config_ID", "n_pcs", "size", "decay")
]

cat("\n========== HYPERPARAMETER GRID ==========\n")
print(parameter_grid)

stopifnot(nrow(parameter_grid) == 18)

n_features <- 5000
max_pcs <- 30
maxit_value <- 1000


###############################################################################
# 5. CREATE STRATIFIED 5-FOLD CV
###############################################################################

cat("\n========== CREATING 5-FOLD CV ==========\n")

set.seed(20261106)

fold_assignment <- integer(
  length(y_train)
)

for (class_name in levels(y_train)) {
  
  class_indices <- which(
    y_train == class_name
  )
  
  class_indices <- sample(class_indices)
  
  fold_assignment[class_indices] <- rep(
    1:5,
    length.out = length(class_indices)
  )
}

names(fold_assignment) <- train_ids

fold_distribution <- do.call(
  rbind,
  lapply(1:5, function(fold) {
    
    idx <- which(
      fold_assignment == fold
    )
    
    data.frame(
      Fold = fold,
      N = length(idx),
      Immune_low = sum(y_train[idx] == "Immune-low"),
      Immune_high = sum(y_train[idx] == "Immune-high")
    )
  })
)

print(fold_distribution)

stopifnot(
  all(fold_assignment %in% 1:5),
  all(
    table(
      factor(
        fold_assignment,
        levels = 1:5
      )
    ) > 0
  ),
  all(
    rowSums(
      fold_distribution[, c(
        "Immune_low",
        "Immune_high"
      )]
    ) == fold_distribution$N
  )
)

###############################################################################
# 6. PREPROCESSING FUNCTION
#
# All preprocessing is fitted on the training portion only:
#   - CpG variance selection
#   - PCA
#   - PCA score standardization
###############################################################################

prepare_mlp_data <- function(
    train_sample_ids,
    valid_sample_ids,
    beta_matrix,
    y_train,
    n_features = 5000,
    max_pcs = 30
) {
  
  train_pos <- match(
    train_sample_ids,
    colnames(beta_matrix)
  )
  
  valid_pos <- match(
    valid_sample_ids,
    colnames(beta_matrix)
  )
  
  if (anyNA(train_pos) || anyNA(valid_pos)) {
    stop("Sample ID alignment failed during preprocessing.")
  }
  
  # 1. CpG variance calculated using training samples only
  beta_train_all <- beta_matrix[
    ,
    train_pos,
    drop = FALSE
  ]
  
  feature_variance <- DelayedMatrixStats::rowVars(
    beta_train_all,
    na.rm = TRUE
  )
  
  valid_features <- which(
    is.finite(feature_variance) &
      feature_variance > 0
  )
  
  if (length(valid_features) < n_features) {
    stop("Insufficient variable CpGs.")
  }
  
  selected_order <- order(
    feature_variance[valid_features],
    decreasing = TRUE
  )
  
  selected_indices <- valid_features[
    selected_order[seq_len(n_features)]
  ]
  
  selected_cpgs <- rownames(beta_matrix)[
    selected_indices
  ]
  
  if (is.null(selected_cpgs)) {
    selected_cpgs <- as.character(selected_indices)
  }
  
  # 2. Extract selected CpGs
  beta_train <- as.matrix(
    beta_matrix[
      selected_indices,
      train_pos,
      drop = FALSE
    ]
  )
  
  beta_valid <- as.matrix(
    beta_matrix[
      selected_indices,
      valid_pos,
      drop = FALSE
    ]
  )
  
  # Samples x CpGs
  x_train_raw <- t(beta_train)
  x_valid_raw <- t(beta_valid)
  
  if (
    any(!is.finite(x_train_raw)) ||
    any(!is.finite(x_valid_raw))
  ) {
    stop("Non-finite methylation values detected.")
  }
  
  # 3. PCA fitted on training samples only
  pca_fit <- prcomp(
    x_train_raw,
    center = TRUE,
    scale. = TRUE,
    rank. = max_pcs
  )
  
  train_scores <- pca_fit$x[
    ,
    seq_len(max_pcs),
    drop = FALSE
  ]
  
  valid_scores <- predict(
    pca_fit,
    newdata = x_valid_raw
  )[
    ,
    seq_len(max_pcs),
    drop = FALSE
  ]
  
  # 4. Standardization fitted on training PCA scores only
  score_means <- colMeans(
    train_scores
  )
  
  score_sds <- apply(
    train_scores,
    2,
    sd
  )
  
  if (
    any(!is.finite(score_sds)) ||
    any(score_sds <= 0)
  ) {
    stop("Invalid PCA score standard deviations.")
  }
  
  train_scaled <- scale(
    train_scores,
    center = score_means,
    scale = score_sds
  )
  
  valid_scaled <- scale(
    valid_scores,
    center = score_means,
    scale = score_sds
  )
  
  train_scaled <- as.matrix(train_scaled)
  valid_scaled <- as.matrix(valid_scaled)
  
  rownames(train_scaled) <- train_sample_ids
  rownames(valid_scaled) <- valid_sample_ids
  
  stopifnot(
    nrow(train_scaled) == length(train_sample_ids),
    nrow(valid_scaled) == length(valid_sample_ids),
    ncol(train_scaled) == max_pcs,
    ncol(valid_scaled) == max_pcs,
    all(is.finite(train_scaled)),
    all(is.finite(valid_scaled)),
    identical(rownames(train_scaled), train_sample_ids),
    identical(rownames(valid_scaled), valid_sample_ids)
  )
  
  list(
    train_x = train_scaled,
    valid_x = valid_scaled,
    train_y = y_train[train_sample_ids],
    valid_y = y_train[valid_sample_ids],
    selected_cpgs = selected_cpgs,
    pca = pca_fit,
    score_means = score_means,
    score_sds = score_sds,
    train_ids = train_sample_ids,
    valid_ids = valid_sample_ids
  )
}


#


# 7. PREPARE DATA FOR EACH CV FOLD
#
# Preprocessing is performed once per fold and reused across configurations.
# No validation sample contributes to feature selection, PCA or scaling.
###############################################################################

cat("\n========== PREPARING CV FOLDS ==========\n")

cv_prepared <- vector(
  "list",
  5
)

for (fold in 1:5) {
  
  cat("\nPreparing fold:", fold, "\n")
  
  valid_idx <- which(
    fold_assignment == fold
  )
  
  train_idx <- which(
    fold_assignment != fold
  )
  
  cv_prepared[[fold]] <- prepare_mlp_data(
    train_sample_ids = train_ids[train_idx],
    valid_sample_ids = train_ids[valid_idx],
    beta_matrix = X_train,
    y_train = y_train,
    n_features = n_features,
    max_pcs = max_pcs
  )
  
  cat(
    "Training:",
    nrow(cv_prepared[[fold]]$train_x),
    "| Validation:",
    nrow(cv_prepared[[fold]]$valid_x),
    "| CpGs:",
    length(cv_prepared[[fold]]$selected_cpgs),
    "\n"
  )
}


###############################################################################
# 8. RUN 5-FOLD CV GRID SEARCH
###############################################################################

cat("\n========== RUNNING MLP GRID SEARCH ==========\n")

cv_results <- list()

result_counter <- 1

for (config_idx in seq_len(nrow(parameter_grid))) {
  
  config <- parameter_grid[config_idx, ]
  
  cat(
    "\nConfiguration",
    config$Config_ID,
    "of",
    nrow(parameter_grid),
    "\n"
  )
  
  for (fold in 1:5) {
    
    prep <- cv_prepared[[fold]]
    
    n_pcs <- config$n_pcs
    
    x_train_fold <- prep$train_x[
      ,
      seq_len(n_pcs),
      drop = FALSE
    ]
    
    x_valid_fold <- prep$valid_x[
      ,
      seq_len(n_pcs),
      drop = FALSE
    ]
    
    y_fold <- prep$train_y
    
    # One-hot response
    y_matrix <- cbind(
      "Immune-low" = as.numeric(
        y_fold == "Immune-low"
      ),
      "Immune-high" = as.numeric(
        y_fold == "Immune-high"
      )
    )
    
    set.seed(
      20261106 +
        config$Config_ID * 100 +
        fold
    )
    
    fit_result <- tryCatch({
      
      model <- nnet::nnet(
        x = x_train_fold,
        y = y_matrix,
        size = config$size,
        decay = config$decay,
        maxit = maxit_value,
        MaxNWts = 5000,
        softmax = TRUE,
        trace = FALSE
      )
      
      probabilities <- predict(
        model,
        newdata = x_valid_fold,
        type = "raw"
      )
      
      probabilities <- as.matrix(
        probabilities
      )
      
      if (is.null(colnames(probabilities))) {
        stop("Prediction probability columns are unnamed.")
      }
      
      if (
        !all(
          c("Immune-low", "Immune-high") %in%
          colnames(probabilities)
        )
      ) {
        stop("Unexpected prediction class names.")
      }
      
      prob_high <- probabilities[
        ,
        "Immune-high"
      ]
      
      if (
        any(!is.finite(prob_high)) ||
        any(prob_high < 0) ||
        any(prob_high > 1)
      ) {
        stop("Invalid predicted probabilities.")
      }
      
      actual_binary <- as.numeric(
        prep$valid_y == "Immune-high"
      )
      
      eps <- 1e-15
      
      fold_logloss <- -mean(
        actual_binary * log(
          pmax(prob_high, eps)
        ) +
          (1 - actual_binary) * log(
            pmax(1 - prob_high, eps)
          )
      )
      
      list(
        logloss = fold_logloss,
        convergence = model$convergence,
        probabilities = prob_high,
        actual = as.character(prep$valid_y),
        sample_ids = prep$valid_ids,
        predicted = ifelse(
          prob_high >= 0.5,
          "Immune-high",
          "Immune-low"
        ),
        error = NA_character_
      )
      
    }, error = function(e) {
      
      list(
        logloss = Inf,
        convergence = NA_integer_,
        probabilities = rep(
          NA_real_,
          nrow(x_valid_fold)
        ),
        actual = as.character(prep$valid_y),
        sample_ids = prep$valid_ids,
        predicted = rep(
          NA_character_,
          nrow(x_valid_fold)
        ),
        error = conditionMessage(e)
      )
    })
    
    cv_results[[result_counter]] <- data.frame(
      Config_ID = config$Config_ID,
      n_pcs = config$n_pcs,
      size = config$size,
      decay = config$decay,
      Fold = fold,
      N_validation = length(fit_result$actual),
      LogLoss = fit_result$logloss,
      Convergence = fit_result$convergence,
      Error = fit_result$error,
      stringsAsFactors = FALSE
    )
    
    attr(
      cv_results[[result_counter]],
      "prediction_result"
    ) <- fit_result
    
    result_counter <- result_counter + 1
    
    cat(
      "Fold:",
      fold,
      "| Log-loss:",
      round(fit_result$logloss, 5),
      "| Convergence:",
      fit_result$convergence,
      "\n"
    )
  }
}

cv_results_df <- do.call(
  rbind,
  cv_results
)




# 9. SUMMARIZE CONFIGURATION PERFORMANCE
###############################################################################

cat("\n========== CV GRID SUMMARY ==========\n")

config_summary <- do.call(
  rbind,
  lapply(
    seq_len(nrow(parameter_grid)),
    function(i) {
      
      config_id <- parameter_grid$Config_ID[i]
      
      idx <- which(
        cv_results_df$Config_ID == config_id
      )
      
      current <- cv_results_df[idx, ]
      
      data.frame(
        Config_ID = config_id,
        n_pcs = parameter_grid$n_pcs[i],
        size = parameter_grid$size[i],
        decay = parameter_grid$decay[i],
        Mean_LogLoss = if (
          all(is.finite(current$LogLoss)) &&
          all(current$Convergence == 0)
        ) {
          mean(current$LogLoss)
        } else {
          Inf
        },
        SD_LogLoss = if (
          all(is.finite(current$LogLoss))
        ) {
          sd(current$LogLoss)
        } else {
          NA_real_
        },
        Successful_Folds = sum(
          is.finite(current$LogLoss) &
            current$Convergence == 0
        )
      )
    }
  )
)

config_summary <- config_summary[
  order(config_summary$Mean_LogLoss),
]

print(
  config_summary,
  row.names = FALSE
)

eligible_configs <- config_summary[
  is.finite(config_summary$Mean_LogLoss) &
    config_summary$Successful_Folds == 5,
]

if (nrow(eligible_configs) == 0) {
  stop("No configuration completed all five folds successfully.")
}

best_config <- eligible_configs[1, ]

cat("\n========== SELECTED FINAL CONFIGURATION ==========\n")
print(best_config)

# 10. COLLECT OOF PREDICTIONS FOR SELECTED CONFIGURATION
###############################################################################

cat("\n========== SELECTED CONFIGURATION CV PREDICTIONS ==========\n")

best_config_id <- best_config$Config_ID

best_prediction_objects <- lapply(
  cv_results,
  function(x) {
    if (x$Config_ID[1] == best_config_id) {
      attr(x, "prediction_result")
    } else {
      NULL
    }
  }
)

best_prediction_objects <- Filter(
  Negate(is.null),
  best_prediction_objects
)

selected_cv_predictions <- do.call(
  rbind,
  lapply(
    best_prediction_objects,
    function(x) {
      data.frame(
        Sample_ID = x$sample_ids,
        Actual = x$actual,
        Predicted = x$predicted,
        Probability_Immune_high = x$probabilities,
        stringsAsFactors = FALSE
      )
    }
  )
)

selected_cv_predictions <- selected_cv_predictions[
  match(
    train_ids,
    selected_cv_predictions$Sample_ID
  ),
]

stopifnot(
  nrow(selected_cv_predictions) == 368,
  !anyNA(selected_cv_predictions$Sample_ID),
  !anyDuplicated(selected_cv_predictions$Sample_ID),
  identical(
    selected_cv_predictions$Sample_ID,
    train_ids
  )
)

# 11. FINAL PREPROCESSING ON ALL 368 TRAINING SAMPLES
###############################################################################

cat("\n========== FINAL PREPROCESSING ==========\n")

final_train_pos <- match(
  train_ids,
  colnames(X_train)
)

if (anyNA(final_train_pos)) {
  stop("Final training sample alignment failed.")
}

final_variance <- DelayedMatrixStats::rowVars(
  X_train,
  na.rm = TRUE
)

final_valid_features <- which(
  is.finite(final_variance) &
    final_variance > 0
)

if (length(final_valid_features) < n_features) {
  stop("Insufficient variable CpGs for final model.")
}

final_feature_order <- order(
  final_variance[final_valid_features],
  decreasing = TRUE
)

final_selected_indices <- final_valid_features[
  final_feature_order[seq_len(n_features)]
]

final_selected_cpgs <- rownames(X_train)[
  final_selected_indices
]

if (is.null(final_selected_cpgs)) {
  final_selected_cpgs <- as.character(
    final_selected_indices
  )
}

final_beta_train <- as.matrix(
  X_train[
    final_selected_indices,
    final_train_pos,
    drop = FALSE
  ]
)

final_x_raw <- t(final_beta_train)

if (any(!is.finite(final_x_raw))) {
  stop("Non-finite values in final training matrix.")
}

final_pca <- prcomp(
  final_x_raw,
  center = TRUE,
  scale. = TRUE,
  rank. = max_pcs
)

final_train_scores <- final_pca$x[
  ,
  seq_len(best_config$n_pcs),
  drop = FALSE
]

final_score_means <- colMeans(
  final_train_scores
)

final_score_sds <- apply(
  final_train_scores,
  2,
  sd
)

if (
  any(!is.finite(final_score_sds)) ||
  any(final_score_sds <= 0)
) {
  stop("Invalid final PCA score standard deviations.")
}

final_train_scaled <- scale(
  final_train_scores,
  center = final_score_means,
  scale = final_score_sds
)

final_train_scaled <- as.matrix(
  final_train_scaled
)

rownames(final_train_scaled) <- train_ids


# 12. FIT FINAL MLP ON ALL 368 TRAINING SAMPLES
###############################################################################

cat("\n========== FITTING FINAL MLP ==========\n")

final_y_matrix <- cbind(
  "Immune-low" = as.numeric(
    y_train == "Immune-low"
  ),
  "Immune-high" = as.numeric(
    y_train == "Immune-high"
  )
)

set.seed(20261106)

final_mlp <- nnet::nnet(
  x = final_train_scaled,
  y = final_y_matrix,
  size = best_config$size,
  decay = best_config$decay,
  maxit = maxit_value,
  MaxNWts = 5000,
  softmax = TRUE,
  trace = FALSE
)

cat("Final model convergence:\n")
print(final_mlp$convergence)

if (final_mlp$convergence != 0) {
  stop("Final MLP did not converge.")
}

# 13. FINAL MODEL OBJECT
###############################################################################

final_model_object <- list(
  
  model_name = "TCGA-SKCM MLP",
  
  model = final_mlp,
  
  selected_configuration = best_config,
  
  selected_cpgs = final_selected_cpgs,
  
  pca = final_pca,
  
  score_means = final_score_means,
  
  score_sds = final_score_sds,
  
  n_features = n_features,
  
  n_pcs = best_config$n_pcs,
  
  training_sample_ids = train_ids,
  
  training_labels = y_train,
  
  training_dimensions = dim(final_train_scaled),
  
  class_levels = levels(y_train),
  
  feature_selection_method = "Top 5000 CpGs by training-set variance",
  
  pca_method = "prcomp; fitted on training samples only",
  
  score_scaling_method = "Training PCA score mean and SD",
  
  hyperparameter_selection = "5-fold stratified CV",
  
  tuning_metric = "Mean validation log-loss",
  
  random_seed = 20261106,
  
  holdout_used = FALSE,
  
  nested_cv_oof_metrics = mlp_oof_metrics,
  
  elastic_net_oof_comparison = comparison_table,
  
  session_info = sessionInfo()
)

# 14. QUALITY CONTROL
###############################################################################

cat("\n========== FINAL MODEL QC ==========\n")

final_qc <- c(
  
  training_sample_count =
    length(train_ids) == 368,
  
  training_IDs_unique =
    !anyDuplicated(train_ids),
  
  training_labels_complete =
    !anyNA(y_train),
  
  selected_CpGs_count =
    length(final_selected_cpgs) == 5000,
  
  PCA_dimensions_correct =
    ncol(final_train_scaled) == best_config$n_pcs,
  
  scaled_training_values_finite =
    all(is.finite(final_train_scaled)),
  
  final_model_converged =
    final_mlp$convergence == 0,
  
  selected_configuration_has_5_successful_folds =
    best_config$Successful_Folds == 5,
  
  CV_predictions_cover_368_samples =
    nrow(selected_cv_predictions) == 368,
  
  CV_prediction_IDs_preserved =
    identical(
      selected_cv_predictions$Sample_ID,
      train_ids
    ),
  
  holdout_not_used =
    identical(
      final_model_object$holdout_used,
      FALSE
    )
)

print(final_qc)

stopifnot(all(final_qc))


#
# 15. SAVE FINAL ARTIFACTS
###############################################################################

cat("\n========== SAVING FINAL MLP ARTIFACTS ==========\n")

dir.create(
  "results/models",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

# Recheck output paths immediately before writing
if (any(file.exists(output_files))) {
  stop("An output file appeared before saving. Nothing was overwritten.")
}

final_tuning_object <- list(
  parameter_grid = parameter_grid,
  fold_distribution = fold_distribution,
  fold_assignment = fold_assignment,
  fold_results = cv_results_df,
  configuration_summary = config_summary,
  selected_configuration = best_config,
  seed = 20261106
)

saveRDS(
  final_model_object,
  output_model
)

saveRDS(
  final_tuning_object,
  output_tuning
)

saveRDS(
  selected_cv_predictions,
  output_cv_predictions
)

write.csv(
  config_summary,
  output_tuning_csv,
  row.names = FALSE
)

write.csv(
  selected_cv_predictions,
  output_cv_csv,
  row.names = FALSE
)


###############################################################################
# 16. OUTPUT VERIFICATION
###############################################################################

cat("\n========== OUTPUT FILE VERIFICATION ==========\n")

saved_files <- c(
  output_model,
  output_tuning,
  output_cv_predictions,
  output_tuning_csv,
  output_cv_csv
)

output_verification <- data.frame(
  File = saved_files,
  Exists = file.exists(saved_files),
  Size_MB = round(
    file.info(saved_files)$size / 1024^2,
    3
  )
)

print(output_verification, row.names = FALSE)

stopifnot(
  all(output_verification$Exists),
  all(output_verification$Size_MB >= 0)
)

cat("\n============================================\n")
cat("SCRIPT 24 COMPLETED SUCCESSFULLY\n")
cat("Final MLP trained on all 368 training samples.\n")
cat("Holdout samples were not used.\n")
cat("============================================\n")

###############################################################################
# END OF SCRIPT 24

# Визначаємо шлях до папки з таблицями
tables_dir <- file.path(getwd(), "results", "tables")

# Стара та нова назви файлів
old_name <- file.path(tables_dir, "TCGA_SKCM_MLP_final_cv_summary.csv")
new_name <- file.path(tables_dir, "TCGA_SKCM_MLP_final_cv_predictions.csv")

# Перевіряємо наявність старого файлу та перейменовуємо
if (file.exists(old_name)) {
  file.rename(old_name, new_name)
  cat("Файл успішно перейменовано:\n", new_name, "\n")
} else {
  cat("Попередження: файл зі старою назвою не знайдено за шляхом:\n", old_name, "\n")
}


