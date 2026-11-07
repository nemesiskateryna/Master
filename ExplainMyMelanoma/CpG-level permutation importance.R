###############################################################################
# SCRIPT 85 — LEAKAGE-FREE 5-FOLD CpG PERMUTATION IMPORTANCE
#
# Model: MLP (nnet)
# Configuration: Config 18
# CV: Original stratified 5-fold assignment
# Metric: Increase in validation log-loss
# Repeats: 5 per CpG per fold
#
# No test-set usage.
# New output files only. Existing files are never overwritten.
###############################################################################

cat("\n========== SCRIPT 85: CpG PERMUTATION IMPORTANCE ==========\n")

# 1. PRE-FLIGHT CHECKS ---------------------------------------------------------

stopifnot(
  exists("beta_81"),
  exists("train_ids_81"),
  exists("y_81"),
  exists("fold_assignment_81"),
  exists("cv_prepared"),
  exists("mlp"),
  exists("oof_84"),
  isTRUE(reproduction_passed_84),
  requireNamespace("nnet", quietly = TRUE)
)

repeats_85 <- 5L
n_folds_85 <- 5L

config_85 <- mlp$selected_configuration

config_id_85 <- config_85$Config_ID
n_pcs_85 <- config_85$n_pcs
size_85 <- config_85$size
decay_85 <- config_85$decay
maxit_85 <- 1000L

final_cpgs_85 <- as.character(mlp$selected_cpgs)

stopifnot(
  length(final_cpgs_85) == 5000L,
  !anyDuplicated(final_cpgs_85),
  all(final_cpgs_85 %in% rownames(beta_81))
)

cat("Final model CpGs:", length(final_cpgs_85), "\n")
cat("Folds:", n_folds_85, "\n")
cat("Permutation repeats:", repeats_85, "\n")
cat("MLP configuration:\n")
print(config_85)

# 2. OUTPUT PATHS AND OVERWRITE PROTECTION ------------------------------------

rds_file_85 <- file.path(
  project_root,
  "results", "models",
  "TCGA_SKCM_MLP_CpG_permutation_importance_5fold.rds"
)

csv_file_85 <- file.path(
  project_root,
  "results", "tables",
  "TCGA_SKCM_MLP_CpG_permutation_importance_5fold.csv"
)

if (file.exists(rds_file_85)) {
  stop("Output RDS already exists. No files overwritten.")
}

if (file.exists(csv_file_85)) {
  stop("Output CSV already exists. No files overwritten.")
}

# 3. LOG-LOSS FUNCTION ---------------------------------------------------------

logloss_85 <- function(probability, actual) {
  
  probability <- pmin(
    pmax(as.numeric(probability), 1e-15),
    1 - 1e-15
  )
  
  actual_binary <- as.numeric(
    as.character(actual) == "Immune-high"
  )
  
  if (
    length(probability) != length(actual_binary) ||
    anyNA(probability) ||
    anyNA(actual_binary)
  ) {
    stop("Invalid probability or outcome vector.")
  }
  
  -mean(
    actual_binary * log(probability) +
      (1 - actual_binary) * log(1 - probability)
  )
}

# 4. STORAGE ------------------------------------------------------------------

fold_importance_85 <- vector("list", n_folds_85)

fold_baselines_85 <- vector("list", n_folds_85)

# 5. MAIN CV LOOP --------------------------------------------------------------

for (fold_85 in seq_len(n_folds_85)) {
  
  cat("\n----------------------------------------\n")
  cat("PROCESSING FOLD:", fold_85, "\n")
  cat("----------------------------------------\n")
  
  prep_85 <- cv_prepared[[fold_85]]
  
  # 5.1. Refit the original fold-specific MLP
  
  x_train_85 <- prep_85$train_x[
    ,
    seq_len(n_pcs_85),
    drop = FALSE
  ]
  
  x_valid_85 <- prep_85$valid_x[
    ,
    seq_len(n_pcs_85),
    drop = FALSE
  ]
  
  y_fold_85 <- prep_85$train_y
  
  y_matrix_85 <- cbind(
    "Immune-low" = as.numeric(
      y_fold_85 == "Immune-low"
    ),
    "Immune-high" = as.numeric(
      y_fold_85 == "Immune-high"
    )
  )
  
  set.seed(
    20261106 + config_id_85 * 100 + fold_85
  )
  
  model_85 <- nnet::nnet(
    x = x_train_85,
    y = y_matrix_85,
    size = size_85,
    decay = decay_85,
    maxit = maxit_85,
    MaxNWts = 5000,
    softmax = TRUE,
    trace = FALSE
  )
  
  if (model_85$convergence != 0) {
    stop("MLP convergence failure in fold ", fold_85)
  }
  
  # 5.2. Reconstruct validation raw CpG matrix
  
  valid_ids_85 <- prep_85$valid_ids
  
  valid_positions_85 <- match(
    valid_ids_85,
    colnames(beta_81)
  )
  
  if (anyNA(valid_positions_85)) {
    stop("Missing validation sample IDs in fold ", fold_85)
  }
  
  fold_cpgs_85 <- as.character(
    prep_85$selected_cpgs
  )
  
  if (
    length(fold_cpgs_85) != 5000L ||
    anyDuplicated(fold_cpgs_85) ||
    !all(fold_cpgs_85 %in% rownames(beta_81))
  ) {
    stop("Invalid fold-specific CpG selection in fold ", fold_85)
  }
  
  beta_valid_85 <- as.matrix(
    beta_81[
      fold_cpgs_85,
      valid_positions_85,
      drop = FALSE
    ]
  )
  
  # Samples x CpGs
  x_valid_raw_85 <- t(beta_valid_85)
  
  colnames(x_valid_raw_85) <- fold_cpgs_85
  rownames(x_valid_raw_85) <- valid_ids_85
  
  if (!all(is.finite(x_valid_raw_85))) {
    stop("Non-finite validation beta values in fold ", fold_85)
  }
  
  # 5.3. Reconstruct PCA projection and verify against saved preprocessing
  
  x_valid_pca_85 <- predict(
    prep_85$pca,
    newdata = x_valid_raw_85
  )
  
  x_valid_pca_85 <- x_valid_pca_85[
    ,
    seq_len(n_pcs_85),
    drop = FALSE
  ]
  
  x_valid_scaled_85 <- scale(
    x_valid_pca_85,
    center = prep_85$score_means[
      seq_len(n_pcs_85)
    ],
    scale = prep_85$score_sds[
      seq_len(n_pcs_85)
    ]
  )
  
  preprocessing_difference_85 <- max(
    abs(
      as.matrix(x_valid_scaled_85) -
        as.matrix(x_valid_85)
    )
  )
  
  cat(
    "PCA reconstruction max difference:",
    preprocessing_difference_85,
    "\n"
  )
  
  if (preprocessing_difference_85 > 1e-8) {
    stop(
      "PCA reconstruction mismatch in fold ",
      fold_85
    )
  }
  
  # 5.4. Baseline predictions
  
  baseline_matrix_85 <- predict(
    model_85,
    newdata = x_valid_scaled_85,
    type = "raw"
  )
  
  baseline_probability_85 <- as.numeric(
    baseline_matrix_85[, "Immune-high"]
  )
  
  baseline_logloss_85 <- logloss_85(
    baseline_probability_85,
    prep_85$valid_y
  )
  
  # Compare against saved OOF predictions
  saved_positions_85 <- match(
    valid_ids_85,
    oof_84$Sample_ID
  )
  
  if (anyNA(saved_positions_85)) {
    stop("Missing saved OOF predictions in fold ", fold_85)
  }
  
  saved_probability_85 <- oof_84$Probability_Immune_high[
    saved_positions_85
  ]
  
  baseline_difference_85 <- max(
    abs(
      baseline_probability_85 -
        saved_probability_85
    )
  )
  
  if (baseline_difference_85 > 1e-8) {
    stop(
      "Baseline prediction mismatch in fold ",
      fold_85
    )
  }
  
  fold_baselines_85[[fold_85]] <- data.frame(
    Fold = fold_85,
    Validation_n = length(valid_ids_85),
    Baseline_LogLoss = baseline_logloss_85,
    OOF_Max_Difference = baseline_difference_85,
    PCA_Max_Difference = preprocessing_difference_85
  )
  
  cat("Baseline log-loss:", baseline_logloss_85, "\n")
  
  # 5.5. Only final-model CpGs that were selected in this fold
  # This ensures that the resulting ranking is compatible with
  # the final model's 5,000-CpG feature set.
  
  cpgs_to_test_85 <- final_cpgs_85[
    final_cpgs_85 %in% fold_cpgs_85
  ]
  
  cpg_positions_85 <- match(
    cpgs_to_test_85,
    fold_cpgs_85
  )
  
  cat(
    "Final-model CpGs tested in this fold:",
    length(cpgs_to_test_85),
    "\n"
  )
  
  if (length(cpgs_to_test_85) == 0L) {
    stop("No overlapping CpGs in fold ", fold_85)
  }
  
  # 5.6. Permutation importance
  
  fold_results_85 <- vector(
    "list",
    length(cpgs_to_test_85)
  )
  
  for (j_85 in seq_along(cpgs_to_test_85)) {
    
    cpg_85 <- cpgs_to_test_85[j_85]
    
    raw_column_85 <- cpg_positions_85[j_85]
    
    delta_values_85 <- numeric(repeats_85)
    
    for (repeat_85 in seq_len(repeats_85)) {
      
      # Deterministic independent permutation seed
      set.seed(
        20261106 +
          fold_85 * 100000 +
          j_85 * 10 +
          repeat_85
      )
      
      permutation_order_85 <- sample(
        seq_len(nrow(x_valid_raw_85))
      )
      
      # Permute one CpG across validation samples
      x_permuted_raw_85 <- x_valid_raw_85
      
      x_permuted_raw_85[, raw_column_85] <-
        x_valid_raw_85[
          permutation_order_85,
          raw_column_85
        ]
      
      # Apply the original fold-specific PCA
      x_permuted_pca_85 <- predict(
        prep_85$pca,
        newdata = x_permuted_raw_85
      )
      
      x_permuted_pca_85 <- x_permuted_pca_85[
        ,
        seq_len(n_pcs_85),
        drop = FALSE
      ]
      
      # Apply original training-fold scaling
      x_permuted_scaled_85 <- scale(
        x_permuted_pca_85,
        center = prep_85$score_means[
          seq_len(n_pcs_85)
        ],
        scale = prep_85$score_sds[
          seq_len(n_pcs_85)
        ]
      )
      
      # Predict with the unchanged fold-specific MLP
      permuted_matrix_85 <- predict(
        model_85,
        newdata = x_permuted_scaled_85,
        type = "raw"
      )
      
      permuted_probability_85 <- as.numeric(
        permuted_matrix_85[, "Immune-high"]
      )
      
      permuted_logloss_85 <- logloss_85(
        permuted_probability_85,
        prep_85$valid_y
      )
      
      delta_values_85[repeat_85] <-
        permuted_logloss_85 - baseline_logloss_85
    }
    
    fold_results_85[[j_85]] <- data.frame(
      CpG = cpg_85,
      Fold = fold_85,
      Baseline_LogLoss = baseline_logloss_85,
      Mean_Delta_LogLoss = mean(delta_values_85),
      SD_Delta_LogLoss = sd(delta_values_85),
      Positive_Repeats = sum(delta_values_85 > 0),
      Repeats = repeats_85
    )
    
    if (
      j_85 %% 500 == 0 ||
      j_85 == length(cpgs_to_test_85)
    ) {
      cat(
        "Fold", fold_85,
        "| CpGs completed:", j_85,
        "/", length(cpgs_to_test_85),
        "\n"
      )
    }
  }
  
  fold_importance_85[[fold_85]] <- do.call(
    rbind,
    fold_results_85
  )
  
  cat(
    "Fold", fold_85,
    "completed. Tested CpGs:",
    nrow(fold_importance_85[[fold_85]]),
    "\n"
  )
  
  rm(
    model_85,
    beta_valid_85,
    x_valid_raw_85,
    x_valid_pca_85,
    x_valid_scaled_85
  )
  
  invisible(gc())
}

# 6. COMBINE FOLD RESULTS -------------------------------------------------------

importance_long_85 <- do.call(
  rbind,
  fold_importance_85
)

baseline_summary_85 <- do.call(
  rbind,
  fold_baselines_85
)

rownames(importance_long_85) <- NULL
rownames(baseline_summary_85) <- NULL

# 7. AGGREGATE CpG IMPORTANCE --------------------------------------------------

importance_summary_85 <- do.call(
  rbind,
  lapply(final_cpgs_85, function(cpg_85) {
    
    z_85 <- importance_long_85[
      importance_long_85$CpG == cpg_85,
      ,
      drop = FALSE
    ]
    
    if (nrow(z_85) == 0L) {
      
      return(data.frame(
        CpG = cpg_85,
        Mean_Delta_LogLoss = NA_real_,
        SD_Fold_Delta_LogLoss = NA_real_,
        Selected_Folds = 0L,
        Positive_Folds = 0L,
        Positive_Repeats_Total = 0L,
        Total_Repeats = 0L
      ))
    }
    
    data.frame(
      CpG = cpg_85,
      Mean_Delta_LogLoss = mean(
        z_85$Mean_Delta_LogLoss
      ),
      SD_Fold_Delta_LogLoss = if (nrow(z_85) > 1L) {
        sd(z_85$Mean_Delta_LogLoss)
      } else {
        NA_real_
      },
      Selected_Folds = nrow(z_85),
      Positive_Folds = sum(
        z_85$Mean_Delta_LogLoss > 0
      ),
      Positive_Repeats_Total = sum(
        z_85$Positive_Repeats
      ),
      Total_Repeats = sum(z_85$Repeats)
    )
  })
)

importance_summary_85 <- importance_summary_85[
  order(
    importance_summary_85$Mean_Delta_LogLoss,
    decreasing = TRUE,
    na.last = TRUE
  ),
]

importance_summary_85$Importance_Rank <- NA_integer_

valid_rank_85 <- which(
  is.finite(importance_summary_85$Mean_Delta_LogLoss)
)

importance_summary_85$Importance_Rank[valid_rank_85] <-
  seq_along(valid_rank_85)

rownames(importance_summary_85) <- NULL

# 8. FINAL QUALITY CHECKS ------------------------------------------------------

cat("\n========== SCRIPT 85 QUALITY CHECKS ==========\n")

cat(
  "Unique final-model CpGs:",
  length(final_cpgs_85),
  "\n"
)

cat(
  "CpGs with at least one fold importance:",
  sum(importance_summary_85$Selected_Folds > 0),
  "\n"
)

cat(
  "CpGs selected in all five folds:",
  sum(importance_summary_85$Selected_Folds == 5),
  "\n"
)

cat(
  "CpGs with positive mean importance:",
  sum(
    importance_summary_85$Mean_Delta_LogLoss > 0,
    na.rm = TRUE
  ),
  "\n"
)

cat("\nFold baseline summary:\n")
print(baseline_summary_85, row.names = FALSE)

cat("\nTop 20 CpGs:\n")
print(
  head(
    importance_summary_85[
      is.finite(importance_summary_85$Mean_Delta_LogLoss),
    ],
    20
  ),
  row.names = FALSE
)

stopifnot(
  nrow(importance_summary_85) == 5000L,
  all(baseline_summary_85$OOF_Max_Difference < 1e-8),
  all(baseline_summary_85$PCA_Max_Difference < 1e-8),
  all(importance_summary_85$Selected_Folds <= 5L),
  all(importance_summary_85$Selected_Folds >= 0L),
  sum(importance_summary_85$Selected_Folds > 0) > 0L
)

# 9. SAVE NEW OUTPUTS ----------------------------------------------------------

importance_object_85 <- list(
  method = "5-fold out-of-fold permutation importance",
  metric = "Increase in validation log-loss",
  repeats_per_cpg_per_fold = repeats_85,
  model = "nnet MLP",
  selected_configuration = config_85,
  fold_assignment = fold_assignment_81,
  baseline_summary = baseline_summary_85,
  fold_importance = importance_long_85,
  importance_summary = importance_summary_85,
  final_model_selected_cpgs = final_cpgs_85,
  test_set_used = FALSE,
  generated_at = Sys.time(),
  session_info = utils::sessionInfo()
)

saveRDS(
  importance_object_85,
  rds_file_85
)

utils::write.csv(
  importance_summary_85,
  csv_file_85,
  row.names = FALSE
)

cat("\n========== SCRIPT 85 COMPLETED ==========\n")

cat("RDS saved:", rds_file_85, "\n")
cat("CSV saved:", csv_file_85, "\n")

cat(
  "Importance records:",
  nrow(importance_long_85),
  "\n"
)

cat(
  "Final ranked CpGs:",
  sum(is.finite(importance_summary_85$Mean_Delta_LogLoss)),
  "\n"
)

cat("Test set used: FALSE\n")
