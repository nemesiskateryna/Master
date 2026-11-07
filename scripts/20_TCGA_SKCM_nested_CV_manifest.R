###############################################################################
# SCRIPT 20. TCGA-SKCM NESTED CV MANIFEST
#
# Outer CV: 5 folds
# Inner CV: 4 folds within each outer-training partition
# Stratification: Immune-low / Immune-high
# No model training in this script
###############################################################################

set.seed(20261105)

# ---------------------------------------------------------------------------
# STEP 1. INPUT VALIDATION
# ---------------------------------------------------------------------------

stopifnot(
  exists("X_train"),
  exists("y_train"),
  exists("cv_manifest"),
  exists("split_manifest")
)

stopifnot(
  length(y_train) == ncol(X_train),
  length(cv_manifest$cv_folds) == 5,
  !anyNA(y_train)
)

if (!identical(
  names(y_train),
  colnames(X_train)
)) {
  stop("Sample ID mismatch between X_train and y_train.")
}

if (any(table(y_train) < 4)) {
  stop("Insufficient samples in at least one class.")
}

# ---------------------------------------------------------------------------
# STEP 2. CREATE STRATIFIED INNER FOLDS
# ---------------------------------------------------------------------------

n_outer <- length(cv_manifest$cv_folds)
n_inner <- 4

nested_cv <- vector("list", n_outer)

for (outer_fold in seq_len(n_outer)) {
  
  cat(
    "\nCreating inner folds for outer fold",
    outer_fold,
    "\n"
  )
  
  # Outer validation indices (relative to X_train)
  outer_valid_idx <- cv_manifest$cv_folds[[outer_fold]]
  
  # Outer training indices
  outer_train_idx <- setdiff(
    seq_len(ncol(X_train)),
    outer_valid_idx
  )
  
  y_outer_train <- y_train[outer_train_idx]
  
  # Assign inner validation folds separately within each class
  inner_fold_assignment <- integer(length(outer_train_idx))
  
  for (class_label in levels(y_train)) {
    
    class_idx <- which(y_outer_train == class_label)
    
    shuffled_idx <- sample(class_idx)
    
    assigned_folds <- rep(
      seq_len(n_inner),
      length.out = length(shuffled_idx)
    )
    
    inner_fold_assignment[shuffled_idx] <- assigned_folds
  }
  
  # Store inner validation indices relative to outer training data
  inner_folds <- lapply(
    seq_len(n_inner),
    function(k) {
      which(inner_fold_assignment == k)
    }
  )
  
  names(inner_folds) <- paste0("Inner_Fold_", seq_len(n_inner))
  
  nested_cv[[outer_fold]] <- list(
    outer_fold = outer_fold,
    outer_train_idx = outer_train_idx,
    outer_valid_idx = outer_valid_idx,
    outer_train_sample_ids = names(y_train)[outer_train_idx],
    outer_valid_sample_ids = names(y_train)[outer_valid_idx],
    inner_folds = inner_folds
  )
}

# ---------------------------------------------------------------------------
# STEP 3. BUILD MANIFEST TABLE
# ---------------------------------------------------------------------------

manifest_rows <- list()
row_counter <- 1

for (outer_fold in seq_len(n_outer)) {
  
  outer_data <- nested_cv[[outer_fold]]
  
  for (inner_fold in seq_len(n_inner)) {
    
    inner_valid_idx <- outer_data$inner_folds[[inner_fold]]
    
    for (idx in inner_valid_idx) {
      
      sample_id <- outer_data$outer_train_sample_ids[idx]
      
      manifest_rows[[row_counter]] <- data.frame(
        Outer_Fold = outer_fold,
        Inner_Fold = inner_fold,
        Sample_ID = sample_id,
        Immune_Group = as.character(y_train[
          match(sample_id, names(y_train))
        ]),
        stringsAsFactors = FALSE
      )
      
      row_counter <- row_counter + 1
    }
  }
}

nested_manifest_table <- do.call(
  rbind,
  manifest_rows
)

rownames(nested_manifest_table) <- NULL

# ---------------------------------------------------------------------------
# STEP 4. NESTED CV QUALITY CONTROL
# ---------------------------------------------------------------------------

nested_qc <- c(
  five_outer_folds = length(nested_cv) == 5,
  
  four_inner_folds_each = all(
    vapply(
      nested_cv,
      function(x) length(x$inner_folds) == 4,
      logical(1)
    )
  ),
  
  all_training_samples_assigned = all(
    vapply(
      nested_cv,
      function(x) {
        assigned <- unlist(x$inner_folds)
        identical(
          sort(assigned),
          seq_along(x$outer_train_idx)
        )
      },
      logical(1)
    )
  ),
  
  no_duplicate_inner_assignments = all(
    vapply(
      nested_cv,
      function(x) {
        assigned <- unlist(x$inner_folds)
        !anyDuplicated(assigned)
      },
      logical(1)
    )
  ),
  
  no_outer_validation_in_inner_cv = all(
    vapply(
      nested_cv,
      function(x) {
        length(intersect(
          x$outer_valid_sample_ids,
          unlist(lapply(
            x$inner_folds,
            function(idx) x$outer_train_sample_ids[idx]
          ))
        )) == 0
      },
      logical(1)
    )
  ),
  
  correct_sample_count = nrow(nested_manifest_table) ==
    sum(vapply(
      nested_cv,
      function(x) length(x$outer_train_idx),
      integer(1)
    )),
  
  no_missing_sample_ids = !anyNA(
    nested_manifest_table$Sample_ID
  )
)

# ---------------------------------------------------------------------------
# STEP 5. INNER FOLD CLASS DISTRIBUTION
# ---------------------------------------------------------------------------

inner_distribution <- aggregate(
  Sample_ID ~ Outer_Fold + Inner_Fold + Immune_Group,
  data = nested_manifest_table,
  FUN = length
)

names(inner_distribution)[4] <- "Samples"

cat("\n--- INNER FOLD CLASS DISTRIBUTION ---\n")
print(inner_distribution, row.names = FALSE)

cat("\n--- NESTED CV QC ---\n")
print(nested_qc)

cat(
  "\nNESTED CV MANIFEST QC:",
  if (all(nested_qc)) "PASS" else "FAIL",
  "\n"
)

stopifnot(all(nested_qc))

# ---------------------------------------------------------------------------
# STEP 6. SAVE OUTPUTS WITHOUT OVERWRITING
# ---------------------------------------------------------------------------

rds_path <- "results/models/TCGA_SKCM_nested_CV_manifest.rds"

csv_path <- "results/tables/TCGA_SKCM_nested_CV_manifest_QC.csv"

if (file.exists(rds_path) || file.exists(csv_path)) {
  stop(
    "An output file already exists. No files were overwritten."
  )
}

nested_cv_manifest <- list(
  seed = 20261105,
  outer_folds = 5,
  inner_folds = 4,
  nested_cv = nested_cv,
  manifest_table = nested_manifest_table,
  inner_distribution = inner_distribution,
  qc = nested_qc
)

saveRDS(
  nested_cv_manifest,
  rds_path
)

write.csv(
  nested_manifest_table,
  csv_path,
  row.names = FALSE
)

cat("\n--- SAVED FILES ---\n")
cat(rds_path, "\n")
cat(csv_path, "\n")

cat("\nNESTED CV MANIFEST CREATION: COMPLETE\n")