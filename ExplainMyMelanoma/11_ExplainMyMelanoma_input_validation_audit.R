###############################################################################
# ExplainMyMelanoma AI
# SCRIPT 11 — INPUT VALIDATION AUDIT
# READ-ONLY — TESTS EXISTING INFERENCE FUNCTION
###############################################################################

project_root <- "/Users/katie/master/Master"

app_r_dir <- file.path(
  project_root,
  "app/ExplainMyMelanoma/R"
)

mlp_file <- file.path(
  project_root,
  "results/models/TCGA_SKCM_MLP_final_model.rds"
)

tcga_file <- file.path(
  project_root,
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune.rds"
)

source(file.path(app_r_dir, "01_mlp_inference.R"))

mlp <- readRDS(mlp_file)
tcga_object <- readRDS(tcga_file)

beta <- tcga_object$beta_matrix

test_id <- "PLACEHOLDER"

###############################################################################
# 1. PREPARE VALID TEST PROFILE
###############################################################################

test_id <- readRDS(
  file.path(
    project_root,
    "results/models/TCGA_SKCM_ML_train_test_split.rds"
  )
)$test_sample_ids[1]

profile <- as.numeric(
  beta[mlp$selected_cpgs, test_id, drop = TRUE]
)

names(profile) <- mlp$selected_cpgs

baseline <- predict_mlp_profile(profile, mlp)

cat("\n========== 1. BASELINE ==========\n")
print(baseline)

###############################################################################
# 2. SHUFFLED CpG ORDER
###############################################################################

cat("\n========== 2. SHUFFLED CpG ORDER ==========\n")

set.seed(2026)

shuffled_profile <- profile[sample(seq_along(profile))]

shuffled_prediction <- predict_mlp_profile(
  shuffled_profile,
  mlp
)

cat(
  "Prediction unchanged after shuffling:",
  identical(
    baseline$predicted_class,
    shuffled_prediction$predicted_class
  ),
  "\n"
)

cat(
  "Probability difference:",
  max(abs(
    c(
      baseline$probability_immune_low,
      baseline$probability_immune_high
    ) -
      c(
        shuffled_prediction$probability_immune_low,
        shuffled_prediction$probability_immune_high
      )
  )),
  "\n"
)

###############################################################################
# 3. DATA.FRAME INPUT
###############################################################################

cat("\n========== 3. DATA.FRAME INPUT ==========\n")

profile_df <- data.frame(
  CpG = names(profile),
  Beta = as.numeric(profile)
)

df_prediction <- predict_mlp_profile(
  profile_df,
  mlp
)

cat(
  "Data.frame prediction matches vector:",
  identical(
    baseline$predicted_class,
    df_prediction$predicted_class
  ),
  "\n"
)

###############################################################################
# 4. MISSING CpG TEST
###############################################################################

cat("\n========== 4. MISSING CpG TEST ==========\n")

missing_profile <- profile[-1]

missing_result <- tryCatch(
  predict_mlp_profile(missing_profile, mlp),
  error = function(e) conditionMessage(e)
)

cat("Expected error:\n")
print(missing_result)

###############################################################################
# 5. DUPLICATE CpG TEST
###############################################################################

cat("\n========== 5. DUPLICATE CpG TEST ==========\n")

duplicate_profile <- c(
  profile,
  setNames(profile[1], names(profile)[1])
)

duplicate_result <- tryCatch(
  predict_mlp_profile(duplicate_profile, mlp),
  error = function(e) conditionMessage(e)
)

cat("Expected error:\n")
print(duplicate_result)

###############################################################################
# 6. INVALID BETA VALUE TEST
###############################################################################

cat("\n========== 6. INVALID BETA TEST ==========\n")

invalid_profile <- profile
invalid_profile[1] <- 1.5

invalid_result <- tryCatch(
  predict_mlp_profile(invalid_profile, mlp),
  error = function(e) conditionMessage(e)
)

cat("Expected error:\n")
print(invalid_result)

###############################################################################
# 7. NA TEST
###############################################################################

cat("\n========== 7. NA TEST ==========\n")

na_profile <- profile
na_profile[1] <- NA_real_

na_result <- tryCatch(
  predict_mlp_profile(na_profile, mlp),
  error = function(e) conditionMessage(e)
)

cat("Expected error:\n")
print(na_result)

cat("\nSCRIPT 11 — COMPLETED\n")
###############################################################################
