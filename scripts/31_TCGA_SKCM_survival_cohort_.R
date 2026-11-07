
###############################################################################
# 30G_TCGA_SKCM_CpG_model_coefficient_verification.R
#
# PURPOSE:
# Verify that the 140-CpG signature coefficients are identical to the
# coefficients extracted directly from the saved final Elastic Net model.
#
# READ-ONLY AUDIT. No model retraining. No source files modified.
###############################################################################

library(glmnet)

cat("\n========== 1. LOAD OBJECTS ==========\n")

signature_object <- readRDS(
  "results/models/TCGA_SKCM_ElasticNet_CpG_signature.rds"
)

final_model_object <- readRDS(
  "results/models/TCGA_SKCM_final_ElasticNet_model.rds"
)

signature <- signature_object$nonzero_signature
model <- final_model_object$model

cat("\nModel class:\n")
print(class(model))

cat("\nModel lambda:\n")
print(final_model_object$lambda)

cat("\nSignature dimensions:\n")
print(dim(signature))

cat("\n========== 2. EXTRACT MODEL COEFFICIENTS ==========\n")

model_coef <- as.matrix(
  stats::coef(
    model,
    s = final_model_object$lambda
  )
)

model_coef_df <- data.frame(
  Feature = rownames(model_coef),
  Model_Coefficient = as.numeric(model_coef[, 1]),
  stringsAsFactors = FALSE
)

nonzero_model_coef <- subset(
  model_coef_df,
  Feature != "(Intercept)" & Model_Coefficient != 0
)

cat("\nTotal model coefficients (excluding intercept):\n")
print(nrow(model_coef_df) - 1)

cat("\nNonzero model coefficients:\n")
print(nrow(nonzero_model_coef))

cat("\nIntercept:\n")
print(
  model_coef_df$Model_Coefficient[
    model_coef_df$Feature == "(Intercept)"
  ]
)

cat("\n========== 3. MATCH SIGNATURE TO MODEL ==========\n")

idx <- match(
  signature$Feature,
  nonzero_model_coef$Feature
)

coefficient_comparison <- data.frame(
  Feature = signature$Feature,
  Signature_Coefficient = signature$Coefficient,
  Model_Coefficient = nonzero_model_coef$Model_Coefficient[idx],
  stringsAsFactors = FALSE
)

coefficient_comparison$Difference <-
  coefficient_comparison$Signature_Coefficient -
  coefficient_comparison$Model_Coefficient

cat("\nMatched CpGs:\n")
print(sum(!is.na(idx)))

cat("\nMissing signature CpGs in model coefficients:\n")
print(signature$Feature[is.na(idx)])

cat("\nCoefficient difference summary:\n")
print(summary(coefficient_comparison$Difference))

cat("\nMaximum absolute coefficient difference:\n")
print(max(
  abs(coefficient_comparison$Difference),
  na.rm = TRUE
))

cat("\nExact equality within tolerance 1e-10:\n")
print(
  all(
    !is.na(idx) &
      abs(coefficient_comparison$Difference) < 1e-10
  )
)

cat("\n========== 4. MODEL FEATURE ORDER ==========\n")

cat("Model selected CpG count:\n")
print(length(final_model_object$selected_cpg_ids))

cat("Model coefficient feature count:\n")
print(nrow(model_coef_df) - 1)

cat("All signature CpGs among model selected CpGs:\n")
print(all(
  signature$Feature %in% final_model_object$selected_cpg_ids
))

###############################################################################
# 5. SAVE AUDIT
###############################################################################

output_csv <- paste0(
  "results/tables/",
  "TCGA_SKCM_survival_CpG_direct_model_coefficient_audit.csv"
)

output_rds <- paste0(
  "results/models/",
  "TCGA_SKCM_survival_CpG_direct_model_coefficient_audit.rds"
)

if (file.exists(output_csv) || file.exists(output_rds)) {
  stop("Output already exists. No files were overwritten.")
}

write.csv(
  coefficient_comparison,
  output_csv,
  row.names = FALSE
)

direct_model_audit <- list(
  model_class = class(model),
  lambda = final_model_object$lambda,
  alpha = final_model_object$alpha,
  intercept = model_coef_df$Model_Coefficient[
    model_coef_df$Feature == "(Intercept)"
  ],
  number_nonzero_model_coefficients = nrow(nonzero_model_coef),
  coefficient_comparison = coefficient_comparison,
  all_signature_features_matched = all(!is.na(idx)),
  coefficients_equal_within_tolerance = all(
    !is.na(idx) &
      abs(coefficient_comparison$Difference) < 1e-10
  ),
  audit_date = Sys.Date()
)

saveRDS(direct_model_audit, output_rds)

cat("\n========== SAVED FILES ==========\n")
print(c(output_csv, output_rds))

cat("\nDIRECT MODEL COEFFICIENT AUDIT COMPLETED\n")



###############################################################################
# 30H_TCGA_SKCM_CpG_signature_score.R
#
# PURPOSE:
# Calculate the fixed Elastic Net CpG signature linear predictor
# for all 460 TCGA-SKCM samples.
#
# No retraining. No survival-based feature selection.
# No modification of existing data.
###############################################################################

library(glmnet)

cat("\n========== 1. LOAD INPUTS ==========\n")

integrated <- readRDS(
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune.rds"
)

final_model_object <- readRDS(
  "results/models/TCGA_SKCM_final_ElasticNet_model.rds"
)

direct_audit <- readRDS(
  "results/models/TCGA_SKCM_survival_CpG_direct_model_coefficient_audit.rds"
)

beta_matrix <- integrated$beta_matrix
immune_scores <- integrated$immune_scores
clinical <- integrated$clinical_data

model <- final_model_object$model

stopifnot(
  isTRUE(direct_audit$all_signature_features_matched),
  isTRUE(direct_audit$coefficients_equal_within_tolerance)
)

cat("Beta dimensions:\n")
print(dim(beta_matrix))

cat("\nImmune scores dimensions:\n")
print(dim(immune_scores))

cat("\nClinical dimensions:\n")
print(dim(clinical))

###############################################################################
# 2. VERIFY FEATURE IDENTIFIERS AND ORDER
###############################################################################

cat("\n========== 2. FEATURE ALIGNMENT ==========\n")

model_feature_ids <- rownames(model$beta)

stopifnot(
  length(model_feature_ids) == final_model_object$n_features,
  !anyDuplicated(model_feature_ids),
  all(model_feature_ids %in% rownames(beta_matrix))
)

cat("Model features:", length(model_feature_ids), "\n")
cat("All model features present in beta matrix:",
    all(model_feature_ids %in% rownames(beta_matrix)), "\n")

###############################################################################
# 3. PREPARE PREDICTION MATRIX
###############################################################################

cat("\n========== 3. PREPARE BETA MATRIX ==========\n")

x_beta <- t(
  as.matrix(
    beta_matrix[
      model_feature_ids,
      ,
      drop = FALSE
    ]
  )
)

stopifnot(
  identical(colnames(x_beta), model_feature_ids),
  all(is.finite(x_beta))
)

cat("Prediction matrix dimensions:\n")
print(dim(x_beta))

###############################################################################
# 4. CALCULATE FIXED MODEL LINEAR PREDICTOR
###############################################################################

cat("\n========== 4. CALCULATE SIGNATURE SCORE ==========\n")

signature_score <- as.numeric(
  predict(
    model,
    newx = x_beta,
    s = final_model_object$lambda,
    type = "link"
  )
)

stopifnot(
  length(signature_score) == nrow(x_beta),
  all(is.finite(signature_score))
)

score_data <- data.frame(
  Sample_ID = rownames(x_beta),
  CpG_signature_score = signature_score,
  stringsAsFactors = FALSE
)

###############################################################################
# 5. VERIFY SAMPLE ALIGNMENT
###############################################################################

cat("\n========== 5. SAMPLE ALIGNMENT ==========\n")

stopifnot(
  !anyDuplicated(score_data$Sample_ID),
  setequal(score_data$Sample_ID, immune_scores$Sample_ID),
  setequal(score_data$Sample_ID, clinical$barcode)
)

cat("Number of scored samples:\n")
print(nrow(score_data))

cat("Unique sample IDs:\n")
print(length(unique(score_data$Sample_ID)))

cat("All samples match immune scores:\n")
print(setequal(score_data$Sample_ID, immune_scores$Sample_ID))

cat("All samples match clinical data:\n")
print(setequal(score_data$Sample_ID, clinical$barcode))

###############################################################################
# 6. DESCRIPTIVE SUMMARY
###############################################################################

cat("\n========== 6. SCORE DISTRIBUTION ==========\n")

print(summary(score_data$CpG_signature_score))

cat("\nScore standard deviation:\n")
print(sd(score_data$CpG_signature_score))

cat("\nNumber of unique scores:\n")
print(length(unique(score_data$CpG_signature_score)))

###############################################################################
# 7. SAVE RESULTS
###############################################################################

output_csv <- paste0(
  "results/tables/",
  "TCGA_SKCM_CpG_signature_scores.csv"
)

output_rds <- paste0(
  "results/models/",
  "TCGA_SKCM_CpG_signature_scores.rds"
)

if (file.exists(output_csv) || file.exists(output_rds)) {
  stop("Output already exists. No files were overwritten.")
}

write.csv(
  score_data,
  output_csv,
  row.names = FALSE
)

score_object <- list(
  scores = score_data,
  model_lambda = final_model_object$lambda,
  model_alpha = final_model_object$alpha,
  model_intercept = direct_audit$intercept,
  number_of_model_features = length(model_feature_ids),
  number_of_nonzero_signature_features = 140,
  number_of_samples = nrow(score_data),
  score_type = "Elastic Net logistic linear predictor",
  source_model = "TCGA_SKCM_final_ElasticNet_model.rds",
  audit_date = Sys.Date()
)

saveRDS(score_object, output_rds)

cat("\n========== SAVED OUTPUTS ==========\n")
print(c(output_csv, output_rds))

cat("\nCpG SIGNATURE SCORE CALCULATION COMPLETED\n")


###############################################################################
# 30I_TCGA_SKCM_survival_cohort_preparation.R
#
# PURPOSE:
# Prepare primary and sensitivity survival cohorts by joining:
# - clinical survival data
# - immune phenotype
# - fixed Elastic Net CpG signature score
#
# No model retraining. No source data modification.
###############################################################################

library(ggplot2)

cat("\n========== 1. LOAD INPUTS ==========\n")

integrated <- readRDS(
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune.rds"
)

score_object <- readRDS(
  "results/models/TCGA_SKCM_CpG_signature_scores.rds"
)

clinical <- integrated$clinical_data
immune_scores <- integrated$immune_scores
score_data <- score_object$scores

cat("Clinical dimensions:\n")
print(dim(clinical))

cat("\nImmune dimensions:\n")
print(dim(immune_scores))

cat("\nScore dimensions:\n")
print(dim(score_data))

###############################################################################
# 2. ALIGN ALL DATA BY SAMPLE ID
###############################################################################

cat("\n========== 2. ALIGN DATA ==========\n")

stopifnot(
  !anyDuplicated(clinical$barcode),
  !anyDuplicated(immune_scores$Sample_ID),
  !anyDuplicated(score_data$Sample_ID),
  setequal(clinical$barcode, immune_scores$Sample_ID),
  setequal(clinical$barcode, score_data$Sample_ID)
)

immune_idx <- match(clinical$barcode, immune_scores$Sample_ID)
score_idx <- match(clinical$barcode, score_data$Sample_ID)

survival_data <- data.frame(
  Sample_ID = clinical$barcode,
  Patient_ID = clinical$patient,
  OS_time = clinical$OS_time_analysis,
  OS_event = clinical$OS_event_analysis,
  Vital_status = clinical$vital_status,
  Immune_group = immune_scores$immune_group[immune_idx],
  TIS_score = immune_scores$TIS_mean_VST[immune_idx],
  CpG_signature_score = score_data$CpG_signature_score[score_idx],
  Age = clinical$age_at_index,
  Sex = clinical$sex_at_birth,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(survival_data) == 460,
  !anyDuplicated(survival_data$Patient_ID),
  all(is.finite(survival_data$CpG_signature_score))
)

cat("Joined samples:", nrow(survival_data), "\n")
cat("Unique patients:", length(unique(survival_data$Patient_ID)), "\n")

###############################################################################
# 3. DEFINE EXCLUSION REASONS
###############################################################################

cat("\n========== 3. SURVIVAL TIME QC ==========\n")

survival_data$Exclusion_reason <- "Included in primary cohort"

survival_data$Exclusion_reason[
  is.na(survival_data$OS_time)
] <- "Missing OS time"

survival_data$Exclusion_reason[
  !is.na(survival_data$OS_time) &
    survival_data$OS_time < 0
] <- "Negative OS time"

survival_data$Exclusion_reason[
  !is.na(survival_data$OS_time) &
    survival_data$OS_time == 0 &
    survival_data$OS_event == 0
] <- "Zero-time censored observation"

survival_data$Exclusion_reason[
  !is.na(survival_data$OS_time) &
    survival_data$OS_time == 0 &
    survival_data$OS_event == 1
] <- "Retained: death at day zero"

cat("\nExclusion/reason counts:\n")
print(table(survival_data$Exclusion_reason, useNA = "ifany"))

###############################################################################
# 4. PRIMARY COHORT
###############################################################################

cat("\n========== 4. PRIMARY SURVIVAL COHORT ==========\n")

primary_cohort <- subset(
  survival_data,
  !is.na(OS_time) &
    OS_time >= 0 &
    !(OS_time == 0 & OS_event == 0)
)

primary_cohort$Exclusion_reason <- NULL

cat("Primary cohort N:\n")
print(nrow(primary_cohort))

cat("Primary cohort events:\n")
print(sum(primary_cohort$OS_event == 1))

cat("Primary cohort censored:\n")
print(sum(primary_cohort$OS_event == 0))

cat("Zero-time events retained:\n")
print(sum(
  primary_cohort$OS_time == 0 &
    primary_cohort$OS_event == 1
))

###############################################################################
# 5. SENSITIVITY COHORT
###############################################################################

cat("\n========== 5. SENSITIVITY COHORT ==========\n")

sensitivity_cohort <- subset(
  survival_data,
  !is.na(OS_time) &
    OS_time > 0
)

sensitivity_cohort$Exclusion_reason <- NULL

cat("Sensitivity cohort N:\n")
print(nrow(sensitivity_cohort))

cat("Sensitivity cohort events:\n")
print(sum(sensitivity_cohort$OS_event == 1))

cat("Sensitivity cohort censored:\n")
print(sum(sensitivity_cohort$OS_event == 0))

###############################################################################
# 6. COHORT AUDIT TABLE
###############################################################################

cat("\n========== 6. COHORT AUDIT ==========\n")

cohort_audit <- data.frame(
  Cohort = c("All integrated samples", "Primary survival cohort",
             "Sensitivity survival cohort"),
  N = c(
    nrow(survival_data),
    nrow(primary_cohort),
    nrow(sensitivity_cohort)
  ),
  Events = c(
    sum(survival_data$OS_event == 1, na.rm = TRUE),
    sum(primary_cohort$OS_event == 1),
    sum(sensitivity_cohort$OS_event == 1)
  ),
  Censored = c(
    sum(survival_data$OS_event == 0, na.rm = TRUE),
    sum(primary_cohort$OS_event == 0),
    sum(sensitivity_cohort$OS_event == 0)
  ),
  stringsAsFactors = FALSE
)

print(cohort_audit, row.names = FALSE)

cat("\nPrimary cohort immune groups:\n")
print(table(primary_cohort$Immune_group))

cat("\nPrimary cohort sex:\n")
print(table(primary_cohort$Sex, useNA = "ifany"))

cat("\nPrimary cohort age summary:\n")
print(summary(primary_cohort$Age))

###############################################################################
# 7. SCORE SUMMARY BY IMMUNE GROUP — CORRECTED
###############################################################################

cat("\n========== 7. CpG SCORE BY IMMUNE GROUP ==========\n")

score_by_immune_table <- do.call(
  rbind,
  lapply(
    split(
      primary_cohort$CpG_signature_score,
      primary_cohort$Immune_group
    ),
    function(x) {
      data.frame(
        N = length(x),
        Mean = mean(x),
        SD = sd(x),
        Median = median(x),
        Q1 = unname(quantile(x, 0.25)),
        Q3 = unname(quantile(x, 0.75))
      )
    }
  )
)

score_by_immune_table$Immune_group <- rownames(score_by_immune_table)
rownames(score_by_immune_table) <- NULL

score_by_immune_table <- score_by_immune_table[
  ,
  c("Immune_group", "N", "Mean", "SD", "Median", "Q1", "Q3")
]

print(score_by_immune_table, row.names = FALSE)


###############################################################################
# 8. FIGURES
###############################################################################

cat("\n========== 8. SAVE FIGURES ==========\n")

p_score_distribution <- ggplot(
  score_data,
  aes(x = CpG_signature_score)
) +
  geom_histogram(bins = 30, color = "black") +
  labs(
    title = "TCGA-SKCM CpG signature score distribution",
    x = "Elastic Net linear predictor",
    y = "Number of samples"
  ) +
  theme_bw()

p_score_immune <- ggplot(
  primary_cohort,
  aes(x = Immune_group, y = CpG_signature_score)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.35) +
  labs(
    title = "CpG signature score by immune phenotype",
    x = "Immune phenotype",
    y = "Elastic Net linear predictor"
  ) +
  theme_bw()

figure_files <- c(
  "results/figures/TCGA_SKCM_CpG_signature_score_distribution.png",
  "results/figures/TCGA_SKCM_CpG_signature_score_by_immune_group.png"
)

ggplot_files <- list(p_score_distribution, p_score_immune)

###############################################################################
# 9. SAVE TABLES AND OBJECTS
###############################################################################

table_files <- c(
  "results/tables/TCGA_SKCM_survival_cohort_audit.csv",
  "results/tables/TCGA_SKCM_survival_cohort_exclusion_reasons.csv",
  "results/tables/TCGA_SKCM_CpG_score_by_immune_group.csv"
)

model_files <- c(
  "results/models/TCGA_SKCM_survival_primary_cohort.rds",
  "results/models/TCGA_SKCM_survival_sensitivity_cohort.rds"
)

all_outputs <- c(table_files, model_files, figure_files)

if (any(file.exists(all_outputs))) {
  stop("At least one output already exists. No files were overwritten.")
}

write.csv(cohort_audit, table_files[1], row.names = FALSE)

write.csv(
  as.data.frame(table(survival_data$Exclusion_reason)),
  table_files[2],
  row.names = FALSE
)

write.csv(score_by_immune_table, table_files[3], row.names = FALSE)

saveRDS(
  list(
    cohort = primary_cohort,
    cohort_audit = cohort_audit,
    definition = paste(
      "OS time >= 0; excludes missing time, negative time,",
      "and zero-time censored observations;",
      "retains death at day zero."
    )
  ),
  model_files[1]
)

saveRDS(
  list(
    cohort = sensitivity_cohort,
    cohort_audit = cohort_audit,
    definition = "OS time > 0; excludes all zero-time and missing/negative observations."
  ),
  model_files[2]
)

for (i in seq_along(ggplot_files)) {
  ggsave(
    filename = figure_files[i],
    plot = ggplot_files[[i]],
    width = 8,
    height = 6,
    dpi = 300
  )
}

cat("\n========== SAVED OUTPUTS ==========\n")
print(all_outputs)

cat("\nSURVIVAL COHORT PREPARATION COMPLETED\n")

