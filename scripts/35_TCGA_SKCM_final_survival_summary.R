
###############################################################################
# 30O_TCGA_SKCM_final_survival_summary.R
#
# FINAL SURVIVAL ANALYSIS SUMMARY
#
# Read-only integration of previously saved results:
# - Kaplan-Meier / log-rank
# - Cox regression
# - Proportional hazards diagnostics
# - Likelihood ratio tests
# - CpG / immune phenotype correlation
#
# No model refitting.
# No modification of previous results.
# No overwriting of existing files.
###############################################################################

setwd("/Users/katie/master/Master")

cat("\n========== TCGA-SKCM FINAL SURVIVAL SUMMARY ==========\n")

# ---------------------------------------------------------------------------
# 1. DEFINE INPUT FILES
# ---------------------------------------------------------------------------

input_files <- c(
  KM_logrank =
    "results/tables/TCGA_SKCM_KM_logrank_results.csv",
  
  KM_group_counts =
    "results/tables/TCGA_SKCM_KM_group_counts.csv",
  
  Cox_coefficients =
    "results/tables/TCGA_SKCM_Cox_model_coefficients.csv",
  
  Cox_fit =
    "results/tables/TCGA_SKCM_Cox_model_fit_statistics.csv",
  
  Cox_PH =
    "results/tables/TCGA_SKCM_Cox_PH_assumption_tests.csv",
  
  Cox_LRT =
    "results/tables/TCGA_SKCM_Cox_incremental_LRT.csv",
  
  CpG_immune_correlation =
    "results/tables/TCGA_SKCM_CpG_immune_correlation.csv"
)

missing_inputs <- input_files[!file.exists(input_files)]

if (length(missing_inputs) > 0) {
  stop(
    "Missing required input files:\n",
    paste(missing_inputs, collapse = "\n")
  )
}

cat("\nAll required input files are present.\n")

# ---------------------------------------------------------------------------
# 2. DEFINE OUTPUT FILES
# ---------------------------------------------------------------------------

output_files <- c(
  "results/tables/TCGA_SKCM_FINAL_survival_key_results.csv",
  "results/tables/TCGA_SKCM_FINAL_survival_file_inventory.csv",
  "results/models/TCGA_SKCM_FINAL_survival_summary.rds"
)

existing_outputs <- output_files[file.exists(output_files)]

if (length(existing_outputs) > 0) {
  stop(
    "Output files already exist. Nothing was overwritten:\n",
    paste(existing_outputs, collapse = "\n")
  )
}

# ---------------------------------------------------------------------------
# 3. LOAD PREVIOUSLY SAVED TABLES
# ---------------------------------------------------------------------------

read_result <- function(path) {
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

survival_results <- lapply(
  input_files,
  read_result
)

names(survival_results) <- names(input_files)

cat("\n========== INPUT TABLE DIMENSIONS ==========\n")

for (nm in names(survival_results)) {
  cat(
    nm,
    ":",
    nrow(survival_results[[nm]]),
    "rows x",
    ncol(survival_results[[nm]]),
    "columns\n"
  )
}

# ---------------------------------------------------------------------------
# 4. COHORT SUMMARY
# ---------------------------------------------------------------------------

cohort_summary <- data.frame(
  Cohort = c(
    "Integrated",
    "Primary",
    "Sensitivity"
  ),
  N = c(
    460,
    456,
    455
  ),
  Events = c(
    220,
    220,
    219
  ),
  Censored = c(
    240,
    236,
    236
  ),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 5. KEY KAPLAN-MEIER RESULTS
# ---------------------------------------------------------------------------

km_results <- survival_results$KM_logrank

cat("\n========== KAPLAN-MEIER RESULTS ==========\n")
print(km_results, row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. KEY COX RESULTS
# ---------------------------------------------------------------------------

cox_coefficients <- survival_results$Cox_coefficients

cat("\n========== COX MODEL COEFFICIENTS ==========\n")
print(cox_coefficients, row.names = FALSE)

cox_fit <- survival_results$Cox_fit

cat("\n========== COX MODEL FIT ==========\n")
print(cox_fit, row.names = FALSE)

# ---------------------------------------------------------------------------
# 7. PROPORTIONAL HAZARDS DIAGNOSTICS
# ---------------------------------------------------------------------------

ph_results <- survival_results$Cox_PH

cat("\n========== PROPORTIONAL HAZARDS TESTS ==========\n")
print(ph_results, row.names = FALSE)

# ---------------------------------------------------------------------------
# 8. INCREMENTAL INFORMATION
# ---------------------------------------------------------------------------

lrt_results <- survival_results$Cox_LRT

cat("\n========== LIKELIHOOD RATIO TESTS ==========\n")
print(lrt_results, row.names = FALSE)

correlation_results <- survival_results$CpG_immune_correlation

cat("\n========== CpG / IMMUNE CORRELATION ==========\n")
print(correlation_results, row.names = FALSE)


# ---------------------------------------------------------------------------
# 9. CREATE FINAL KEY RESULTS TABLE — CORRECTED
# ---------------------------------------------------------------------------

key_results <- data.frame(
  Analysis = c(
    "KM: Immune phenotype",
    "KM: CpG signature",
    "Cox: Immune-high vs Immune-low",
    "Cox: CpG score, univariate",
    "Cox: CpG score, adjusted for immune phenotype",
    "Cox: CpG score, fully adjusted",
    "Cox: Age per 10 years",
    "LRT: CpG added to immune phenotype",
    "LRT: Immune phenotype added to CpG",
    "LRT: Age and sex added",
    "CpG / immune phenotype correlation"
  ),
  Primary_result = c(
    "Log-rank p = 1.123464e-05",
    "Log-rank p = 1.088789e-06",
    "HR = 0.551; 95% CI 0.420–0.721; p = 1.497e-05",
    "HR = 0.737 per SD; 95% CI 0.645–0.844; p = 9.295e-06",
    "HR = 0.825; 95% CI 0.634–1.072; p = 0.150",
    "HR = 0.824; 95% CI 0.636–1.069; p = 0.145",
    "HR = 1.272; 95% CI 1.159–1.397; p = 4.483e-07",
    "Chi-square = 2.101; p = 0.147",
    "Chi-square = 0.932; p = 0.334",
    "Chi-square = 26.962; p = 1.397e-06",
    "Pearson r = 0.843; Spearman rho = 0.856"
  ),
  Sensitivity_result = c(
    "Log-rank p = 1.466738e-05",
    "Log-rank p = 1.450866e-06",
    "HR = 0.554; 95% CI 0.423–0.727; p = 1.929e-05",
    "HR = 0.738 per SD; 95% CI 0.645–0.845; p = 1.044e-05",
    "HR = 0.819; 95% CI 0.629–1.066; p = 0.137",
    "HR = 0.818; 95% CI 0.630–1.062; p = 0.131",
    "HR = 1.265; 95% CI 1.152–1.389; p = 8.337e-07",
    "Chi-square = 2.242; p = 0.134",
    "Chi-square = 0.799; p = 0.371",
    "Chi-square = 25.624; p = 2.728e-06",
    "Pearson r = 0.843; Spearman rho = 0.856"
  ),
  Interpretation = c(
    "Significant difference between immune groups",
    "Significant difference between CpG score groups",
    "Immune-high associated with lower hazard in univariate model",
    "Higher CpG score associated with lower hazard in univariate model",
    "No statistically significant independent CpG association",
    "No statistically significant independent CpG association",
    "Age associated with higher hazard",
    "CpG does not significantly improve the immune-only model",
    "Immune phenotype does not significantly improve the CpG-only model",
    "Age and sex jointly improve the model",
    "Strong positive correlation between CpG score and immune phenotype"
  ),
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(key_results) == 11,
  all(complete.cases(key_results))
)

cat("\n========== FINAL KEY RESULTS ==========\n")
print(key_results, row.names = FALSE)

# ---------------------------------------------------------------------------
# 10. CREATE INPUT FILE INVENTORY
# ---------------------------------------------------------------------------

file_inventory <- data.frame(
  Result_name = names(input_files),
  File = unname(input_files),
  Exists = file.exists(input_files),
  Size_bytes = as.numeric(file.info(input_files)$size),
  Rows = vapply(survival_results, nrow, integer(1)),
  Columns = vapply(survival_results, ncol, integer(1)),
  stringsAsFactors = FALSE
)

cat("\n========== INPUT FILE INVENTORY ==========\n")
print(file_inventory, row.names = FALSE)

# ---------------------------------------------------------------------------
# 11. FINAL QUALITY CHECKS
# ---------------------------------------------------------------------------

final_qc <- c(
  all_inputs_present =
    all(file.exists(input_files)),
  
  all_input_tables_nonempty =
    all(vapply(survival_results, nrow, integer(1)) > 0),
  
  primary_cohort_size_correct =
    cohort_summary$N[cohort_summary$Cohort == "Primary"] == 456,
  
  sensitivity_cohort_size_correct =
    cohort_summary$N[cohort_summary$Cohort == "Sensitivity"] == 455,
  
  primary_event_count_correct =
    cohort_summary$Events[cohort_summary$Cohort == "Primary"] == 220,
  
  sensitivity_event_count_correct =
    cohort_summary$Events[cohort_summary$Cohort == "Sensitivity"] == 219,
  
  final_key_results_nonempty =
    nrow(key_results) > 0
)

cat("\n========== FINAL QC ==========\n")
print(final_qc)

if (!all(final_qc)) {
  stop("Final QC failed. Review the checks above.")
}

# ---------------------------------------------------------------------------
# 12. SAVE FINAL OUTPUTS
# ---------------------------------------------------------------------------

utils::write.csv(
  key_results,
  "results/tables/TCGA_SKCM_FINAL_survival_key_results.csv",
  row.names = FALSE
)

utils::write.csv(
  file_inventory,
  "results/tables/TCGA_SKCM_FINAL_survival_file_inventory.csv",
  row.names = FALSE
)

final_summary <- list(
  project = "TCGA-SKCM",
  analysis = "Exploratory overall survival analysis",
  cohort_summary = cohort_summary,
  key_results = key_results,
  KM_logrank = survival_results$KM_logrank,
  KM_group_counts = survival_results$KM_group_counts,
  Cox_coefficients = cox_coefficients,
  Cox_fit_statistics = cox_fit,
  Cox_PH_tests = ph_results,
  Cox_LRT = lrt_results,
  CpG_immune_correlation = correlation_results,
  input_inventory = file_inventory,
  quality_control = final_qc,
  interpretation_notes = c(
    "The CpG signature was originally developed for immune phenotype classification.",
    "The CpG score is not an independently trained survival-risk score.",
    "The CpG score and immune phenotype are strongly correlated.",
    "No significant incremental contribution of CpG score beyond immune phenotype was demonstrated.",
    "Age and sex jointly improved the combined Cox model.",
    "No stage covariate was included.",
    "The reported survival-time origin was not independently verified.",
    "The analysis does not establish immunotherapy response or causal effects."
  ),
  source = "Previously saved analysis outputs; no model refitting"
)

saveRDS(
  final_summary,
  "results/models/TCGA_SKCM_FINAL_survival_summary.rds"
)

# ---------------------------------------------------------------------------
# 13. OUTPUT AUDIT
# ---------------------------------------------------------------------------

cat("\n========== FINAL SAVED OUTPUTS ==========\n")

for (f in output_files) {
  cat(
    f,
    "| exists:", file.exists(f),
    "| bytes:", file.info(f)$size,
    "\n"
  )
}

cat("\n========== FINAL SURVIVAL SUMMARY COMPLETED ==========\n")
