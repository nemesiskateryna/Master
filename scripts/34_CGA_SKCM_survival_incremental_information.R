###############################################################################
# 30N_TCGA_SKCM_survival_incremental_information.R
#
# Read-only analysis of saved Cox models:
# - CpG score / immune phenotype correlation
# - Likelihood ratio tests for nested models
# - Incremental contribution of covariates
#
# No model refitting and no modification of previous results.
###############################################################################

setwd("/Users/katie/master/Master")

# ---------------------------------------------------------------------------
# 1. LOAD SAVED RESULTS
# ---------------------------------------------------------------------------

cox_analysis <- readRDS(
  "results/models/TCGA_SKCM_Cox_survival_analysis.rds"
)

primary <- readRDS(
  "results/models/TCGA_SKCM_survival_primary_cohort.rds"
)$cohort

sensitivity <- readRDS(
  "results/models/TCGA_SKCM_survival_sensitivity_cohort.rds"
)$cohort

# ---------------------------------------------------------------------------
# 2. OUTPUT PROTECTION
# ---------------------------------------------------------------------------

output_files <- c(
  "results/tables/TCGA_SKCM_Cox_incremental_LRT.csv",
  "results/tables/TCGA_SKCM_CpG_immune_correlation.csv",
  "results/models/TCGA_SKCM_survival_incremental_information.rds"
)

existing <- output_files[file.exists(output_files)]

if (length(existing) > 0) {
  stop(
    "Files already exist. Nothing was overwritten:\n",
    paste(existing, collapse = "\n")
  )
}

# ---------------------------------------------------------------------------
# 3. CORRELATION BETWEEN CpG SCORE AND IMMUNE PHENOTYPE
# ---------------------------------------------------------------------------

calculate_correlation <- function(data, cohort_name) {
  
  immune_binary <- ifelse(
    data$Immune_group == "Immune-high",
    1,
    0
  )
  
  pearson <- stats::cor.test(
    data$CpG_signature_score,
    immune_binary,
    method = "pearson"
  )
  
  spearman <- stats::cor.test(
    data$CpG_signature_score,
    immune_binary,
    method = "spearman",
    exact = FALSE
  )
  
  data.frame(
    Cohort = cohort_name,
    N = nrow(data),
    Pearson_r = unname(pearson$estimate),
    Pearson_p = pearson$p.value,
    Spearman_rho = unname(spearman$estimate),
    Spearman_p = spearman$p.value
  )
}

correlation_results <- rbind(
  calculate_correlation(primary, "Primary"),
  calculate_correlation(sensitivity, "Sensitivity")
)

cat("\n========== CpG / IMMUNE CORRELATION ==========\n")
print(correlation_results, row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. LIKELIHOOD RATIO TESTS
# ---------------------------------------------------------------------------

run_lrt <- function(models, cohort_name) {
  
  comparisons <- list(
    CpG_added_to_immune = c(
      "Model_1_Immune",
      "Model_3_CpG_Immune"
    ),
    
    Immune_added_to_CpG = c(
      "Model_2_CpG",
      "Model_3_CpG_Immune"
    ),
    
    Age_sex_added = c(
      "Model_3_CpG_Immune",
      "Model_4_Adjusted"
    )
  )
  
  results <- lapply(names(comparisons), function(comparison_name) {
    
    model_names <- comparisons[[comparison_name]]
    
    comparison <- stats::anova(
      models[[model_names[1]]],
      models[[model_names[2]]],
      test = "LRT"
    )
    
    comparison <- as.data.frame(comparison)
    
    data.frame(
      Cohort = cohort_name,
      Comparison = comparison_name,
      Reduced_model = model_names[1],
      Full_model = model_names[2],
      Chi_square = comparison$Chisq[2],
      Degrees_freedom = comparison$Df[2],
      P_value = comparison$`Pr(>|Chi|)`[2]
    )
  })
  
  do.call(rbind, results)
}

lrt_results <- rbind(
  run_lrt(
    cox_analysis$primary$models,
    "Primary"
  ),
  run_lrt(
    cox_analysis$sensitivity$models,
    "Sensitivity"
  )
)

cat("\n========== LIKELIHOOD RATIO TESTS ==========\n")
print(lrt_results, row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. SAVE RESULTS
# ---------------------------------------------------------------------------

utils::write.csv(
  correlation_results,
  "results/tables/TCGA_SKCM_CpG_immune_correlation.csv",
  row.names = FALSE
)

utils::write.csv(
  lrt_results,
  "results/tables/TCGA_SKCM_Cox_incremental_LRT.csv",
  row.names = FALSE
)

incremental_analysis <- list(
  correlation = correlation_results,
  LRT = lrt_results,
  source = "Previously saved Cox models; no refitting"
)

saveRDS(
  incremental_analysis,
  "results/models/TCGA_SKCM_survival_incremental_information.rds"
)

# ---------------------------------------------------------------------------
# 6. OUTPUT AUDIT
# ---------------------------------------------------------------------------

cat("\n========== SAVED OUTPUTS ==========\n")

for (f in output_files) {
  cat(
    f,
    "| exists:", file.exists(f),
    "| bytes:", file.info(f)$size,
    "\n"
  )
}

cat("\nINCREMENTAL INFORMATION ANALYSIS COMPLETED\n")
