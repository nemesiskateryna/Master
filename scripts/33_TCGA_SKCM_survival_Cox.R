###############################################################################
# 30M_TCGA_SKCM_survival_Cox_models.R
#
# Cox proportional hazards models:
# Primary and sensitivity cohorts.
#
# Model 1: Immune phenotype
# Model 2: CpG signature score
# Model 3: CpG score + immune phenotype
# Model 4: CpG score + immune phenotype + age + sex
#
# PH assumption: Schoenfeld residuals
# Royal colour palette
# No overwriting of existing results
###############################################################################

setwd("/Users/katie/master/Master")

# ---------------------------------------------------------------------------
# 1. PACKAGES
# ---------------------------------------------------------------------------

required_packages <- c("survival", "ggplot2")

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

# ---------------------------------------------------------------------------
# 2. OUTPUT FILES AND OVERWRITE PROTECTION
# ---------------------------------------------------------------------------

output_files <- c(
  "results/tables/TCGA_SKCM_Cox_model_coefficients.csv",
  "results/tables/TCGA_SKCM_Cox_model_fit_statistics.csv",
  "results/tables/TCGA_SKCM_Cox_PH_assumption_tests.csv",
  "results/models/TCGA_SKCM_Cox_survival_analysis.rds",
  "results/figures/TCGA_SKCM_Cox_forest_primary.png",
  "results/figures/TCGA_SKCM_Cox_forest_sensitivity.png",
  "results/figures/TCGA_SKCM_Cox_PH_diagnostics_primary.png",
  "results/figures/TCGA_SKCM_Cox_PH_diagnostics_sensitivity.png"
)

existing_files <- output_files[file.exists(output_files)]

if (length(existing_files) > 0) {
  stop(
    "Output files already exist. Nothing was overwritten:\n",
    paste(existing_files, collapse = "\n")
  )
}

# ---------------------------------------------------------------------------
# 3. LOAD COHORTS
# ---------------------------------------------------------------------------

primary <- readRDS(
  "results/models/TCGA_SKCM_survival_primary_cohort.rds"
)$cohort

sensitivity <- readRDS(
  "results/models/TCGA_SKCM_survival_sensitivity_cohort.rds"
)$cohort

# ---------------------------------------------------------------------------
# 4. ROYAL PALETTE
# ---------------------------------------------------------------------------

royal_blue <- "#4169E1"
royal_gold <- "#D4AF37"
deep_purple <- "#6A0DAD"
navy <- "#191970"

# ---------------------------------------------------------------------------
# 5. PREPARE COVARIATES
# ---------------------------------------------------------------------------

# CpG score is standardized using PRIMARY cohort SD.
# The same SD is applied to the sensitivity cohort.

cpg_sd_primary <- sd(
  primary$CpG_signature_score,
  na.rm = TRUE
)

if (!is.finite(cpg_sd_primary) || cpg_sd_primary <= 0) {
  stop("Invalid primary-cohort CpG score SD.")
}

prepare_cohort <- function(data, cpg_sd) {
  
  data$Immune_group <- factor(
    data$Immune_group,
    levels = c("Immune-low", "Immune-high")
  )
  
  data$Sex <- factor(
    tolower(data$Sex),
    levels = c("female", "male")
  )
  
  data$CpG_score_SD <- data$CpG_signature_score / cpg_sd
  
  data$Age_per_10y <- data$Age / 10
  
  data
}

primary <- prepare_cohort(primary, cpg_sd_primary)
sensitivity <- prepare_cohort(sensitivity, cpg_sd_primary)

# Verify all model variables.

model_variables <- c(
  "OS_time",
  "OS_event",
  "Immune_group",
  "CpG_score_SD",
  "Age_per_10y",
  "Sex"
)

stopifnot(
  all(model_variables %in% names(primary)),
  all(model_variables %in% names(sensitivity)),
  all(complete.cases(primary[, model_variables])),
  all(complete.cases(sensitivity[, model_variables])),
  all(primary$OS_time >= 0),
  all(sensitivity$OS_time > 0),
  all(primary$OS_event %in% c(0, 1)),
  all(sensitivity$OS_event %in% c(0, 1))
)

cat("\n========== MODEL INPUT SUMMARY ==========\n")

cat("Primary N:", nrow(primary), "\n")
cat("Primary events:", sum(primary$OS_event), "\n")

cat("Sensitivity N:", nrow(sensitivity), "\n")
cat("Sensitivity events:", sum(sensitivity$OS_event), "\n")

cat("Primary CpG SD:", cpg_sd_primary, "\n")

cat("\nImmune group reference: Immune-low\n")
cat("Sex reference: female\n")
cat("CpG HR unit: one primary-cohort SD\n")
cat("Age HR unit: 10 years\n")

# ---------------------------------------------------------------------------
# 6. DEFINE COX MODELS
# ---------------------------------------------------------------------------

model_formulas <- list(
  
  Model_1_Immune = survival::Surv(
    OS_time, OS_event
  ) ~ Immune_group,
  
  Model_2_CpG = survival::Surv(
    OS_time, OS_event
  ) ~ CpG_score_SD,
  
  Model_3_CpG_Immune = survival::Surv(
    OS_time, OS_event
  ) ~ CpG_score_SD + Immune_group,
  
  Model_4_Adjusted = survival::Surv(
    OS_time, OS_event
  ) ~ CpG_score_SD + Immune_group + Age_per_10y + Sex
)

# ---------------------------------------------------------------------------
# 7. FIT MODELS
# ---------------------------------------------------------------------------

fit_cox_models <- function(data, cohort_name) {
  
  fits <- lapply(
    model_formulas,
    function(formula) {
      survival::coxph(
        formula,
        data = data,
        ties = "efron",
        x = TRUE,
        y = TRUE,
        model = TRUE
      )
    }
  )
  
  names(fits) <- names(model_formulas)
  
  fits
}

primary_models <- fit_cox_models(primary, "Primary")

sensitivity_models <- fit_cox_models(
  sensitivity,
  "Sensitivity"
)

cat("\n========== COX MODELS FITTED ==========\n")
cat("Primary models:", length(primary_models), "\n")
cat("Sensitivity models:", length(sensitivity_models), "\n")

# ---------------------------------------------------------------------------
# 8. EXTRACT COEFFICIENTS, HR, CI AND P-VALUES
# ---------------------------------------------------------------------------

extract_cox_coefficients <- function(models, cohort_name) {
  
  output <- list()
  
  for (model_name in names(models)) {
    
    fit <- models[[model_name]]
    s <- summary(fit)
    
    coefficient_table <- as.data.frame(
      s$coefficients
    )
    
    confidence_intervals <- as.data.frame(
      s$conf.int
    )
    
    result <- data.frame(
      Cohort = cohort_name,
      Model = model_name,
      Term = rownames(coefficient_table),
      Beta = coefficient_table$coef,
      HR = coefficient_table$`exp(coef)`,
      CI_lower = confidence_intervals$`lower .95`,
      CI_upper = confidence_intervals$`upper .95`,
      SE = coefficient_table$`se(coef)`,
      Z = coefficient_table$z,
      P_value = coefficient_table$`Pr(>|z|)`,
      row.names = NULL
    )
    
    output[[model_name]] <- result
  }
  
  do.call(rbind, output)
}

cox_coefficients <- rbind(
  extract_cox_coefficients(primary_models, "Primary"),
  extract_cox_coefficients(sensitivity_models, "Sensitivity")
)

cat("\n========== COX COEFFICIENT RESULTS ==========\n")
print(cox_coefficients, row.names = FALSE)

# ---------------------------------------------------------------------------
# 9. MODEL FIT STATISTICS
# ---------------------------------------------------------------------------

extract_fit_statistics <- function(models, data, cohort_name) {
  
  output <- lapply(names(models), function(model_name) {
    
    fit <- models[[model_name]]
    s <- summary(fit)
    
    data.frame(
      Cohort = cohort_name,
      Model = model_name,
      N = nrow(data),
      Events = sum(data$OS_event == 1),
      Parameters = length(stats::coef(fit)),
      AIC = stats::AIC(fit),
      Concordance = unname(s$concordance[1]),
      Concordance_SE = unname(s$concordance[2])
    )
  })
  
  do.call(rbind, output)
}

cox_fit_statistics <- rbind(
  extract_fit_statistics(
    primary_models, primary, "Primary"
  ),
  extract_fit_statistics(
    sensitivity_models, sensitivity, "Sensitivity"
  )
)

cat("\n========== MODEL FIT STATISTICS ==========\n")
print(cox_fit_statistics, row.names = FALSE)

# ---------------------------------------------------------------------------
# 10. SCHOENFELD PROPORTIONAL HAZARDS TESTS
# ---------------------------------------------------------------------------

extract_ph_tests <- function(models, cohort_name) {
  
  output <- list()
  zph_objects <- list()
  
  for (model_name in names(models)) {
    
    zph <- survival::cox.zph(
      models[[model_name]],
      transform = "km"
    )
    
    zph_objects[[model_name]] <- zph
    
    ph_table <- as.data.frame(zph$table)
    ph_table$Term <- rownames(ph_table)
    
    ph_table$Cohort <- cohort_name
    ph_table$Model <- model_name
    
    output[[model_name]] <- ph_table
  }
  
  list(
    table = do.call(rbind, output),
    zph = zph_objects
  )
}

primary_ph <- extract_ph_tests(
  primary_models,
  "Primary"
)

sensitivity_ph <- extract_ph_tests(
  sensitivity_models,
  "Sensitivity"
)

ph_results <- rbind(
  primary_ph$table,
  sensitivity_ph$table
)

rownames(ph_results) <- NULL

cat("\n========== SCHOENFELD PH TESTS ==========\n")
print(ph_results, row.names = FALSE)

# ---------------------------------------------------------------------------
# 11. FOREST PLOT — FULLY ADJUSTED MODEL
# ---------------------------------------------------------------------------

make_forest_plot <- function(coefficient_table, cohort_name, output_file) {
  
  forest_data <- coefficient_table[
    coefficient_table$Cohort == cohort_name &
      coefficient_table$Model == "Model_4_Adjusted",
    ,
    drop = FALSE
  ]
  
  term_labels <- c(
    "CpG_score_SD" = "CpG score (per 1 SD)",
    "Immune_groupImmune-high" = "Immune-high vs Immune-low",
    "Age_per_10y" = "Age (per 10 years)",
    "Sexmale" = "Male vs female"
  )
  
  forest_data$Label <- unname(
    term_labels[forest_data$Term]
  )
  
  if (anyNA(forest_data$Label)) {
    stop("Could not map all Cox terms to forest plot labels.")
  }
  
  forest_data$Label <- factor(
    forest_data$Label,
    levels = rev(forest_data$Label)
  )
  
  p <- ggplot2::ggplot(
    forest_data,
    ggplot2::aes(
      x = HR,
      y = Label
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = deep_purple,
      linewidth = 0.7
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = CI_lower,
        xmax = CI_upper
      ),
      orientation = "y",
      width = 0.18,
      colour = royal_blue,
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      size = 3.2,
      colour = royal_gold
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = paste(
        "Adjusted Cox model —",
        cohort_name,
        "cohort"
      ),
      subtitle = "Hazard ratios with 95% confidence intervals",
      x = "Hazard ratio (log scale)",
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        colour = navy
      ),
      plot.subtitle = ggplot2::element_text(
        colour = deep_purple
      )
    )
  
  ggplot2::ggsave(
    output_file,
    p,
    width = 8,
    height = 5.5,
    dpi = 300
  )
  
  p
}

forest_primary <- make_forest_plot(
  cox_coefficients,
  "Primary",
  "results/figures/TCGA_SKCM_Cox_forest_primary.png"
)

forest_sensitivity <- make_forest_plot(
  cox_coefficients,
  "Sensitivity",
  "results/figures/TCGA_SKCM_Cox_forest_sensitivity.png"
)

# ---------------------------------------------------------------------------
# 12. PH DIAGNOSTIC PLOTS — FULLY ADJUSTED MODEL
# ---------------------------------------------------------------------------

save_ph_plot <- function(zph_object, title, output_file) {
  
  grDevices::png(
    filename = output_file,
    width = 2400,
    height = 1800,
    res = 300
  )
  
  old_par <- graphics::par(
    mfrow = c(2, 2),
    mar = c(4, 4, 3, 1)
  )
  
  graphics::plot(
    zph_object,
    resid = TRUE,
    se = TRUE,
    col = royal_blue,
    lwd = 2,
    main = title
  )
  
  graphics::abline(
    h = 0,
    col = royal_gold,
    lty = 2
  )
  
  graphics::par(old_par)
  grDevices::dev.off()
}

save_ph_plot(
  primary_ph$zph$Model_4_Adjusted,
  "Schoenfeld diagnostics — Primary cohort",
  "results/figures/TCGA_SKCM_Cox_PH_diagnostics_primary.png"
)

save_ph_plot(
  sensitivity_ph$zph$Model_4_Adjusted,
  "Schoenfeld diagnostics — Sensitivity cohort",
  "results/figures/TCGA_SKCM_Cox_PH_diagnostics_sensitivity.png"
)

# ---------------------------------------------------------------------------
# 13. SAVE TABLES
# ---------------------------------------------------------------------------

utils::write.csv(
  cox_coefficients,
  "results/tables/TCGA_SKCM_Cox_model_coefficients.csv",
  row.names = FALSE
)

utils::write.csv(
  cox_fit_statistics,
  "results/tables/TCGA_SKCM_Cox_model_fit_statistics.csv",
  row.names = FALSE
)

utils::write.csv(
  ph_results,
  "results/tables/TCGA_SKCM_Cox_PH_assumption_tests.csv",
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# 14. SAVE MODEL OBJECTS
# ---------------------------------------------------------------------------

cox_analysis <- list(
  CpG_SD_primary = cpg_sd_primary,
  primary = list(
    models = primary_models,
    PH_tests = primary_ph
  ),
  sensitivity = list(
    models = sensitivity_models,
    PH_tests = sensitivity_ph
  ),
  coefficients = cox_coefficients,
  fit_statistics = cox_fit_statistics,
  PH_results = ph_results,
  royal_palette = c(
    royal_blue = royal_blue,
    royal_gold = royal_gold,
    deep_purple = deep_purple,
    navy = navy
  )
)

saveRDS(
  cox_analysis,
  "results/models/TCGA_SKCM_Cox_survival_analysis.rds"
)

# ---------------------------------------------------------------------------
# 15. FINAL OUTPUT AUDIT
# ---------------------------------------------------------------------------

cat("\n========== SAVED OUTPUT AUDIT ==========\n")

for (f in output_files) {
  cat(
    f,
    "| exists:", file.exists(f),
    "| bytes:",
    if (file.exists(f)) file.info(f)$size else NA,
    "\n"
  )
}

cat("\n========== COX ANALYSIS COMPLETED ==========\n")
