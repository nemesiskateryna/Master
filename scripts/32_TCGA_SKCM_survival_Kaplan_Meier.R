###############################################################################
# 30K_TCGA_SKCM_survival_Kaplan_Meier.R
#
# Kaplan-Meier analysis:
# 1. Immune-high vs Immune-low
# 2. High vs Low CpG signature score
#
# Primary and sensitivity cohorts.
# Royal colour palette.
# No overwrite of existing results.
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
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

# ---------------------------------------------------------------------------
# 2. OUTPUT DIRECTORIES AND OVERWRITE PROTECTION
# ---------------------------------------------------------------------------

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/models", recursive = TRUE, showWarnings = FALSE)

output_files <- c(
  "results/tables/TCGA_SKCM_KM_logrank_results.csv",
  "results/tables/TCGA_SKCM_KM_group_counts.csv",
  "results/models/TCGA_SKCM_KM_analysis.rds",
  "results/figures/TCGA_SKCM_KM_immune_primary.png",
  "results/figures/TCGA_SKCM_KM_CpG_primary.png",
  "results/figures/TCGA_SKCM_KM_immune_sensitivity.png",
  "results/figures/TCGA_SKCM_KM_CpG_sensitivity.png"
)

existing_files <- output_files[file.exists(output_files)]

if (length(existing_files) > 0) {
  stop(
    "The following output files already exist. ",
    "No files were overwritten:\n",
    paste(existing_files, collapse = "\n")
  )
}

# ---------------------------------------------------------------------------
# 3. LOAD SAVED COHORTS
# ---------------------------------------------------------------------------

primary <- readRDS(
  "results/models/TCGA_SKCM_survival_primary_cohort.rds"
)$cohort

sensitivity <- readRDS(
  "results/models/TCGA_SKCM_survival_sensitivity_cohort.rds"
)$cohort

# ---------------------------------------------------------------------------
# 4. DEFINE ROYAL COLOUR PALETTE
# ---------------------------------------------------------------------------

royal_palette <- c(
  "Immune-high" = "#D4AF37",
  "Immune-low"  = "#4169E1",
  "High"        = "#D4AF37",
  "Low"         = "#4169E1"
)

# ---------------------------------------------------------------------------
# 5. DEFINE CpG SCORE GROUPS
# ---------------------------------------------------------------------------

cpg_cutoff <- median(
  primary$CpG_signature_score,
  na.rm = TRUE
)

primary$CpG_score_group <- factor(
  ifelse(
    primary$CpG_signature_score >= cpg_cutoff,
    "High",
    "Low"
  ),
  levels = c("Low", "High")
)

sensitivity$CpG_score_group <- factor(
  ifelse(
    sensitivity$CpG_signature_score >= cpg_cutoff,
    "High",
    "Low"
  ),
  levels = c("Low", "High")
)

cat("\n========== CpG SCORE CUTOFF ==========\n")
print(cpg_cutoff)

cat("\nPrimary CpG groups:\n")
print(table(primary$CpG_score_group))

cat("\nSensitivity CpG groups:\n")
print(table(sensitivity$CpG_score_group))

# ---------------------------------------------------------------------------
# 6. KAPLAN-MEIER AND LOG-RANK FUNCTION
# ---------------------------------------------------------------------------

run_km <- function(data, group_variable, cohort_name, analysis_name) {
  
  group <- data[[group_variable]]
  
  formula <- survival::Surv(OS_time, OS_event) ~ group
  
  fit <- survival::survfit(
    formula,
    data = data,
    conf.type = "log"
  )
  
  logrank <- survival::survdiff(
    formula,
    data = data
  )
  
  p_value <- stats::pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1,
    lower.tail = FALSE
  )
  
  km_summary <- summary(fit)
  
  plot_data <- data.frame(
    time = km_summary$time,
    survival = km_summary$surv,
    lower = km_summary$lower,
    upper = km_summary$upper,
    strata = as.character(km_summary$strata)
  )
  
  plot_data$Group <- sub(
    "^group=",
    "",
    plot_data$strata
  )
  
  plot_data$Group <- factor(
    plot_data$Group,
    levels = levels(group)
  )
  
  plot_data <- plot_data[
    order(plot_data$Group, plot_data$time),
  ]
  
  result <- data.frame(
    Cohort = cohort_name,
    Analysis = analysis_name,
    N = nrow(data),
    Events = sum(data$OS_event == 1),
    Log_rank_chisq = unname(logrank$chisq),
    Log_rank_df = length(logrank$n) - 1,
    Log_rank_p = p_value
  )
  
  list(
    fit = fit,
    plot_data = plot_data,
    logrank = logrank,
    result = result
  )
}

# ---------------------------------------------------------------------------
# 7. RUN FOUR KM ANALYSES
# ---------------------------------------------------------------------------

km_primary_immune <- run_km(
  primary,
  "Immune_group",
  "Primary",
  "Immune phenotype"
)

km_primary_cpg <- run_km(
  primary,
  "CpG_score_group",
  "Primary",
  "CpG signature score"
)

km_sensitivity_immune <- run_km(
  sensitivity,
  "Immune_group",
  "Sensitivity",
  "Immune phenotype"
)

km_sensitivity_cpg <- run_km(
  sensitivity,
  "CpG_score_group",
  "Sensitivity",
  "CpG signature score"
)

km_results <- rbind(
  km_primary_immune$result,
  km_primary_cpg$result,
  km_sensitivity_immune$result,
  km_sensitivity_cpg$result
)

cat("\n========== LOG-RANK RESULTS ==========\n")
print(km_results, row.names = FALSE)

# ---------------------------------------------------------------------------
# 8. GROUP COUNTS AND EVENTS
# ---------------------------------------------------------------------------

make_group_counts <- function(data, group_variable, cohort_name,
                              analysis_name) {
  
  counts <- as.data.frame(
    table(
      Group = data[[group_variable]],
      Event = data$OS_event
    )
  )
  
  names(counts)[3] <- "N"
  
  counts$Cohort <- cohort_name
  counts$Analysis <- analysis_name
  
  counts[, c("Cohort", "Analysis", "Group", "Event", "N")]
}

group_counts <- rbind(
  make_group_counts(
    primary, "Immune_group",
    "Primary", "Immune phenotype"
  ),
  make_group_counts(
    primary, "CpG_score_group",
    "Primary", "CpG signature score"
  ),
  make_group_counts(
    sensitivity, "Immune_group",
    "Sensitivity", "Immune phenotype"
  ),
  make_group_counts(
    sensitivity, "CpG_score_group",
    "Sensitivity", "CpG signature score"
  )
)

cat("\n========== GROUP COUNTS ==========\n")
print(group_counts, row.names = FALSE)

# ---------------------------------------------------------------------------
# 9. PLOT FUNCTION
# ---------------------------------------------------------------------------

make_km_plot <- function(km_object, title, output_file) {
  
  p_value <- km_object$result$Log_rank_p
  
  p_label <- if (p_value < 0.001) {
    "Log-rank p < 0.001"
  } else {
    paste0("Log-rank p = ", signif(p_value, 3))
  }
  
  p <- ggplot2::ggplot(
    km_object$plot_data,
    ggplot2::aes(
      x = time,
      y = survival,
      colour = Group
    )
  ) +
    ggplot2::geom_step(linewidth = 1.05) +
    ggplot2::geom_step(
      ggplot2::aes(y = lower),
      linewidth = 0.35,
      linetype = "dashed",
      alpha = 0.65
    ) +
    ggplot2::geom_step(
      ggplot2::aes(y = upper),
      linewidth = 0.35,
      linetype = "dashed",
      alpha = 0.65
    ) +
    ggplot2::scale_colour_manual(
      values = royal_palette,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = p_label,
      x = "Reported survival time (days)",
      y = "Overall survival probability",
      colour = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        colour = "#191970"
      ),
      plot.subtitle = ggplot2::element_text(
        colour = "#6A0DAD"
      ),
      legend.position = "bottom"
    )
  
  ggplot2::ggsave(
    filename = output_file,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  p
}

# ---------------------------------------------------------------------------
# 10. SAVE FOUR FIGURES
# ---------------------------------------------------------------------------

p1 <- make_km_plot(
  km_primary_immune,
  "Overall survival by immune phenotype — Primary cohort",
  "results/figures/TCGA_SKCM_KM_immune_primary.png"
)

p2 <- make_km_plot(
  km_primary_cpg,
  "Overall survival by CpG signature score — Primary cohort",
  "results/figures/TCGA_SKCM_KM_CpG_primary.png"
)

p3 <- make_km_plot(
  km_sensitivity_immune,
  "Overall survival by immune phenotype — Sensitivity cohort",
  "results/figures/TCGA_SKCM_KM_immune_sensitivity.png"
)

p4 <- make_km_plot(
  km_sensitivity_cpg,
  "Overall survival by CpG signature score — Sensitivity cohort",
  "results/figures/TCGA_SKCM_KM_CpG_sensitivity.png"
)

# ---------------------------------------------------------------------------
# 11. SAVE TABLES AND MODEL OBJECTS
# ---------------------------------------------------------------------------

utils::write.csv(
  km_results,
  "results/tables/TCGA_SKCM_KM_logrank_results.csv",
  row.names = FALSE
)

utils::write.csv(
  group_counts,
  "results/tables/TCGA_SKCM_KM_group_counts.csv",
  row.names = FALSE
)

km_analysis <- list(
  cutoff_primary_median = cpg_cutoff,
  royal_palette = royal_palette,
  primary = list(
    immune = km_primary_immune,
    cpg = km_primary_cpg
  ),
  sensitivity = list(
    immune = km_sensitivity_immune,
    cpg = km_sensitivity_cpg
  ),
  logrank_results = km_results,
  group_counts = group_counts
)

saveRDS(
  km_analysis,
  "results/models/TCGA_SKCM_KM_analysis.rds"
)

# ---------------------------------------------------------------------------
# 12. FINAL OUTPUT AUDIT
# ---------------------------------------------------------------------------

cat("\n========== SAVED OUTPUTS ==========\n")

for (f in output_files) {
  cat(
    f,
    " | exists:", file.exists(f),
    " | bytes:", if (file.exists(f)) file.info(f)$size else NA,
    "\n"
  )
}

cat("\n========== KAPLAN-MEIER ANALYSIS COMPLETED ==========\n")
