
###############################################################################
# 30Q_TCGA_SKCM_TAP_gene_expression_mapping.R
#
# Map HLA-F and TAP1 symbols to Ensembl IDs and extract expression.
# No model fitting. No modification of previous results.
###############################################################################

setwd("/Users/katie/master/Master")

# ---------------------------------------------------------------------------
# 1. CHECK ANNOTATION PACKAGE
# ---------------------------------------------------------------------------

if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  stop("AnnotationDbi is not installed.")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stop("org.Hs.eg.db is not installed.")
}

# ---------------------------------------------------------------------------
# 2. LOAD EXPRESSION AND SURVIVAL DATA
# ---------------------------------------------------------------------------

integrated <- readRDS(
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune.rds"
)

primary <- readRDS(
  "results/models/TCGA_SKCM_survival_primary_cohort.rds"
)$cohort

sensitivity <- readRDS(
  "results/models/TCGA_SKCM_survival_sensitivity_cohort.rds"
)$cohort

expression_matrix <- integrated$expression_matrix

# ---------------------------------------------------------------------------
# 3. REMOVE ENSEMBL VERSION SUFFIXES
# ---------------------------------------------------------------------------

ensembl_original <- rownames(expression_matrix)

ensembl_clean <- sub(
  "\\.[0-9]+$",
  "",
  ensembl_original
)

cat("\n========== ENSEMBL ID QC ==========\n")

cat("Total expression rows:", length(ensembl_original), "\n")
cat("Unique versionless IDs:", length(unique(ensembl_clean)), "\n")
cat("Duplicated versionless IDs:", sum(duplicated(ensembl_clean)), "\n")

if (anyDuplicated(ensembl_clean)) {
  stop(
    "Duplicated Ensembl IDs after removing version suffixes. ",
    "Resolve duplicates before extracting expression."
  )
}
# ---------------------------------------------------------------------------
# 4. SELECT ONLY ENSEMBL IDS PRESENT IN RNA-SEQ
# ---------------------------------------------------------------------------

gene_mapping_available <- gene_mapping[
  !is.na(gene_mapping$Expression_row),
  ,
  drop = FALSE
]

# Confirm exactly one RNA-seq row per target gene
mapping_counts <- table(gene_mapping_available$SYMBOL)

print(mapping_counts)

if (!all(target_genes %in% names(mapping_counts))) {
  stop("At least one target gene is missing from RNA-seq.")
}

if (any(mapping_counts[target_genes] != 1)) {
  stop("Expected exactly one RNA-seq row per target gene.")
}

gene_mapping_final <- gene_mapping_available[
  match(target_genes, gene_mapping_available$SYMBOL),
  ,
  drop = FALSE
]

cat("\n========== FINAL GENE MAPPING ==========\n")
print(
  gene_mapping_final[
    ,
    c("SYMBOL", "ENSEMBL", "Expression_ID", "Expression_row"),
    drop = FALSE
  ],
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# 5. EXTRACT EXPRESSION
# ---------------------------------------------------------------------------

selected_rows <- gene_mapping_final$Expression_row

gene_expression <- as.data.frame(
  t(
    as.matrix(
      expression_matrix[selected_rows, , drop = FALSE]
    )
  )
)

colnames(gene_expression) <- target_genes

gene_expression$Sample_ID <- rownames(gene_expression)

gene_expression <- gene_expression[
  ,
  c("Sample_ID", "HLA-F", "TAP1"),
  drop = FALSE
]

cat("\n========== EXTRACTED EXPRESSION ==========\n")
cat("Dimensions:", dim(gene_expression), "\n")
print(head(gene_expression))

# ---------------------------------------------------------------------------
# 6. EXPRESSION SUMMARY
# ---------------------------------------------------------------------------

cat("\n========== EXPRESSION SUMMARY ==========\n")

print(
  summary(
    gene_expression[, target_genes, drop = FALSE]
  )
)

cat("\n========== EXPRESSION CORRELATION ==========\n")

print(
  cor.test(
    gene_expression$`HLA-F`,
    gene_expression$TAP1,
    method = "spearman",
    exact = FALSE
  )
)

# ---------------------------------------------------------------------------
# 7. VERIFY SURVIVAL COHORT ALIGNMENT
# ---------------------------------------------------------------------------

primary_match <- primary$Sample_ID %in% gene_expression$Sample_ID
sensitivity_match <- sensitivity$Sample_ID %in% gene_expression$Sample_ID

alignment_summary <- data.frame(
  Cohort = c("Primary", "Sensitivity"),
  N = c(nrow(primary), nrow(sensitivity)),
  Matched = c(
    sum(primary_match),
    sum(sensitivity_match)
  ),
  Missing = c(
    sum(!primary_match),
    sum(!sensitivity_match)
  )
)

cat("\n========== SURVIVAL ALIGNMENT ==========\n")
print(alignment_summary, row.names = FALSE)

stopifnot(
  all(primary_match),
  all(sensitivity_match),
  !anyDuplicated(gene_expression$Sample_ID),
  all(is.finite(as.matrix(
    gene_expression[, target_genes, drop = FALSE]
  )))
)

# ---------------------------------------------------------------------------
# 8. SAVE NEW AUDIT OUTPUTS
# ---------------------------------------------------------------------------

output_files <- c(
  "results/tables/TCGA_SKCM_TAP_gene_expression_ID_mapping.csv",
  "results/tables/TCGA_SKCM_TAP_gene_expression_alignment.csv",
  "results/models/TCGA_SKCM_TAP_gene_expression_data.rds"
)

existing <- output_files[file.exists(output_files)]

if (length(existing) > 0) {
  stop(
    "Output files already exist. Nothing was overwritten:\n",
    paste(existing, collapse = "\n")
  )
}

utils::write.csv(
  gene_mapping,
  output_files[1],
  row.names = FALSE
)

utils::write.csv(
  alignment_summary,
  output_files[2],
  row.names = FALSE
)

saveRDS(
  list(
    expression = gene_expression,
    mapping = gene_mapping,
    alignment = alignment_summary,
    genes = target_genes,
    source = "TCGA-SKCM integrated RNA-seq expression matrix"
  ),
  output_files[3]
)

cat("\n========== SAVED OUTPUTS ==========\n")

for (f in output_files) {
  cat(
    f,
    "| exists:", file.exists(f),
    "| bytes:", file.info(f)$size,
    "\n"
  )
}

cat("\nGENE EXPRESSION MAPPING COMPLETED\n")


###############################################################################
# 30R_TCGA_SKCM_TAP_gene_survival_analysis.R
#
# Exploratory survival analysis of HLA-F and TAP1 expression.
# Primary and sensitivity cohorts.
# No modification of previous results.
###############################################################################

setwd("/Users/katie/master/Master")

# ---------------------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------------------

gene_data <- readRDS(
  "results/models/TCGA_SKCM_TAP_gene_expression_data.rds"
)$expression

primary <- readRDS(
  "results/models/TCGA_SKCM_survival_primary_cohort.rds"
)$cohort

sensitivity <- readRDS(
  "results/models/TCGA_SKCM_survival_sensitivity_cohort.rds"
)$cohort

target_genes <- c("HLA-F", "TAP1")

# ---------------------------------------------------------------------------
# 2. ALIGN EXPRESSION WITH SURVIVAL DATA
# ---------------------------------------------------------------------------

prepare_survival_data <- function(cohort, expression_data) {
  
  idx <- match(cohort$Sample_ID, expression_data$Sample_ID)
  
  if (anyNA(idx)) {
    stop("Some survival samples have no gene expression data.")
  }
  
  expression_subset <- expression_data[
    idx,
    c("HLA-F", "TAP1"),
    drop = FALSE
  ]
  
  result <- cbind(
    cohort,
    expression_subset
  )
  
  if (!identical(result$Sample_ID, cohort$Sample_ID)) {
    stop("Sample alignment failed.")
  }
  
  result
}

primary_data <- prepare_survival_data(
  primary,
  gene_data
)

sensitivity_data <- prepare_survival_data(
  sensitivity,
  gene_data
)


# ---------------------------------------------------------------------------
# 3. STANDARDIZE GENE EXPRESSION
# ---------------------------------------------------------------------------

gene_scaling <- data.frame(
  Gene = target_genes,
  Mean = vapply(
    target_genes,
    function(g) mean(primary_data[[g]]),
    numeric(1)
  ),
  SD = vapply(
    target_genes,
    function(g) sd(primary_data[[g]]),
    numeric(1)
  )
)

# Safe variable names for R formulas
gene_variables <- c(
  "HLA-F" = "HLA_F_SD",
  "TAP1" = "TAP1_SD"
)

for (gene in target_genes) {
  
  gene_mean <- gene_scaling$Mean[
    gene_scaling$Gene == gene
  ]
  
  gene_sd <- gene_scaling$SD[
    gene_scaling$Gene == gene
  ]
  
  safe_name <- gene_variables[[gene]]
  
  primary_data[[safe_name]] <-
    (primary_data[[gene]] - gene_mean) / gene_sd
  
  sensitivity_data[[safe_name]] <-
    (sensitivity_data[[gene]] - gene_mean) / gene_sd
}

cat("\n========== GENE SCALING ==========\n")
print(gene_scaling, row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. PREPARE COVARIATES
# ---------------------------------------------------------------------------

prepare_factors <- function(data) {
  
  data$Immune_group <- factor(
    data$Immune_group,
    levels = c("Immune-low", "Immune-high")
  )
  
  data$Sex <- factor(
    data$Sex,
    levels = c("female", "male")
  )
  
  data$Age_per_10y <- data$Age / 10
  
  data
}

primary_data <- prepare_factors(primary_data)
sensitivity_data <- prepare_factors(sensitivity_data)

# ---------------------------------------------------------------------------
# 5. FIT COX MODELS
# ---------------------------------------------------------------------------

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' is required.")
}

fit_gene_models <- function(data, cohort_name) {
  
  results <- list()
  models <- list()
  
  for (gene in target_genes) {
    
    expression_variable <- gene_variables[[gene]]
    
    formula_univariate <- stats::as.formula(
      paste(
        "survival::Surv(OS_time, OS_event) ~",
        expression_variable
      )
    )
    
    formula_adjusted <- stats::as.formula(
      paste(
        "survival::Surv(OS_time, OS_event) ~",
        expression_variable,
        "+ Immune_group + Age_per_10y + Sex"
      )
    )
    
    model_univariate <- survival::coxph(
      formula_univariate,
      data = data,
      ties = "efron",
      x = TRUE
    )
    
    model_adjusted <- survival::coxph(
      formula_adjusted,
      data = data,
      ties = "efron",
      x = TRUE
    )
    
    model_names <- c(
      paste0(cohort_name, "_", gene, "_Univariate"),
      paste0(cohort_name, "_", gene, "_Adjusted")
    )
    
    models[[model_names[1]]] <- model_univariate
    models[[model_names[2]]] <- model_adjusted
    
    extract_gene_result <- function(model, model_name) {
      
      coefficient_table <- summary(model)$coefficients
      confidence_table <- summary(model)$conf.int
      
      row_index <- match(
        expression_variable,
        rownames(coefficient_table)
      )
      
      if (is.na(row_index)) {
        stop("Gene coefficient not found in Cox model.")
      }
      
      data.frame(
        Cohort = cohort_name,
        Gene = gene,
        Model = model_name,
        N = model$n,
        Events = model$nevent,
        HR = confidence_table[row_index, "exp(coef)"],
        CI_lower = confidence_table[row_index, "lower .95"],
        CI_upper = confidence_table[row_index, "upper .95"],
        P_value = coefficient_table[row_index, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
    }
    
    results[[length(results) + 1]] <- extract_gene_result(
      model_univariate,
      "Univariate"
    )
    
    results[[length(results) + 1]] <- extract_gene_result(
      model_adjusted,
      "Fully_adjusted"
    )
  }
  
  list(
    results = do.call(rbind, results),
    models = models
  )
}

# ---------------------------------------------------------------------------
# 6. RUN PRIMARY AND SENSITIVITY ANALYSES
# ---------------------------------------------------------------------------

primary_fit <- fit_gene_models(
  primary_data,
  "Primary"
)

sensitivity_fit <- fit_gene_models(
  sensitivity_data,
  "Sensitivity"
)

cox_results <- rbind(
  primary_fit$results,
  sensitivity_fit$results
)

cox_results$FDR_BH <- NA_real_

correction_groups <- unique(
  cox_results[, c("Cohort", "Model")]
)

for (i in seq_len(nrow(correction_groups))) {
  
  selected <- (
    cox_results$Cohort == correction_groups$Cohort[i] &
      cox_results$Model == correction_groups$Model[i]
  )
  
  cox_results$FDR_BH[selected] <- p.adjust(
    cox_results$P_value[selected],
    method = "BH"
  )
}

cat("\n========== GENE COX RESULTS ==========\n")
print(cox_results, row.names = FALSE)




###############################################################################
# SCRIPT 30R — BLOCK 7–9
# HLA-F / TAP1: model diagnostics, PH tests, forest plots and saving
###############################################################################

# ---------------------------------------------------------------------------
# 7. MODEL FIT STATISTICS AND PROPORTIONAL HAZARDS TESTS
# ---------------------------------------------------------------------------

all_gene_models <- c(
  primary_fit$models,
  sensitivity_fit$models
)

model_fit_results <- list()
ph_results <- list()

for (model_name in names(all_gene_models)) {
  
  model <- all_gene_models[[model_name]]
  model_summary <- summary(model)
  
  model_fit_results[[model_name]] <- data.frame(
    Model_ID = model_name,
    N = model$n,
    Events = model$nevent,
    AIC = AIC(model),
    Concordance = unname(model_summary$concordance[1]),
    Concordance_SE = unname(model_summary$concordance[2]),
    stringsAsFactors = FALSE
  )
  
  ph_test <- survival::cox.zph(
    model,
    transform = "km"
  )
  
  ph_table <- as.data.frame(ph_test$table)
  ph_table$Variable <- rownames(ph_table)
  rownames(ph_table) <- NULL
  ph_table$Model_ID <- model_name
  
  ph_results[[model_name]] <- ph_table
}

model_fit_results <- do.call(rbind, model_fit_results)
rownames(model_fit_results) <- NULL

ph_results <- do.call(rbind, ph_results)
rownames(ph_results) <- NULL

cat("\n========== GENE COX MODEL FIT ==========\n")
print(model_fit_results, row.names = FALSE)

cat("\n========== GENE COX PH TESTS ==========\n")
print(ph_results, row.names = FALSE)


# ---------------------------------------------------------------------------
# 8. PREPARE ROYAL-STYLE FOREST PLOTS
# ---------------------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}

plot_data <- cox_results

plot_data$Model <- factor(
  plot_data$Model,
  levels = c("Univariate", "Fully_adjusted")
)

plot_data$Gene <- factor(
  plot_data$Gene,
  levels = c("HLA-F", "TAP1")
)

royal_colors <- c(
  "Univariate" = "#4169E1",
  "Fully_adjusted" = "#D4AF37"
)

# Output directories
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/models", recursive = TRUE, showWarnings = FALSE)

# Unique output filenames
output_files <- c(
  "results/tables/TCGA_SKCM_TAP_gene_Cox_results.csv",
  "results/tables/TCGA_SKCM_TAP_gene_Cox_model_fit.csv",
  "results/tables/TCGA_SKCM_TAP_gene_Cox_PH_tests.csv",
  "results/tables/TCGA_SKCM_TAP_gene_scaling.csv",
  "results/models/TCGA_SKCM_TAP_gene_survival_analysis.rds",
  "results/figures/TCGA_SKCM_TAP_gene_Cox_forest_primary.png",
  "results/figures/TCGA_SKCM_TAP_gene_Cox_forest_sensitivity.png"
)

existing_files <- output_files[file.exists(output_files)]

if (length(existing_files) > 0) {
  stop(
    "To avoid overwriting existing files, rename the outputs first:\n",
    paste(existing_files, collapse = "\n")
  )
}

# Create one forest plot per cohort
for (cohort_name in c("Primary", "Sensitivity")) {
  
  cohort_plot_data <- plot_data[
    plot_data$Cohort == cohort_name,
    ,
    drop = FALSE
  ]
  
  p <- ggplot2::ggplot(
    cohort_plot_data,
    ggplot2::aes(
      x = HR,
      y = Gene,
      xmin = CI_lower,
      xmax = CI_upper,
      color = Model
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "#191970",
      linewidth = 0.7
    ) +
    ggplot2::geom_errorbarh(
      height = 0.18,
      position = ggplot2::position_dodge(width = 0.5),
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.5),
      size = 3
    ) +
    ggplot2::scale_color_manual(values = royal_colors) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = paste(
        "HLA-F and TAP1 —",
        cohort_name,
        "cohort"
      ),
      subtitle = "Cox proportional hazards models",
      x = "Hazard ratio (log scale), per 1 SD increase",
      y = NULL,
      color = "Model"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        color = "#191970"
      ),
      axis.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "bottom"
    )
  
  output_plot <- if (cohort_name == "Primary") {
    "results/figures/TCGA_SKCM_TAP_gene_Cox_forest_primary.png"
  } else {
    "results/figures/TCGA_SKCM_TAP_gene_Cox_forest_sensitivity.png"
  }
  
  ggplot2::ggsave(
    filename = output_plot,
    plot = p,
    width = 8,
    height = 5.5,
    dpi = 300
  )
}


# ---------------------------------------------------------------------------
# 9. SAVE ALL RESULTS
# ---------------------------------------------------------------------------

utils::write.csv(
  cox_results,
  "results/tables/TCGA_SKCM_TAP_gene_Cox_results.csv",
  row.names = FALSE
)

utils::write.csv(
  model_fit_results,
  "results/tables/TCGA_SKCM_TAP_gene_Cox_model_fit.csv",
  row.names = FALSE
)

utils::write.csv(
  ph_results,
  "results/tables/TCGA_SKCM_TAP_gene_Cox_PH_tests.csv",
  row.names = FALSE
)

utils::write.csv(
  gene_scaling,
  "results/tables/TCGA_SKCM_TAP_gene_scaling.csv",
  row.names = FALSE
)

gene_survival_analysis <- list(
  cox_results = cox_results,
  model_fit = model_fit_results,
  PH_tests = ph_results,
  gene_scaling = gene_scaling,
  primary_models = primary_fit$models,
  sensitivity_models = sensitivity_fit$models,
  analysis_parameters = list(
    genes = target_genes,
    adjustment_variables = c(
      "Immune_group",
      "Age_per_10y",
      "Sex"
    ),
    reference_immune_group = "Immune-low",
    reference_sex = "female",
    gene_expression_scaling = "Primary cohort mean and SD",
    PH_test = "Schoenfeld residuals",
    ties = "Efron"
  )
)

saveRDS(
  gene_survival_analysis,
  "results/models/TCGA_SKCM_TAP_gene_survival_analysis.rds"
)

# ---------------------------------------------------------------------------
# FINAL OUTPUT AUDIT
# ---------------------------------------------------------------------------

cat("\n========== SCRIPT 30R OUTPUT AUDIT ==========\n")

final_inventory <- data.frame(
  File = output_files,
  Exists = file.exists(output_files),
  Size_bytes = ifelse(
    file.exists(output_files),
    file.info(output_files)$size,
    NA_real_
  )
)

print(final_inventory, row.names = FALSE)

stopifnot(
  all(final_inventory$Exists),
  all(final_inventory$Size_bytes > 0),
  nrow(cox_results) == 8,
  nrow(model_fit_results) == 8,
  nrow(ph_results) >= 8
)

cat("\nSCRIPT 30R COMPLETED SUCCESSFULLY\n")

