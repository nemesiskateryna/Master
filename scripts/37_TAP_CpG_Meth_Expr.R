
###############################################################################
# SCRIPT 30S — CpG–GENE EXPRESSION ANALYSIS
# STEP 1: READ-ONLY INPUT AND SAMPLE ALIGNMENT AUDIT
###############################################################################

# 1. Load integrated multiomics data
integrated_data <- readRDS(
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN_immune.rds"
)

# 2. Load previously extracted HLA-F and TAP1 expression
tap_gene_data <- readRDS(
  "results/models/TCGA_SKCM_TAP_gene_expression_data.rds"
)

# 3. Load the established survival cohorts
primary_cohort <- readRDS(
  "results/models/TCGA_SKCM_survival_primary_cohort.rds"
)

sensitivity_cohort <- readRDS(
  "results/models/TCGA_SKCM_survival_sensitivity_cohort.rds"
)

# ---------------------------------------------------------------------------
# 4. INSPECT OBJECT STRUCTURES
# ---------------------------------------------------------------------------

cat("\n========== INPUT OBJECT STRUCTURES ==========\n")

cat("\nIntegrated object:\n")
print(names(integrated_data))

cat("\nTAP gene expression object:\n")
print(names(tap_gene_data))
str(tap_gene_data, max.level = 1)

cat("\nPrimary cohort:\n")
print(dim(primary_cohort))
print(names(primary_cohort))

cat("\nSensitivity cohort:\n")
print(dim(sensitivity_cohort))
print(names(sensitivity_cohort))

# ---------------------------------------------------------------------------
# 5. INSPECT METHYLATION MATRIX
# ---------------------------------------------------------------------------

beta_matrix <- integrated_data$beta_matrix

target_cpgs <- c(
  "cg16890093",
  "cg07016276"
)

cat("\n========== METHYLATION MATRIX ==========\n")
cat("Dimensions:", paste(dim(beta_matrix), collapse = " x "), "\n")

cat("\nTarget CpG availability:\n")
print(data.frame(
  CpG = target_cpgs,
  Available = target_cpgs %in% rownames(beta_matrix)
))


# ---------------------------------------------------------------------------
# 6. EXTRACT GENE EXPRESSION
# ---------------------------------------------------------------------------

gene_expression <- as.data.frame(
  tap_gene_data$expression
)

cat("\n========== GENE EXPRESSION ==========\n")
print(dim(gene_expression))
print(names(gene_expression))
print(head(gene_expression))

# ---------------------------------------------------------------------------
# 7. EXTRACT SAMPLE IDs FROM SAVED COHORT OBJECTS
# ---------------------------------------------------------------------------

primary_data <- as.data.frame(
  primary_cohort$cohort
)

sensitivity_data <- as.data.frame(
  sensitivity_cohort$cohort
)

cat("\n========== COHORT STRUCTURES ==========\n")
print(dim(primary_data))
print(names(primary_data))

print(dim(sensitivity_data))
print(names(sensitivity_data))

# ---------------------------------------------------------------------------
# 8. CHECK SAMPLE ID ALIGNMENT
# ---------------------------------------------------------------------------

required_ids <- list(
  Primary = primary_data$Sample_ID,
  Sensitivity = sensitivity_data$Sample_ID
)

expression_ids <- gene_expression$Sample_ID
methylation_ids <- colnames(beta_matrix)

alignment_results <- do.call(
  rbind,
  lapply(names(required_ids), function(cohort_name) {
    
    ids <- required_ids[[cohort_name]]
    
    data.frame(
      Cohort = cohort_name,
      N = length(ids),
      Unique_IDs = length(unique(ids)),
      All_in_expression = all(ids %in% expression_ids),
      All_in_methylation = all(ids %in% methylation_ids),
      No_duplicate_IDs = !anyDuplicated(ids),
      stringsAsFactors = FALSE
    )
  })
)

cat("\n========== SAMPLE ALIGNMENT AUDIT ==========\n")
print(alignment_results, row.names = FALSE)

# ---------------------------------------------------------------------------
# 9. CHECK TARGET GENES
# ---------------------------------------------------------------------------

cat("\n========== TARGET GENE AVAILABILITY ==========\n")

print(data.frame(
  Gene = c("HLA-F", "TAP1"),
  Available = c(
    "HLA-F" %in% names(gene_expression),
    "TAP1" %in% names(gene_expression)
  )
))

cat("\nSTEP 1 AUDIT COMPLETED\n")


###############################################################################
# SCRIPT 30S — STEP 2
# CpG METHYLATION vs GENE EXPRESSION CORRELATION
# Primary and sensitivity cohorts
###############################################################################

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.")
}

# ---------------------------------------------------------------------------
# 1. DEFINE TARGET PAIRS
# ---------------------------------------------------------------------------

target_pairs <- data.frame(
  CpG = c("cg16890093", "cg07016276"),
  Gene = c("TAP1", "HLA-F"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 2. EXTRACT MATCHED CpG AND GENE VALUES
# ---------------------------------------------------------------------------

build_pair_data <- function(cohort_data, cohort_name, cpg, gene) {
  
  sample_ids <- cohort_data$Sample_ID
  
  # Exact sample matching
  methylation_values <- as.numeric(
    beta_matrix[cpg, sample_ids]
  )
  
  expression_index <- match(
    sample_ids,
    gene_expression$Sample_ID
  )
  
  expression_values <- gene_expression[[gene]][
    expression_index
  ]
  
  result <- data.frame(
    Sample_ID = sample_ids,
    CpG = cpg,
    Gene = gene,
    Methylation_beta = methylation_values,
    Gene_expression = as.numeric(expression_values),
    Cohort = cohort_name,
    stringsAsFactors = FALSE
  )
  
  # Remove only incomplete pairs, if any
  result <- result[
    is.finite(result$Methylation_beta) &
      is.finite(result$Gene_expression),
    ,
    drop = FALSE
  ]
  
  result
}

correlation_data <- list()

for (cohort_name in c("Primary", "Sensitivity")) {
  
  cohort_data <- if (cohort_name == "Primary") {
    primary_data
  } else {
    sensitivity_data
  }
  
  for (i in seq_len(nrow(target_pairs))) {
    
    pair_data <- build_pair_data(
      cohort_data = cohort_data,
      cohort_name = cohort_name,
      cpg = target_pairs$CpG[i],
      gene = target_pairs$Gene[i]
    )
    
    object_name <- paste(
      cohort_name,
      target_pairs$CpG[i],
      sep = "_"
    )
    
    correlation_data[[object_name]] <- pair_data
  }
}

# ---------------------------------------------------------------------------
# 3. CALCULATE SPEARMAN CORRELATIONS
# ---------------------------------------------------------------------------

correlation_results <- list()

for (object_name in names(correlation_data)) {
  
  dat <- correlation_data[[object_name]]
  
  test <- stats::cor.test(
    dat$Methylation_beta,
    dat$Gene_expression,
    method = "spearman",
    exact = FALSE
  )
  
  correlation_results[[object_name]] <- data.frame(
    Cohort = unique(dat$Cohort),
    CpG = unique(dat$CpG),
    Gene = unique(dat$Gene),
    N = nrow(dat),
    Spearman_rho = unname(test$estimate),
    P_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

correlation_results <- do.call(
  rbind,
  correlation_results
)

rownames(correlation_results) <- NULL

# BH correction separately within each cohort
correlation_results$FDR_BH <- NA_real_

for (cohort_name in unique(correlation_results$Cohort)) {
  
  selected <- correlation_results$Cohort == cohort_name
  
  correlation_results$FDR_BH[selected] <- stats::p.adjust(
    correlation_results$P_value[selected],
    method = "BH"
  )
}

cat("\n========== CpG–GENE SPEARMAN CORRELATIONS ==========\n")
print(correlation_results, row.names = FALSE)


# ---------------------------------------------------------------------------
# 4. PREPARE UNIQUE OUTPUT FILENAMES — VERSION 2
# ---------------------------------------------------------------------------

plot_files <- c(
  "results/figures/TCGA_SKCM_cg16890093_TAP1_correlation_primary_v2.png",
  "results/figures/TCGA_SKCM_cg07016276_HLA-F_correlation_primary_v2.png",
  "results/figures/TCGA_SKCM_cg16890093_TAP1_correlation_sensitivity_v2.png",
  "results/figures/TCGA_SKCM_cg07016276_HLA-F_correlation_sensitivity_v2.png"
)

table_file <- "results/tables/TCGA_SKCM_CpG_gene_expression_correlations_v2.csv"

model_file <- "results/models/TCGA_SKCM_CpG_gene_expression_correlation_v2.rds"

output_files <- c(
  plot_files,
  table_file,
  model_file
)

existing_files <- output_files[file.exists(output_files)]

if (length(existing_files) > 0) {
  stop(
    "Some v2 outputs already exist. No files were overwritten:\n",
    paste(existing_files, collapse = "\n")
  )
}

cat("\nNew output filenames confirmed. No existing files will be overwritten.\n")

# ---------------------------------------------------------------------------
# 5. CREATE AND SAVE FOUR PLOTS
# ---------------------------------------------------------------------------

for (object_name in names(correlation_data)) {
  
  dat <- correlation_data[[object_name]]
  
  result_row <- correlation_results[
    correlation_results$Cohort == unique(dat$Cohort) &
      correlation_results$CpG == unique(dat$CpG),
    ,
    drop = FALSE
  ]
  
  rho_label <- sprintf(
    "Spearman rho = %.3f",
    result_row$Spearman_rho
  )
  
  p_label <- paste0(
    "p = ",
    format.pval(result_row$P_value, digits = 3, eps = 0.001)
  )
  
  fdr_label <- paste0(
    "FDR = ",
    format.pval(result_row$FDR_BH, digits = 3, eps = 0.001)
  )
  
  annotation_text <- paste(
    rho_label,
    p_label,
    fdr_label,
    paste0("n = ", nrow(dat)),
    sep = "\n"
  )
  
  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = Methylation_beta,
      y = Gene_expression
    )
  ) +
    ggplot2::geom_point(
      color = royal_blue,
      alpha = 0.72,
      size = 2.3
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = royal_gold,
      fill = royal_purple,
      alpha = 0.13,
      linewidth = 1
    ) +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = Inf,
      label = annotation_text,
      hjust = 1.08,
      vjust = 1.15,
      size = 3.8,
      color = royal_navy
    ) +
    ggplot2::labs(
      title = paste(
        unique(dat$CpG),
        "and",
        unique(dat$Gene)
      ),
      subtitle = paste(
        unique(dat$Cohort),
        "cohort — Spearman correlation"
      ),
      x = paste0(
        "DNA methylation beta value (",
        unique(dat$CpG),
        ")"
      ),
      y = paste0(
        "Gene expression (VST): ",
        unique(dat$Gene)
      )
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        color = royal_navy
      ),
      plot.subtitle = ggplot2::element_text(
        color = royal_purple
      ),
      axis.title = ggplot2::element_text(
        face = "bold"
      )
    )
  
  if (unique(dat$Cohort) == "Primary") {
    output_plot <- if (unique(dat$CpG) == "cg16890093") {
      plot_files[1]
    } else {
      plot_files[2]
    }
  } else {
    output_plot <- if (unique(dat$CpG) == "cg16890093") {
      plot_files[3]
    } else {
      plot_files[4]
    }
  }
  
  ggplot2::ggsave(
    filename = output_plot,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# ---------------------------------------------------------------------------
# 6. SAVE TABLE AND COMPLETE ANALYSIS OBJECT
# ---------------------------------------------------------------------------

utils::write.csv(
  correlation_results,
  table_file,
  row.names = FALSE
)

correlation_analysis <- list(
  results = correlation_results,
  paired_data = correlation_data,
  target_pairs = target_pairs,
  method = "Spearman rank correlation",
  multiple_testing = "Benjamini-Hochberg within each cohort",
  plot_trend = "Linear regression line for visualization only",
  cohorts = c(
    "Primary",
    "Sensitivity"
  )
)

saveRDS(
  correlation_analysis,
  model_file
)

# ---------------------------------------------------------------------------
# 7. FINAL OUTPUT AUDIT
# ---------------------------------------------------------------------------

final_inventory <- data.frame(
  File = output_files,
  Exists = file.exists(output_files),
  Size_bytes = ifelse(
    file.exists(output_files),
    file.info(output_files)$size,
    NA_real_
  )
)

cat("\n========== CORRELATION OUTPUT AUDIT ==========\n")
print(final_inventory, row.names = FALSE)

stopifnot(
  nrow(correlation_results) == 4,
  all(is.finite(correlation_results$Spearman_rho)),
  all(final_inventory$Exists),
  all(final_inventory$Size_bytes > 0)
)

cat("\nSCRIPT 30S STEP 2 COMPLETED SUCCESSFULLY\n")

