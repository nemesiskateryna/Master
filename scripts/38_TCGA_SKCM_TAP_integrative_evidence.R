###############################################################################
# SCRIPT 30T
# INTEGRATIVE EVIDENCE: TAP1 / HLA-F / CpG / SURVIVAL
#
# STEP 1 — READ-ONLY INPUT INVENTORY
###############################################################################

setwd("/Users/katie/master/Master")

cat("\n========== SCRIPT 30T: STEP 1 ==========\n")

# Directories to inspect
directories <- c(
  "results/tables",
  "results/models"
)

# Relevant filename patterns
patterns <- paste(
  c(
    "GO",
    "enrichment",
    "TAP",
    "HLA",
    "CpG",
    "correlation",
    "Cox",
    "survival"
  ),
  collapse = "|"
)

# Collect matching files
input_files <- unlist(
  lapply(directories, function(directory) {
    
    files <- list.files(
      directory,
      pattern = patterns,
      full.names = TRUE,
      recursive = FALSE,
      ignore.case = TRUE
    )
    
    files
  }),
  use.names = FALSE
)

# Inventory
input_inventory <- data.frame(
  File = input_files,
  Exists = file.exists(input_files),
  Size_bytes = ifelse(
    file.exists(input_files),
    file.info(input_files)$size,
    NA
  ),
  stringsAsFactors = FALSE
)

input_inventory <- input_inventory[
  order(input_inventory$File),
  ,
  drop = FALSE
]

cat("\nMatching input files:", nrow(input_inventory), "\n\n")

print(input_inventory, row.names = FALSE)

cat("\n========== FILE CATEGORY COUNTS ==========\n")

cat(
  "GO / enrichment:",
  sum(grepl("GO|enrichment", input_inventory$File, ignore.case = TRUE)),
  "\n"
)

cat(
  "CpG:",
  sum(grepl("CpG", input_inventory$File, ignore.case = TRUE)),
  "\n"
)

cat(
  "TAP / HLA:",
  sum(grepl("TAP|HLA", input_inventory$File, ignore.case = TRUE)),
  "\n"
)

cat(
  "Correlation:",
  sum(grepl("correlation", input_inventory$File, ignore.case = TRUE)),
  "\n"
)

cat(
  "Cox / survival:",
  sum(grepl("Cox|survival", input_inventory$File, ignore.case = TRUE)),
  "\n"
)

cat("\nSCRIPT 30T STEP 1 COMPLETED\n")



###############################################################################
# SCRIPT 30T
# STEP 3 — VERIFY TAP BINDING ENRICHMENT
###############################################################################

cat("\n========== SCRIPT 30T: STEP 3 ==========\n")

# Exact source files
go_mf_expanded_file <-
  "results/tables/TCGA_SKCM_GO_MF_expanded_background.csv"

go_expanded_global_file <-
  "results/tables/TCGA_SKCM_GO_expanded_global_BH.csv"

# Check availability
go_files <- c(
  GO_MF_expanded = go_mf_expanded_file,
  GO_global_expanded = go_expanded_global_file
)

print(data.frame(
  Source = names(go_files),
  File = unname(go_files),
  Exists = file.exists(go_files)
), row.names = FALSE)

stopifnot(all(file.exists(go_files)))

# Read without modifying source files
go_mf_expanded <- read.csv(
  go_mf_expanded_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

go_global_expanded <- read.csv(
  go_expanded_global_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n========== EXPANDED MF STRUCTURE ==========\n")
cat("Dimensions:", paste(dim(go_mf_expanded), collapse = " x "), "\n")
print(names(go_mf_expanded))

cat("\n========== TAP BINDING: MF RESULT ==========\n")

tap_mf <- go_mf_expanded[
  grepl(
    "TAP binding",
    go_mf_expanded$Description,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]

print(tap_mf, row.names = FALSE)

cat("\n========== EXPANDED GLOBAL STRUCTURE ==========\n")
cat("Dimensions:", paste(dim(go_global_expanded), collapse = " x "), "\n")
print(names(go_global_expanded))

cat("\n========== TAP BINDING: GLOBAL BH RESULT ==========\n")

tap_global <- go_global_expanded[
  grepl(
    "TAP binding",
    go_global_expanded$Description,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]

print(tap_global, row.names = FALSE)

cat("\n========== STEP 3 QC ==========\n")

cat("MF TAP binding rows:", nrow(tap_mf), "\n")
cat("Global TAP binding rows:", nrow(tap_global), "\n")

stopifnot(
  nrow(tap_mf) >= 1,
  nrow(tap_global) >= 1
)

cat("\nSCRIPT 30T STEP 3 COMPLETED\n")


###############################################################################
# SCRIPT 30T
# STEP 4 — AUDIT HLA-F / TAP1 CORRELATION AND COX RESULTS
###############################################################################

cat("\n========== SCRIPT 30T: STEP 4 ==========\n")

# Read exact source tables
cpg_gene_cor <- read.csv(
  "results/tables/TCGA_SKCM_CpG_gene_expression_correlations_v2.csv",
  stringsAsFactors = FALSE
)

tap_cox <- read.csv(
  "results/tables/TCGA_SKCM_TAP_gene_Cox_results.csv",
  stringsAsFactors = FALSE
)

# 1. CpG-gene correlations
cat("\n========== CpG–GENE CORRELATIONS ==========\n")

correlation_subset <- cpg_gene_cor[
  cpg_gene_cor$Gene %in% c("HLA-F", "TAP1"),
  ,
  drop = FALSE
]

correlation_subset <- correlation_subset[
  order(correlation_subset$Gene, correlation_subset$Cohort),
  ,
  drop = FALSE
]

print(correlation_subset, row.names = FALSE)

# 2. Gene-level Cox results
cat("\n========== HLA-F / TAP1 COX RESULTS ==========\n")

cox_subset <- tap_cox[
  tap_cox$Gene %in% c("HLA-F", "TAP1"),
  ,
  drop = FALSE
]

cox_subset <- cox_subset[
  order(cox_subset$Gene, cox_subset$Cohort, cox_subset$Model),
  ,
  drop = FALSE
]

print(cox_subset, row.names = FALSE)

# 3. Check expected groups and models
cat("\n========== GROUP AND MODEL COUNTS ==========\n")

cat("\nCorrelation rows by cohort:\n")
print(table(correlation_subset$Cohort))

cat("\nCox rows by cohort:\n")
print(table(cox_subset$Cohort))

cat("\nCox rows by gene and model:\n")
print(with(cox_subset, table(Gene, Model)))

# 4. Numerical QC
correlation_numeric_ok <- all(
  is.finite(correlation_subset$Spearman_rho) &
    is.finite(correlation_subset$P_value) &
    is.finite(correlation_subset$FDR_BH)
)

cox_numeric_ok <- all(
  is.finite(cox_subset$HR) &
    is.finite(cox_subset$CI_lower) &
    is.finite(cox_subset$CI_upper) &
    is.finite(cox_subset$P_value) &
    is.finite(cox_subset$FDR_BH)
)

cox_ci_ok <- all(
  cox_subset$CI_lower < cox_subset$HR &
    cox_subset$HR < cox_subset$CI_upper
)

cat("\n========== NUMERICAL QC ==========\n")
cat("Correlation values finite:", correlation_numeric_ok, "\n")
cat("Cox values finite:", cox_numeric_ok, "\n")
cat("Cox confidence intervals valid:", cox_ci_ok, "\n")

# 5. Expected cohort sizes
expected_correlation_n <- c(
  Primary = 456,
  Sensitivity = 455
)

correlation_n_ok <- all(
  vapply(
    names(expected_correlation_n),
    function(cohort_name) {
      all(
        correlation_subset$N[
          correlation_subset$Cohort == cohort_name
        ] == expected_correlation_n[[cohort_name]]
      )
    },
    logical(1)
  )
)

cox_n_ok <- all(
  vapply(
    names(expected_correlation_n),
    function(cohort_name) {
      all(
        cox_subset$N[
          cox_subset$Cohort == cohort_name
        ] == expected_correlation_n[[cohort_name]]
      )
    },
    logical(1)
  )
)

cat("Correlation cohort sizes correct:", correlation_n_ok, "\n")
cat("Cox cohort sizes correct:", cox_n_ok, "\n")

# 6. Final assertions
stopifnot(
  nrow(correlation_subset) == 4,
  nrow(cox_subset) == 8,
  setequal(correlation_subset$Gene, c("HLA-F", "TAP1")),
  setequal(cox_subset$Gene, c("HLA-F", "TAP1")),
  setequal(correlation_subset$Cohort, c("Primary", "Sensitivity")),
  setequal(cox_subset$Cohort, c("Primary", "Sensitivity")),
  correlation_numeric_ok,
  cox_numeric_ok,
  cox_ci_ok,
  correlation_n_ok,
  cox_n_ok
)

cat("\nSCRIPT 30T STEP 4 COMPLETED SUCCESSFULLY\n")



###############################################################################
# SCRIPT 30T
# STEP 5 — BUILD INTEGRATIVE EVIDENCE TABLES
###############################################################################

cat("\n========== SCRIPT 30T: STEP 5 ==========\n")

# ---------------------------------------------------------------------------
# 1. Load exact source tables
# ---------------------------------------------------------------------------

cpg_annotation <- read.csv(
  "results/tables/TCGA_SKCM_ElasticNet_140_CpG_annotated_hg19.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cpg_gene_cor <- read.csv(
  "results/tables/TCGA_SKCM_CpG_gene_expression_correlations_v2.csv",
  stringsAsFactors = FALSE
)

tap_cox <- read.csv(
  "results/tables/TCGA_SKCM_TAP_gene_Cox_results.csv",
  stringsAsFactors = FALSE
)

go_mf_expanded <- read.csv(
  "results/tables/TCGA_SKCM_GO_MF_expanded_background.csv",
  stringsAsFactors = FALSE
)

go_global_expanded <- read.csv(
  "results/tables/TCGA_SKCM_GO_expanded_global_BH.csv",
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 2. Define target CpG-gene pairs
# ---------------------------------------------------------------------------

target_pairs <- data.frame(
  CpG = c("cg16890093", "cg07016276"),
  Gene = c("TAP1", "HLA-F"),
  stringsAsFactors = FALSE
)

# Select corresponding annotation records
annotation_subset <- cpg_annotation[
  cpg_annotation$Feature %in% target_pairs$CpG,
  ,
  drop = FALSE
]

annotation_subset <- annotation_subset[
  match(target_pairs$CpG, annotation_subset$Feature),
  ,
  drop = FALSE
]

stopifnot(
  nrow(annotation_subset) == 2,
  identical(annotation_subset$Feature, target_pairs$CpG)
)

# Confirm that each CpG annotation contains its expected gene
annotation_gene_check <- vapply(
  seq_len(nrow(target_pairs)),
  function(i) {
    
    annotated_genes <- trimws(
      unlist(
        strsplit(
          annotation_subset$UCSC_RefGene_Name[i],
          ";",
          fixed = TRUE
        )
      )
    )
    
    target_pairs$Gene[i] %in% annotated_genes
  },
  logical(1)
)

stopifnot(all(annotation_gene_check))

# ---------------------------------------------------------------------------
# 3. Build annotation table
# ---------------------------------------------------------------------------

annotation_evidence <- data.frame(
  CpG = annotation_subset$Feature,
  Gene = target_pairs$Gene,
  Coefficient = annotation_subset$Coefficient,
  ElasticNet_Direction = annotation_subset$Direction,
  Chromosome = annotation_subset$chr,
  Position_hg19 = annotation_subset$pos,
  Strand = annotation_subset$strand,
  Gene_Region = annotation_subset$UCSC_RefGene_Group,
  Island_Relation = annotation_subset$Relation_to_Island,
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 4. Add CpG-gene correlations
# ---------------------------------------------------------------------------

correlation_evidence <- cpg_gene_cor[
  cpg_gene_cor$Gene %in% target_pairs$Gene,
  ,
  drop = FALSE
]

correlation_primary <- correlation_evidence[
  correlation_evidence$Cohort == "Primary",
  c("Gene", "Spearman_rho", "P_value", "FDR_BH"),
  drop = FALSE
]

names(correlation_primary)[-1] <- c(
  "Correlation_Primary_Rho",
  "Correlation_Primary_P",
  "Correlation_Primary_FDR"
)

correlation_sensitivity <- correlation_evidence[
  correlation_evidence$Cohort == "Sensitivity",
  c("Gene", "Spearman_rho", "P_value", "FDR_BH"),
  drop = FALSE
]

names(correlation_sensitivity)[-1] <- c(
  "Correlation_Sensitivity_Rho",
  "Correlation_Sensitivity_P",
  "Correlation_Sensitivity_FDR"
)

# ---------------------------------------------------------------------------
# 5. Reshape gene-level Cox results
# ---------------------------------------------------------------------------

cox_evidence <- tap_cox[
  tap_cox$Gene %in% target_pairs$Gene,
  ,
  drop = FALSE
]

make_cox_wide <- function(data, cohort_name, model_name, suffix) {
  
  selected <- data[
    data$Cohort == cohort_name &
      data$Model == model_name,
    c(
      "Gene", "HR", "CI_lower", "CI_upper",
      "P_value", "FDR_BH"
    ),
    drop = FALSE
  ]
  
  names(selected)[-1] <- paste0(
    c("HR", "CI_lower", "CI_upper", "P", "FDR"),
    "_", suffix
  )
  
  selected
}

cox_primary_univariate <- make_cox_wide(
  cox_evidence, "Primary", "Univariate",
  "Primary_Univariate"
)

cox_primary_adjusted <- make_cox_wide(
  cox_evidence, "Primary", "Fully_adjusted",
  "Primary_Adjusted"
)

cox_sensitivity_univariate <- make_cox_wide(
  cox_evidence, "Sensitivity", "Univariate",
  "Sensitivity_Univariate"
)

cox_sensitivity_adjusted <- make_cox_wide(
  cox_evidence, "Sensitivity", "Fully_adjusted",
  "Sensitivity_Adjusted"
)

# ---------------------------------------------------------------------------
# 6. Merge all gene-level evidence
# ---------------------------------------------------------------------------

integrative_evidence <- annotation_evidence

tables_to_merge <- list(
  correlation_primary,
  correlation_sensitivity,
  cox_primary_univariate,
  cox_primary_adjusted,
  cox_sensitivity_univariate,
  cox_sensitivity_adjusted
)

for (current_table in tables_to_merge) {
  
  integrative_evidence <- merge(
    integrative_evidence,
    current_table,
    by = "Gene",
    all.x = TRUE,
    sort = FALSE
  )
}

# Restore target order
integrative_evidence <- integrative_evidence[
  match(target_pairs$Gene, integrative_evidence$Gene),
  ,
  drop = FALSE
]

# ---------------------------------------------------------------------------
# 7. Add explicit statistical interpretation
# ---------------------------------------------------------------------------

integrative_evidence$Correlation_Primary_Significant <-
  integrative_evidence$Correlation_Primary_FDR < 0.05

integrative_evidence$Correlation_Sensitivity_Significant <-
  integrative_evidence$Correlation_Sensitivity_FDR < 0.05

integrative_evidence$Cox_Primary_Univariate_Significant <-
  integrative_evidence$FDR_Primary_Univariate < 0.05

integrative_evidence$Cox_Primary_Adjusted_Significant <-
  integrative_evidence$FDR_Primary_Adjusted < 0.05

integrative_evidence$Cox_Sensitivity_Univariate_Significant <-
  integrative_evidence$FDR_Sensitivity_Univariate < 0.05

integrative_evidence$Cox_Sensitivity_Adjusted_Significant <-
  integrative_evidence$FDR_Sensitivity_Adjusted < 0.05

# ---------------------------------------------------------------------------
# 8. Build separate GO TAP binding context table
# ---------------------------------------------------------------------------

tap_mf <- go_mf_expanded[
  grepl(
    "TAP binding",
    go_mf_expanded$Description,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]

tap_global <- go_global_expanded[
  grepl(
    "TAP binding",
    go_global_expanded$Description,
    ignore.case = TRUE
  ),
  ,
  drop = FALSE
]

stopifnot(
  nrow(tap_mf) == 1,
  nrow(tap_global) == 1,
  tap_mf$ID == tap_global$ID
)

go_tap_context <- data.frame(
  GO_ID = tap_mf$ID,
  Description = tap_mf$Description,
  Ontology = tap_global$Ontology,
  Genes = tap_mf$geneID,
  GeneRatio = tap_mf$GeneRatio,
  BackgroundRatio = tap_mf$BgRatio,
  FoldEnrichment = tap_mf$FoldEnrichment,
  Nominal_P = tap_mf$pvalue,
  MF_adjusted_FDR = tap_mf$p.adjust,
  MF_qvalue = tap_mf$qvalue,
  Global_BH_FDR = tap_global$Global_FDR_BH,
  Global_FDR_Significant = tap_global$Global_FDR_Significant,
  Interpretation = ifelse(
    tap_mf$p.adjust < 0.05 &
      !tap_global$Global_FDR_Significant,
    "Significant within MF; not significant after global BH",
    "Requires interpretation using adjusted significance"
  ),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 9. Output paths — unique names, no overwriting
# ---------------------------------------------------------------------------

evidence_file <-
  "results/tables/TCGA_SKCM_TAP_integrative_gene_evidence_v1.csv"

go_context_file <-
  "results/tables/TCGA_SKCM_TAP_integrative_GO_context_v1.csv"

model_file <-
  "results/models/TCGA_SKCM_TAP_integrative_evidence_v1.rds"

output_files <- c(
  evidence_file,
  go_context_file,
  model_file
)

existing_files <- output_files[file.exists(output_files)]

if (length(existing_files) > 0) {
  stop(
    "Overwrite protection: these files already exist:\n",
    paste(existing_files, collapse = "\n"),
    "\nChoose new filenames before saving."
  )
}

# ---------------------------------------------------------------------------
# 10. Save outputs
# ---------------------------------------------------------------------------

write.csv(
  integrative_evidence,
  evidence_file,
  row.names = FALSE
)

write.csv(
  go_tap_context,
  go_context_file,
  row.names = FALSE
)

saveRDS(
  list(
    gene_evidence = integrative_evidence,
    GO_TAP_context = go_tap_context,
    source_files = c(
      "TCGA_SKCM_ElasticNet_140_CpG_annotated_hg19.csv",
      "TCGA_SKCM_CpG_gene_expression_correlations_v2.csv",
      "TCGA_SKCM_TAP_gene_Cox_results.csv",
      "TCGA_SKCM_GO_MF_expanded_background.csv",
      "TCGA_SKCM_GO_expanded_global_BH.csv"
    ),
    cohort_definitions = c(
      Primary = "N=456; 220 events",
      Sensitivity = "N=455; 219 events"
    )
  ),
  model_file
)

# ---------------------------------------------------------------------------
# 11. Final output audit
# ---------------------------------------------------------------------------

final_files <- data.frame(
  File = output_files,
  Exists = file.exists(output_files),
  Size_bytes = file.info(output_files)$size,
  stringsAsFactors = FALSE
)

cat("\n========== INTEGRATIVE GENE EVIDENCE ==========\n")
print(integrative_evidence, row.names = FALSE)

cat("\n========== GO TAP BINDING CONTEXT ==========\n")
print(go_tap_context, row.names = FALSE)

cat("\n========== FINAL OUTPUT AUDIT ==========\n")
print(final_files, row.names = FALSE)

stopifnot(
  nrow(integrative_evidence) == 2,
  nrow(go_tap_context) == 1,
  all(final_files$Exists),
  all(final_files$Size_bytes > 0),
  all(integrative_evidence$Correlation_Primary_Significant),
  all(integrative_evidence$Correlation_Sensitivity_Significant)
)

cat("\nSCRIPT 30T STEP 5 COMPLETED SUCCESSFULLY\n")
