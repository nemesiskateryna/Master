#06.1_GEO_expression_annotation.R
#
#
# Purpose:
# Annotate GSE65904 expression probes using GPL10558.




# 1. Project setup --------------------------------------------

cat("\n====================================================\n")
cat("06. GEO EXPRESSION PROBE ANNOTATION\n")
cat("====================================================\n")

# Run from the project root directory
if (!file.exists("Master.Rproj")) {
  stop(
    "Master.Rproj not found. ",
    "Please set the working directory to the project root."
  )
}

project_dir <- getwd()

data_dir <- file.path(project_dir, "data")
geo_processed_dir <- file.path(
  data_dir, "GEO", "processed"
)

results_dir <- file.path(
  project_dir, "results"
)

tables_dir <- file.path(
  results_dir, "tables"
)

logs_dir <- file.path(
  results_dir, "logs"
)

dir.create(
  geo_processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  logs_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\nProject directory:", project_dir, "\n")


# 2. Input files ----------------------------------------------

expression_file <- file.path(
  geo_processed_dir,
  "GSE65904_expression_matrix_high_confidence_log2.rds"
)

annotation_file <- file.path(
  geo_processed_dir,
  "GPL10558_annotation_full.rds"
)

if (!file.exists(expression_file)) {
  stop("Expression matrix file not found: ", expression_file)
}

if (!file.exists(annotation_file)) {
  stop("Full GPL10558 annotation file not found: ", annotation_file)
}


# 3. Load input data ------------------------------------------

cat("\n[1/7] Loading expression matrix and GPL annotation...\n")

expr_matrix_log2 <- readRDS(expression_file)

gpl_table_full <- readRDS(annotation_file)

if (!is.matrix(expr_matrix_log2)) {
  stop("Expression data must be a matrix.")
}

if (is.null(rownames(expr_matrix_log2))) {
  stop("Expression matrix has no probe IDs in rownames.")
}

if (anyDuplicated(rownames(expr_matrix_log2))) {
  stop("Duplicated probe IDs found in expression matrix.")
}

if (!"ID" %in% colnames(gpl_table_full)) {
  stop("GPL annotation does not contain the ID column.")
}

cat(
  "Expression matrix:",
  nrow(expr_matrix_log2), "probes x",
  ncol(expr_matrix_log2), "samples\n"
)

cat(
  "GPL annotation:",
  nrow(gpl_table_full), "probes x",
  ncol(gpl_table_full), "columns\n"
)


# 4. Clean annotation fields ----------------------------------

cat("\n[2/7] Cleaning annotation fields...\n")

clean_field <- function(x) {
  
  x <- trimws(as.character(x))
  
  x[
    is.na(x) |
      x == "" |
      toupper(x) %in% c("NA", "N/A", "NULL")
  ] <- NA_character_
  
  x
}

required_columns <- c(
  "ID",
  "Symbol",
  "Entrez_Gene_ID",
  "ILMN_Gene",
  "Source",
  "Probe_Type",
  "Species"
)

missing_columns <- setdiff(
  required_columns,
  colnames(gpl_table_full)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing expected GPL annotation columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

gpl_table_full$ID_clean <- clean_field(
  gpl_table_full$ID
)

gpl_table_full$Symbol_clean <- clean_field(
  gpl_table_full$Symbol
)

gpl_table_full$Entrez_clean <- clean_field(
  gpl_table_full$Entrez_Gene_ID
)

gpl_table_full$ILMN_Gene_clean <- clean_field(
  gpl_table_full$ILMN_Gene
)

gpl_table_full$Source_clean <- clean_field(
  gpl_table_full$Source
)

gpl_table_full$Probe_Type_clean <- clean_field(
  gpl_table_full$Probe_Type
)

gpl_table_full$Species_clean <- clean_field(
  gpl_table_full$Species
)

# Keep only strictly numeric Entrez identifiers
gpl_table_full$Entrez_numeric <- ifelse(
  !is.na(gpl_table_full$Entrez_clean) &
    grepl("^[0-9]+$", gpl_table_full$Entrez_clean),
  gpl_table_full$Entrez_clean,
  NA_character_
)


# 5. Match expression probes to GPL --------------------------

cat("\n[3/7] Matching expression probes to GPL annotation...\n")

expression_probe_ids <- rownames(expr_matrix_log2)

annotation_index <- match(
  expression_probe_ids,
  gpl_table_full$ID_clean
)

matched_annotation <- gpl_table_full[
  annotation_index,
  ,
  drop = FALSE
]

# Restore expression probe IDs as the reference
matched_annotation$ID <- expression_probe_ids

cat(
  "Expression probes:", length(expression_probe_ids), "\n"
)

cat(
  "Successfully matched:",
  sum(!is.na(annotation_index)), "\n"
)

cat(
  "Not found in GPL:",
  sum(is.na(annotation_index)), "\n"
)


# 6. Build final probe annotation -----------------------------

cat("\n[4/7] Building probe annotation table...\n")

probe_annotation <- data.frame(
  
  Probe_ID = expression_probe_ids,
  
  Symbol = matched_annotation$Symbol_clean,
  
  Entrez_Gene_ID = matched_annotation$Entrez_numeric,
  
  ILMN_Gene = matched_annotation$ILMN_Gene_clean,
  
  Source = matched_annotation$Source_clean,
  
  Probe_Type = matched_annotation$Probe_Type_clean,
  
  Species = matched_annotation$Species_clean,
  
  RefSeq_ID = clean_field(matched_annotation$RefSeq_ID),
  
  Unigene_ID = clean_field(matched_annotation$Unigene_ID),
  
  Transcript = clean_field(matched_annotation$Transcript),
  
  Definition = clean_field(matched_annotation$Definition),
  
  stringsAsFactors = FALSE
)

# Annotation status
probe_annotation$Annotation_Status <- ifelse(
  is.na(annotation_index),
  "probe_not_found",
  ifelse(
    !is.na(probe_annotation$Entrez_Gene_ID) &
      !is.na(probe_annotation$Symbol),
    "entrez_and_symbol",
    ifelse(
      !is.na(probe_annotation$Entrez_Gene_ID),
      "entrez_without_symbol",
      ifelse(
        !is.na(probe_annotation$Symbol),
        "symbol_without_entrez",
        ifelse(
          !is.na(probe_annotation$Unigene_ID),
          "unigene_without_entrez_symbol",
          "no_gene_identifier"
        )
      )
    )
  )
)

# Species status
probe_annotation$Species_Status <- ifelse(
  is.na(probe_annotation$Species),
  "species_missing",
  ifelse(
    probe_annotation$Species == "Homo sapiens",
    "human",
    "non_human_or_other"
  )
)

# Flag symbols that may contain multiple mappings
probe_annotation$Potential_Multiple_Symbols <- (
  !is.na(probe_annotation$Symbol) &
    grepl(
      "///|;|\\|",
      probe_annotation$Symbol
    )
)

# Flag probes suitable for initial gene-level mapping
probe_annotation$Eligible_Entrez_Mapping <- (
  !is.na(probe_annotation$Entrez_Gene_ID) &
    !is.na(probe_annotation$Symbol) &
    probe_annotation$Species_Status == "human" &
    !probe_annotation$Potential_Multiple_Symbols
)


# 7. QC summary -----------------------------------------------

cat("\n[5/7] Generating annotation QC...\n")

qc_summary <- data.frame(
  
  Metric = c(
    "Total expression probes",
    "Matched to GPL",
    "Not found in GPL",
    "Probes with Symbol",
    "Probes with numeric Entrez ID",
    "Probes with both Symbol and Entrez ID",
    "Probes with UniGene ID",
    "Probes with neither Symbol nor Entrez ID",
    "Human probes",
    "Non-human or other species",
    "Potential multiple-symbol probes",
    "Eligible for initial Entrez mapping"
  ),
  
  Value = c(
    nrow(probe_annotation),
    
    sum(!is.na(annotation_index)),
    
    sum(is.na(annotation_index)),
    
    sum(!is.na(probe_annotation$Symbol)),
    
    sum(!is.na(probe_annotation$Entrez_Gene_ID)),
    
    sum(
      !is.na(probe_annotation$Symbol) &
        !is.na(probe_annotation$Entrez_Gene_ID)
    ),
    
    sum(!is.na(probe_annotation$Unigene_ID)),
    
    sum(
      is.na(probe_annotation$Symbol) &
        is.na(probe_annotation$Entrez_Gene_ID)
    ),
    
    sum(probe_annotation$Species_Status == "human"),
    
    sum(probe_annotation$Species_Status == "non_human_or_other"),
    
    sum(probe_annotation$Potential_Multiple_Symbols),
    
    sum(probe_annotation$Eligible_Entrez_Mapping)
  ),
  
  stringsAsFactors = FALSE
)

cat("\n--- ANNOTATION QC SUMMARY ---\n")
print(qc_summary, row.names = FALSE)

cat("\n--- ANNOTATION STATUS ---\n")
print(
  sort(
    table(probe_annotation$Annotation_Status),
    decreasing = TRUE
  )
)

cat("\n--- PROBE TYPE DISTRIBUTION ---\n")
print(
  table(
    probe_annotation$Probe_Type,
    useNA = "ifany"
  )
)

cat("\n--- SOURCE DISTRIBUTION ---\n")
print(
  sort(
    table(
      probe_annotation$Source,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)


# 8. Entrez mapping audit -------------------------------------

cat("\n[6/7] Auditing Entrez-to-Symbol relationships...\n")

valid_mapping <- probe_annotation[
  probe_annotation$Eligible_Entrez_Mapping,
  c(
    "Probe_ID",
    "Entrez_Gene_ID",
    "Symbol",
    "Probe_Type"
  ),
  drop = FALSE
]

# Number of unique symbols per Entrez ID
symbols_per_entrez <- aggregate(
  Symbol ~ Entrez_Gene_ID,
  data = valid_mapping,
  FUN = function(x) length(unique(x))
)

# Number of unique Entrez IDs per symbol
entrez_per_symbol <- aggregate(
  Entrez_Gene_ID ~ Symbol,
  data = valid_mapping,
  FUN = function(x) length(unique(x))
)

multi_symbol_entrez <- symbols_per_entrez[
  symbols_per_entrez$Symbol > 1,
  ,
  drop = FALSE
]

multi_entrez_symbols <- entrez_per_symbol[
  entrez_per_symbol$Entrez_Gene_ID > 1,
  ,
  drop = FALSE
]

cat(
  "Entrez IDs associated with multiple symbols:",
  nrow(multi_symbol_entrez), "\n"
)

cat(
  "Symbols associated with multiple Entrez IDs:",
  nrow(multi_entrez_symbols), "\n"
)

# Save probe-level mapping audit
entrez_mapping_audit <- unique(
  valid_mapping[
    ,
    c("Entrez_Gene_ID", "Symbol"),
    drop = FALSE
  ]
)

entrez_mapping_audit <- entrez_mapping_audit[
  order(
    entrez_mapping_audit$Entrez_Gene_ID,
    entrez_mapping_audit$Symbol
  ),
  ,
  drop = FALSE
]


# 9. Save outputs ---------------------------------------------

cat("\n[7/7] Saving annotation outputs...\n")

saveRDS(
  probe_annotation,
  file.path(
    geo_processed_dir,
    "GSE65904_probe_annotation_high_confidence.rds"
  )
)

write.csv(
  probe_annotation,
  file.path(
    geo_processed_dir,
    "GSE65904_probe_annotation_high_confidence.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  qc_summary,
  file.path(
    tables_dir,
    "GSE65904_probe_annotation_QC.csv"
  ),
  row.names = FALSE
)

write.csv(
  probe_annotation[
    probe_annotation$Annotation_Status != "entrez_and_symbol",
    ,
    drop = FALSE
  ],
  file.path(
    tables_dir,
    "GSE65904_unannotated_probes.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  entrez_mapping_audit,
  file.path(
    tables_dir,
    "GSE65904_Entrez_mapping_audit.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  multi_symbol_entrez,
  file.path(
    tables_dir,
    "GSE65904_Entrez_with_multiple_symbols.csv"
  ),
  row.names = FALSE
)

write.csv(
  multi_entrez_symbols,
  file.path(
    tables_dir,
    "GSE65904_symbols_with_multiple_Entrez.csv"
  ),
  row.names = FALSE
)

cat("\nAll annotation outputs saved successfully.\n")

cat("\nMain annotation file:\n")
cat(
  file.path(
    geo_processed_dir,
    "GSE65904_probe_annotation_high_confidence.rds"
  ),
  "\n"
)

cat("\nQC table:\n")
cat(
  file.path(
    tables_dir,
    "GSE65904_probe_annotation_QC.csv"
  ),
  "\n"
)

cat("\n====================================================\n")
cat("SCRIPT 06 COMPLETED\n")
cat("====================================================\n")

# 10. Final gene mapping audit -------------------------------

cat("\n--- FINAL GENE MAPPING AUDIT ---\n")

valid_mapping <- probe_annotation[
  probe_annotation$Eligible_Entrez_Mapping,
  ,
  drop = FALSE
]

# Number of unique genes
unique_entrez <- unique(
  valid_mapping$Entrez_Gene_ID
)

cat(
  "Eligible probes:",
  nrow(valid_mapping), "\n"
)

cat(
  "Unique Entrez IDs:",
  length(unique_entrez), "\n"
)

cat(
  "Mean probes per gene:",
  round(
    nrow(valid_mapping) / length(unique_entrez),
    2
  ), "\n"
)

# Number of probes per Entrez ID
probes_per_gene <- as.data.frame(
  table(valid_mapping$Entrez_Gene_ID),
  stringsAsFactors = FALSE
)

colnames(probes_per_gene) <- c(
  "Entrez_Gene_ID",
  "Number_of_probes"
)

cat("\n--- PROBES PER GENE ---\n")

print(
  summary(probes_per_gene$Number_of_probes)
)

cat(
  "Genes with one probe:",
  sum(probes_per_gene$Number_of_probes == 1), "\n"
)

cat(
  "Genes with multiple probes:",
  sum(probes_per_gene$Number_of_probes > 1), "\n"
)

# Check Entrez IDs associated with multiple symbols
symbols_per_entrez <- aggregate(
  Symbol ~ Entrez_Gene_ID,
  data = valid_mapping,
  FUN = function(x) paste(sort(unique(x)), collapse = " | ")
)

symbol_counts <- aggregate(
  Symbol ~ Entrez_Gene_ID,
  data = valid_mapping,
  FUN = function(x) length(unique(x))
)

ambiguous_entrez <- symbol_counts[
  symbol_counts$Symbol > 1,
  ,
  drop = FALSE
]

cat(
  "\nEntrez IDs with multiple symbols:",
  nrow(ambiguous_entrez), "\n"
)

if (nrow(ambiguous_entrez) > 0) {
  print(
    merge(
      ambiguous_entrez,
      symbols_per_entrez,
      by = "Entrez_Gene_ID"
    )
  )
}

# Save probe count per gene
write.csv(
  probes_per_gene,
  file.path(
    tables_dir,
    "GSE65904_probes_per_gene_QC.csv"
  ),
  row.names = FALSE
)

cat("\nFinal gene mapping audit completed.\n")





