# 09_GEO_CpG_genomic_distribution.R

# Purpose:
# Describe the genomic distribution of GSE144487 CpG probes.
#
# Analyses:
#   1. Gene-region distribution
#   2. Promoter-associated CpGs
#   3. CpG island distribution
#   4. Chromosome distribution
#   5. Gene-region x CpG-island distribution

# 1. Project directories --------------------------------------

project_dir <- getwd()

processed_dir <- file.path(
  project_dir,
  "data", "GEO", "processed"
)

tables_dir <- file.path(
  project_dir,
  "results", "tables"
)

figures_dir <- file.path(
  project_dir,
  "results", "figures"
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. LOAD & VALIDATE INPUT ANNOTATION -------------------------

annotation_file <- file.path(processed_dir, "GSE144487_CpG_annotation.rds")

if (!file.exists(annotation_file)) {
  stop("CpG annotation file not found: ", annotation_file)
}

cpg_annotation <- readRDS(annotation_file)

required_cols <- c("CpG_ID", "UCSC_RefGene_Group", "Relation_to_Island", "chr")
missing_cols  <- setdiff(required_cols, colnames(cpg_annotation))

if (length(missing_cols) > 0) {
  stop("Missing required columns in annotation: ", paste(missing_cols, collapse = ", "))
}

cat("\n--- INPUT ANNOTATION LOADED ---\n")
cat("Total probes:", nrow(cpg_annotation), "\n")


# 3. HELPER FUNCTION: EXPAND SEMICOLON-SEPARATED VALUES -------

expand_cpg_annotations <- function(cpg_ids, annotation_values, value_name = "Category") {
  annotation_values <- as.character(annotation_values)
  annotation_values[is.na(annotation_values)] <- ""
  
  split_values <- strsplit(annotation_values, split = ";", fixed = TRUE)
  n_values     <- lengths(split_values)
  
  expanded_ids    <- rep(cpg_ids, times = n_values)
  expanded_values <- trimws(unlist(split_values, use.names = FALSE))
  
  valid <- !is.na(expanded_values) & nzchar(expanded_values)
  
  result <- data.frame(
    CpG_ID   = expanded_ids[valid],
    Category = expanded_values[valid],
    stringsAsFactors = FALSE
  )
  
  result <- unique(result)
  colnames(result)[2] <- value_name
  return(result)
}


# 4. DISTRIBUTIONS & SUMMARIES --------------------------------

# А. Розподіл за генними регіонами (RefGene Group)
gene_region_data <- expand_cpg_annotations(
  cpg_ids = cpg_annotation$CpG_ID,
  annotation_values = cpg_annotation$UCSC_RefGene_Group,
  value_name = "Gene_region"
)

gene_region_summary <- as.data.frame(table(gene_region_data$Gene_region), stringsAsFactors = FALSE)
colnames(gene_region_summary) <- c("Gene_region", "Number_of_CpGs")
gene_region_summary <- gene_region_summary[order(-gene_region_summary$Number_of_CpGs), , drop = FALSE]
gene_region_summary$Percentage_of_all_CpGs <- round(100 * gene_region_summary$Number_of_CpGs / nrow(cpg_annotation), 4)

# Б. Промоторні CpG (TSS200 та TSS1500)
promoter_regions <- c("TSS200", "TSS1500")
promoter_data    <- gene_region_data[gene_region_data$Gene_region %in% promoter_regions, , drop = FALSE]
promoter_cpg_ids <- unique(promoter_data$CpG_ID)

promoter_summary <- data.frame(
  Metric = c(
    "CpGs annotated to TSS200",
    "CpGs annotated to TSS1500",
    "Unique CpGs annotated to TSS200 or TSS1500",
    "Percentage of all input CpGs"
  ),
  Value = c(
    length(unique(promoter_data$CpG_ID[promoter_data$Gene_region == "TSS200"])),
    length(unique(promoter_data$CpG_ID[promoter_data$Gene_region == "TSS1500"])),
    length(promoter_cpg_ids),
    round(100 * length(promoter_cpg_ids) / nrow(cpg_annotation), 4)
  ),
  stringsAsFactors = FALSE
)

# В. Розподіл відносно CpG-острівців (Relation to Island)
island_data <- cpg_annotation[
  !is.na(cpg_annotation$Relation_to_Island) & nzchar(trimws(cpg_annotation$Relation_to_Island)),
  c("CpG_ID", "Relation_to_Island"), drop = FALSE
]
island_data <- unique(island_data)

island_summary <- as.data.frame(table(island_data$Relation_to_Island), stringsAsFactors = FALSE)
colnames(island_summary) <- c("Island_relation", "Number_of_CpGs")
island_summary <- island_summary[order(-island_summary$Number_of_CpGs), , drop = FALSE]
island_summary$Percentage_of_annotated_CpGs <- round(100 * island_summary$Number_of_CpGs / length(unique(island_data$CpG_ID)), 4)

# Г. Хромосомний розподіл
chromosome_data <- cpg_annotation[
  !is.na(cpg_annotation$chr) & nzchar(trimws(cpg_annotation$chr)),
  c("CpG_ID", "chr"), drop = FALSE
]
chromosome_data <- unique(chromosome_data)

chromosome_summary <- as.data.frame(table(chromosome_data$chr), stringsAsFactors = FALSE)
colnames(chromosome_summary) <- c("Chromosome", "Number_of_CpGs")

chromosome_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
chromosome_summary$sort_order <- match(chromosome_summary$Chromosome, chromosome_order)
chromosome_summary$sort_order[is.na(chromosome_summary$sort_order)] <- 999
chromosome_summary <- chromosome_summary[order(chromosome_summary$sort_order), , drop = FALSE]
chromosome_summary$Percentage_of_annotated_CpGs <- round(100 * chromosome_summary$Number_of_CpGs / length(unique(chromosome_data$CpG_ID)), 4)
chromosome_summary$sort_order <- NULL

# Д. Перетин: Генні регіони х Острівці
region_island_data <- unique(merge(gene_region_data, island_data, by = "CpG_ID"))
region_island_summary <- as.data.frame(xtabs(~ Gene_region + Relation_to_Island, data = region_island_data), stringsAsFactors = FALSE)
colnames(region_island_summary) <- c("Gene_region", "Island_relation", "Number_of_CpGs")
region_island_summary <- region_island_summary[order(region_island_summary$Gene_region, -region_island_summary$Number_of_CpGs), , drop = FALSE]


# 5. SAVE TABLES TO DISK --------------------------------------

write.csv(gene_region_summary,   file.path(tables_dir, "GSE144487_CpG_gene_region_distribution.csv"), row.names = FALSE)
write.csv(promoter_summary,      file.path(tables_dir, "GSE144487_CpG_promoter_distribution.csv"), row.names = FALSE)
write.csv(island_summary,        file.path(tables_dir, "GSE144487_CpG_island_distribution.csv"), row.names = FALSE)
write.csv(chromosome_summary,    file.path(tables_dir, "GSE144487_CpG_chromosome_distribution.csv"), row.names = FALSE)
write.csv(region_island_summary, file.path(tables_dir, "GSE144487_CpG_gene_region_by_island.csv"), row.names = FALSE)

cat("\n[OK] All summary tables successfully saved to results/tables/\n")


# 6. VISUALIZATION (GGPLOT2) ----------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# Графік 1: Розподіл за генними регіонами
p_gene_region <- ggplot(gene_region_summary, aes(x = reorder(Gene_region, Number_of_CpGs), y = Number_of_CpGs)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "GSE144487: CpG distribution by gene region", x = "Gene region", y = "Number of CpGs") +
  theme_minimal()

ggsave(filename = file.path(figures_dir, "GSE144487_CpG_gene_region_distribution.png"), plot = p_gene_region, width = 9, height = 6, dpi = 300)

# Графік 2: Розподіл відносно острівців
p_island <- ggplot(island_summary, aes(x = reorder(Island_relation, Number_of_CpGs), y = Number_of_CpGs)) +
  geom_col(fill = "darkseagreen") +
  coord_flip() +
  labs(title = "GSE144487: CpG distribution by island relation", x = "CpG island relation", y = "Number of CpGs") +
  theme_minimal()

ggsave(filename = file.path(figures_dir, "GSE144487_CpG_island_distribution.png"), plot = p_island, width = 9, height = 6, dpi = 300)

# Графік 3: Хромосомний розподіл
p_chromosome <- ggplot(chromosome_summary, aes(x = factor(Chromosome, levels = Chromosome), y = Number_of_CpGs)) +
  geom_col(fill = "indianred") +
  coord_flip() +
  labs(title = "GSE144487: CpG distribution across chromosomes", x = "Chromosome", y = "Number of CpGs") +
  theme_minimal()

ggsave(filename = file.path(figures_dir, "GSE144487_CpG_chromosome_distribution.png"), plot = p_chromosome, width = 9, height = 7, dpi = 300)

cat("[OK] All visualization plots successfully generated and saved to results/figures/\n")


# 7. FINAL QC AUDIT LOG ---------------------------------------

cat("\n============================================\n")
cat("CPG ANNOTATION ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("Input CpGs:", nrow(cpg_annotation), "\n")
cat("CpGs with gene-region annotation:", length(unique(gene_region_data$CpG_ID)), "\n")
cat("CpGs with island annotation:", length(unique(island_data$CpG_ID)), "\n")
cat("Unique promoter-associated CpGs:", length(promoter_cpg_ids), "\n")
cat("============================================\n")


