
# ============================================================
# STEP 5.1
# FUNCTIONAL ENRICHMENT: GENE SET PREPARATION
# Project: TCGA-SKCM
# ============================================================

cat("\n========================================\n")
cat("STEP 5.1: ENRICHMENT PREPARATION\n")
cat("========================================\n")

# ------------------------------------------------------------
# 1. Перевірка необхідних об'єктів
# ------------------------------------------------------------

required_objects <- c(
  "robust_pairs_to_save",
  "comparison_correlations",
  "cpg_gene_expression_mapping",
  "gene_mapping_final"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = .GlobalEnv,
    inherits = FALSE
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Відсутні об'єкти: ",
    paste(missing_objects, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 2. Перевірка пакетів
# ------------------------------------------------------------

required_packages <- c(
  "AnnotationDbi",
  "org.Hs.eg.db",
  "clusterProfiler",
  "enrichplot",
  "ReactomePA"
)

package_status <- data.frame(
  Package = required_packages,
  Installed = vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
)

cat("\n--- ENRICHMENT PACKAGE STATUS ---\n")
print(package_status, row.names = FALSE)

# ------------------------------------------------------------
# 3. Формування основного набору генів
# ------------------------------------------------------------

robust_genes <- unique(
  as.character(robust_pairs_to_save$Gene)
)

robust_genes <- sort(
  robust_genes[
    !is.na(robust_genes) &
      nzchar(robust_genes)
  ]
)

# ------------------------------------------------------------
# 4. Формування background gene set
# ------------------------------------------------------------

background_genes <- unique(
  as.character(
    cpg_gene_expression_mapping$Gene[
      !is.na(cpg_gene_expression_mapping$ENSEMBL) &
        !is.na(cpg_gene_expression_mapping$Gene)
    ]
  )
)

background_genes <- sort(
  background_genes[nzchar(background_genes)]
)

# ------------------------------------------------------------
# 5. Основна інформація
# ------------------------------------------------------------

cat("\n--- GENE SET SUMMARY ---\n")

cat(
  "Robust CpG-Gene pairs:",
  nrow(robust_pairs_to_save), "\n"
)

cat(
  "Unique robust genes:",
  length(robust_genes), "\n"
)

cat(
  "Background genes:",
  length(background_genes), "\n"
)

cat("\n--- ROBUST GENES ---\n")
print(robust_genes)

# ------------------------------------------------------------
# 6. Перевірка gene mapping
# ------------------------------------------------------------

cat("\n--- GENE MAPPING COLUMNS ---\n")
print(colnames(gene_mapping_final))

robust_mapping_check <- gene_mapping_final[
  gene_mapping_final$Original_Gene %in% robust_genes,
  ,
  drop = FALSE
]

background_mapping_check <- gene_mapping_final[
  gene_mapping_final$Original_Gene %in% background_genes,
  ,
  drop = FALSE
]

cat("\n--- MAPPING COUNTS ---\n")

cat(
  "Robust genes with mapping:",
  length(unique(robust_mapping_check$Original_Gene)),
  "\n"
)

cat(
  "Background genes with mapping:",
  length(unique(background_mapping_check$Original_Gene)),
  "\n"
)

# ------------------------------------------------------------
# 7. Перевірка alias-derived відповідностей
# ------------------------------------------------------------

cat("\n--- ROBUST GENE MAPPING TYPES ---\n")

if ("Mapping_Type" %in% colnames(robust_mapping_check)) {
  print(
    table(
      robust_mapping_check$Mapping_Type,
      useNA = "ifany"
    )
  )
}

# ------------------------------------------------------------
# 8. Підготовчий QC
# ------------------------------------------------------------

enrichment_prep_qc <- c(
  
  robust_gene_set_nonempty =
    length(robust_genes) > 0,
  
  background_nonempty =
    length(background_genes) > 0,
  
  robust_genes_subset_of_background =
    all(robust_genes %in% background_genes),
  
  robust_mapping_unique =
    !anyDuplicated(
      robust_mapping_check$Original_Gene
    ),
  
  background_mapping_unique =
    !anyDuplicated(
      background_mapping_check$Original_Gene
    ),
  
  robust_gene_count_matches =
    length(unique(robust_mapping_check$Original_Gene)) ==
    length(robust_genes)
  
)

cat("\n--- ENRICHMENT PREPARATION QC ---\n")
print(enrichment_prep_qc)

stopifnot(all(enrichment_prep_qc))

cat("\nSTEP 5.1: PREPARATION QC PASS\n")




# ============================================================
# STEP 5.2
# CONVERT GENE SYMBOLS TO ENTREZ ID
# Project: TCGA-SKCM
# ============================================================

cat("\n========================================\n")
cat("STEP 5.2: ENTREZ ID MAPPING\n")
cat("========================================\n")

# 1. Перевірка необхідних пакетів
stopifnot(
  requireNamespace("AnnotationDbi", quietly = TRUE),
  requireNamespace("org.Hs.eg.db", quietly = TRUE)
)

# 2. Формуємо таблицю Original Gene -> Ensembl ID
gene_ensembl_map <- unique(
  gene_mapping_final[
    ,
    c("Original_Gene", "ENSEMBL", "Mapping_Type"),
    drop = FALSE
  ]
)

gene_ensembl_map <- gene_ensembl_map[
  !is.na(gene_ensembl_map$ENSEMBL) &
    nzchar(gene_ensembl_map$ENSEMBL),
  ,
  drop = FALSE
]

# 3. Перевіряємо унікальність Ensembl ID
cat("\n--- ENSEMBL MAPPING CHECK ---\n")
cat("Mapping rows:", nrow(gene_ensembl_map), "\n")
cat(
  "Unique Ensembl IDs:",
  length(unique(gene_ensembl_map$ENSEMBL)),
  "\n"
)

stopifnot(
  !anyDuplicated(gene_ensembl_map$Original_Gene),
  !anyDuplicated(gene_ensembl_map$ENSEMBL)
)

# 4. Отримуємо Entrez ID через Ensembl ID
ensembl_to_entrez <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = unique(gene_ensembl_map$ENSEMBL),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME")
)

# 5. Видаляємо записи без Entrez ID
ensembl_to_entrez <- ensembl_to_entrez[
  !is.na(ensembl_to_entrez$ENTREZID) &
    nzchar(as.character(ensembl_to_entrez$ENTREZID)),
  ,
  drop = FALSE
]

# 6. Об'єднуємо початкові назви генів із Entrez ID
gene_entrez_mapping <- merge(
  gene_ensembl_map,
  ensembl_to_entrez,
  by = "ENSEMBL",
  all.x = TRUE,
  sort = FALSE
)

gene_entrez_mapping <- gene_entrez_mapping[
  order(gene_entrez_mapping$Original_Gene),
  ,
  drop = FALSE
]

# 7. Формуємо Entrez ID для основного набору
robust_entrez_mapping <- gene_entrez_mapping[
  gene_entrez_mapping$Original_Gene %in% robust_genes &
    !is.na(gene_entrez_mapping$ENTREZID),
  ,
  drop = FALSE
]

background_entrez_mapping <- gene_entrez_mapping[
  gene_entrez_mapping$Original_Gene %in% background_genes &
    !is.na(gene_entrez_mapping$ENTREZID),
  ,
  drop = FALSE
]

robust_entrez <- unique(
  as.character(robust_entrez_mapping$ENTREZID)
)

background_entrez <- unique(
  as.character(background_entrez_mapping$ENTREZID)
)

# 8. Гени без Entrez ID
unmapped_robust_genes <- setdiff(
  robust_genes,
  unique(robust_entrez_mapping$Original_Gene)
)

unmapped_background_genes <- setdiff(
  background_genes,
  unique(background_entrez_mapping$Original_Gene)
)

# 9. Перевірка багатозначних відповідностей
entrez_per_gene <- aggregate(
  ENTREZID ~ Original_Gene,
  data = gene_entrez_mapping[
    !is.na(gene_entrez_mapping$ENTREZID),
    ,
    drop = FALSE
  ],
  FUN = function(x) length(unique(x))
)

multiple_entrez_genes <- entrez_per_gene[
  entrez_per_gene$ENTREZID > 1,
  ,
  drop = FALSE
]

# 10. Виведення результатів
cat("\n--- ENTREZ MAPPING SUMMARY ---\n")
cat("Robust genes:", length(robust_genes), "\n")
cat(
  "Robust genes mapped:",
  length(unique(robust_entrez_mapping$Original_Gene)),
  "\n"
)
cat("Unique robust Entrez IDs:", length(robust_entrez), "\n")

cat("\nBackground genes:", length(background_genes), "\n")
cat(
  "Background genes mapped:",
  length(unique(background_entrez_mapping$Original_Gene)),
  "\n"
)
cat("Unique background Entrez IDs:", length(background_entrez), "\n")

cat("\nUnmapped robust genes:\n")
print(unmapped_robust_genes)

cat("\nUnmapped background genes:\n")
print(unmapped_background_genes)

cat("\nGenes with multiple Entrez IDs:\n")
print(multiple_entrez_genes, row.names = FALSE)

cat("\n--- MAPPING TYPE SUMMARY ---\n")
print(
  table(
    robust_entrez_mapping$Mapping_Type,
    useNA = "ifany"
  )
)

# 11. QC
entrez_mapping_qc <- c(
  robust_entrez_nonempty =
    length(robust_entrez) > 0,
  
  background_entrez_nonempty =
    length(background_entrez) > 0,
  
  robust_entrez_subset_of_background =
    all(robust_entrez %in% background_entrez),
  
  robust_gene_mapping_available =
    length(unique(robust_entrez_mapping$Original_Gene)) > 0,
  
  no_missing_entrez_in_robust_mapping =
    all(!is.na(robust_entrez_mapping$ENTREZID)),
  
  no_duplicate_entrez_in_gene_sets =
    !anyDuplicated(robust_entrez) &&
    !anyDuplicated(background_entrez)
)

cat("\n--- ENTREZ MAPPING QC ---\n")
print(entrez_mapping_qc)

stopifnot(all(entrez_mapping_qc))

cat("\nSTEP 5.2: ENTREZ MAPPING QC PASS\n")




# ============================================================
# STEP 5.3
# ENRICHMENT PACKAGES AND GENE SET CHECKPOINT
# Project: TCGA-SKCM
# ============================================================

cat("\n========================================\n")
cat("STEP 5.3: ENRICHMENT CHECKPOINT\n")
cat("========================================\n")

# 1. Перевірка пакетів
enrichment_packages <- c(
  "AnnotationDbi",
  "org.Hs.eg.db",
  "clusterProfiler",
  "enrichplot",
  "ReactomePA",
  "DOSE"
)

enrichment_package_status <- data.frame(
  Package = enrichment_packages,
  Installed = vapply(
    enrichment_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
)

cat("\n--- ENRICHMENT PACKAGE STATUS ---\n")
print(enrichment_package_status, row.names = FALSE)

# 2. Перевірка Entrez ID
stopifnot(
  length(robust_entrez) == 42,
  length(background_entrez) == 85,
  !anyDuplicated(robust_entrez),
  !anyDuplicated(background_entrez),
  all(robust_entrez %in% background_entrez)
)

# 3. Перевірка, що всі Entrez ID числові
robust_entrez_numeric <- suppressWarnings(
  as.numeric(robust_entrez)
)

background_entrez_numeric <- suppressWarnings(
  as.numeric(background_entrez)
)

cat("\n--- ENTREZ ID VALIDATION ---\n")
cat(
  "Robust IDs numeric:",
  all(is.finite(robust_entrez_numeric)),
  "\n"
)
cat(
  "Background IDs numeric:",
  all(is.finite(background_entrez_numeric)),
  "\n"
)

# 4. Підготовка фінальних наборів
enrichment_gene_sets <- list(
  robust_gene_symbols = robust_genes,
  background_gene_symbols = background_genes,
  robust_entrez = robust_entrez,
  background_entrez = background_entrez,
  robust_mapping = robust_entrez_mapping,
  background_mapping = background_entrez_mapping,
  unmapped_robust_genes = unmapped_robust_genes,
  unmapped_background_genes = unmapped_background_genes
)

# 5. Створення директорій
dir.create(
  "results/models",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

# 6. Збереження наборів
saveRDS(
  enrichment_gene_sets,
  "results/models/TCGA_SKCM_enrichment_gene_sets.rds"
)

write.csv(
  robust_entrez_mapping,
  "results/tables/TCGA_SKCM_enrichment_robust_gene_mapping.csv",
  row.names = FALSE
)

write.csv(
  background_entrez_mapping,
  "results/tables/TCGA_SKCM_enrichment_background_gene_mapping.csv",
  row.names = FALSE
)

# 7. Перевірка збережених файлів
enrichment_checkpoint_files <- c(
  "results/models/TCGA_SKCM_enrichment_gene_sets.rds",
  "results/tables/TCGA_SKCM_enrichment_robust_gene_mapping.csv",
  "results/tables/TCGA_SKCM_enrichment_background_gene_mapping.csv"
)

checkpoint_qc <- c(
  all_packages_available =
    all(enrichment_package_status$Installed),
  
  robust_gene_count_42 =
    length(robust_entrez) == 42,
  
  background_gene_count_85 =
    length(background_entrez) == 85,
  
  robust_subset_of_background =
    all(robust_entrez %in% background_entrez),
  
  all_ids_numeric =
    all(is.finite(robust_entrez_numeric)) &&
    all(is.finite(background_entrez_numeric)),
  
  all_checkpoint_files_exist =
    all(file.exists(enrichment_checkpoint_files))
)

cat("\n--- STEP 5.3 QC ---\n")
print(checkpoint_qc)

cat("\n--- SAVED FILES ---\n")
print(
  data.frame(
    File = enrichment_checkpoint_files,
    Exists = file.exists(enrichment_checkpoint_files)
  ),
  row.names = FALSE
)

cat("\nSTEP 5.3 CHECKPOINT COMPLETED\n")

cat("\n--- MISSING ENRICHMENT PACKAGES ---\n")

print(
  enrichment_package_status[
    !enrichment_package_status$Installed,
  ],
  row.names = FALSE
)

cat("\n--- INSTALLED ENRICHMENT PACKAGES ---\n")

print(
  enrichment_package_status[
    enrichment_package_status$Installed,
  ],
  row.names = FALSE
)

# Перевірка BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Встановлення лише відсутнього пакета
BiocManager::install(
  "ReactomePA",
  ask = FALSE,
  update = FALSE
)

# Перевірка встановлення
reactomepa_installed <- requireNamespace(
  "ReactomePA",
  quietly = TRUE
)

cat(
  "\nReactomePA installed:",
  reactomepa_installed,
  "\n"
)

if (reactomepa_installed) {
  cat(
    "ReactomePA version:",
    as.character(packageVersion("ReactomePA")),
    "\n"
  )
}


# ============================================================
# STEP 5.5
# GO ENRICHMENT ANALYSIS: BP, MF, CC
# Project: TCGA-SKCM
# ============================================================

cat("\n========================================\n")
cat("STEP 5.5: GO ENRICHMENT ANALYSIS\n")
cat("========================================\n")

# 1. Завантаження пакетів
library(clusterProfiler)
library(org.Hs.eg.db)

# 2. Завантаження збережених наборів генів
enrichment_gene_sets <- readRDS(
  "results/models/TCGA_SKCM_enrichment_gene_sets.rds"
)

robust_entrez <- as.character(
  enrichment_gene_sets$robust_entrez
)

background_entrez <- as.character(
  enrichment_gene_sets$background_entrez
)

# 3. Перевірка вхідних даних
stopifnot(
  length(robust_entrez) == 42,
  length(background_entrez) == 85,
  !anyDuplicated(robust_entrez),
  !anyDuplicated(background_entrez),
  all(robust_entrez %in% background_entrez)
)

cat("\nInput gene set:", length(robust_entrez), "\n")
cat("Background:", length(background_entrez), "\n")

# 4. Функція GO enrichment
run_go_enrichment <- function(ontology) {
  
  cat("\nRunning GO:", ontology, "\n")
  
  clusterProfiler::enrichGO(
    gene = robust_entrez,
    universe = background_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 500,
    readable = TRUE
  )
}

# 5. Запуск трьох категорій GO
go_bp <- run_go_enrichment("BP")
go_mf <- run_go_enrichment("MF")
go_cc <- run_go_enrichment("CC")

go_results <- list(
  BP = go_bp,
  MF = go_mf,
  CC = go_cc
)

# 6. Перетворення результатів у таблиці
go_tables <- lapply(
  go_results,
  function(x) as.data.frame(x)
)

# 7. Підсумкова статистика
go_summary <- data.frame(
  Ontology = names(go_tables),
  Total_terms = vapply(
    go_tables,
    nrow,
    integer(1)
  ),
  Significant_FDR_05 = vapply(
    go_tables,
    function(x) {
      if (nrow(x) == 0) return(0L)
      as.integer(sum(x$p.adjust < 0.05, na.rm = TRUE))
    },
    integer(1)
  )
)

cat("\n--- GO ENRICHMENT SUMMARY ---\n")
print(go_summary, row.names = FALSE)

# 8. Перегляд значущих термінів
for (ontology in names(go_tables)) {
  
  cat("\n---", ontology, ": FDR < 0.05 ---\n")
  
  current_table <- go_tables[[ontology]]
  
  if (nrow(current_table) == 0) {
    cat("No GO terms returned.\n")
  } else {
    
    significant_table <- current_table[
      !is.na(current_table$p.adjust) &
        current_table$p.adjust < 0.05,
      ,
      drop = FALSE
    ]
    
    if (nrow(significant_table) == 0) {
      cat("No terms passed FDR < 0.05.\n")
    } else {
      significant_table <- significant_table[
        order(significant_table$p.adjust),
        ,
        drop = FALSE
      ]
      
      print(
        head(
          significant_table[
            ,
            intersect(
              c(
                "ID",
                "Description",
                "GeneRatio",
                "BgRatio",
                "Count",
                "pvalue",
                "p.adjust",
                "geneID"
              ),
              colnames(significant_table)
            ),
            drop = FALSE
          ],
          20
        ),
        row.names = FALSE
      )
    }
  }
}

# 9. Збереження результатів
dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/models",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  go_tables$BP,
  "results/tables/TCGA_SKCM_GO_BP_enrichment.csv",
  row.names = FALSE
)

write.csv(
  go_tables$MF,
  "results/tables/TCGA_SKCM_GO_MF_enrichment.csv",
  row.names = FALSE
)

write.csv(
  go_tables$CC,
  "results/tables/TCGA_SKCM_GO_CC_enrichment.csv",
  row.names = FALSE
)

write.csv(
  go_summary,
  "results/tables/TCGA_SKCM_GO_enrichment_summary.csv",
  row.names = FALSE
)

saveRDS(
  go_results,
  "results/models/TCGA_SKCM_GO_enrichment_results.rds"
)

# 10. Фінальний QC
go_files <- c(
  "results/tables/TCGA_SKCM_GO_BP_enrichment.csv",
  "results/tables/TCGA_SKCM_GO_MF_enrichment.csv",
  "results/tables/TCGA_SKCM_GO_CC_enrichment.csv",
  "results/tables/TCGA_SKCM_GO_enrichment_summary.csv",
  "results/models/TCGA_SKCM_GO_enrichment_results.rds"
)

go_qc <- c(
  BP_result_created = !is.null(go_bp),
  MF_result_created = !is.null(go_mf),
  CC_result_created = !is.null(go_cc),
  BP_table_created = is.data.frame(go_tables$BP),
  MF_table_created = is.data.frame(go_tables$MF),
  CC_table_created = is.data.frame(go_tables$CC),
  all_output_files_exist = all(file.exists(go_files)),
  input_gene_set_valid =
    length(robust_entrez) == 42 &&
    all(robust_entrez %in% background_entrez)
)

cat("\n--- GO ENRICHMENT QC ---\n")
print(go_qc)

stopifnot(all(go_qc))

cat("\n--- SAVED FILES ---\n")
print(
  data.frame(
    File = go_files,
    Exists = file.exists(go_files)
  ),
  row.names = FALSE
)

cat("\nSTEP 5.5: GO ENRICHMENT COMPLETED\n")



# ============================================================
# STEP 5.6
# REACTOME PATHWAY ENRICHMENT
# Project: TCGA-SKCM
# ============================================================

cat("\n========================================\n")
cat("STEP 5.6: REACTOME ENRICHMENT\n")
cat("========================================\n")

# 1. Завантаження пакетів
library(ReactomePA)
library(org.Hs.eg.db)

# 2. Завантаження наборів генів
enrichment_gene_sets <- readRDS(
  "results/models/TCGA_SKCM_enrichment_gene_sets.rds"
)

robust_entrez <- as.character(
  enrichment_gene_sets$robust_entrez
)

background_entrez <- as.character(
  enrichment_gene_sets$background_entrez
)

# 3. Перевірка вхідних даних
stopifnot(
  length(robust_entrez) == 42,
  length(background_entrez) == 85,
  !anyDuplicated(robust_entrez),
  !anyDuplicated(background_entrez),
  all(robust_entrez %in% background_entrez)
)

cat("\nInput genes:", length(robust_entrez), "\n")
cat("Background genes:", length(background_entrez), "\n")

# 4. Reactome enrichment
reactome_result <- ReactomePA::enrichPathway(
  gene = robust_entrez,
  organism = "human",
  universe = background_entrez,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  minGSSize = 5,
  maxGSSize = 500,
  readable = TRUE
)

# 5. Перетворення результатів у таблицю
reactome_table <- as.data.frame(reactome_result)

cat("\nTotal Reactome pathways returned:",
    nrow(reactome_table), "\n")

# 6. Визначення значущих шляхів
if (nrow(reactome_table) > 0) {
  
  reactome_significant <- reactome_table[
    !is.na(reactome_table$p.adjust) &
      reactome_table$p.adjust < 0.05,
    ,
    drop = FALSE
  ]
  
  reactome_significant <- reactome_significant[
    order(reactome_significant$p.adjust),
    ,
    drop = FALSE
  ]
  
} else {
  
  reactome_significant <- reactome_table
}

cat(
  "Significant pathways (FDR < 0.05):",
  nrow(reactome_significant),
  "\n"
)

# 7. Перегляд значущих шляхів
cat("\n--- SIGNIFICANT REACTOME PATHWAYS ---\n")

if (nrow(reactome_significant) == 0) {
  
  cat("No pathways passed FDR < 0.05.\n")
  
} else {
  
  print(
    head(
      reactome_significant[
        ,
        intersect(
          c(
            "ID",
            "Description",
            "GeneRatio",
            "BgRatio",
            "Count",
            "pvalue",
            "p.adjust",
            "geneID"
          ),
          colnames(reactome_significant)
        ),
        drop = FALSE
      ],
      20
    ),
    row.names = FALSE
  )
}

# 8. Перегляд 20 найкращих шляхів за FDR
cat("\n--- TOP 20 REACTOME PATHWAYS ---\n")

if (nrow(reactome_table) > 0) {
  
  reactome_top20 <- reactome_table[
    order(reactome_table$p.adjust, reactome_table$pvalue),
    ,
    drop = FALSE
  ]
  
  print(
    head(
      reactome_top20[
        ,
        intersect(
          c(
            "ID",
            "Description",
            "GeneRatio",
            "BgRatio",
            "Count",
            "pvalue",
            "p.adjust",
            "geneID"
          ),
          colnames(reactome_top20)
        ),
        drop = FALSE
      ],
      20
    ),
    row.names = FALSE
  )
  
} else {
  
  reactome_top20 <- reactome_table
  cat("No Reactome pathways returned.\n")
}

# 9. Підсумкова таблиця
reactome_summary <- data.frame(
  Analysis = "Reactome",
  Input_genes = length(robust_entrez),
  Background_genes = length(background_entrez),
  Total_pathways = nrow(reactome_table),
  Significant_FDR_05 = nrow(reactome_significant)
)

cat("\n--- REACTOME SUMMARY ---\n")
print(reactome_summary, row.names = FALSE)

# 10. Збереження результатів
write.csv(
  reactome_table,
  "results/tables/TCGA_SKCM_Reactome_enrichment.csv",
  row.names = FALSE
)

write.csv(
  reactome_significant,
  "results/tables/TCGA_SKCM_Reactome_significant_pathways.csv",
  row.names = FALSE
)

write.csv(
  reactome_summary,
  "results/tables/TCGA_SKCM_Reactome_enrichment_summary.csv",
  row.names = FALSE
)

saveRDS(
  reactome_result,
  "results/models/TCGA_SKCM_Reactome_enrichment_results.rds"
)

# 11. Фінальний QC
reactome_files <- c(
  "results/tables/TCGA_SKCM_Reactome_enrichment.csv",
  "results/tables/TCGA_SKCM_Reactome_significant_pathways.csv",
  "results/tables/TCGA_SKCM_Reactome_enrichment_summary.csv",
  "results/models/TCGA_SKCM_Reactome_enrichment_results.rds"
)

reactome_qc <- c(
  result_created = !is.null(reactome_result),
  result_is_table = is.data.frame(reactome_table),
  significant_table_is_valid =
    is.data.frame(reactome_significant),
  input_gene_count_correct =
    length(robust_entrez) == 42,
  background_count_correct =
    length(background_entrez) == 85,
  robust_genes_in_background =
    all(robust_entrez %in% background_entrez),
  all_files_saved = all(file.exists(reactome_files))
)

cat("\n--- REACTOME QC ---\n")
print(reactome_qc)

stopifnot(all(reactome_qc))

cat("\n--- SAVED FILES ---\n")
print(
  data.frame(
    File = reactome_files,
    Exists = file.exists(reactome_files)
  ),
  row.names = FALSE
)

cat("\nSTEP 5.6: REACTOME ENRICHMENT COMPLETED\n")




# ============================================================
# STEP 5.7A
# PREPARE EXPANDED BACKGROUND GENE UNIVERSE
# TCGA-SKCM | All 409,556 selected CpGs
# ============================================================

cat("\n========================================\n")
cat("STEP 5.7A: EXPANDED BACKGROUND PREPARATION\n")
cat("========================================\n")

# 1. Перевірка необхідних пакетів
required_packages <- c(
  "minfi",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "AnnotationDbi",
  "org.Hs.eg.db"
)

package_status <- vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

cat("\n--- PACKAGE CHECK ---\n")
print(package_status)

stopifnot(all(package_status))

# 2. Завантаження повного списку відібраних CpG
cpg_keep_check <- readRDS(
  "data/TCGA/processed/TCGA_SKCM_selected_CpG_list.rds"
)

stopifnot(
  is.data.frame(cpg_keep_check),
  "CpG_ID" %in% colnames(cpg_keep_check)
)

selected_cpg_ids <- unique(
  as.character(cpg_keep_check$CpG_ID)
)

cat("\nSelected CpGs:", length(selected_cpg_ids), "\n")

stopifnot(length(selected_cpg_ids) == 409556)

# 3. Завантаження анотації Illumina 450K
suppressPackageStartupMessages(
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
)

annotation_object <- IlluminaHumanMethylation450kanno.ilmn12.hg19

annotation_450k <- minfi::getAnnotation(
  annotation_object
)

annotation_450k_df <- as.data.frame(
  annotation_450k
)

# 4. Зіставлення CpG з анотацією
annotation_idx <- match(
  selected_cpg_ids,
  rownames(annotation_450k_df)
)

annotation_match_qc <- c(
  all_selected_cpgs_annotated =
    all(!is.na(annotation_idx)),
  annotation_rows_unique =
    !anyDuplicated(rownames(annotation_450k_df))
)

cat("\n--- ANNOTATION MATCH QC ---\n")
print(annotation_match_qc)

stopifnot(all(annotation_match_qc))

selected_annotation <- annotation_450k_df[
  annotation_idx,
  ,
  drop = FALSE
]

rownames(selected_annotation) <- selected_cpg_ids

# 5. Отримання символів генів з UCSC RefGene
stopifnot(
  "UCSC_RefGene_Name" %in% colnames(selected_annotation)
)

gene_symbols_by_cpg <- as.character(
  selected_annotation$UCSC_RefGene_Name
)

gene_symbols_list <- strsplit(
  gene_symbols_by_cpg,
  split = ";",
  fixed = TRUE
)

all_annotation_symbols <- unique(
  trimws(unlist(gene_symbols_list))
)

all_annotation_symbols <- all_annotation_symbols[
  !is.na(all_annotation_symbols) &
    nzchar(all_annotation_symbols)
]

cat(
  "\nUnique RefGene symbols from selected CpGs:",
  length(all_annotation_symbols),
  "\n"
)

# 6. Завантаження експресійної матриці
integrated_check <- readRDS(
  "data/TCGA/processed/TCGA_SKCM_integrated_multiomics_KNN.rds"
)

expression_matrix_check <- integrated_check$expression_matrix

expression_ensembl_check <- sub(
  "\\..*$",
  "",
  rownames(expression_matrix_check)
)

expression_ensembl_check <- unique(
  expression_ensembl_check
)

cat(
  "Unique Ensembl IDs in expression matrix:",
  length(expression_ensembl_check),
  "\n"
)

# 7. Пряме зіставлення SYMBOL -> ENSEMBL
direct_symbol_mapping <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = all_annotation_symbols,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENSEMBL", "GENENAME")
)

direct_symbol_mapping <- direct_symbol_mapping[
  !is.na(direct_symbol_mapping$SYMBOL) &
    !is.na(direct_symbol_mapping$ENSEMBL),
  ,
  drop = FALSE
]

direct_symbol_mapping <- direct_symbol_mapping[
  direct_symbol_mapping$ENSEMBL %in% expression_ensembl_check,
  ,
  drop = FALSE
]

direct_symbol_mapping <- unique(
  direct_symbol_mapping
)

direct_symbol_mapping$Mapping_Type <- "Direct"

colnames(direct_symbol_mapping)[
  colnames(direct_symbol_mapping) == "SYMBOL"
] <- "Original_Gene"

# 8. Додавання раніше перевірених alias-derived зіставлень
# Використовуємо лише alias mapping, який уже був перевірений
# на попередньому етапі проєкту.

saved_gene_sets <- readRDS(
  "results/models/TCGA_SKCM_enrichment_gene_sets.rds"
)

validated_background_mapping <- unique(
  saved_gene_sets$background_mapping[
    saved_gene_sets$background_mapping$Mapping_Type ==
      "Alias-derived",
    c("Original_Gene", "ENSEMBL", "Mapping_Type"),
    drop = FALSE
  ]
)

validated_background_mapping <- validated_background_mapping[
  validated_background_mapping$ENSEMBL %in%
    expression_ensembl_check,
  ,
  drop = FALSE
]

# 9. Об'єднання прямих та перевірених alias mapping
wide_gene_mapping <- rbind(
  direct_symbol_mapping[
    ,
    c("Original_Gene", "ENSEMBL", "Mapping_Type"),
    drop = FALSE
  ],
  validated_background_mapping
)

wide_gene_mapping <- unique(wide_gene_mapping)

# 10. Перевірка зіставлення досліджуваних 42 генів
robust_genes_check <- unique(
  saved_gene_sets$robust_gene_symbols
)

mapped_robust_genes <- intersect(
  robust_genes_check,
  unique(wide_gene_mapping$Original_Gene)
)

unmapped_robust_wide <- setdiff(
  robust_genes_check,
  unique(wide_gene_mapping$Original_Gene)
)

# 11. Підсумкові показники
wide_background_symbols <- unique(
  wide_gene_mapping$Original_Gene
)

cat("\n--- EXPANDED BACKGROUND SUMMARY ---\n")

cat(
  "Selected CpGs:",
  length(selected_cpg_ids),
  "\n"
)

cat(
  "Unique RefGene symbols:",
  length(all_annotation_symbols),
  "\n"
)

cat(
  "Direct mapped gene symbols:",
  length(unique(
    direct_symbol_mapping$Original_Gene
  )),
  "\n"
)

cat(
  "Validated alias-derived symbols:",
  length(unique(
    validated_background_mapping$Original_Gene
  )),
  "\n"
)

cat(
  "Expanded background gene symbols:",
  length(wide_background_symbols),
  "\n"
)

cat(
  "Original robust genes:",
  length(robust_genes_check),
  "\n"
)

cat(
  "Robust genes present in expanded background:",
  length(mapped_robust_genes),
  "\n"
)

cat("\nUnmapped robust genes:\n")
print(unmapped_robust_wide)

cat("\n--- MAPPING TYPE SUMMARY ---\n")
print(
  table(wide_gene_mapping$Mapping_Type)
)

# 12. QC
wide_background_qc <- c(
  selected_cpg_count_correct =
    length(selected_cpg_ids) == 409556,
  
  all_selected_cpgs_annotated =
    all(!is.na(annotation_idx)),
  
  annotation_symbols_nonempty =
    length(all_annotation_symbols) > 0,
  
  direct_mapping_nonempty =
    nrow(direct_symbol_mapping) > 0,
  
  wide_background_nonempty =
    length(wide_background_symbols) > 0,
  
  robust_gene_count_correct =
    length(robust_genes_check) == 42,
  
  all_robust_genes_in_wide_background =
    length(unmapped_robust_wide) == 0,
  
  all_mapped_ensembl_in_expression =
    all(
      wide_gene_mapping$ENSEMBL %in%
        expression_ensembl_check
    )
)

cat("\n--- EXPANDED BACKGROUND QC ---\n")
print(wide_background_qc)

stopifnot(all(wide_background_qc))

cat("\nSTEP 5.7A: EXPANDED BACKGROUND PREPARATION PASS\n")



# ============================================================
# STEP 5.7B
# AUDIT EXPANDED BACKGROUND GENE MAPPING
# ============================================================

cat("\n========================================\n")
cat("STEP 5.7B: MAPPING AMBIGUITY AUDIT\n")
cat("========================================\n")

# 1. Кількість Ensembl ID на один символ гена
ensembl_per_symbol <- aggregate(
  ENSEMBL ~ Original_Gene,
  data = wide_gene_mapping,
  FUN = function(x) length(unique(x))
)

ambiguous_symbols <- ensembl_per_symbol[
  ensembl_per_symbol$ENSEMBL > 1,
  ,
  drop = FALSE
]

cat("\nSymbols with multiple Ensembl IDs:",
    nrow(ambiguous_symbols), "\n")

cat("\n--- AMBIGUOUS SYMBOLS ---\n")
print(ambiguous_symbols, row.names = FALSE)

# 2. Перевірка повторних Ensembl ID для різних символів
symbols_per_ensembl <- aggregate(
  Original_Gene ~ ENSEMBL,
  data = wide_gene_mapping,
  FUN = function(x) length(unique(x))
)

shared_ensembl <- symbols_per_ensembl[
  symbols_per_ensembl$Original_Gene > 1,
  ,
  drop = FALSE
]

cat("\nEnsembl IDs assigned to multiple symbols:",
    nrow(shared_ensembl), "\n")

cat("\n--- SHARED ENSEMBL IDs ---\n")
print(shared_ensembl, row.names = FALSE)

# 3. Перевірка mapping для 42 robust genes
old_robust_mapping <- unique(
  saved_gene_sets$robust_mapping[
    ,
    c("Original_Gene", "ENSEMBL", "Mapping_Type"),
    drop = FALSE
  ]
)

new_robust_mapping <- wide_gene_mapping[
  wide_gene_mapping$Original_Gene %in% robust_genes_check,
  ,
  drop = FALSE
]

old_robust_mapping <- unique(old_robust_mapping)
new_robust_mapping <- unique(new_robust_mapping)

robust_mapping_comparison <- merge(
  old_robust_mapping,
  new_robust_mapping,
  by = "Original_Gene",
  suffixes = c("_old", "_new"),
  all = TRUE
)

robust_mapping_comparison$Same_ENSEMBL <-
  robust_mapping_comparison$ENSEMBL_old ==
  robust_mapping_comparison$ENSEMBL_new

robust_mapping_comparison$Same_Mapping_Type <-
  robust_mapping_comparison$Mapping_Type_old ==
  robust_mapping_comparison$Mapping_Type_new

cat("\n--- ROBUST MAPPING COMPARISON ---\n")
print(robust_mapping_comparison, row.names = FALSE)

# 4. Перевірка унікальності Ensembl для robust genes
robust_ensembl_counts <- aggregate(
  ENSEMBL ~ Original_Gene,
  data = new_robust_mapping,
  FUN = function(x) length(unique(x))
)

cat("\n--- ROBUST GENES WITH MULTIPLE ENSEMBL IDs ---\n")
print(
  robust_ensembl_counts[
    robust_ensembl_counts$ENSEMBL > 1,
    ,
    drop = FALSE
  ],
  row.names = FALSE
)

# 5. QC
mapping_audit_qc <- c(
  all_robust_genes_present =
    setequal(
      unique(new_robust_mapping$Original_Gene),
      robust_genes_check
    ),
  
  no_robust_gene_has_multiple_ensembl =
    all(robust_ensembl_counts$ENSEMBL == 1),
  
  robust_ensembl_matches_previous =
    all(robust_mapping_comparison$Same_ENSEMBL),
  
  robust_mapping_type_matches =
    all(robust_mapping_comparison$Same_Mapping_Type),
  
  no_missing_ensembl =
    !anyNA(wide_gene_mapping$ENSEMBL),
  
  all_ensembl_present_in_expression =
    all(
      wide_gene_mapping$ENSEMBL %in%
        expression_ensembl_check
    )
)

cat("\n--- MAPPING AUDIT QC ---\n")
print(mapping_audit_qc)

cat("\nSTEP 5.7B AUDIT FINISHED\n")



# ============================================================
# STEP 5.7C
# DETAILED SHARED ENSEMBL ID AUDIT
# ============================================================

cat("\n========================================\n")
cat("STEP 5.7C: SHARED ENSEMBL AUDIT\n")
cat("========================================\n")

# 1. Вивести всі пари SYMBOL - ENSEMBL
shared_ensembl_details <- wide_gene_mapping[
  wide_gene_mapping$ENSEMBL %in% shared_ensembl$ENSEMBL,
  c("ENSEMBL", "Original_Gene", "Mapping_Type"),
  drop = FALSE
]

shared_ensembl_details <- shared_ensembl_details[
  order(
    shared_ensembl_details$ENSEMBL,
    shared_ensembl_details$Original_Gene
  ),
  ,
  drop = FALSE
]

cat("\n--- SHARED ENSEMBL DETAILS ---\n")
print(shared_ensembl_details, row.names = FALSE)

# 2. Перевірити неоднозначні SYMBOL -> ENSEMBL
cat("\n--- SYMBOLS WITH MULTIPLE ENSEMBL IDs ---\n")
print(ambiguous_symbols, row.names = FALSE)

# 3. Підрахувати унікальні ідентифікатори
cat("\n--- BACKGROUND MAPPING COUNTS ---\n")

cat(
  "Mapping rows:",
  nrow(wide_gene_mapping),
  "\n"
)

cat(
  "Unique gene symbols:",
  length(unique(wide_gene_mapping$Original_Gene)),
  "\n"
)

cat(
  "Unique Ensembl IDs:",
  length(unique(wide_gene_mapping$ENSEMBL)),
  "\n"
)

cat(
  "Ensembl IDs shared by multiple symbols:",
  length(unique(shared_ensembl$ENSEMBL)),
  "\n"
)

# 4. QC
shared_mapping_qc <- c(
  shared_ensembl_count_correct =
    length(unique(shared_ensembl$ENSEMBL)) == 12,
  
  shared_mapping_details_complete =
    nrow(shared_ensembl_details) ==
    sum(shared_ensembl$Original_Gene),
  
  no_missing_symbols =
    !anyNA(shared_ensembl_details$Original_Gene),
  
  no_missing_ensembl =
    !anyNA(shared_ensembl_details$ENSEMBL),
  
  robust_mapping_still_unchanged =
    all(robust_mapping_comparison$Same_ENSEMBL) &&
    all(robust_mapping_comparison$Same_Mapping_Type)
)

cat("\n--- SHARED MAPPING QC ---\n")
print(shared_mapping_qc)

stopifnot(all(shared_mapping_qc))

cat("\nSTEP 5.7C: AUDIT PASS\n")




# ============================================================
# STEP 5.7D
# ENSEMBL TO ENTREZ MAPPING: EXPANDED BACKGROUND
# ============================================================

cat("\n========================================\n")
cat("STEP 5.7D: ENSEMBL TO ENTREZ MAPPING\n")
cat("========================================\n")

# 1. Унікальні Ensembl ID
wide_ensembl_ids <- unique(
  as.character(wide_gene_mapping$ENSEMBL)
)

wide_ensembl_ids <- wide_ensembl_ids[
  !is.na(wide_ensembl_ids) &
    nzchar(wide_ensembl_ids)
]

cat(
  "\nUnique Ensembl IDs:",
  length(wide_ensembl_ids),
  "\n"
)

# 2. Перетворення Ensembl -> Entrez
wide_ensembl_entrez_raw <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = wide_ensembl_ids,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "ENTREZID", "SYMBOL")
)

wide_ensembl_entrez_raw <- unique(
  wide_ensembl_entrez_raw
)

wide_ensembl_entrez_raw$ENTREZID <- as.character(
  wide_ensembl_entrez_raw$ENTREZID
)

# 3. Відокремлення зіставлених і незіставлених ID
wide_ensembl_entrez <- wide_ensembl_entrez_raw[
  !is.na(wide_ensembl_entrez_raw$ENTREZID) &
    nzchar(wide_ensembl_entrez_raw$ENTREZID),
  ,
  drop = FALSE
]

unmapped_ensembl_ids <- setdiff(
  wide_ensembl_ids,
  unique(wide_ensembl_entrez$ENSEMBL)
)

# 4. Перевірка Ensembl -> multiple Entrez
entrez_per_ensembl <- aggregate(
  ENTREZID ~ ENSEMBL,
  data = wide_ensembl_entrez,
  FUN = function(x) length(unique(x))
)

multiple_entrez_per_ensembl <- entrez_per_ensembl[
  entrez_per_ensembl$ENTREZID > 1,
  ,
  drop = FALSE
]

# 5. Перевірка Entrez -> multiple Ensembl
ensembl_per_entrez <- aggregate(
  ENSEMBL ~ ENTREZID,
  data = wide_ensembl_entrez,
  FUN = function(x) length(unique(x))
)

multiple_ensembl_per_entrez <- ensembl_per_entrez[
  ensembl_per_entrez$ENSEMBL > 1,
  ,
  drop = FALSE
]

# 6. Формування розширеного Entrez background
wide_background_entrez <- unique(
  wide_ensembl_entrez$ENTREZID
)

# 7. Перевірка 42 досліджуваних генів
robust_ensembl_ids <- unique(
  new_robust_mapping$ENSEMBL
)

robust_entrez_mapping_wide <- wide_ensembl_entrez[
  wide_ensembl_entrez$ENSEMBL %in%
    robust_ensembl_ids,
  ,
  drop = FALSE
]

robust_entrez_wide <- unique(
  robust_entrez_mapping_wide$ENTREZID
)

# 8. Порівняння з попереднім Entrez mapping
old_robust_entrez <- unique(
  as.character(robust_entrez)
)

robust_entrez_comparison <- list(
  old_count = length(old_robust_entrez),
  new_count = length(robust_entrez_wide),
  old_not_in_new = setdiff(
    old_robust_entrez,
    robust_entrez_wide
  ),
  new_not_in_old = setdiff(
    robust_entrez_wide,
    old_robust_entrez
  )
)

# 9. Підсумок
cat("\n--- EXPANDED ENTREZ SUMMARY ---\n")

cat(
  "Unique Ensembl IDs:",
  length(wide_ensembl_ids),
  "\n"
)

cat(
  "Mapped Ensembl IDs:",
  length(unique(wide_ensembl_entrez$ENSEMBL)),
  "\n"
)

cat(
  "Unmapped Ensembl IDs:",
  length(unmapped_ensembl_ids),
  "\n"
)

cat(
  "Unique Entrez IDs:",
  length(wide_background_entrez),
  "\n"
)

cat(
  "Ensembl IDs with multiple Entrez IDs:",
  nrow(multiple_entrez_per_ensembl),
  "\n"
)

cat(
  "Entrez IDs with multiple Ensembl IDs:",
  nrow(multiple_ensembl_per_entrez),
  "\n"
)

cat(
  "Robust genes: previous Entrez count:",
  robust_entrez_comparison$old_count,
  "\n"
)

cat(
  "Robust genes: expanded Entrez count:",
  robust_entrez_comparison$new_count,
  "\n"
)

# 10. Детальні таблиці неоднозначностей
cat("\n--- MULTIPLE ENTREZ PER ENSEMBL ---\n")
print(multiple_entrez_per_ensembl, row.names = FALSE)

cat("\n--- MULTIPLE ENSEMBL PER ENTREZ ---\n")
print(multiple_ensembl_per_entrez, row.names = FALSE)

cat("\n--- UNMAPPED ENSEMBL IDs: FIRST 30 ---\n")
print(head(unmapped_ensembl_ids, 30))

cat("\n--- ROBUST ENTREZ COMPARISON ---\n")
print(robust_entrez_comparison)

# 11. QC
wide_entrez_qc <- c(
  background_entrez_nonempty =
    length(wide_background_entrez) > 0,
  
  all_42_robust_genes_have_entrez =
    length(robust_entrez_wide) == 42,
  
  robust_entrez_subset_of_wide_background =
    all(robust_entrez_wide %in% wide_background_entrez),
  
  no_missing_entrez =
    !anyNA(wide_ensembl_entrez$ENTREZID),
  
  no_duplicate_background_entrez =
    !anyDuplicated(wide_background_entrez),
  
  all_robust_entrez_preserved =
    length(robust_entrez_comparison$old_not_in_new) == 0 &&
    length(robust_entrez_comparison$new_not_in_old) == 0
)

cat("\n--- WIDE ENTREZ QC ---\n")
print(wide_entrez_qc)

cat("\nSTEP 5.7D AUDIT FINISHED\n")



# ============================================================
# STEP 5.7E
# AUDIT OF AMBIGUOUS ENSEMBL -> ENTREZ MAPPINGS
# ============================================================

cat("\n========================================\n")
cat("STEP 5.7E: AMBIGUOUS MAPPING AUDIT\n")
cat("========================================\n")

# 1. Повна таблиця неоднозначних Ensembl ID
ambiguous_ensembl_ids <- multiple_entrez_per_ensembl$ENSEMBL

ambiguous_ensembl_details <- wide_ensembl_entrez[
  wide_ensembl_entrez$ENSEMBL %in% ambiguous_ensembl_ids,
  c("ENSEMBL", "ENTREZID", "SYMBOL"),
  drop = FALSE
]

ambiguous_ensembl_details <- unique(
  ambiguous_ensembl_details
)

ambiguous_ensembl_details <- ambiguous_ensembl_details[
  order(
    ambiguous_ensembl_details$ENSEMBL,
    ambiguous_ensembl_details$ENTREZID
  ),
  ,
  drop = FALSE
]

cat("\n--- AMBIGUOUS ENSEMBL MAPPING DETAILS ---\n")
cat(
  "Rows:",
  nrow(ambiguous_ensembl_details),
  "\n"
)

print(
  head(ambiguous_ensembl_details, 40),
  row.names = FALSE
)

# 2. Розподіл кількості Entrez на один Ensembl
cat("\n--- NUMBER OF ENTREZ IDS PER ENSEMBL ---\n")

print(
  table(multiple_entrez_per_ensembl$ENTREZID)
)

# 3. Чи є неоднозначні Ensembl серед 42 robust genes?
robust_ambiguous_ensembl <- intersect(
  ambiguous_ensembl_ids,
  robust_ensembl_ids
)

cat("\n--- ROBUST SIGNATURE AMBIGUOUS ENSEMBL ---\n")
cat(
  "Count:",
  length(robust_ambiguous_ensembl),
  "\n"
)

if (length(robust_ambiguous_ensembl) > 0) {
  print(
    ambiguous_ensembl_details[
      ambiguous_ensembl_details$ENSEMBL %in%
        robust_ambiguous_ensembl,
      ,
      drop = FALSE
    ],
    row.names = FALSE
  )
} else {
  cat("None\n")
}

# 4. Перевірка множинних Ensembl на один Entrez
multiple_ensembl_details <- wide_ensembl_entrez[
  wide_ensembl_entrez$ENTREZID %in%
    multiple_ensembl_per_entrez$ENTREZID,
  c("ENTREZID", "ENSEMBL", "SYMBOL"),
  drop = FALSE
]

multiple_ensembl_details <- unique(
  multiple_ensembl_details
)

multiple_ensembl_details <- multiple_ensembl_details[
  order(
    multiple_ensembl_details$ENTREZID,
    multiple_ensembl_details$ENSEMBL
  ),
  ,
  drop = FALSE
]

cat("\n--- MULTIPLE ENSEMBL PER ENTREZ: DETAILS ---\n")
print(
  multiple_ensembl_details,
  row.names = FALSE
)

# 5. Перевірка повторних пар Ensembl-Entrez
duplicate_ensembl_entrez_pairs <- duplicated(
  wide_ensembl_entrez[
    ,
    c("ENSEMBL", "ENTREZID"),
    drop = FALSE
  ]
)

cat("\n--- DUPLICATE ENSEMBL-ENTREZ PAIRS ---\n")
cat(
  "Duplicate rows:",
  sum(duplicate_ensembl_entrez_pairs),
  "\n"
)

# 6. QC
ambiguity_audit_qc <- c(
  ambiguous_details_cover_all_102 =
    length(unique(ambiguous_ensembl_details$ENSEMBL)) ==
    length(ambiguous_ensembl_ids),
  
  robust_signature_has_no_ambiguous_ensembl =
    length(robust_ambiguous_ensembl) == 0,
  
  multiple_ensembl_entrez_details_complete =
    all(
      multiple_ensembl_per_entrez$ENTREZID %in%
        multiple_ensembl_details$ENTREZID
    ),
  
  no_missing_ensembl =
    !anyNA(ambiguous_ensembl_details$ENSEMBL),
  
  no_missing_entrez =
    !anyNA(ambiguous_ensembl_details$ENTREZID),
  
  no_missing_symbol =
    !anyNA(ambiguous_ensembl_details$SYMBOL)
)

cat("\n--- AMBIGUITY AUDIT QC ---\n")
print(ambiguity_audit_qc)

stopifnot(all(ambiguity_audit_qc))

cat("\nSTEP 5.7E: AUDIT PASS\n")


cat("\n--- WIDE GENE MAPPING COLUMNS ---\n")
print(names(wide_gene_mapping))

cat("\n--- DIRECT MAPPING EXAMPLES ---\n")
print(
  head(
    wide_gene_mapping[
      wide_gene_mapping$Mapping_Type == "Direct",
      ,
      drop = FALSE
    ],
    5
  )
)

cat("\n--- ALIAS-DERIVED MAPPING EXAMPLES ---\n")
print(
  head(
    wide_gene_mapping[
      wide_gene_mapping$Mapping_Type == "Alias-derived",
      ,
      drop = FALSE
    ],
    15
  )
)

cat("\n--- PREVIOUS ALIAS MAPPING TABLE COLUMNS ---\n")
print(names(gene_mapping_final))

cat("\n--- PREVIOUS ALIAS MAPPING EXAMPLES ---\n")
print(
  gene_mapping_final[
    gene_mapping_final$Mapping_Type == "Alias-derived",
    ,
    drop = FALSE
  ]
)



# ============================================================
# STEP 5.7G
# STRICT SYMBOL-BASED ENSEMBL -> ENTREZ MAPPING
# ============================================================

cat("\n========================================\n")
cat("STEP 5.7G: STRICT BACKGROUND MAPPING\n")
cat("========================================\n")

# 1. Додаємо очікуваний символ для кожного зіставлення
wide_gene_mapping_checked <- wide_gene_mapping

wide_gene_mapping_checked$Expression_Gene <-
  wide_gene_mapping_checked$Original_Gene

# Для Alias-derived використовуємо перевірений сучасний символ
alias_rows <- wide_gene_mapping_checked$Mapping_Type ==
  "Alias-derived"

alias_lookup <- gene_mapping_final[
  gene_mapping_final$Mapping_Type == "Alias-derived",
  c("Original_Gene", "Expression_Gene"),
  drop = FALSE
]

alias_match <- match(
  wide_gene_mapping_checked$Original_Gene[alias_rows],
  alias_lookup$Original_Gene
)

wide_gene_mapping_checked$Expression_Gene[alias_rows] <-
  alias_lookup$Expression_Gene[alias_match]

# 2. Перевірка, що всі alias отримали сучасний символ
alias_symbol_qc <- (
  !anyNA(
    wide_gene_mapping_checked$Expression_Gene[alias_rows]
  ) &&
    all(
      nzchar(
        wide_gene_mapping_checked$Expression_Gene[alias_rows]
      )
    )
)

cat(
  "\nAlias rows with valid expression symbols:",
  sum(
    !is.na(
      wide_gene_mapping_checked$Expression_Gene[alias_rows]
    )
  ),
  "of",
  sum(alias_rows),
  "\n"
)

stopifnot(alias_symbol_qc)

# 3. Приєднуємо Entrez mapping за Ensembl ID
wide_gene_entrez_candidates <- merge(
  wide_gene_mapping_checked,
  wide_ensembl_entrez,
  by = "ENSEMBL",
  all.x = TRUE,
  sort = FALSE
)

# 4. Залишаємо тільки точний збіг символу
wide_gene_entrez_strict <- wide_gene_entrez_candidates[
  !is.na(wide_gene_entrez_candidates$ENTREZID) &
    !is.na(wide_gene_entrez_candidates$SYMBOL) &
    wide_gene_entrez_candidates$Expression_Gene ==
    wide_gene_entrez_candidates$SYMBOL,
  ,
  drop = FALSE
]

wide_gene_entrez_strict <- unique(
  wide_gene_entrez_strict
)

# 5. Перевіряємо гени, для яких не залишилося Entrez
mapped_original_genes <- unique(
  wide_gene_entrez_strict$Original_Gene
)

unmapped_original_genes <- setdiff(
  unique(wide_gene_mapping_checked$Original_Gene),
  mapped_original_genes
)

# 6. Перевіряємо, чи один символ має кілька Entrez ID
entrez_per_original_gene <- aggregate(
  ENTREZID ~ Original_Gene,
  data = wide_gene_entrez_strict,
  FUN = function(x) length(unique(x))
)

multiple_entrez_per_gene <- entrez_per_original_gene[
  entrez_per_original_gene$ENTREZID > 1,
  ,
  drop = FALSE
]

# 7. Формуємо унікальний background
wide_background_entrez_strict <- unique(
  as.character(wide_gene_entrez_strict$ENTREZID)
)

# 8. Перевіряємо основні 42 гени
robust_strict_mapping <- wide_gene_entrez_strict[
  wide_gene_entrez_strict$Original_Gene %in% robust_genes,
  ,
  drop = FALSE
]

robust_entrez_strict <- unique(
  as.character(robust_strict_mapping$ENTREZID)
)

# 9. Порівнюємо з попереднім набором 42 Entrez ID
robust_strict_comparison <- list(
  previous_count = length(old_robust_entrez),
  strict_count = length(robust_entrez_strict),
  previous_not_in_strict = setdiff(
    old_robust_entrez,
    robust_entrez_strict
  ),
  strict_not_in_previous = setdiff(
    robust_entrez_strict,
    old_robust_entrez
  )
)

# 10. Підсумок
cat("\n--- STRICT MAPPING SUMMARY ---\n")

cat(
  "Original gene symbols:",
  length(unique(wide_gene_mapping_checked$Original_Gene)),
  "\n"
)

cat(
  "Genes mapped after strict matching:",
  length(mapped_original_genes),
  "\n"
)

cat(
  "Genes without strict Entrez mapping:",
  length(unmapped_original_genes),
  "\n"
)

cat(
  "Unique strict Entrez background:",
  length(wide_background_entrez_strict),
  "\n"
)

cat(
  "Genes with multiple Entrez IDs:",
  nrow(multiple_entrez_per_gene),
  "\n"
)

cat(
  "Robust signature Entrez IDs:",
  length(robust_entrez_strict),
  "\n"
)

cat("\n--- UNMAPPED ORIGINAL GENES: FIRST 40 ---\n")
print(head(unmapped_original_genes, 40))

cat("\n--- GENES WITH MULTIPLE ENTREZ IDs ---\n")
print(multiple_entrez_per_gene, row.names = FALSE)

cat("\n--- ROBUST SIGNATURE COMPARISON ---\n")
print(robust_strict_comparison)

# 11. QC
strict_background_qc <- c(
  alias_symbols_valid = alias_symbol_qc,
  
  strict_background_nonempty =
    length(wide_background_entrez_strict) > 0,
  
  all_42_robust_genes_mapped =
    length(unique(robust_strict_mapping$Original_Gene)) == 42,
  
  all_42_robust_entrez_preserved =
    length(robust_strict_comparison$previous_not_in_strict) == 0 &&
    length(robust_strict_comparison$strict_not_in_previous) == 0,
  
  no_missing_entrez =
    !anyNA(wide_gene_entrez_strict$ENTREZID),
  
  no_duplicate_background_entrez =
    !anyDuplicated(wide_background_entrez_strict),
  
  exact_symbol_matching =
    all(
      wide_gene_entrez_strict$Expression_Gene ==
        wide_gene_entrez_strict$SYMBOL
    )
)

cat("\n--- STRICT BACKGROUND QC ---\n")
print(strict_background_qc)

stopifnot(all(strict_background_qc))

cat("\nSTEP 5.7G: STRICT MAPPING QC PASS\n")





# ============================================================
# STEP 5.7H
# SAVE EXPANDED ENRICHMENT GENE SETS
# ============================================================

cat("\n========================================\n")
cat("STEP 5.7H: SAVE EXPANDED GENE SETS\n")
cat("========================================\n")

# 1. Каталоги
dir.create(
  "results/models",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. Формуємо окремий об'єкт
enrichment_gene_sets_expanded <- list(
  
  robust_genes = robust_genes,
  
  robust_entrez = sort(
    unique(as.character(robust_entrez_strict))
  ),
  
  background_genes = unique(
    wide_gene_mapping_checked$Original_Gene
  ),
  
  background_entrez = sort(
    unique(as.character(wide_background_entrez_strict))
  ),
  
  robust_mapping = robust_strict_mapping,
  
  background_mapping = wide_gene_entrez_strict,
  
  mapping_metadata = list(
    background_type =
      "Genes annotated to selected TCGA-SKCM 450K CpGs and represented in expression matrix",
    
    annotation =
      "Illumina HumanMethylation450K hg19",
    
    gene_mapping =
      "org.Hs.eg.db SYMBOL/ALIAS to Ensembl, followed by exact SYMBOL matching to Entrez",
    
    robust_gene_count =
      length(unique(robust_genes)),
    
    robust_entrez_count =
      length(unique(robust_entrez_strict)),
    
    background_gene_count =
      length(unique(wide_gene_mapping_checked$Original_Gene)),
    
    background_entrez_count =
      length(unique(wide_background_entrez_strict)),
    
    ambiguous_mapping_policy =
      "Retained only exact matches between expected expression symbol and org.Hs.eg.db SYMBOL",
    
    alias_mapping_policy =
      "Used previously validated Expression_Gene mappings",
    
    date_created =
      as.character(Sys.Date())
  ),
  
  qc = list(
    strict_background_qc = strict_background_qc,
    ambiguity_audit_qc = ambiguity_audit_qc,
    robust_mapping_comparison = robust_strict_comparison
  )
)

# 3. Шляхи
expanded_rds_file <-
  "results/models/TCGA_SKCM_enrichment_gene_sets_expanded.rds"

expanded_background_file <-
  "results/tables/TCGA_SKCM_enrichment_expanded_background_mapping.csv"

expanded_robust_file <-
  "results/tables/TCGA_SKCM_enrichment_expanded_robust_mapping.csv"

expanded_entrez_file <-
  "results/tables/TCGA_SKCM_enrichment_expanded_background_entrez.csv"

expanded_qc_file <-
  "results/tables/TCGA_SKCM_enrichment_expanded_background_QC.csv"

# 4. Збереження
saveRDS(
  enrichment_gene_sets_expanded,
  expanded_rds_file
)

write.csv(
  wide_gene_entrez_strict,
  expanded_background_file,
  row.names = FALSE
)

write.csv(
  robust_strict_mapping,
  expanded_robust_file,
  row.names = FALSE
)

write.csv(
  data.frame(
    ENTREZID = enrichment_gene_sets_expanded$background_entrez
  ),
  expanded_entrez_file,
  row.names = FALSE
)

expanded_save_qc <- data.frame(
  Check = c(
    "RDS exists",
    "Background mapping CSV exists",
    "Robust mapping CSV exists",
    "Entrez background CSV exists",
    "Background contains 14206 Entrez IDs",
    "Robust set contains 42 Entrez IDs",
    "All robust IDs belong to background"
  ),
  
  Pass = c(
    file.exists(expanded_rds_file),
    file.exists(expanded_background_file),
    file.exists(expanded_robust_file),
    file.exists(expanded_entrez_file),
    
    length(enrichment_gene_sets_expanded$background_entrez) == 14206,
    
    length(enrichment_gene_sets_expanded$robust_entrez) == 42,
    
    all(
      enrichment_gene_sets_expanded$robust_entrez %in%
        enrichment_gene_sets_expanded$background_entrez
    )
  )
)

write.csv(
  expanded_save_qc,
  expanded_qc_file,
  row.names = FALSE
)

# 5. Перевірка повторного завантаження
expanded_check <- readRDS(expanded_rds_file)

reload_qc <- c(
  correct_background_count =
    length(expanded_check$background_entrez) == 14206,
  
  correct_robust_count =
    length(expanded_check$robust_entrez) == 42,
  
  robust_subset_of_background =
    all(
      expanded_check$robust_entrez %in%
        expanded_check$background_entrez
    ),
  
  no_duplicate_background =
    !anyDuplicated(expanded_check$background_entrez),
  
  no_duplicate_robust =
    !anyDuplicated(expanded_check$robust_entrez)
)

cat("\n--- SAVE QC ---\n")
print(expanded_save_qc)

cat("\n--- RELOAD QC ---\n")
print(reload_qc)

stopifnot(
  all(expanded_save_qc$Pass),
  all(reload_qc)
)

cat("\n--- SAVED FILES ---\n")
print(
  c(
    expanded_rds_file,
    expanded_background_file,
    expanded_robust_file,
    expanded_entrez_file,
    expanded_qc_file
  )
)

cat("\nSTEP 5.7H: SAVE AND RELOAD QC PASS\n")





# ============================================================
# STEP 5.8
# GO ENRICHMENT WITH EXPANDED BACKGROUND
# ============================================================

cat("\n========================================\n")
cat("STEP 5.8: EXPANDED GO ENRICHMENT\n")
cat("========================================\n")

# 1. Завантаження наборів генів
enrichment_gene_sets_expanded <- readRDS(
  "results/models/TCGA_SKCM_enrichment_gene_sets_expanded.rds"
)

go_genes <- unique(
  as.character(enrichment_gene_sets_expanded$robust_entrez)
)

go_background <- unique(
  as.character(enrichment_gene_sets_expanded$background_entrez)
)

cat("\n--- INPUT GENE SETS ---\n")
cat("Robust genes:", length(go_genes), "\n")
cat("Expanded background:", length(go_background), "\n")

# 2. Перевірка пакетів
required_go_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi"
)

go_package_status <- vapply(
  required_go_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

cat("\n--- PACKAGE CHECK ---\n")
print(go_package_status)

stopifnot(all(go_package_status))

# 3. Вхідний QC
go_input_qc <- c(
  robust_gene_count_correct =
    length(go_genes) == 42,
  
  background_gene_count_correct =
    length(go_background) == 14206,
  
  robust_genes_unique =
    !anyDuplicated(go_genes),
  
  background_unique =
    !anyDuplicated(go_background),
  
  robust_subset_of_background =
    all(go_genes %in% go_background),
  
  no_missing_robust_IDs =
    !anyNA(go_genes),
  
  no_missing_background_IDs =
    !anyNA(go_background)
)

cat("\n--- GO INPUT QC ---\n")
print(go_input_qc)

stopifnot(all(go_input_qc))

# 4. Функція GO enrichment
run_expanded_go <- function(ontology) {
  
  clusterProfiler::enrichGO(
    gene = go_genes,
    universe = go_background,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 500,
    readable = TRUE
  )
}

# 5. Запуск BP, MF, CC
cat("\nRunning GO Biological Process...\n")
go_bp_expanded <- run_expanded_go("BP")

cat("Running GO Molecular Function...\n")
go_mf_expanded <- run_expanded_go("MF")

cat("Running GO Cellular Component...\n")
go_cc_expanded <- run_expanded_go("CC")

# 6. Перетворення результатів у таблиці
go_bp_expanded_df <- as.data.frame(go_bp_expanded)
go_mf_expanded_df <- as.data.frame(go_mf_expanded)
go_cc_expanded_df <- as.data.frame(go_cc_expanded)

# 7. Кількість термінів і значущих результатів
go_expanded_summary <- data.frame(
  Ontology = c("BP", "MF", "CC"),
  
  Total_Terms = c(
    nrow(go_bp_expanded_df),
    nrow(go_mf_expanded_df),
    nrow(go_cc_expanded_df)
  ),
  
  Significant_FDR_005 = c(
    sum(go_bp_expanded_df$p.adjust < 0.05, na.rm = TRUE),
    sum(go_mf_expanded_df$p.adjust < 0.05, na.rm = TRUE),
    sum(go_cc_expanded_df$p.adjust < 0.05, na.rm = TRUE)
  ),
  
  stringsAsFactors = FALSE
)

cat("\n--- EXPANDED GO SUMMARY ---\n")
print(go_expanded_summary, row.names = FALSE)

# 8. Топ-10 термінів за FDR
print_top_go <- function(df, ontology) {
  
  cat("\n--- TOP 10", ontology, "---\n")
  
  if (nrow(df) == 0) {
    cat("No GO terms returned.\n")
    return(invisible(NULL))
  }
  
  df <- df[order(df$p.adjust, df$pvalue), , drop = FALSE]
  
  columns_to_print <- intersect(
    c(
      "ID",
      "Description",
      "GeneRatio",
      "BgRatio",
      "Count",
      "pvalue",
      "p.adjust",
      "geneID"
    ),
    names(df)
  )
  
  print(
    head(df[, columns_to_print, drop = FALSE], 10),
    row.names = FALSE
  )
}

print_top_go(go_bp_expanded_df, "BP")
print_top_go(go_mf_expanded_df, "MF")
print_top_go(go_cc_expanded_df, "CC")

# 9. Збереження результатів
write.csv(
  go_bp_expanded_df,
  "results/tables/TCGA_SKCM_GO_BP_expanded_background.csv",
  row.names = FALSE
)

write.csv(
  go_mf_expanded_df,
  "results/tables/TCGA_SKCM_GO_MF_expanded_background.csv",
  row.names = FALSE
)

write.csv(
  go_cc_expanded_df,
  "results/tables/TCGA_SKCM_GO_CC_expanded_background.csv",
  row.names = FALSE
)

write.csv(
  go_expanded_summary,
  "results/tables/TCGA_SKCM_GO_expanded_background_summary.csv",
  row.names = FALSE
)

# 10. Збереження повного об'єкта
go_expanded_results <- list(
  BP = go_bp_expanded,
  MF = go_mf_expanded,
  CC = go_cc_expanded,
  summary = go_expanded_summary,
  robust_entrez = go_genes,
  background_entrez = go_background,
  parameters = list(
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 500,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  ),
  input_qc = go_input_qc
)

saveRDS(
  go_expanded_results,
  "results/models/TCGA_SKCM_GO_enrichment_expanded_background.rds"
)

# 11. Фінальний QC
go_expanded_files <- c(
  "results/tables/TCGA_SKCM_GO_BP_expanded_background.csv",
  "results/tables/TCGA_SKCM_GO_MF_expanded_background.csv",
  "results/tables/TCGA_SKCM_GO_CC_expanded_background.csv",
  "results/tables/TCGA_SKCM_GO_expanded_background_summary.csv",
  "results/models/TCGA_SKCM_GO_enrichment_expanded_background.rds"
)

go_expanded_save_qc <- c(
  all_files_exist = all(file.exists(go_expanded_files)),
  BP_result_created = !is.null(go_bp_expanded),
  MF_result_created = !is.null(go_mf_expanded),
  CC_result_created = !is.null(go_cc_expanded),
  correct_robust_gene_count = length(go_genes) == 42,
  correct_background_count = length(go_background) == 14206,
  input_qc_pass = all(go_input_qc)
)

cat("\n--- GO EXPANDED SAVE QC ---\n")
print(go_expanded_save_qc)

stopifnot(all(go_expanded_save_qc))

cat("\nSTEP 5.8: EXPANDED GO ENRICHMENT PASS\n")




# ============================================================
# STEP 5.9
# REACTOME ENRICHMENT WITH EXPANDED BACKGROUND
# ============================================================

cat("\n========================================\n")
cat("STEP 5.9: EXPANDED REACTOME ENRICHMENT\n")
cat("========================================\n")

# 1. Завантаження наборів генів
enrichment_gene_sets_expanded <- readRDS(
  "results/models/TCGA_SKCM_enrichment_gene_sets_expanded.rds"
)

reactome_genes <- unique(
  as.character(enrichment_gene_sets_expanded$robust_entrez)
)

reactome_background <- unique(
  as.character(enrichment_gene_sets_expanded$background_entrez)
)

cat("\n--- INPUT GENE SETS ---\n")
cat("Robust genes:", length(reactome_genes), "\n")
cat("Expanded background:", length(reactome_background), "\n")

# 2. Перевірка пакетів
required_reactome_packages <- c(
  "ReactomePA",
  "org.Hs.eg.db",
  "AnnotationDbi"
)

reactome_package_status <- vapply(
  required_reactome_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

cat("\n--- PACKAGE CHECK ---\n")
print(reactome_package_status)

stopifnot(all(reactome_package_status))

# 3. Вхідний QC
reactome_input_qc <- c(
  robust_gene_count_correct =
    length(reactome_genes) == 42,
  
  background_gene_count_correct =
    length(reactome_background) == 14206,
  
  robust_genes_unique =
    !anyDuplicated(reactome_genes),
  
  background_unique =
    !anyDuplicated(reactome_background),
  
  robust_subset_of_background =
    all(reactome_genes %in% reactome_background),
  
  no_missing_robust_IDs =
    !anyNA(reactome_genes),
  
  no_missing_background_IDs =
    !anyNA(reactome_background)
)

cat("\n--- REACTOME INPUT QC ---\n")
print(reactome_input_qc)

stopifnot(all(reactome_input_qc))

# 4. Reactome enrichment
cat("\nRunning Reactome pathway enrichment...\n")

reactome_expanded <- ReactomePA::enrichPathway(
  gene = reactome_genes,
  universe = reactome_background,
  organism = "human",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  minGSSize = 5,
  maxGSSize = 500,
  readable = TRUE
)

# 5. Перетворення результату в таблицю
reactome_expanded_df <- as.data.frame(
  reactome_expanded
)

cat("\nTotal Reactome pathways returned:",
    nrow(reactome_expanded_df), "\n")

# 6. Підсумок
reactome_expanded_summary <- data.frame(
  Analysis = "Expanded background",
  Robust_Genes = length(reactome_genes),
  Background_Genes = length(reactome_background),
  Total_Pathways = nrow(reactome_expanded_df),
  Significant_FDR_005 = sum(
    reactome_expanded_df$p.adjust < 0.05,
    na.rm = TRUE
  ),
  stringsAsFactors = FALSE
)

cat("\n--- EXPANDED REACTOME SUMMARY ---\n")
print(reactome_expanded_summary, row.names = FALSE)

# 7. Топ-10 шляхів
if (nrow(reactome_expanded_df) > 0) {
  
  reactome_expanded_df <- reactome_expanded_df[
    order(
      reactome_expanded_df$p.adjust,
      reactome_expanded_df$pvalue
    ),
    ,
    drop = FALSE
  ]
  
  reactome_top_columns <- intersect(
    c(
      "ID",
      "Description",
      "GeneRatio",
      "BgRatio",
      "Count",
      "pvalue",
      "p.adjust",
      "geneID"
    ),
    names(reactome_expanded_df)
  )
  
  cat("\n--- TOP 10 REACTOME PATHWAYS ---\n")
  
  print(
    head(
      reactome_expanded_df[
        ,
        reactome_top_columns,
        drop = FALSE
      ],
      10
    ),
    row.names = FALSE
  )
  
} else {
  cat("\nNo Reactome pathways returned.\n")
}

# 8. Окрема таблиця значущих шляхів
reactome_significant_expanded <- reactome_expanded_df[
  !is.na(reactome_expanded_df$p.adjust) &
    reactome_expanded_df$p.adjust < 0.05,
  ,
  drop = FALSE
]

cat("\n--- SIGNIFICANT REACTOME PATHWAYS ---\n")
cat("Number of significant pathways:",
    nrow(reactome_significant_expanded), "\n")

if (nrow(reactome_significant_expanded) > 0) {
  print(
    reactome_significant_expanded[
      ,
      reactome_top_columns,
      drop = FALSE
    ],
    row.names = FALSE
  )
}

# 9. Збереження таблиць
write.csv(
  reactome_expanded_df,
  "results/tables/TCGA_SKCM_Reactome_expanded_background.csv",
  row.names = FALSE
)

write.csv(
  reactome_significant_expanded,
  "results/tables/TCGA_SKCM_Reactome_significant_expanded_background.csv",
  row.names = FALSE
)

write.csv(
  reactome_expanded_summary,
  "results/tables/TCGA_SKCM_Reactome_expanded_background_summary.csv",
  row.names = FALSE
)

# 10. Збереження повного об'єкта
reactome_expanded_results <- list(
  results = reactome_expanded,
  results_table = reactome_expanded_df,
  significant_pathways = reactome_significant_expanded,
  summary = reactome_expanded_summary,
  robust_entrez = reactome_genes,
  background_entrez = reactome_background,
  parameters = list(
    organism = "human",
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 500,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  ),
  input_qc = reactome_input_qc
)

saveRDS(
  reactome_expanded_results,
  "results/models/TCGA_SKCM_Reactome_enrichment_expanded_background.rds"
)

# 11. Фінальний QC
reactome_expanded_files <- c(
  "results/tables/TCGA_SKCM_Reactome_expanded_background.csv",
  "results/tables/TCGA_SKCM_Reactome_significant_expanded_background.csv",
  "results/tables/TCGA_SKCM_Reactome_expanded_background_summary.csv",
  "results/models/TCGA_SKCM_Reactome_enrichment_expanded_background.rds"
)

reactome_expanded_save_qc <- c(
  all_files_exist =
    all(file.exists(reactome_expanded_files)),
  
  result_created =
    !is.null(reactome_expanded),
  
  correct_robust_gene_count =
    length(reactome_genes) == 42,
  
  correct_background_count =
    length(reactome_background) == 14206,
  
  input_qc_pass =
    all(reactome_input_qc),
  
  significant_table_consistent =
    nrow(reactome_significant_expanded) ==
    sum(reactome_expanded_df$p.adjust < 0.05, na.rm = TRUE)
)

cat("\n--- REACTOME EXPANDED SAVE QC ---\n")
print(reactome_expanded_save_qc)

stopifnot(all(reactome_expanded_save_qc))

cat("\nSTEP 5.9: EXPANDED REACTOME ENRICHMENT PASS\n")


# ============================================================
# STEP 5.10A
# REACTOME BACKGROUND COMPARISON AUDIT
# ============================================================

cat("\n========================================\n")
cat("STEP 5.10A: REACTOME BACKGROUND AUDIT\n")
cat("========================================\n")

# 1. Завантаження розширеного результату
reactome_expanded_loaded <- readRDS(
  "results/models/TCGA_SKCM_Reactome_enrichment_expanded_background.rds"
)

# 2. Перевірка структури вузького об'єкта
cat("\n--- NARROW OBJECT ---\n")

cat("Class:", class(reactome_narrow), "\n")
cat("Input genes:", length(reactome_narrow@gene), "\n")
cat("Universe genes:", length(reactome_narrow@universe), "\n")
cat("Pathways:", nrow(as.data.frame(reactome_narrow)), "\n")

# 3. Перевірка структури розширеного об'єкта
cat("\n--- EXPANDED OBJECT ---\n")

cat("Class:", class(reactome_expanded_loaded$results), "\n")
cat("Input genes:", length(reactome_expanded_loaded$results@gene), "\n")
cat("Universe genes:", length(reactome_expanded_loaded$results@universe), "\n")
cat("Pathways:", nrow(as.data.frame(
  reactome_expanded_loaded$results
)), "\n")

# 4. Порівняння вхідних наборів
cat("\n--- INPUT GENE COUNTS ---\n")

cat(
  "Narrow input genes:",
  length(reactome_narrow@gene),
  "\n"
)

cat(
  "Expanded input genes:",
  length(reactome_expanded_loaded$results@gene),
  "\n"
)

cat(
  "Narrow universe:",
  length(reactome_narrow@universe),
  "\n"
)

cat(
  "Expanded universe:",
  length(reactome_expanded_loaded$results@universe),
  "\n"
)

# 5. Перевірка перетину effective universes
narrow_universe <- unique(
  as.character(reactome_narrow@universe)
)

expanded_universe <- unique(
  as.character(reactome_expanded_loaded$results@universe)
)

cat("\n--- UNIVERSE OVERLAP ---\n")

cat(
  "Common genes:",
  length(intersect(narrow_universe, expanded_universe)),
  "\n"
)

cat(
  "Narrow-only genes:",
  length(setdiff(narrow_universe, expanded_universe)),
  "\n"
)

cat(
  "Expanded-only genes:",
  length(setdiff(expanded_universe, narrow_universe)),
  "\n"
)

# 6. Перевірка GeneRatio та BgRatio
cat("\n--- NARROW TOP PATHWAYS ---\n")

narrow_df <- as.data.frame(reactome_narrow)

if (nrow(narrow_df) > 0) {
  print(
    head(
      narrow_df[
        order(narrow_df$p.adjust, narrow_df$pvalue),
        intersect(
          c(
            "ID", "Description", "GeneRatio",
            "BgRatio", "Count", "pvalue",
            "p.adjust", "geneID"
          ),
          names(narrow_df)
        ),
        drop = FALSE
      ],
      10
    ),
    row.names = FALSE
  )
}

cat("\n--- EXPANDED TOP PATHWAYS ---\n")

expanded_df <- as.data.frame(
  reactome_expanded_loaded$results
)

if (nrow(expanded_df) > 0) {
  print(
    head(
      expanded_df[
        order(expanded_df$p.adjust, expanded_df$pvalue),
        intersect(
          c(
            "ID", "Description", "GeneRatio",
            "BgRatio", "Count", "pvalue",
            "p.adjust", "geneID"
          ),
          names(expanded_df)
        ),
        drop = FALSE
      ],
      10
    ),
    row.names = FALSE
  )
}

# 7. Перевірка значущих шляхів
cat("\n--- SIGNIFICANT PATHWAYS ---\n")

cat(
  "Narrow FDR < 0.05:",
  sum(narrow_df$p.adjust < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Expanded FDR < 0.05:",
  sum(expanded_df$p.adjust < 0.05, na.rm = TRUE),
  "\n"
)




# ============================================================
# STEP 5.10B
# FINAL GO + REACTOME BACKGROUND COMPARISON
# ============================================================

cat("\n========================================\n")
cat("STEP 5.10B: FINAL ENRICHMENT COMPARISON\n")
cat("========================================\n")

# 1. Завантаження результатів
go_narrow_results <- readRDS(
  "results/models/TCGA_SKCM_GO_enrichment_results.rds"
)

go_expanded_results <- readRDS(
  "results/models/TCGA_SKCM_GO_enrichment_expanded_background.rds"
)

reactome_narrow <- readRDS(
  "results/models/TCGA_SKCM_Reactome_enrichment_results.rds"
)

reactome_expanded_results <- readRDS(
  "results/models/TCGA_SKCM_Reactome_enrichment_expanded_background.rds"
)

# 2. Функція отримання метрик
get_enrichment_metrics <- function(object, analysis, ontology) {
  
  result_df <- as.data.frame(object)
  
  data.frame(
    Analysis = analysis,
    Ontology = ontology,
    Total_Terms = nrow(result_df),
    Significant_FDR_005 = sum(
      result_df$p.adjust < 0.05,
      na.rm = TRUE
    ),
    Best_P_value = if (nrow(result_df) > 0) {
      min(result_df$pvalue, na.rm = TRUE)
    } else {
      NA_real_
    },
    Best_FDR = if (nrow(result_df) > 0) {
      min(result_df$p.adjust, na.rm = TRUE)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

# 3. GO comparison
go_comparison <- do.call(
  rbind,
  list(
    get_enrichment_metrics(
      go_narrow_results$BP,
      "Narrow",
      "GO_BP"
    ),
    get_enrichment_metrics(
      go_narrow_results$MF,
      "Narrow",
      "GO_MF"
    ),
    get_enrichment_metrics(
      go_narrow_results$CC,
      "Narrow",
      "GO_CC"
    ),
    get_enrichment_metrics(
      go_expanded_results$BP,
      "Expanded",
      "GO_BP"
    ),
    get_enrichment_metrics(
      go_expanded_results$MF,
      "Expanded",
      "GO_MF"
    ),
    get_enrichment_metrics(
      go_expanded_results$CC,
      "Expanded",
      "GO_CC"
    )
  )
)

cat("\n--- GO COMPARISON ---\n")
print(go_comparison, row.names = FALSE)

# 4. Reactome comparison
reactome_comparison <- rbind(
  get_enrichment_metrics(
    reactome_narrow,
    "Narrow",
    "Reactome"
  ),
  get_enrichment_metrics(
    reactome_expanded_results$results,
    "Expanded",
    "Reactome"
  )
)

cat("\n--- REACTOME COMPARISON ---\n")
print(reactome_comparison, row.names = FALSE)

# 5. Effective universe comparison
effective_universe_comparison <- data.frame(
  Analysis = c("Narrow", "Expanded"),
  Input_Background = c(85, 14206),
  Reactome_Effective_Universe = c(
    length(reactome_narrow@universe),
    length(reactome_expanded_results$results@universe)
  ),
  Reactome_Input_Genes = c(
    length(reactome_narrow@gene),
    length(reactome_expanded_results$results@gene)
  ),
  stringsAsFactors = FALSE
)

cat("\n--- EFFECTIVE REACTOME UNIVERSE ---\n")
print(effective_universe_comparison, row.names = FALSE)

# 6. Перевірка однаковості вхідних robust genes
narrow_robust <- sort(
  as.character(reactome_narrow@gene)
)

expanded_robust <- sort(
  as.character(reactome_expanded_results$results@gene)
)

robust_gene_overlap <- c(
  same_robust_gene_set =
    identical(narrow_robust, expanded_robust),
  
  narrow_robust_count =
    length(narrow_robust) == 42,
  
  expanded_robust_count =
    length(expanded_robust) == 42
)

cat("\n--- ROBUST GENE SET QC ---\n")
print(robust_gene_overlap)

# 7. Загальний підсумок значущих результатів
final_enrichment_summary <- data.frame(
  Analysis = c(
    "GO narrow",
    "GO expanded",
    "Reactome narrow",
    "Reactome expanded"
  ),
  Significant_Results = c(
    sum(go_comparison$Significant_FDR_005[
      go_comparison$Analysis == "Narrow"
    ]),
    sum(go_comparison$Significant_FDR_005[
      go_comparison$Analysis == "Expanded"
    ]),
    reactome_comparison$Significant_FDR_005[1],
    reactome_comparison$Significant_FDR_005[2]
  )
)

cat("\n--- FINAL ENRICHMENT SUMMARY ---\n")
print(final_enrichment_summary, row.names = FALSE)

# 8. Збереження таблиць
write.csv(
  go_comparison,
  "results/tables/TCGA_SKCM_GO_narrow_vs_expanded_comparison.csv",
  row.names = FALSE
)

write.csv(
  reactome_comparison,
  "results/tables/TCGA_SKCM_Reactome_narrow_vs_expanded_comparison.csv",
  row.names = FALSE
)

write.csv(
  effective_universe_comparison,
  "results/tables/TCGA_SKCM_Reactome_effective_universe_comparison.csv",
  row.names = FALSE
)

write.csv(
  final_enrichment_summary,
  "results/tables/TCGA_SKCM_final_enrichment_summary.csv",
  row.names = FALSE
)

# 9. Фінальний QC
comparison_files <- c(
  "results/tables/TCGA_SKCM_GO_narrow_vs_expanded_comparison.csv",
  "results/tables/TCGA_SKCM_Reactome_narrow_vs_expanded_comparison.csv",
  "results/tables/TCGA_SKCM_Reactome_effective_universe_comparison.csv",
  "results/tables/TCGA_SKCM_final_enrichment_summary.csv"
)

final_enrichment_qc <- c(
  all_comparison_files_exist =
    all(file.exists(comparison_files)),
  
  same_robust_gene_set =
    robust_gene_overlap["same_robust_gene_set"],
  
  correct_narrow_robust_count =
    robust_gene_overlap["narrow_robust_count"],
  
  correct_expanded_robust_count =
    robust_gene_overlap["expanded_robust_count"],
  
  no_significant_narrow_reactome =
    reactome_comparison$Significant_FDR_005[1] == 0,
  
  no_significant_expanded_reactome =
    reactome_comparison$Significant_FDR_005[2] == 0
)

cat("\n--- FINAL ENRICHMENT QC ---\n")
print(final_enrichment_qc)

stopifnot(all(final_enrichment_qc))

cat("\nSTEP 5.10B: FINAL ENRICHMENT COMPARISON PASS\n")




# ============================================================
# STEP 5.11 — GLOBAL BH CORRECTION FOR ALL GO ONTOLOGIES
# TCGA-SKCM | Narrow and expanded backgrounds
# ============================================================

cat("\n===== STEP 5.11: GLOBAL GO BH CORRECTION =====\n")

# 1. Перевірка необхідних пакетів і каталогів
required_packages <- c("clusterProfiler")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Відсутні пакети: ",
    paste(missing_packages, collapse = ", ")
  )
}

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/models", recursive = TRUE, showWarnings = FALSE)

# 2. Завантаження збережених результатів GO
narrow_file <- "results/models/TCGA_SKCM_GO_enrichment_results.rds"

expanded_file <- "results/models/TCGA_SKCM_GO_enrichment_expanded_background.rds"

if (!file.exists(narrow_file)) {
  stop("Не знайдено файл вузького GO-аналізу: ", narrow_file)
}

if (!file.exists(expanded_file)) {
  stop("Не знайдено файл розширеного GO-аналізу: ", expanded_file)
}

go_narrow_results <- readRDS(narrow_file)
go_expanded_results <- readRDS(expanded_file)

cat("\nЗавантаження результатів:\n")
cat("Narrow:", class(go_narrow_results), "\n")
cat("Expanded:", class(go_expanded_results), "\n")

# 3. Функція для об'єднання BP, MF, CC
#    та глобальної BH-корекції
perform_global_go_bh <- function(go_results, analysis_name) {
  
  required_ontologies <- c("BP", "MF", "CC")
  
  if (!all(required_ontologies %in% names(go_results))) {
    stop(
      analysis_name,
      ": у списку мають бути компоненти BP, MF і CC. Наявні: ",
      paste(names(go_results), collapse = ", ")
    )
  }
  
  ontology_tables <- lapply(required_ontologies, function(ontology) {
    
    result_object <- go_results[[ontology]]
    
    if (is.null(result_object)) {
      stop(analysis_name, ": відсутній результат для ", ontology)
    }
    
    result_df <- as.data.frame(result_object)
    
    if (nrow(result_df) == 0) {
      return(NULL)
    }
    
    required_columns <- c("ID", "Description", "pvalue")
    
    if (!all(required_columns %in% colnames(result_df))) {
      stop(
        analysis_name,
        ": у результаті ",
        ontology,
        " відсутні необхідні колонки."
      )
    }
    
    result_df$Ontology <- ontology
    
    result_df
  })
  
  # Об'єднання всіх трьох онтологій
  combined_df <- do.call(rbind, ontology_tables)
  
  if (is.null(combined_df) || nrow(combined_df) == 0) {
    stop(analysis_name, ": немає GO-термінів для корекції.")
  }
  
  rownames(combined_df) <- NULL
  
  # Перевірка p-value
  valid_p <- !is.na(combined_df$pvalue) &
    is.finite(combined_df$pvalue)
  
  if (any(combined_df$pvalue[valid_p] < 0 |
          combined_df$pvalue[valid_p] > 1)) {
    stop(analysis_name, ": виявлено p-value поза межами [0, 1].")
  }
  
  # Глобальна BH-корекція:
  # усі BP + MF + CC в одному наборі
  combined_df$Global_FDR_BH <- NA_real_
  
  combined_df$Global_FDR_BH[valid_p] <- p.adjust(
    combined_df$pvalue[valid_p],
    method = "BH"
  )
  
  combined_df$Global_FDR_Significant <- (
    !is.na(combined_df$Global_FDR_BH) &
      combined_df$Global_FDR_BH < 0.05
  )
  
  combined_df$Original_FDR_Significant <- (
    !is.na(combined_df$p.adjust) &
      combined_df$p.adjust < 0.05
  )
  
  # Сортування за глобальним FDR
  combined_df <- combined_df[
    order(combined_df$Global_FDR_BH, combined_df$pvalue),
    ,
    drop = FALSE
  ]
  
  rownames(combined_df) <- NULL
  
  # 4. Підсумок окремо для кожної онтології
  summary_df <- do.call(
    rbind,
    lapply(required_ontologies, function(ontology) {
      
      subset_df <- combined_df[
        combined_df$Ontology == ontology,
        ,
        drop = FALSE
      ]
      
      data.frame(
        Analysis = analysis_name,
        Ontology = ontology,
        Total_Terms = nrow(subset_df),
        Original_FDR_Significant = sum(
          subset_df$Original_FDR_Significant,
          na.rm = TRUE
        ),
        Global_FDR_Significant = sum(
          subset_df$Global_FDR_Significant,
          na.rm = TRUE
        ),
        Best_P_value = if (nrow(subset_df) > 0) {
          min(subset_df$pvalue, na.rm = TRUE)
        } else {
          NA_real_
        },
        Best_Global_FDR = if (any(!is.na(subset_df$Global_FDR_BH))) {
          min(subset_df$Global_FDR_BH, na.rm = TRUE)
        } else {
          NA_real_
        }
      )
    })
  )
  
  rownames(summary_df) <- NULL
  
  list(
    analysis = analysis_name,
    results = combined_df,
    summary = summary_df,
    total_terms = nrow(combined_df),
    valid_p_values = sum(valid_p),
    missing_p_values = sum(!valid_p)
  )
}

# 5. Глобальна корекція для вузького фону
global_go_narrow <- perform_global_go_bh(
  go_results = go_narrow_results,
  analysis_name = "Narrow_background"
)

# 6. Глобальна корекція для розширеного фону
global_go_expanded <- perform_global_go_bh(
  go_results = go_expanded_results,
  analysis_name = "Expanded_background"
)

# 7. Перевірка TAP binding у розширеному аналізі
tap_binding <- global_go_expanded$results[
  global_go_expanded$results$ID == "GO:0046977",
  ,
  drop = FALSE
]

cat("\n===== TAP BINDING — GLOBAL BH RESULT =====\n")

if (nrow(tap_binding) == 0) {
  
  cat("GO:0046977 TAP binding не знайдено.\n")
  
} else {
  
  print(
    tap_binding[
      ,
      c(
        "ID",
        "Description",
        "Ontology",
        "GeneRatio",
        "BgRatio",
        "pvalue",
        "p.adjust",
        "Global_FDR_BH",
        "Original_FDR_Significant",
        "Global_FDR_Significant",
        "geneID"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )
}

# 8. Загальні підсумки
global_go_summary <- rbind(
  global_go_narrow$summary,
  global_go_expanded$summary
)

cat("\n===== GLOBAL GO SUMMARY =====\n")
print(global_go_summary, row.names = FALSE)

cat("\n===== TOTAL NUMBER OF TERMS =====\n")
cat(
  "Narrow:",
  global_go_narrow$total_terms,
  "\n"
)

cat(
  "Expanded:",
  global_go_expanded$total_terms,
  "\n"
)

# 9. Збереження повних таблиць
write.csv(
  global_go_narrow$results,
  "results/tables/TCGA_SKCM_GO_narrow_global_BH.csv",
  row.names = FALSE
)

write.csv(
  global_go_expanded$results,
  "results/tables/TCGA_SKCM_GO_expanded_global_BH.csv",
  row.names = FALSE
)

write.csv(
  global_go_summary,
  "results/tables/TCGA_SKCM_GO_global_BH_summary.csv",
  row.names = FALSE
)

# 10. Збереження RDS
global_go_bh_results <- list(
  narrow = global_go_narrow,
  expanded = global_go_expanded,
  summary = global_go_summary,
  method = list(
    correction = "Benjamini-Hochberg",
    scope = "All GO BP, MF and CC terms combined within each background",
    significance_threshold = 0.05,
    interpretation = paste(
      "Global correction is a sensitivity analysis.",
      "Original ontology-specific adjusted p-values are retained."
    )
  )
)

saveRDS(
  global_go_bh_results,
  "results/models/TCGA_SKCM_GO_global_BH_results.rds"
)

# 11. Фінальний QC збережених файлів
output_files <- c(
  "results/tables/TCGA_SKCM_GO_narrow_global_BH.csv",
  "results/tables/TCGA_SKCM_GO_expanded_global_BH.csv",
  "results/tables/TCGA_SKCM_GO_global_BH_summary.csv",
  "results/models/TCGA_SKCM_GO_global_BH_results.rds"
)

global_bh_save_qc <- file.exists(output_files)

cat("\n===== SAVE QC =====\n")

print(
  data.frame(
    File = output_files,
    Exists = global_bh_save_qc
  ),
  row.names = FALSE
)

# 12. Контрольні перевірки
global_bh_qc <- c(
  narrow_terms_present =
    global_go_narrow$total_terms > 0,
  
  expanded_terms_present =
    global_go_expanded$total_terms > 0,
  
  all_three_ontologies_present =
    all(c("BP", "MF", "CC") %in%
          global_go_expanded$results$Ontology),
  
  global_fdr_values_valid =
    all(
      is.na(global_go_expanded$results$Global_FDR_BH) |
        (
          global_go_expanded$results$Global_FDR_BH >= 0 &
            global_go_expanded$results$Global_FDR_BH <= 1
        )
    ),
  
  tap_binding_found =
    nrow(tap_binding) == 1,
  
  all_outputs_saved =
    all(global_bh_save_qc)
)

cat("\n===== GLOBAL BH QC =====\n")
print(global_bh_qc)

if (all(global_bh_qc)) {
  cat("\nGLOBAL GO BH ANALYSIS: PASS\n")
} else {
  cat("\nGLOBAL GO BH ANALYSIS: CHECK FAILED ITEMS\n")
}


