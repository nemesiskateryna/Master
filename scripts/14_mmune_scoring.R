# 1. Match Ensembl IDs to RNA-seq matrix

idx <- match(
  unique_gene_rows$ENSEMBL_original,
  rownames(expr.mat)
)

# 2. Verify matching

cat("Total annotation rows:", nrow(unique_gene_rows), "\n")
cat("Unmatched Ensembl IDs:", sum(is.na(idx)), "\n")

stopifnot(!anyNA(idx))

# 3. Create Gene Symbol expression matrix

expr_symbol <- expr.mat[idx, , drop = FALSE]

rownames(expr_symbol) <- unique_gene_rows$Gene_Symbol

# 4. Verify matrix

cat("\n--- GENE SYMBOL MATRIX QC ---\n")

cat("Genes:", nrow(expr_symbol), "\n")
cat("Samples:", ncol(expr_symbol), "\n")
cat("Duplicated Gene Symbols:", sum(duplicated(rownames(expr_symbol))), "\n")

cat(
  "Sample order preserved:",
  identical(colnames(expr_symbol), colnames(expr.mat)),
  "\n"
)

stopifnot(!anyDuplicated(rownames(expr_symbol)))

stopifnot(
  identical(
    colnames(expr_symbol),
    colnames(expr.mat)
  )
)
# Define immune signatures

cyt_genes <- c("GZMA", "PRF1")

cd8_genes <- c("CD8A", "CD8B")

ifng_genes <- c(
  "IFNG", "CXCL9", "CXCL10", "STAT1", "IDO1"
)

tis_genes <- c(
  "CCL5", "CD27", "CD274", "CD276",
  "CD8A", "CMKLR1", "CXCL9", "CXCR6",
  "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1",
  "LAG3", "NKG7", "PDCD1LG2", "PSMB10",
  "STAT1", "TIGIT"
)

# Combine signatures

signatures <- list(
  CYT = cyt_genes,
  CD8 = cd8_genes,
  IFNG = ifng_genes,
  TIS = tis_genes
)

# Check gene coverage

signature_coverage <- do.call(
  rbind,
  lapply(names(signatures), function(nm) {
    
    genes <- signatures[[nm]]
    
    present <- intersect(
      genes,
      rownames(expr_symbol)
    )
    
    missing <- setdiff(
      genes,
      rownames(expr_symbol)
    )
    
    data.frame(
      Signature = nm,
      Total_genes = length(genes),
      Present_genes = length(present),
      Missing_genes = paste(missing, collapse = "; "),
      Coverage_percent = round(
        100 * length(present) / length(genes),
        2
      ),
      Complete = length(missing) == 0,
      stringsAsFactors = FALSE
    )
  })
)

# Display results

cat("\n--- IMMUNE SIGNATURE COVERAGE ---\n")

print(
  signature_coverage,
  row.names = FALSE
)




# 1. Check expression matrix

stopifnot(
  all(unlist(signatures) %in% rownames(expr_symbol))
)

# 2. Check gene variability

gene_sd <- apply(
  expr_symbol,
  1,
  sd
)

if (any(!is.finite(gene_sd) | gene_sd == 0)) {
  stop("Some genes have zero or non-finite variance.")
}

# Save immune scores

saveRDS(
  immune_scores,
  "data/TCGA/processed/TCGA_SKCM_immune_scores.rds"
)

write.csv(
  immune_scores,
  "results/tables/TCGA_SKCM_immune_scores.csv",
  row.names = FALSE
)



