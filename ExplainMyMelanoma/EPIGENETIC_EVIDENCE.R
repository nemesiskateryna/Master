
###############################################################################
# SCRIPT 34 — CREATE EPIGENETIC EVIDENCE MODULE
#
# Creates:
#   app/ExplainMyMelanoma/R/05_epigenetic_evidence_module.R
#
# Reads existing CSV and PNG results only.
# Does NOT modify app.R or existing analysis outputs.
###############################################################################

cat("\n========== SCRIPT 34 ==========\n")

# 1. CHECK REQUIRED PATHS ------------------------------------------------------

stopifnot(exists("app_dir"))
stopifnot(exists("tables_dir"))
stopifnot(exists("figures_dir"))

module_dir <- file.path(app_dir, "R")
module_file <- file.path(
  module_dir,
  "05_epigenetic_evidence_module.R"
)

if (!dir.exists(module_dir)) {
  stop("Module directory does not exist: ", module_dir)
}

if (file.exists(module_file)) {
  stop(
    "STOP: module already exists. ",
    "No file was overwritten: ",
    module_file
  )
}

# 2. DEFINE MODULE -------------------------------------------------------------

module_code <- r"(
epigeneticEvidenceUI <- function(id) {

  ns <- shiny::NS(id)

  shiny::tagList(

    # GEO EXTERNAL VALIDATION
    shiny::h4("GEO external CpG validation"),

    shiny::p(
      "GSE144487: external assessment of the TCGA-derived ",
      "DNA methylation signature."
    ),

    shiny::tableOutput(ns("geo_validation_table")),

    shiny::div(
      class = "scientific-disclaimer",
      shiny::strong("Interpretation"),
      shiny::p(
        "The CpG signature shows substantial overlap and ",
        "partial directional concordance in GEO. However, ",
        "the permutation results do not reach the conventional ",
        "0.05 significance threshold. This is evidence of ",
        "a biological association, not definitive global replication."
      )
    ),

    shiny::hr(),

    # CpG-GENE ASSOCIATIONS
    shiny::h4("CpG–gene expression associations"),

    shiny::p(
      "Selected direct CpG–gene associations from the ",
      "TCGA-SKCM analysis."
    ),

    shiny::tableOutput(ns("cpg_gene_table")),

    shiny::div(
      class = "evidence-figure",
      shiny::h5("cg16890093 — TAP1"),
      shiny::imageOutput(ns("tap1_plot"), height = "auto")
    ),

    shiny::div(
      class = "evidence-figure",
      shiny::h5("cg07016276 — HLA-F"),
      shiny::imageOutput(ns("hlaf_plot"), height = "auto")
    ),

    shiny::hr(),

    # GEO IMMUNE SIGNATURES
    shiny::h4("GSE65904: immune signature analysis"),

    shiny::p(
      "Association of the TIS-17 gene signature with ",
      "auxiliary immune-related expression signatures."
    ),

    shiny::tableOutput(ns("tis_table")),

    shiny::div(
      class = "scientific-disclaimer",
      shiny::strong("Cohort matching limitation"),
      shiny::p(
        "The available matching audit labels all 74 records ",
        "as assumed matches based on numeric sample suffixes. ",
        "These matches should not be interpreted as independently ",
        "verified patient-level matches."
      )
    ),

    shiny::hr(),

    # GEO EPIGENETIC DISTRIBUTIONS
    shiny::h4("GSE144487: CpG annotation overview"),

    shiny::p(
      "Genomic-region and CpG-island distributions ",
      "of the overlapping CpG signature."
    ),

    shiny::div(
      class = "evidence-figure",
      shiny::h5("Genomic region distribution"),
      shiny::imageOutput(ns("region_plot"), height = "auto")
    ),

    shiny::div(
      class = "evidence-figure",
      shiny::h5("CpG-island distribution"),
      shiny::imageOutput(ns("island_plot"), height = "auto")
    ),

    shiny::p(
      class = "scientific-disclaimer",
      "These analyses describe associations and annotation patterns. ",
      "They do not establish that DNA methylation directly causes ",
      "changes in gene expression or immune phenotype."
    )
  )
}


epigeneticEvidenceServer <- function(
  id,
  tables_dir,
  figures_dir
) {

  shiny::moduleServer(id, function(input, output, session) {

    # Safe CSV reader
    read_result <- function(filename) {

      path <- file.path(tables_dir, filename)

      if (!file.exists(path)) {
        warning("Evidence table not found: ", path)
        return(NULL)
      }

      tryCatch(
        read.csv(
          path,
          stringsAsFactors = FALSE,
          check.names = FALSE
        ),
        error = function(e) {
          warning(
            "Could not read ", filename, ": ",
            conditionMessage(e)
          )
          NULL
        }
      )
    }


    # GEO VALIDATION TABLE
    output$geo_validation_table <- shiny::renderTable({

      dat <- read_result(
        "GEO_external_validation_final_summary.csv"
      )

      if (is.null(dat)) return(NULL)

      dat <- dat[
        dat$Section %in% c(
          "CpG overlap",
          "Directional validation",
          "Permutation validation"
        ),
        c("Section", "Metric", "Value"),
        drop = FALSE
      ]

      dat

    }, striped = TRUE, bordered = FALSE, spacing = "s")


    # CpG-GENE ASSOCIATION TABLE
    output$cpg_gene_table <- shiny::renderTable({

      dat <- read_result(
        "TCGA_SKCM_CpG_gene_correlation_robust_pairs.csv"
      )

      if (is.null(dat)) return(NULL)

      required <- c(
        "CpG_ID", "Gene", "Pooled_rho",
        "Pooled_FDR", "Primary_rho", "Metastatic_rho"
      )

      if (!all(required %in% names(dat))) {
        warning("Unexpected CpG-gene table structure.")
        return(NULL)
      }

      selected <- dat[
        (dat$CpG_ID == "cg16890093" & dat$Gene == "TAP1") |
          (dat$CpG_ID == "cg07016276" & dat$Gene == "HLA-F"),
        required,
        drop = FALSE
      ]

      selected

    }, digits = 4, striped = TRUE, bordered = FALSE, spacing = "s")


    # TIS-17 CORRELATIONS
    output$tis_table <- shiny::renderTable({

      dat <- read_result(
        "GEO_TIS17_auxiliary_correlations.csv"
      )

      if (is.null(dat)) return(NULL)

      dat

    }, digits = 4, striped = TRUE, bordered = FALSE, spacing = "s")


    # SAFE IMAGE RENDERING
    render_evidence_image <- function(filename, description) {

      image_path <- file.path(figures_dir, filename)

      if (!file.exists(image_path)) {
        warning("Evidence figure not found: ", image_path)
        return(NULL)
      }

      list(
        src = normalizePath(image_path),
        contentType = "image/png",
        alt = description,
        width = "100%"
      )
    }


    output$tap1_plot <- shiny::renderImage({
      render_evidence_image(
        "TCGA_SKCM_cg16890093_TAP1_correlation_primary_v2.png",
        "Correlation between cg16890093 methylation and TAP1 expression"
      )
    }, deleteFile = FALSE)


    output$hlaf_plot <- shiny::renderImage({
      render_evidence_image(
        "TCGA_SKCM_cg07016276_HLA-F_correlation_primary_v2.png",
        "Correlation between cg07016276 methylation and HLA-F expression"
      )
    }, deleteFile = FALSE)


    output$region_plot <- shiny::renderImage({
      render_evidence_image(
        "GSE144487_CpG_gene_region_distribution.png",
        "Genomic region distribution of GSE144487 CpGs"
      )
    }, deleteFile = FALSE)


    output$island_plot <- shiny::renderImage({
      render_evidence_image(
        "GSE144487_CpG_island_distribution.png",
        "CpG island distribution in GSE144487"
      )
    }, deleteFile = FALSE)

  })
}
)"

# 3. WRITE NEW MODULE ---------------------------------------------------------

writeLines(module_code, module_file)

cat("\nModule created:\n", module_file, "\n")

# 4. VERIFY SYNTAX ------------------------------------------------------------

syntax_ok <- tryCatch({
  parse(module_file)
  TRUE
}, error = function(e) {
  cat("Parse error:", conditionMessage(e), "\n")
  FALSE
})

cat("Module syntax valid:", syntax_ok, "\n")

# 5. VERIFY REQUIRED SOURCE TABLES AND FIGURES --------------------------------

required_tables <- c(
  "GEO_external_validation_final_summary.csv",
  "TCGA_SKCM_CpG_gene_correlation_robust_pairs.csv",
  "GEO_TIS17_auxiliary_correlations.csv"
)

required_figures <- c(
  "TCGA_SKCM_cg16890093_TAP1_correlation_primary_v2.png",
  "TCGA_SKCM_cg07016276_HLA-F_correlation_primary_v2.png",
  "GSE144487_CpG_gene_region_distribution.png",
  "GSE144487_CpG_island_distribution.png"
)

cat("\n--- SOURCE TABLE CHECK ---\n")

for (f in required_tables) {
  cat(
    if (file.exists(file.path(tables_dir, f))) "[FOUND] " else "[MISSING] ",
    f, "\n"
  )
}

cat("\n--- SOURCE FIGURE CHECK ---\n")

for (f in required_figures) {
  cat(
    if (file.exists(file.path(figures_dir, f))) "[FOUND] " else "[MISSING] ",
    f, "\n"
  )
}

# 6. FINAL STATUS -------------------------------------------------------------

cat("\n========== SCRIPT 34 COMPLETED ==========\n")
cat("New module exists:", file.exists(module_file), "\n")
cat("Syntax valid:", syntax_ok, "\n")
cat("Existing app.R modified: NO\n")
cat("Existing analysis results modified: NO\n")

if (!syntax_ok) {
  stop("Module was created but has a syntax error. Do not integrate it yet.")
}

cat("\nSCRIPT 34 PASSED.\n")
