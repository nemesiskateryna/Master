###############################################################################
# SCRIPT 43 — CREATE DYNAMIC EPIGENETIC VISUALIZATION MODULE
###############################################################################

cat("\n========== SCRIPT 43: CREATE DYNAMIC VISUALIZATION MODULE ==========\n")

# 1. CHECK ENVIRONMENT --------------------------------------------------------

if (!exists("app_dir") || !dir.exists(app_dir)) {
  stop("app_dir is missing or does not exist.")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package ggplot2 is required.")
}

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Package shiny is required.")
}

module_file <- file.path(
  app_dir,
  "R",
  "06_dynamic_epigenetic_visualizations.R"
)

if (!dir.exists(dirname(module_file))) {
  stop("Application R directory does not exist.")
}

# Never overwrite an existing module
if (file.exists(module_file)) {
  stop(
    "Module already exists. No file was overwritten:\n",
    module_file
  )
}

# 2. DEFINE MODULE CODE -------------------------------------------------------

module_code <- r"(###############################################################################)"
# MODULE 06 — DYNAMIC EPIGENETIC VISUALIZATIONS
#
# Three visualizations based on the uploaded methylation profile:
# 1. MLP-selected CpG methylation heatmap
# 2. MLP-selected CpG beta-value distribution
# 3. Uploaded sample projected onto the fixed TCGA training PCA reference
#
# Does not retrain the model or modify training data.
###############################################################################

dynamicEpigeneticUI <- function(id) {

  ns <- shiny::NS(id)

  shiny::tagList(

    shiny::h4("MLP-selected CpG methylation heatmap"),

    shiny::p(
      "Methylation beta values for the 5,000 CpG features ",
      "used by the MLP model. CpGs are displayed in model feature order."
    ),

    shiny::plotOutput(
      ns("heatmap"),
      height = "430px"
    ),

    shiny::hr(),

    shiny::h4("Beta-value distribution"),

    shiny::p(
      "Distribution of beta values across the 5,000 CpG features ",
      "used by the MLP model."
    ),

    shiny::plotOutput(
      ns("distribution"),
      height = "350px"
    ),

    shiny::hr(),

    shiny::h4("TCGA-SKCM reference PCA"),

    shiny::p(
      "The uploaded sample is projected into the original PCA space ",
      "fitted on the 368 TCGA-SKCM training samples. ",
      "Training coordinates remain unchanged."
    ),

    shiny::plotOutput(
      ns("pca"),
      height = "480px"
    ),

    shiny::p(
      class = "scientific-disclaimer",
      "PCA position describes similarity in the model's selected ",
      "methylation feature space. It is not a clinical classification ",
      "or an independent prediction."
    )
  )
}


dynamicEpigeneticServer <- function(id, profile, mlp) {

  shiny::moduleServer(id, function(input, output, session) {

    # Validate and extract the 5,000 MLP features -----------------------------

    selected_profile <- shiny::reactive({

      uploaded <- profile()

      shiny::validate(
        shiny::need(
          !is.null(uploaded),
          "Upload a methylation CSV profile to display the visualizations."
        )
      )

      shiny::validate(
        shiny::need(
          is.data.frame(uploaded),
          "Uploaded profile must be a data frame."
        )
      )

      shiny::validate(
        shiny::need(
          all(c("CpG", "Beta") %in% colnames(uploaded)),
          "Uploaded profile must contain CpG and Beta columns."
        )
      )

      selected_cpgs <- mlp$selected_cpgs

      shiny::validate(
        shiny::need(
          length(selected_cpgs) == 5000,
          "The MLP model does not contain the expected 5,000 CpGs."
        )
      )

      positions <- match(selected_cpgs, uploaded$CpG)

      shiny::validate(
        shiny::need(
          !anyNA(positions),
          "Some required MLP CpGs are missing from the uploaded profile."
        )
      )

      beta <- as.numeric(uploaded$Beta[positions])

      shiny::validate(
        shiny::need(
          all(is.finite(beta)) && all(beta >= 0 & beta <= 1),
          "Invalid beta values detected in the selected CpG profile."
        )
      )

      data.frame(
        feature_index = seq_along(selected_cpgs),
        CpG = selected_cpgs,
        Beta = beta,
        stringsAsFactors = FALSE
      )
    })


    # 1. METHYLATION HEATMAP ---------------------------------------------------

    output$heatmap <- shiny::renderPlot({

      dat <- selected_profile()

      dat$grid_x <- ((dat$feature_index - 1) %% 100) + 1
      dat$grid_y <- ((dat$feature_index - 1) %/% 100) + 1

      ggplot2::ggplot(
        dat,
        ggplot2::aes(
          x = grid_x,
          y = grid_y,
          fill = Beta
        )
      ) +

        ggplot2::geom_tile() +

        ggplot2::scale_fill_gradientn(
          colours = c(
            "#191970",
            "#4169E1",
            "#D4AF37"
          ),
          limits = c(0, 1),
          breaks = c(0, 0.25, 0.5, 0.75, 1),
          name = "Beta value"
        ) +

        ggplot2::scale_y_reverse() +

        ggplot2::coord_fixed() +

        ggplot2::labs(
          x = "MLP feature index",
          y = NULL,
          title = "Methylation profile of the uploaded sample",
          subtitle = "5,000 model-selected CpG features"
        ) +

        ggplot2::theme_minimal(base_size = 12) +

        ggplot2::theme(
          plot.title = ggplot2::element_text(
            colour = "#191970",
            face = "bold"
          ),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
          legend.position = "right"
        )

    })


    # 2. BETA-VALUE DISTRIBUTION -----------------------------------------------

    output$distribution <- shiny::renderPlot({

      dat <- selected_profile()

      ggplot2::ggplot(
        dat,
        ggplot2::aes(x = Beta)
      ) +

        ggplot2::geom_histogram(
          bins = 40,
          fill = "#4169E1",
          colour = "white",
          linewidth = 0.25
        ) +

        ggplot2::scale_x_continuous(
          limits = c(0, 1),
          breaks = seq(0, 1, by = 0.1)
        ) +

        ggplot2::labs(
          x = "DNA methylation beta value",
          y = "Number of CpG features",
          title = "Distribution of MLP-selected CpG beta values",
          subtitle = paste(
            "Features:", format(nrow(dat), big.mark = ",")
          )
        ) +

        ggplot2::theme_minimal(base_size = 13) +

        ggplot2::theme(
          plot.title = ggplot2::element_text(
            colour = "#191970",
            face = "bold"
          ),
          panel.grid.minor = ggplot2::element_blank()
        )

    })


    # 3. FIXED TCGA REFERENCE PCA ----------------------------------------------

    output$pca <- shiny::renderPlot({

      dat <- selected_profile()

      pca_ref <- mlp$pca

      shiny::validate(
        shiny::need(
          !is.null(pca_ref$x) &&
            !is.null(pca_ref$rotation) &&
            !is.null(pca_ref$center) &&
            !is.null(pca_ref$scale),
          "Saved PCA reference is incomplete."
        )
      )

      shiny::validate(
        shiny::need(
          nrow(pca_ref$x) == length(mlp$training_labels) &&
            nrow(pca_ref$x) == length(mlp$training_sample_ids),
          "PCA reference and training metadata are inconsistent."
        )
      )

      shiny::validate(
        shiny::need(
          identical(
            rownames(pca_ref$rotation),
            mlp$selected_cpgs
          ),
          "PCA feature order does not match the MLP feature order."
        )
      )

      shiny::validate(
        shiny::need(
          all(c("PC1", "PC2") %in% colnames(pca_ref$x)),
          "The saved PCA reference does not contain PC1 and PC2."
        )
      )

      # Preserve the exact original PCA transformation
      beta <- dat$Beta

      X <- matrix(
        beta,
        nrow = 1,
        dimnames = list(
          "Uploaded sample",
          mlp$selected_cpgs
        )
      )

      X_scaled <- scale(
        X,
        center = pca_ref$center,
        scale = pca_ref$scale
      )

      new_scores <- X_scaled %*% pca_ref$rotation

      shiny::validate(
        shiny::need(
          all(is.finite(new_scores)),
          "The uploaded sample could not be projected into PCA space."
        )
      )

      training_scores <- data.frame(
        PC1 = pca_ref$x[, "PC1"],
        PC2 = pca_ref$x[, "PC2"],
        Immune_phenotype = factor(
          mlp$training_labels,
          levels = c("Immune-low", "Immune-high")
        ),
        Sample = "TCGA training sample"
      )

      new_sample <- data.frame(
        PC1 = as.numeric(new_scores[1, "PC1"]),
        PC2 = as.numeric(new_scores[1, "PC2"])
      )

      ggplot2::ggplot(
        training_scores,
        ggplot2::aes(
          x = PC1,
          y = PC2,
          colour = Immune_phenotype
        )
      ) +

        ggplot2::geom_point(
          alpha = 0.65,
          size = 2.1
        ) +

        ggplot2::geom_point(
          data = new_sample,
          ggplot2::aes(x = PC1, y = PC2),
          inherit.aes = FALSE,
          shape = 23,
          size = 5,
          stroke = 1.2,
          fill = "#D4AF37",
          colour = "#191970"
        ) +

        ggplot2::scale_colour_manual(
          values = c(
            "Immune-low" = "#6A0DAD",
            "Immune-high" = "#4169E1"
          )
        ) +

        ggplot2::labs(
          x = "PC1",
          y = "PC2",
          colour = "Training phenotype",
          title = "Uploaded sample in TCGA-SKCM PCA space",
          subtitle = "Reference: 368 original TCGA-SKCM training samples"
        ) +

        ggplot2::theme_minimal(base_size = 13) +

        ggplot2::theme(
          plot.title = ggplot2::element_text(
            colour = "#191970",
            face = "bold"
          ),
          legend.position = "bottom",
          panel.grid.minor = ggplot2::element_blank()
        ) +

        ggplot2::annotate(
          "label",
          x = new_sample$PC1,
          y = new_sample$PC2,
          label = "Uploaded sample",
          nudge_y = 0.06 * diff(range(training_scores$PC2)),
          size = 3,
          colour = "#191970",
          fill = "#FFF9E8"
        )

    })

    invisible(NULL)
  })
}

###############################################################################
# END MODULE 06
###############################################################################"

# 3. VALIDATE CODE BEFORE WRITING ---------------------------------------------

cat("\n--- VALIDATING MODULE SYNTAX BEFORE FILE CREATION ---\n")

parsed_candidate <- tryCatch(
    parse(text = module_code),
    error = function(e) {
        stop(
            "Module syntax validation failed. No file created.\n",
            conditionMessage(e)
        )
    }
)

cat("Candidate syntax: VALID\n")

# 4. WRITE NEW FILE -----------------------------------------------------------

writeLines(
    module_code,
    con = module_file,
    useBytes = TRUE
)

cat("\nNew module created:\n")
cat(module_file, "\n")

# 5. VALIDATE SAVED FILE -------------------------------------------------------

parse(
    file = module_file,
    keep.source = FALSE
)

cat("\nSaved module syntax: VALID\n")

# 6. VERIFY REQUIRED FUNCTIONS ------------------------------------------------

module_env <- new.env(parent = globalenv())

sys.source(
    module_file,
    envir = module_env
)

required_functions <- c(
    "dynamicEpigeneticUI",
    "dynamicEpigeneticServer"
)

function_check <- vapply(
    required_functions,
    exists,
    logical(1),
    envir = module_env,
    inherits = FALSE
)

print(function_check)

if (!all(function_check)) {
    stop("One or more required module functions are missing.")
}

cat("\n========== SCRIPT 43 COMPLETED ==========\n")
cat("New module created and validated.\n")
cat("No existing files modified.\n")



###############################################################################
# SCRIPT 44 — RESTORE DYNAMIC MODULE
###############################################################################

cat("\n========== SCRIPT 44: RESTORE MODULE ==========\n")

# 1. PATHS --------------------------------------------------------------------

script_file <- path.expand(
  "~/master/Master/ExplainMyMelanoma/DYNAMIC_module.R"
)

module_file <- file.path(
  app_dir,
  "R",
  "06_dynamic_epigenetic_visualizations.R"
)

stopifnot(file.exists(script_file))
stopifnot(file.exists(module_file))

# 2. READ ORIGINAL SCRIPT -----------------------------------------------------

script_lines <- readLines(
  script_file,
  warn = FALSE
)

# Locate the incorrectly opened raw string
start_idx <- grep(
  'module_code <- r"\\(',
  script_lines
)

if (length(start_idx) != 1) {
  stop("Could not uniquely identify module_code declaration.")
}

# Locate the beginning of validation section
end_idx <- grep(
  "^# 3\\. VALIDATE CODE BEFORE WRITING",
  script_lines
)

if (length(end_idx) != 1 || end_idx <= start_idx) {
  stop("Could not identify the end of the module source block.")
}

# Extract module source, excluding the broken declaration
module_lines <- script_lines[
  (start_idx + 1):(end_idx - 1)
]

# Remove trailing blank lines
while (
  length(module_lines) > 0 &&
  !nzchar(trimws(tail(module_lines, 1)))
) {
  module_lines <- head(
    module_lines,
    -1
  )
}

cat("Extracted module lines:", length(module_lines), "\n")

# 3. VALIDATE EXTRACTED SOURCE -----------------------------------------------

module_code_recovered <- paste(
  module_lines,
  collapse = "\n"
)

parsed_module <- tryCatch(
  parse(text = module_code_recovered),
  error = function(e) {
    stop(
      "Recovered module has a syntax error:\n",
      conditionMessage(e)
    )
  }
)

cat("Recovered module syntax: VALID\n")

# 4. VERIFY FUNCTION DEFINITIONS ---------------------------------------------

required_patterns <- c(
  "dynamicEpigeneticUI <- function",
  "dynamicEpigeneticServer <- function"
)

for (pattern in required_patterns) {
  if (!any(grepl(pattern, module_lines, fixed = TRUE))) {
    stop(
      "Missing function definition: ",
      pattern
    )
  }
}

cat("Both module functions found.\n")

# 5. PRESERVE THE EXISTING SHORT FILE ----------------------------------------

backup_file <- paste0(
  module_file,
  ".incomplete_backup"
)

if (file.exists(backup_file)) {
  stop(
    "Recovery backup already exists. No files changed:\n",
    backup_file
  )
}

if (!file.copy(
  module_file,
  backup_file,
  overwrite = FALSE
)) {
  stop("Could not preserve the incomplete module file.")
}

cat("Incomplete file preserved:\n", backup_file, "\n")

# 6. WRITE RECOVERED MODULE ---------------------------------------------------

writeLines(
  module_code_recovered,
  con = module_file,
  useBytes = TRUE
)

# 7. FINAL VALIDATION ---------------------------------------------------------

parse(
  file = module_file,
  keep.source = FALSE
)

module_env <- new.env(parent = globalenv())

sys.source(
  module_file,
  envir = module_env
)

function_check <- vapply(
  c(
    "dynamicEpigeneticUI",
    "dynamicEpigeneticServer"
  ),
  exists,
  logical(1),
  envir = module_env,
  inherits = FALSE
)

cat("\n--- FUNCTION CHECK ---\n")
print(function_check)

cat("\nFILE SIZE:", file.info(module_file)$size, "bytes\n")

if (!all(function_check)) {
  stop("Module recovery verification failed.")
}

cat("\n========== SCRIPT 44 COMPLETED ==========\n")
cat("Dynamic visualization module successfully recovered.\n")
cat("Original script preserved.\n")

###############################################################################
# SCRIPT 47 — TEST PROFILE AND MODULE UI AUDIT
###############################################################################

cat("\n========== SCRIPT 47: MODULE INPUT AUDIT ==========\n")

# 1. DEFINE PATHS --------------------------------------------------------------

module_file <- file.path(
  app_dir,
  "R",
  "06_dynamic_epigenetic_visualizations.R"
)

project_root <- normalizePath(
  file.path(app_dir, "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

model_file <- file.path(
  project_root,
  "results",
  "models",
  "TCGA_SKCM_MLP_final_model.rds"
)

test_file <- file.path(
  app_dir,
  "data",
  "reference",
  "ExplainMyMelanoma_test_profile.csv"
)

cat("Module exists:", file.exists(module_file), "\n")
cat("Model exists:", file.exists(model_file), "\n")
cat("Test profile exists:", file.exists(test_file), "\n")

stopifnot(
  file.exists(module_file),
  file.exists(model_file),
  file.exists(test_file)
)

# 2. LOAD MODEL AND TEST PROFILE ----------------------------------------------

mlp_test <- readRDS(model_file)

profile_test <- read.csv(
  test_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n--- MODEL ---\n")
cat("Selected CpGs:", length(mlp_test$selected_cpgs), "\n")
cat("PCA dimensions:", paste(dim(mlp_test$pca$x), collapse = " x "), "\n")

cat("\n--- TEST PROFILE ---\n")
cat("Dimensions:", paste(dim(profile_test), collapse = " x "), "\n")
cat("Columns:", paste(colnames(profile_test), collapse = ", "), "\n")

# 3. CHECK INPUT STRUCTURE -----------------------------------------------------

required_columns <- c("CpG", "Beta")

if (!all(required_columns %in% colnames(profile_test))) {
  stop("Test profile does not contain CpG and Beta columns.")
}

if (anyDuplicated(profile_test$CpG)) {
  stop("Duplicate CpG identifiers detected in test profile.")
}

positions <- match(
  mlp_test$selected_cpgs,
  profile_test$CpG
)

cat("\n--- FEATURE MATCHING ---\n")
cat("Required CpGs:", length(mlp_test$selected_cpgs), "\n")
cat("Matched CpGs:", sum(!is.na(positions)), "\n")
cat("Missing CpGs:", sum(is.na(positions)), "\n")

if (anyNA(positions)) {
  stop("Some MLP-selected CpGs are missing from the test profile.")
}

beta_test <- as.numeric(profile_test$Beta[positions])

cat("All beta values finite:", all(is.finite(beta_test)), "\n")
cat(
  "All beta values in [0,1]:",
  all(beta_test >= 0 & beta_test <= 1),
  "\n"
)

if (!all(is.finite(beta_test)) ||
    !all(beta_test >= 0 & beta_test <= 1)) {
  stop("Invalid beta values detected.")
}

# 4. CHECK FEATURE ORDER -------------------------------------------------------

cat("\n--- FEATURE ORDER ---\n")

cat(
  "Profile aligned to model order:",
  identical(
    profile_test$CpG[positions],
    mlp_test$selected_cpgs
  ),
  "\n"
)

cat("First 5 matched CpGs:\n")
print(head(profile_test$CpG[positions], 5))

# 5. CHECK MODULE UI -----------------------------------------------------------

module_env <- new.env(parent = globalenv())

sys.source(
  module_file,
  envir = module_env
)

ui_test <- module_env$dynamicEpigeneticUI("dynamic_test")

cat("\n--- MODULE UI ---\n")
cat("UI object class:", paste(class(ui_test), collapse = ", "), "\n")
cat("UI object created:", !is.null(ui_test), "\n")

# 6. FINAL SUMMARY -------------------------------------------------------------

cat("\n========== SCRIPT 47 SUMMARY ==========\n")
cat("Selected CpGs matched:", sum(!is.na(positions)), "/ 5000\n")
cat(
  "Beta values valid:",
  all(is.finite(beta_test) & beta_test >= 0 & beta_test <= 1),
  "\n"
)
cat("UI created:", !is.null(ui_test), "\n")

cat("\n========== SCRIPT 47 COMPLETED ==========\n")
