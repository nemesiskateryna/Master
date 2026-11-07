# Title: 02_GEO_cohort_audit.R
# Project: Master

# GSE144487: DNA methylation
# GSE65904: Gene expression

#1. Packages
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop("Install GEOquery using BiocManager::install('GEOquery')")
}
library(GEOquery)
library(Biobase)
library(dplyr)
library(tibble)
library(stringr)

# GEO accessions
geo_accessions <- c(
  "GSE144487",
  "GSE65904"
)
geo_raw_dir <- "data/GEO/raw"
geo_processed_dir <- "data/GEO/processed"
dir.create(geo_raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(geo_processed_dir, recursive = TRUE, showWarnings = FALSE)

#3. Function to download GEO metadata
if (inherits(geo_object, "ExpressionSet")) {
  geo_list <- list(geo_object)
} else {
  geo_list <- geo_object
}
if (length(geo_list) == 0) {
  stop("No ExpressionSet returned for ", accession)
}
names(geo_list) <- paste0(
  accession,
  "platform",
  seq_along(geo_list)
)
all_metadata <- list()
for (i in seq_along(geo_list)) {
  
  eset <- geo_list[[i]]
  
  pheno <- Biobase::pData(eset)
  
  if (nrow(pheno) == 0) {
    warning("No sample metadata for ", names(geo_list)[i])
    next
  }
  
  pheno <- as.data.frame(pheno)
  
  pheno$geo_accession <- rownames(pheno)
  pheno$GSE <- accession
  pheno$platform_id <- Biobase::annotation(eset)
  
  # Save metadata
  output_name <- paste0(
    accession,
    "_metadata_platform_",
    i
  )
  
  saveRDS(
    pheno,
    file.path(
      geo_processed_dir,
      paste0(output_name, ".rds")
    )
  )
  
  write.csv(
    pheno,
    file.path(
      geo_processed_dir,
      paste0(output_name, ".csv")
    ),
    row.names = FALSE
  )
  
  all_metadata[[i]] <- pheno
  
  message(
    "Platform: ", Biobase::annotation(eset),
    " | Samples: ", nrow(pheno)
  )
}
combined_metadata <- bind_rows(all_metadata)
saveRDS(
  combined_metadata,
  file.path(
    geo_processed_dir,
    paste0(accession, "_all_metadata.rds")
  )
)
write.csv(
  combined_metadata,
  file.path(
    geo_processed_dir,
    paste0(accession, "_all_metadata.csv")
  ),
  row.names = FALSE
)
return(
  list(
    accession = accession,
    geo_object = geo_list,
    metadata = combined_metadata
  )
)
}

#4. Download both datasets
------------------------------------------------------------
  geo_results <- list()
  for (accession in geo_accessions) {
    geo_results[[accession]] <- download_geo_metadata(
      accession
    )
  }
 # 5. Переглядаємо датасети
  dataset_summary <- bind_rows(
    lapply(names(geo_results), function(accession) {
      
      metadata <- geo_results[[accession]]$metadata
      
      data.frame(
        GSE = accession,
        n_samples = nrow(metadata),
        n_platforms = length(
          geo_results[[accession]]$geo_object
        ),
        platform_ids = paste(
          unique(metadata$platform_id),
          collapse = "; "
        )
      )
    })
  )
  print(dataset_summary)
  write.csv(
    dataset_summary,
    "results/tables/GEO_dataset_summary.csv",
    row.names = FALSE
  )
  
  # 6. Перевірити назви колонок
  metadata_columns <- bind_rows(
    lapply(names(geo_results), function(accession) {
      
      data.frame(
        GSE = accession,
        column_name = names(
          geo_results[[accession]]$metadata
        )
      )
      
    })
  )
  
  write.csv(
    metadata_columns,
    "results/tables/GEO_metadata_columns.csv",
    row.names = FALSE
  )
  
 # 7. Search for survival-related metadata fields
  survival_pattern <- paste(
    c(
      "survival",
      "overall",
      "os.days",
      "os.events",
      "os_days",
      "os_event",
      "follow.up",
      "follow_up",
      "vital",
      "death",
      "event",
      "status",
      "days",
      "months",
      "time"
    ),
    collapse = "|"
  )
  
  survival_fields <- metadata_columns %>%
    filter(
      str_detect(
        tolower(column_name),
        survival_pattern
      )
    )
  
  print(survival_fields)
  
  write.csv(
    survival_fields,
    "results/tables/GEO_survival_fields.csv",
    row.names = FALSE
  )
  
  # 8. Search for clinical and sample identifiers
  
  clinical_pattern <- paste(
    c(
      "characteristics",
      "clinical",
      "patient",
      "sample",
      "title",
      "source_name",
      "tissue",
      "tumor",
      "stage",
      "treatment",
      "response"
    ),
    collapse = "|"
  )
  
  clinical_fields <- metadata_columns %>%
    filter(
      str_detect(
        tolower(column_name),
        clinical_pattern
      )
    )
  
  write.csv(
    clinical_fields,
    "results/tables/GEO_clinical_fields.csv",
    row.names = FALSE
  )
  
  # 9. Extract sample identifiers and titles
  # ------------------------------------------------------------
  
  sample_overview <- bind_rows(
    lapply(names(geo_results), function(accession) {
      
      metadata <- geo_results[[accession]]$metadata
      
      selected_columns <- intersect(
        c(
          "geo_accession",
          "title",
          "source_name_ch1",
          "platform_id",
          "GSE"
        ),
        names(metadata)
      )
      
      metadata[, selected_columns, drop = FALSE]
      
    })
  )
  
  write.csv(
    sample_overview,
    "results/tables/GEO_sample_overview.csv",
    row.names = FALSE
  )
  
  # ------------------------------------------------------------
  # 10. Inspect survival-related values
  # ------------------------------------------------------------
  
  for (accession in names(geo_results)) {
    
    metadata <- geo_results[[accession]]$metadata
    
    relevant_columns <- names(metadata)[
      str_detect(
        tolower(names(metadata)),
        survival_pattern
      )
    ]
    
    cat("\n====================================\n")
    cat("DATASET:", accession, "\n")
    cat("====================================\n")
    
    if (length(relevant_columns) == 0) {
      
      cat("No obvious survival-related fields found.\n")
      
    } else {
      
      for (column in relevant_columns) {
        
        cat("\nColumn:", column, "\n")
        
        print(
          head(
            unique(metadata[[column]]),
            15
          )
        )
        
      }
      
    }
    
  }
  
  # ------------------------------------------------------------
  # 11. Save audit objects
  # ------------------------------------------------------------
  
  saveRDS(
    geo_results,
    file.path(
      geo_processed_dir,
      "GEO_metadata_audit_objects.rds"
    )
  )
  
  writeLines(
    capture.output(sessionInfo()),
    "results/logs/GEO_audit_sessionInfo.txt"
  )
  
  message("\nGEO metadata audit completed.")
  
  
  
  
  
  
  
  
  
  
  