
###############################################################################
###############################################################################
# SCRIPT 16 — INSPECT SHINY UPLOAD SERVER MODULE
###############################################################################

cat("\n========== SCRIPT 16: MODULE SOURCE AUDIT ==========\n")

cat("\n--- SERVER FUNCTION SOURCE ---\n")

cat(
  paste(
    deparse(melanomaUploadServer),
    collapse = "\n"
  )
)

cat("\n\n--- REACTIVE OBJECTS / FUNCTIONS ---\n")

cat(
  "\nresult exists:",
  exists("result"),
  "\n"
)

cat(
  "error_message exists:",
  exists("error_message"),
  "\n"
)

cat("\n--- MODULE ENVIRONMENT ---\n")

print(environment(melanomaUploadServer))

cat("\n--- SHINY VERSION ---\n")

print(packageVersion("shiny"))

cat("\nSCRIPT 16 — COMPLETED\n")

###############################################################################
###############################################################################
# SCRIPT 17 — INSPECT MODULE RETURN VALUE
###############################################################################

cat("\n========== SCRIPT 17: MODULE RETURN VALUE TEST ==========\n")

# 1. Create a temporary valid CSV
test_id <- split$test_sample_ids[1]

test_profile <- as.numeric(
  beta[mlp$selected_cpgs, test_id, drop = TRUE]
)

valid_csv_data <- data.frame(
  CpG = mlp$selected_cpgs,
  Beta = test_profile
)

temp_csv <- tempfile(fileext = ".csv")

write.csv(
  valid_csv_data,
  temp_csv,
  row.names = FALSE
)

diagnostic17 <- new.env(parent = emptyenv())

# 2. Test the module
shiny::testServer(
  melanomaUploadServer,
  args = list(mlp = mlp),
  {
    uploaded_file <- data.frame(
      name = "test_profile.csv",
      size = as.numeric(file.info(temp_csv)$size),
      type = "text/csv",
      datapath = temp_csv,
      stringsAsFactors = FALSE
    )
    
    session$setInputs(profile_file = uploaded_file)
    session$flushReact()
    
    # Inspect the actual return value of the module
    returned_value <- tryCatch(
      session$returned(),
      error = function(e) {
        diagnostic17$return_error <<- conditionMessage(e)
        NULL
      }
    )
    
    diagnostic17$returned <- returned_value
    
    if (!is.null(returned_value)) {
      diagnostic17$returned_names <- names(returned_value)
      
      if ("result" %in% names(returned_value)) {
        diagnostic17$result <- tryCatch(
          shiny::isolate(returned_value$result()),
          error = function(e) {
            diagnostic17$result_error <<- conditionMessage(e)
            NULL
          }
        )
      }
      
      if ("error" %in% names(returned_value)) {
        diagnostic17$error <- tryCatch(
          shiny::isolate(returned_value$error()),
          error = function(e) {
            diagnostic17$error_read_error <<- conditionMessage(e)
            NULL
          }
        )
      }
    }
  }
)

# 3. Display module return information
cat("\n========== RETURN VALUE AUDIT ==========\n")

cat("\nReturned value available:",
    !is.null(diagnostic17$returned), "\n")

if (!is.null(diagnostic17$return_error)) {
  cat("\nReturn access error:\n")
  print(diagnostic17$return_error)
}

if (!is.null(diagnostic17$returned_names)) {
  cat("\nReturned components:\n")
  print(diagnostic17$returned_names)
}

cat("\nPrediction result available:",
    !is.null(diagnostic17$result), "\n")

if (!is.null(diagnostic17$result)) {
  cat("\nPrediction:\n")
  print(diagnostic17$result$prediction)
  
  cat("\nInput summary:\n")
  print(diagnostic17$result$input_summary)
}

cat("\nError value:\n")
print(diagnostic17$error)

if (!is.null(diagnostic17$result_error)) {
  cat("\nResult read error:\n")
  print(diagnostic17$result_error)
}

if (!is.null(diagnostic17$error_read_error)) {
  cat("\nError reactive read error:\n")
  print(diagnostic17$error_read_error)
}

# 4. Independently verify the CSV validator
cat("\n========== INDEPENDENT VALIDATOR CHECK ==========\n")

independent_validation <- validate_methylation_csv(
  temp_csv,
  mlp
)

cat("Validator status:", independent_validation$status, "\n")
cat(
  "Validator predicted class:",
  independent_validation$prediction$predicted_class,
  "\n"
)

# 5. Cleanup
unlink(temp_csv)

cat("\nSCRIPT 17 — COMPLETED\n")

###############################################################################





###############################################################################
# SCRIPT 18 — SHINY REACTIVE ENVIRONMENT AUDIT
###############################################################################

cat("\n========== SCRIPT 18: REACTIVE ENVIRONMENT AUDIT ==========\n")

# 1. Create a fresh valid CSV
test_id <- split$test_sample_ids[1]

valid_csv_data <- data.frame(
  CpG = mlp$selected_cpgs,
  Beta = as.numeric(
    beta[mlp$selected_cpgs, test_id, drop = TRUE]
  )
)

temp_csv <- tempfile(fileext = ".csv")
write.csv(valid_csv_data, temp_csv, row.names = FALSE)

diagnostic18 <- new.env(parent = emptyenv())

# 2. Run module test
shiny::testServer(
  melanomaUploadServer,
  args = list(mlp = mlp),
  {
    cat("\n--- SESSION API ---\n")
    
    cat("session$returned type: ")
    print(typeof(session$returned))
    
    cat("session$returned class: ")
    print(class(session$returned))
    
    cat("\n--- TEST ENVIRONMENT ---\n")
    
    cat("Objects visible in test environment:\n")
    print(ls(all.names = TRUE))
    
    # Find result and error_message in the environment chain
    current_env <- environment()
    found_result <- FALSE
    found_error <- FALSE
    
    while (!identical(current_env, emptyenv())) {
      
      if (
        exists("result", envir = current_env, inherits = FALSE)
      ) {
        found_result <- TRUE
        diagnostic18$result_binding <- get(
          "result",
          envir = current_env,
          inherits = FALSE
        )
        diagnostic18$result_environment <- environmentName(current_env)
      }
      
      if (
        exists("error_message", envir = current_env, inherits = FALSE)
      ) {
        found_error <- TRUE
        diagnostic18$error_binding <- get(
          "error_message",
          envir = current_env,
          inherits = FALSE
        )
        diagnostic18$error_environment <- environmentName(current_env)
      }
      
      current_env <- parent.env(current_env)
    }
    
    cat("\nFound result binding:", found_result, "\n")
    cat("Found error_message binding:", found_error, "\n")
    
    if (found_result) {
      cat("Result binding class:\n")
      print(class(diagnostic18$result_binding))
    }
    
    if (found_error) {
      cat("Error binding class:\n")
      print(class(diagnostic18$error_binding))
    }
    
    # 3. Upload CSV
    uploaded_file <- data.frame(
      name = "test_profile.csv",
      size = as.numeric(file.info(temp_csv)$size),
      type = "text/csv",
      datapath = temp_csv,
      stringsAsFactors = FALSE
    )
    
    session$setInputs(profile_file = uploaded_file)
    session$flushReact()
    
    # 4. Inspect output
    cat("\n--- OUTPUT STATUS ---\n")
    
    diagnostic18$status_output <- tryCatch(
      output$status,
      error = function(e) conditionMessage(e)
    )
    
    print(diagnostic18$status_output)
    
    cat("\n--- OUTPUT PREDICTION ---\n")
    
    diagnostic18$prediction_output <- tryCatch(
      output$prediction,
      error = function(e) conditionMessage(e)
    )
    
    print(diagnostic18$prediction_output)
    
    cat("\n--- OUTPUT SUMMARY ---\n")
    
    diagnostic18$summary_output <- tryCatch(
      output$summary,
      error = function(e) conditionMessage(e)
    )
    
    print(diagnostic18$summary_output)
  }
)

# 5. Display final diagnostics
cat("\n========== SCRIPT 18 RESULTS ==========\n")

cat("\nResult binding found:",
    !is.null(diagnostic18$result_binding), "\n")

cat("Error binding found:",
    !is.null(diagnostic18$error_binding), "\n")

cat("\nStatus output:\n")
print(diagnostic18$status_output)

cat("\nPrediction output:\n")
print(diagnostic18$prediction_output)

cat("\nSummary output:\n")
print(diagnostic18$summary_output)

# 6. Cleanup
unlink(temp_csv)

cat("\nSCRIPT 18 — COMPLETED\n")

###############################################################################


###############################################################################
# SCRIPT 19 — READ SHINY MODULE RETURNED REACTIVES
###############################################################################

cat("\n========== SCRIPT 19: MODULE RETURNED REACTIVES ==========\n")

# 1. Prepare test CSV
test_id <- split$test_sample_ids[1]

valid_csv_data <- data.frame(
  CpG = mlp$selected_cpgs,
  Beta = as.numeric(
    beta[mlp$selected_cpgs, test_id, drop = TRUE]
  )
)

temp_csv <- tempfile(fileext = ".csv")
write.csv(valid_csv_data, temp_csv, row.names = FALSE)

diagnostic19 <- new.env(parent = emptyenv())

# 2. Run test
shiny::testServer(
  melanomaUploadServer,
  args = list(mlp = mlp),
  {
    cat("\n--- MODULE RETURNED OBJECT ---\n")
    
    cat("Returned object class:\n")
    print(class(session$returned))
    
    cat("\nReturned object names:\n")
    print(names(session$returned))
    
    cat("\nResult component class:\n")
    print(class(session$returned$result))
    
    cat("\nError component class:\n")
    print(class(session$returned$error))
    
    # Read initial reactive values
    diagnostic19$initial_result <- session$returned$result()
    diagnostic19$initial_error <- session$returned$error()
    
    # Upload file
    uploaded_file <- data.frame(
      name = "test_profile.csv",
      size = as.numeric(file.info(temp_csv)$size),
      type = "text/csv",
      datapath = temp_csv,
      stringsAsFactors = FALSE
    )
    
    session$setInputs(profile_file = uploaded_file)
    
    # Flush pending reactive work
    session$flushReact()
    
    # Read reactive values after upload
    diagnostic19$final_result <- session$returned$result()
    diagnostic19$final_error <- session$returned$error()
    
    # Inspect status output
    diagnostic19$status <- output$status
  }
)

# 3. Display results
cat("\n========== INITIAL STATE ==========\n")

cat(
  "Initial result is NULL:",
  is.null(diagnostic19$initial_result),
  "\n"
)

cat(
  "Initial error is NULL:",
  is.null(diagnostic19$initial_error),
  "\n"
)

cat("\n========== AFTER UPLOAD ==========\n")

cat(
  "Final result is NULL:",
  is.null(diagnostic19$final_result),
  "\n"
)

cat(
  "Final error is NULL:",
  is.null(diagnostic19$final_error),
  "\n"
)

if (!is.null(diagnostic19$final_error)) {
  cat("\nError message:\n")
  print(diagnostic19$final_error)
}

if (!is.null(diagnostic19$final_result)) {
  cat("\nPrediction:\n")
  print(diagnostic19$final_result$prediction)
  
  cat("\nInput summary:\n")
  print(diagnostic19$final_result$input_summary)
}

cat("\nStatus output:\n")
print(diagnostic19$status)

# 4. Cleanup
unlink(temp_csv)

cat("\nSCRIPT 19 — COMPLETED\n")

###############################################################################
###############################################################################
# SCRIPT 20 — TEST OBSERVER INITIALIZATION
###############################################################################

cat("\n========== SCRIPT 20: OBSERVER INITIALIZATION TEST ==========\n")

# 1. Prepare test CSV
test_id <- split$test_sample_ids[1]

valid_csv_data <- data.frame(
  CpG = mlp$selected_cpgs,
  Beta = as.numeric(
    beta[mlp$selected_cpgs, test_id, drop = TRUE]
  )
)

temp_csv <- tempfile(fileext = ".csv")
write.csv(valid_csv_data, temp_csv, row.names = FALSE)

diagnostic20 <- new.env(parent = emptyenv())

# 2. Run test
shiny::testServer(
  melanomaUploadServer,
  args = list(mlp = mlp),
  {
    cat("\n--- INITIAL REACTIVE FLUSH ---\n")
    
    # Allow observers to initialize before uploading
    session$flushReact()
    
    diagnostic20$initial_result <-
      session$returned$result()
    
    diagnostic20$initial_error <-
      session$returned$error()
    
    cat(
      "Initial result is NULL:",
      is.null(diagnostic20$initial_result),
      "\n"
    )
    
    cat(
      "Initial error is NULL:",
      is.null(diagnostic20$initial_error),
      "\n"
    )
    
    # 3. Upload CSV after initialization
    uploaded_file <- data.frame(
      name = "test_profile.csv",
      size = as.numeric(file.info(temp_csv)$size),
      type = "text/csv",
      datapath = temp_csv,
      stringsAsFactors = FALSE
    )
    
    cat("\n--- UPLOAD EVENT ---\n")
    
    session$setInputs(profile_file = uploaded_file)
    
    session$flushReact()
    
    # 4. Inspect returned reactives
    diagnostic20$final_result <-
      session$returned$result()
    
    diagnostic20$final_error <-
      session$returned$error()
    
    diagnostic20$status <- output$status
    diagnostic20$prediction <- output$prediction
    diagnostic20$summary <- output$summary
  }
)

# 5. Results
cat("\n========== SCRIPT 20 RESULTS ==========\n")

cat(
  "Result returned:",
  !is.null(diagnostic20$final_result),
  "\n"
)

cat(
  "Error returned:",
  !is.null(diagnostic20$final_error),
  "\n"
)

if (!is.null(diagnostic20$final_error)) {
  cat("\nError message:\n")
  print(diagnostic20$final_error)
}

if (!is.null(diagnostic20$final_result)) {
  cat("\nPrediction:\n")
  print(diagnostic20$final_result$prediction)
  
  cat("\nInput summary:\n")
  print(diagnostic20$final_result$input_summary)
}

cat("\nStatus output:\n")
print(diagnostic20$status)

cat("\nPrediction output:\n")
print(diagnostic20$prediction)

cat("\nSummary output:\n")
print(diagnostic20$summary)

# 6. Cleanup
unlink(temp_csv)

cat("\nSCRIPT 20 — COMPLETED\n")

###############################################################################



