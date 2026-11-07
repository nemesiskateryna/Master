###############################################################################
# SCRIPT 39 — LAUNCH EXPLAINMYMELANOMA AI
###############################################################################

cat("\n========== SCRIPT 39: LAUNCH APP ==========\n")

if (!exists("app_env") || !exists("app_object")) {
  stop("Run SCRIPT 38 first.")
}

if (!inherits(app_object, "shiny.appobj")) {
  stop("Invalid Shiny application object.")
}

cat("Launching ExplainMyMelanoma AI...\n")
cat("The browser should open automatically.\n")
cat("Stop the app with the red Stop button in RStudio.\n\n")

shiny::runApp(
  app_object,
  launch.browser = TRUE
)
