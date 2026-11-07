
###############################################################################
# SCRIPT 23 — LAUNCH EXPLAINMYMELANOMA AI
###############################################################################

app_dir <- file.path(
  "/Users/katie/master/Master",
  "app/ExplainMyMelanoma"
)

cat("\n========== LAUNCHING APPLICATION ==========\n")
cat("Application directory:", app_dir, "\n")

shiny::runApp(
  appDir = app_dir,
  launch.browser = TRUE
)

###############################################################################

