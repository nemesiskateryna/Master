###############################################################################
# SCRIPT 70 — LAUNCH EXPLAINMYMELANOMA AI
###############################################################################

stopifnot(
  exists("app_object_69"),
  inherits(app_object_69, "shiny.appobj")
)

shiny::runApp(
  app_object_69,
  launch.browser = TRUE
)


