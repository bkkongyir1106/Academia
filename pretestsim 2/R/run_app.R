#' Launch the normality pretesting dashboard
#'
#' The bundled framework loads on startup. Simulations and their progress
#' indicators start only after the Run simulation button is clicked.
#' @param launch.browser Whether to open the dashboard in a browser.
#' @param ... Additional arguments passed to \code{shiny::runApp()}.
#' @return Invisibly \code{NULL} when the app stops.
#' @export
#' @examples
#' if (interactive()) run_app()
run_app <- function(launch.browser = interactive(), ...) {
  app <- shiny::shinyApp(ui = app_ui(), server = app_server)
  shiny::runApp(app, launch.browser = launch.browser, ...)
  invisible(NULL)
}
