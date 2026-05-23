#' Launch the Pretest Simulation Studio dashboard
#'
#' Starts the Shiny app that drives the bundled, unmodified pretest-simulation
#' framework through its six phases (normality AUROC comparison, trade-off /
#' optimal-threshold selection, power vs Type I error ROC-like curves, power and
#' Type I error across sample sizes with AUC tables, and power vs effect size).
#'
#' @param framework_path Optional path to a framework \code{.R} file to use
#'   instead of the bundled copy. Sets the \code{pretestsim.framework_path}
#'   option for the session.
#' @param max_upload_mb Maximum browser upload size in MB (default 2048). Large
#'   trained-model RData is better loaded by path inside the ML panel.
#' @param launch.browser Passed to \code{\link[shiny]{runApp}}.
#' @param ... Further arguments passed to \code{\link[shiny]{runApp}} (e.g.
#'   \code{port}, \code{host}).
#'
#' @return Invisibly \code{NULL}; called for its side effect of running the app.
#' @examples
#' \dontrun{
#'   library(pretestsim)
#'   run_app()                                   # uses the bundled framework
#'   run_app(framework_path = "~/my_framework.R")
#' }
#' @export
run_app <- function(framework_path = NULL,
                    max_upload_mb = 2048,
                    launch.browser = interactive(),
                    ...) {
  for (p in c("shiny", "bslib", "DT"))
    if (!requireNamespace(p, quietly = TRUE))
      stop("Package '", p, "' is required. Install it with install.packages('", p, "').")

  if (!is.null(framework_path)) {
    if (!file.exists(framework_path))
      stop("framework_path does not exist: ", framework_path)
    options(pretestsim.framework_path = normalizePath(framework_path))
  }
  old <- options(shiny.maxRequestSize = max_upload_mb * 1024^2)
  on.exit(options(old), add = TRUE)

  app <- shiny::shinyApp(ui = app_ui(), server = app_server)
  shiny::runApp(app, launch.browser = launch.browser, ...)
  invisible(NULL)
}
