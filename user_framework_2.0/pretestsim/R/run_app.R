#' Launch the Pretest Simulation Studio dashboard
#'
#' @param raw_base Base raw GitHub directory URL containing the framework assets.
#' @param refresh_assets If \code{TRUE}, refresh cached GitHub files on launch.
#' @param max_upload_mb Maximum Shiny upload size in MB.
#' @param launch.browser Passed to \code{\link[shiny]{runApp}}.
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#' @return Invisibly \code{NULL}.
#' @export
run_app <- function(raw_base = .default_raw_base(),
                    refresh_assets = FALSE,
                    max_upload_mb = 2048,
                    launch.browser = interactive(),
                    ...) {
  for (pkg in c("shiny", "bslib", "DT")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required. Install it with install.packages('", pkg, "').")
    }
  }

  options(pretestsim.github_raw_base = raw_base)
  old <- options(shiny.maxRequestSize = max_upload_mb * 1024^2)
  on.exit(options(old), add = TRUE)

  if (isTRUE(refresh_assets)) {
    cache_asset("framework", refresh = TRUE, raw_base = raw_base)
  }

  app <- shiny::shinyApp(
    ui = app_ui(raw_base = raw_base),
    server = function(input, output, session) app_server(input, output, session, raw_base = raw_base)
  )
  shiny::runApp(app, launch.browser = launch.browser, ...)
  invisible(NULL)
}
