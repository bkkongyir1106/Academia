#' pretestsim: Interactive Studio for the Pretest-Simulation Framework
#'
#' A Shiny dashboard that drives the original, unmodified pretest-simulation
#' framework (bundled in \code{inst/framework}) through the six phases of its
#' \code{run_simulation()} engine. Start it with \code{\link{run_app}}; load the
#' framework programmatically with \code{\link{load_framework}}.
#'
#' @section Main entry points:
#' \itemize{
#'   \item \code{\link{run_app}} - launch the dashboard.
#'   \item \code{\link{load_framework}} - source the framework into an environment.
#'   \item \code{\link{ensure_framework_packages}} - install the framework's deps.
#' }
#'
#' @keywords internal
"_PACKAGE"

# Quiet R CMD check notes for symbols used only inside the Shiny reactive
# context or pulled from the framework environment / uploaded RData.
utils::globalVariables(c("trained_models"))
