#' pretestsim: GitHub-backed adaptive normality-pretesting studio
#'
#' The package launches a Shiny dashboard and programmatic loaders for the
#' adaptive normality-pretesting framework. Required framework and ML assets are
#' downloaded from the package GitHub repository and cached locally.
#'
#' @section Main entry points:
#' \itemize{
#'   \item \code{\link{run_app}} launches the dashboard.
#'   \item \code{\link{load_framework}} downloads and sources the framework.
#'   \item \code{\link{load_ml_resources_from_github}} downloads ML assets.
#' }
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c("trained_models"))
