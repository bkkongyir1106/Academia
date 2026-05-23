#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Ensure the framework's runtime packages are available
#'
#' Attempts to install (quietly) the packages the bundled framework reaches for
#' via \code{::} in its runnable paths. Anything missing simply disables the
#' feature that needs it (for example the Pareto distribution needs
#' \pkg{VGAM}). Never errors.
#'
#' @return (invisibly) the character vector of packages that are available.
#' @export
ensure_framework_packages <- function() {
  pkgs <- c("MASS", "nortest", "moments", "tseries", "DescTools",
            "coin", "Rfit", "LaplacesDemon", "VGAM", "evd")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    suppressWarnings(try(
      utils::install.packages(miss, repos = "https://cloud.r-project.org", quiet = TRUE),
      silent = TRUE))
  invisible(pkgs[vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)])
}

# Decide whether a top-level expression should be skipped during loading.
.fw_should_skip <- function(expr) {
  txt <- paste(deparse(expr), collapse = " ")

  # the run_simulation DEFINITION must be kept; its example CALLS must be skipped
  is_def <- grepl("^\\s*run_simulation\\s*<-\\s*function", txt)
  if (grepl("run_simulation", txt) && !is_def) return(TRUE)

  if (grepl("^\\s*setwd\\(", txt))                      return(TRUE)
  if (grepl("^\\s*source\\(", txt))                     return(TRUE)
  if (grepl("p_load\\(|install\\.packages\\(", txt))    return(TRUE)
  if (grepl("requireNamespace\\(\\s*[\"']pacman", txt)) return(TRUE)
  if (grepl("^\\s*load\\(", txt))                       return(TRUE)
  FALSE
}

#' Load the original framework into an environment, unmodified
#'
#' Parses the framework source file and evaluates every function definition,
#' every \code{norm_config_*} preset and every \code{*_fns} / \code{*_params}
#' bundle exactly as written. Only the machine-specific side-effects
#' (\code{setwd}, the external ML \code{source()}, the \code{pacman} install
#' block, and the \code{load()} of trained models) and the example
#' \code{run_simulation()} calls at the end of the file are skipped. The
#' \code{save.image()} call inside Phase 6 is redirected to a no-op so a long
#' session is never dumped; \code{save()} of the compact results object still
#' runs.
#'
#' @param path Path to the framework \code{.R} file. Defaults to the copy
#'   bundled in the package (\code{inst/framework/user_framework_rpk_v3_0.R}).
#' @return An environment containing the framework's objects (notably
#'   \code{run_simulation}), with a \code{"load_report"} attribute.
#' @examples
#' \dontrun{
#' fw <- load_framework()
#' args <- c(fw$onesample_fns, list(
#'   Nsim = 100, N_tradeoff = 200, test_type = "demo",
#'   distributions = c("exponential", "normal"),
#'   norm_config = list(method = "classical", config = list(norm_test = "SW")),
#'   threshold_grid = seq(0, 1, 0.05), tol_pos = 0.01, loss_tol = 0.01,
#'   test_alpha = 0.05, center_by = "median", effect_size = 0.5,
#'   sample_sizes = c(10, 20, 30), single_n = 20,
#'   effect_sizes_plot = c(0, 0.3, 0.6), sig_levels = seq(0, 1, 0.05),
#'   norm_test = "SW", selected_tests = "SW", phases = c(1, 2, 4)))
#' res <- do.call(fw$run_simulation, args)
#' }
#' @export
load_framework <- function(path = default_framework_path()) {
  if (is.null(path) || is.na(path) || !file.exists(path))
    stop("Framework file not found at: ", path %||% "(NULL)")

  ensure_framework_packages()

  env <- new.env(parent = globalenv())
  assign("save.image", function(...) invisible(NULL), envir = env)

  exprs   <- parse(path)
  kept    <- skipped <- 0L
  errors  <- character(0)

  for (e in exprs) {
    if (.fw_should_skip(e)) { skipped <- skipped + 1L; next }
    tryCatch({ eval(e, envir = env); kept <- kept + 1L },
             error = function(err)
               errors[[length(errors) + 1L]] <<-
                 paste0(substr(paste(deparse(e), collapse = " "), 1, 60),
                        " ... : ", conditionMessage(err)))
  }

  attr(env, "load_report") <- list(kept = kept, skipped = skipped, errors = errors)
  if (!is.function(env$run_simulation))
    stop("Framework loaded but run_simulation() was not found.")
  env
}

#' Path to the framework file shipped with the package
#'
#' Resolution order: the \code{pretestsim.framework_path} option, then the
#' bundled \code{inst/framework} copy.
#' @return A file path, or \code{NA_character_} if none is found.
#' @keywords internal
default_framework_path <- function() {
  opt <- getOption("pretestsim.framework_path", NULL)
  if (!is.null(opt) && nzchar(opt) && file.exists(opt)) return(normalizePath(opt))
  bundled <- system.file("framework", "user_framework_rpk_v3_0.R",
                         package = "pretestsim")
  if (nzchar(bundled) && file.exists(bundled)) return(normalizePath(bundled))
  NA_character_
}
