# =============================================================================
# framework_loader.R
# -----------------------------------------------------------------------------
# Loads the user's ORIGINAL framework file (user_framework_rpk_v3_0.R) WITHOUT
# editing it. It parses the file into top-level expressions and evaluates them
# into a dedicated environment, skipping only:
#   * setwd()                       - machine-specific working directory
#   * source(...)                   - the external ML framework (path-specific)
#   * pacman::p_load / install.*    - package bootstrapping (handled separately)
#   * load("...trained_models...")  - ML model file (user supplies via the app)
#   * the example `run_simulation()` invocations at the bottom of the file
#
# Everything else - every function definition, every norm_config_* preset,
# every *_fns / *_params bundle - is evaluated exactly as written.
#
# save.image() (called inside run_simulation Phase 6) is redirected to a no-op
# in the load environment so a long-running app session is never dumped to disk;
# save() of the compact `results` object is left intact.
# =============================================================================

ensure_framework_packages <- function() {
  # Packages the framework reaches for via :: in its runnable (non-ML-training)
  # paths. Attempted quietly; anything missing simply disables the feature that
  # needs it (e.g. the Pareto distribution needs VGAM).
  pkgs <- c("MASS", "nortest", "moments", "tseries", "DescTools",
            "coin", "Rfit", "LaplacesDemon", "VGAM", "evd")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    suppressWarnings(try(
      install.packages(miss, repos = "https://cloud.r-project.org", quiet = TRUE),
      silent = TRUE))
  invisible(pkgs[vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)])
}

# Decide whether a top-level expression should be skipped during loading.
.fw_should_skip <- function(expr) {
  txt <- paste(deparse(expr), collapse = " ")

  # the run_simulation DEFINITION must be kept; its example CALLS must be skipped
  is_def <- grepl("^\\s*run_simulation\\s*<-\\s*function", txt)
  if (grepl("run_simulation", txt) && !is_def) return(TRUE)

  # machine-specific side effects / package bootstrap / model load
  if (grepl("^\\s*setwd\\(", txt))                      return(TRUE)
  if (grepl("^\\s*source\\(", txt))                     return(TRUE)
  if (grepl("p_load\\(|install\\.packages\\(", txt))    return(TRUE)
  if (grepl("requireNamespace\\(\\s*[\"']pacman", txt)) return(TRUE)
  if (grepl("^\\s*load\\(", txt))                       return(TRUE)
  FALSE
}

load_framework <- function(path) {
  if (!file.exists(path))
    stop("Framework file not found at: ", path)

  ensure_framework_packages()

  env <- new.env(parent = globalenv())
  # Redirect Phase-6 whole-session dump to a no-op inside the framework's scope.
  assign("save.image", function(...) invisible(NULL), envir = env)

  exprs <- parse(path)
  skipped <- kept <- 0L
  errors  <- character(0)

  for (e in exprs) {
    if (.fw_should_skip(e)) { skipped <- skipped + 1L; next }
    res <- tryCatch({ eval(e, envir = env); kept <- kept + 1L; TRUE },
                    error = function(err) {
                      errors[[length(errors) + 1L]] <<-
                        paste0(substr(paste(deparse(e), collapse = " "), 1, 60),
                               " ... : ", conditionMessage(err))
                      FALSE
                    })
  }

  attr(env, "load_report") <- list(kept = kept, skipped = skipped, errors = errors)
  if (!is.function(env$run_simulation))
    stop("Framework loaded but run_simulation() was not found.")
  env
}
