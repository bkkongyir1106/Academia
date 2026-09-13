#' Load the framework from GitHub or a local path
#'
#' @param path Optional local framework R file. When \code{NULL}, the file is
#'   downloaded from GitHub and cached.
#' @param refresh Download a fresh GitHub copy before loading.
#' @param raw_base Base raw GitHub directory URL.
#' @return An environment containing the framework functions.
#' @export
load_framework <- function(path = NULL,
                           refresh = FALSE,
                           raw_base = .default_raw_base()) {
  if (is.null(path)) {
    path <- cache_asset("framework", refresh = refresh, raw_base = raw_base)
  }
  if (!file.exists(path)) {
    stop("Framework file not found: ", path, call. = FALSE)
  }

  ensure_framework_packages()

  env <- new.env(parent = globalenv())
  sys.source(path, envir = env, chdir = TRUE)

  required <- c("run_simulation", "load_ml_resources", "make_ml_normality_configs")
  missing <- required[!vapply(required, exists, logical(1), envir = env, mode = "function", inherits = FALSE)]
  if (length(missing) > 0L) {
    stop(
      "The framework loaded, but required function(s) are missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  attr(env, "framework_file") <- normalizePath(path, winslash = "/", mustWork = TRUE)
  env
}

#' Load ML support functions and trained models from GitHub
#'
#' @param framework_env Environment returned by \code{load_framework()}.
#' @param refresh Download fresh GitHub copies before loading.
#' @param raw_base Base raw GitHub directory URL.
#' @param required_sample_sizes Optional sample sizes to map to available model
#'   bundles.
#' @return The ML resources object produced by the user framework.
#' @export
load_ml_resources_from_github <- function(framework_env,
                                          refresh = FALSE,
                                          raw_base = .default_raw_base(),
                                          required_sample_sizes = NULL) {
  if (!is.environment(framework_env)) {
    stop("framework_env must be an environment returned by load_framework().", call. = FALSE)
  }
  ml_file <- cache_asset("ml_framework", refresh = refresh, raw_base = raw_base)
  model_file <- cache_asset("trained_models", refresh = refresh, raw_base = raw_base)

  framework_env$load_ml_resources(
    framework_file = ml_file,
    model_file = model_file,
    project_dir = dirname(ml_file),
    required_sample_sizes = required_sample_sizes
  )
}
