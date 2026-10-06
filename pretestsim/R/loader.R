#' Load the bundled simulation framework
#'
#' No network access or package installation occurs when loading the framework.
#' Optional dependencies are checked by the selected procedures when needed.
#' @param path Optional path to a compatible local framework source file.
#' @return An environment containing the framework functions.
#' @export
#' @examples
#' fw <- load_framework()
#' names(formals(fw$run_simulation))
load_framework <- function(path = NULL) {
  if (is.null(path)) {
    path <- system.file("framework", "framework.R", package = "pretestsim")
  }
  if (!length(path) || !nzchar(path) || !file.exists(path)) {
    stop("Framework source not found: ", path, call. = FALSE)
  }
  env <- new.env(parent = baseenv())
  # Resolve ordinary statistical/graphics calls independently of user objects.
  for (pkg in c("stats", "utils", "graphics", "grDevices")) {
    for (name in getNamespaceExports(pkg)) {
      if (!exists(name, envir = env, inherits = FALSE)) {
        assign(name, getExportedValue(pkg, name), envir = env)
      }
    }
  }
  sys.source(path, envir = env)
  required <- c("run_simulation", "generate_pval", "fisher_combined")
  stopifnot(all(vapply(required, exists, logical(1), envir = env,
                       mode = "function", inherits = FALSE)))
  attr(env, "framework_file") <- normalizePath(path, winslash = "/")
  env
}

`%||%` <- function(a, b) if (is.null(a)) b else a
