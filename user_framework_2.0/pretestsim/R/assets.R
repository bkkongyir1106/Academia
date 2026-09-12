`%||%` <- function(a, b) if (!is.null(a)) a else b

.default_raw_base <- function() {
  getOption(
    "pretestsim.github_raw_base",
    "https://raw.githubusercontent.com/bkkongyir1106/Academia/main/R_pkg/user_framework_2.0"
  )
}

.asset_registry <- function() {
  c(
    framework = "user_framework_rpkg.R",
    ml_framework = "ML_framework_fun.R",
    trained_models = "trained_models_multiple_sample_sizes.RData",
    run_framework = "Run_framework.R",
    run_ml_framework = "run_ml_framework.R"
  )
}

#' List GitHub-backed package assets
#'
#' @return A named character vector of known asset file names.
#' @export
available_assets <- function() {
  .asset_registry()
}

#' Build a raw GitHub URL for a package asset
#'
#' @param asset Asset key from \code{available_assets()} or a file name.
#' @param raw_base Base raw GitHub directory URL.
#' @return A raw GitHub URL.
#' @export
asset_url <- function(asset, raw_base = .default_raw_base()) {
  registry <- .asset_registry()
  file_name <- unname(registry[[asset]] %||% asset)
  if (is.null(file_name) || !nzchar(file_name)) {
    stop("Unknown or empty asset: ", asset, call. = FALSE)
  }
  paste0(sub("/+$", "", raw_base), "/", file_name)
}

.cache_dir <- function() {
  path <- tools::R_user_dir("pretestsim", "cache")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

#' Download and cache a package asset from GitHub
#'
#' @param asset Asset key from \code{available_assets()} or a file name.
#' @param refresh If \code{TRUE}, download again even when cached.
#' @param raw_base Base raw GitHub directory URL.
#' @param cache_dir Cache directory. Defaults to the user's R cache directory.
#' @return Normalized local file path.
#' @export
cache_asset <- function(asset,
                        refresh = FALSE,
                        raw_base = .default_raw_base(),
                        cache_dir = .cache_dir()) {
  registry <- .asset_registry()
  file_name <- unname(registry[[asset]] %||% asset)
  if (is.null(file_name) || !nzchar(file_name)) {
    stop("Unknown or empty asset: ", asset, call. = FALSE)
  }

  destination <- file.path(cache_dir, file_name)
  if (isTRUE(refresh) || !file.exists(destination) || file.info(destination)$size == 0) {
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    url <- asset_url(asset, raw_base = raw_base)
    status <- tryCatch(
      utils::download.file(url, destination, mode = "wb", quiet = TRUE),
      error = function(e) e
    )
    if (inherits(status, "error") || !file.exists(destination) || file.info(destination)$size == 0) {
      if (file.exists(destination)) unlink(destination)
      stop(
        "Could not download asset '", asset, "' from:\n", url, "\n",
        if (inherits(status, "error")) conditionMessage(status) else "",
        call. = FALSE
      )
    }
  }

  normalizePath(destination, winslash = "/", mustWork = TRUE)
}

#' Ensure runtime packages used by the framework are available
#'
#' Missing packages are installed from CRAN when possible. The function returns
#' the packages that are available after the attempt.
#'
#' @return Invisibly, the available package names.
#' @export
ensure_framework_packages <- function() {
  packages <- c(
    "MASS", "nortest", "moments", "tseries", "DescTools", "coin", "Rfit",
    "LaplacesDemon", "VGAM", "evd", "Lmoments", "robustbase", "pracma",
    "ineq", "caret", "glmnet", "randomForest", "gbm", "nnet", "kernlab",
    "e1071", "foreach", "doParallel", "pROC", "pacman"
  )
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    suppressWarnings(try(
      utils::install.packages(missing, repos = "https://cloud.r-project.org", quiet = TRUE),
      silent = TRUE
    ))
  }
  invisible(packages[vapply(packages, requireNamespace, logical(1), quietly = TRUE)])
}
