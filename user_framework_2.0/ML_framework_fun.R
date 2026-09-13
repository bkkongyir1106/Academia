#' ============================================================
#' ML-BASED NORMALITY TESTING FRAMEWORK
#' Part 1: Setup, Helpers, Distribution Matching, Data Generation
#' ============================================================
#'
#' This script defines the complete machine-learning normality-testing framework.
#' It installs and loads required packages, generates moment-matched samples,
#' calculates normality features, trains and evaluates classifiers, produces
#' ROC and variable-importance results, and saves the requested outputs.
#'
#' Main steps:
#'   1. configure packages, random-number generation, and parallel processing;
#'   2. define training, validation, and evaluation distributions;
#'   3. generate samples and calculate selected normality features;
#'   4. train the requested ML models with cross-validation;
#'   5. evaluate classification performance and ROC behavior;
#'   6. save plots, tables, trained models, and reusable results.
#'
#' Performance revisions in this version:
#'   1. calculate only the requested feature subset;
#'   2. use one configurable cross-validation run by default;
#'   3. parallelize caret tuning without nested sample-size workers;
#'   4. reuse evaluation predictions in pair-specific ML ROC curves;
#'   5. use common Monte Carlo samples across classical ROC tests;
#'   6. report progress and elapsed time during long stages.
#'
#' ============================================================
#' SECTION 0: SETUP
#' ============================================================

RNGversion("4.2.0")

#' Continue the calculation.
set.seed(
  12345,
  kind        = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(mc_doScale_quiet = TRUE)

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pkgs <- c(
  #' Distributions
  "LaplacesDemon", "VGAM", "evd",

  #' Normality tests
  "nortest", "DescTools", "moments", "tseries",

  #' Feature engineering
  "Lmoments", "robustbase", "pracma", "ineq",

  #' Machine learning
  "caret", "glmnet", "randomForest", "gbm",
  "nnet", "kernlab", "e1071", "foreach", "doParallel",

  #' Evaluation
  "pROC"
)

pacman::p_load(char = pkgs)

#' Packages that should be loaded on parallel workers.
options(normality_parallel_packages = pkgs)


#' ============================================================
#' SECTION 1: SMALL HELPER FUNCTIONS
#' ============================================================

#' ------------------------------------------------------------
#' Evaluate an expression safely and return a fallback value on failure.
#' ------------------------------------------------------------
safe_calc <- function(expr, default = NA_real_) {
  #' Return the fallback value when the calculation fails.
  out <- tryCatch(expr, error = function(e) default)

  #' Treat empty results as failed calculations.
  if (is.null(out) || length(out) == 0L) {
    return(default)
  }

  out
}


#' ------------------------------------------------------------
#' Create the output directory when it does not already exist.
#' ------------------------------------------------------------
make_output_dir <- function(output_dir) {
  #' Create nested directories when needed.
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  invisible(output_dir)
}


#' Create a base-R progress bar when progress reporting is enabled.
#'
#' The helper returns NULL when show_progress = FALSE, allowing the calling
#' function to use the same code in sequential and parallel runs.
start_progress_bar <- function(total, label = NULL, show_progress = TRUE) {
  if (!isTRUE(show_progress) || total < 1L) {
    return(NULL)
  }

  if (!is.null(label) && nzchar(label)) {
    cat(label, "\n")
  }

  utils::txtProgressBar(min = 0, max = total, style = 3)
}


#' Advance a progress bar without requiring repeated NULL checks.
update_progress_bar <- function(progress_bar, value) {
  if (!is.null(progress_bar)) {
    utils::setTxtProgressBar(progress_bar, value)
  }

  invisible(NULL)
}


#' Close a progress bar safely.
close_progress_bar <- function(progress_bar) {
  if (!is.null(progress_bar)) {
    close(progress_bar)
    cat("\n")
  }

  invisible(NULL)
}


#' Format elapsed time for readable progress messages.
format_elapsed_time <- function(start_time) {
  #' Convert elapsed seconds into a readable time label.
  elapsed_seconds <- as.numeric(
    difftime(Sys.time(), start_time, units = "secs")
  )
  if (elapsed_seconds < 60) {
    return(sprintf("%.1f seconds", elapsed_seconds))
  }
  if (elapsed_seconds < 3600) {
    return(sprintf("%.1f minutes", elapsed_seconds / 60))
  }
  sprintf("%.2f hours", elapsed_seconds / 3600)
}


#' ------------------------------------------------------------
#' Impute missing numeric values by column.
#' ------------------------------------------------------------
impute_na <- function(df, method = c("median", "mean")) {
  #' Match the requested imputation method.
  method <- match.arg(method)
  #' Process each column independently.
  for (nm in names(df)) {
    #' Skip nonnumeric columns and complete numeric columns.
    if (!is.numeric(df[[nm]]) || !anyNA(df[[nm]])) {
      next
    }
    #' Calculate the replacement value from observed entries.
    fill_value <- switch(
      method,
      median = median(df[[nm]], na.rm = TRUE),
      mean = mean(df[[nm]], na.rm = TRUE)
    )
    #' Use zero when an entire column is missing or invalid.
    if (!is.finite(fill_value)) {
      fill_value <- 0
    }
    #' Replace only the missing entries.
    df[[nm]][is.na(df[[nm]])] <- fill_value
  }
  df
}


#' ============================================================
#' CROSS-PLATFORM PARALLEL HELPER
#' ============================================================
#' - macOS/Linux uses parallel::mclapply();
#' - Windows uses a PSOCK cluster with parallel::parLapplyLB();
#' - n_cores = 1 always falls back to ordinary lapply().
#'
#' The PSOCK branch loads the required packages and exports functions
#' defined in this script to each worker.
#' ============================================================

get_available_cores <- function() {
  cores <- parallel::detectCores(logical = TRUE)
  if (!is.finite(cores) || is.na(cores) || cores < 1L) {
    return(1L)
  }
  as.integer(cores)
}


#' ------------------------------------------------------------
#' Validate and limit the requested number of processing cores.
#' ------------------------------------------------------------
normalize_n_cores <- function(n_cores, max_tasks = NULL) {
  #' Convert the request to an integer core count.
  n_cores <- suppressWarnings(as.integer(n_cores))
  if (!is.finite(n_cores) || is.na(n_cores) || n_cores < 1L) {
    n_cores <- 1L
  }
  #' Do not create more workers than available tasks.
  if (!is.null(max_tasks)) {
    max_tasks <- suppressWarnings(as.integer(max_tasks))
    if (is.finite(max_tasks) && !is.na(max_tasks) && max_tasks > 0L) {
      n_cores <- min(n_cores, max_tasks)
    }
  }

  max(1L, n_cores)
}


#' ------------------------------------------------------------
#' Choose a conservative default number of processing cores.
#' ------------------------------------------------------------
default_parallel_cores <- function(max_tasks = NULL, reserve_cores = 1L) {
  #' Detect available logical cores.
  available <- get_available_cores()
  #' Leave requested cores free for the operating system.
  reserve_cores <- max(0L, as.integer(reserve_cores))
  normalize_n_cores(
    n_cores   = max(1L, available - reserve_cores),
    max_tasks = max_tasks
  )
}


#' ------------------------------------------------------------
#' Identify the parallel backend used on the current operating system.
#' ------------------------------------------------------------
get_parallel_backend <- function(n_cores = 1L) {
  #' Validate the requested worker count.
  n_cores <- normalize_n_cores(n_cores)
  #' Use ordinary sequential processing for one core.
  if (n_cores == 1L) {
    return("sequential")
  }
  #' Windows requires a PSOCK cluster.
  if (.Platform$OS.type == "windows") {
    return("PSOCK")
  }
  "fork"
}


#' ------------------------------------------------------------
#' Identify functions and small constants to export to PSOCK workers.
#' ------------------------------------------------------------
get_parallel_export_names <- function(envir = environment(safe_mclapply)) {
  object_names <- ls(envir = envir, all.names = TRUE)
  #' Export functions and small atomic constants only.
  #' Large fitted objects should be passed through function arguments or
  #' closures, not exported globally to every worker.
  keep <- vapply(
    object_names,
    function(nm) {
      obj <- get(nm, envir = envir)
      is.function(obj) || (is.atomic(obj) && length(obj) <= 1000L)
    },
    logical(1L)
  )
  object_names[keep]
}


#' ------------------------------------------------------------
#' Load the required packages on a parallel worker.
#' ------------------------------------------------------------
load_parallel_packages <- function(packages = getOption("normality_parallel_packages", character(0))) {
  #' Remove duplicate package names.
  packages <- unique(as.character(packages))
  packages <- packages[nzchar(packages)]
  if (length(packages) == 0L) {
    return(invisible(character(0)))
  }
  #' Calculate loaded.
  loaded <- vapply(
    packages,
    function(pkg) {
      suppressPackageStartupMessages( require(pkg, character.only = TRUE, quietly = TRUE))
    },
    logical(1L)
  )
  #' Identify packages that failed to load.
  missing_pkgs <- packages[!loaded]

  if (length(missing_pkgs) > 0L) {
    stop("The following packages could not be loaded on a parallel worker: ",
      paste(missing_pkgs, collapse = ", ")
    )
  }
  invisible(packages)
}


#' ------------------------------------------------------------
#' Apply a function sequentially or in parallel across operating systems.
#' ------------------------------------------------------------
safe_mclapply <- function(X, FUN, n_cores = 1, ...,
                          .packages   = getOption("normality_parallel_packages", character(0)),
                          .export_env = environment(safe_mclapply),
                          .seed       = 12345) {

  n_tasks <- length(X)
  if (n_tasks == 0L) {
    return(list())
  }
  n_cores <- normalize_n_cores(
    n_cores   = n_cores,
    max_tasks = n_tasks
  )
  #' Sequential fallback. This is used whenever n_cores = 1.
  if (n_cores == 1L) {
    return(lapply(X, FUN, ...))
  }
  #' ------------------------------------------------------------
  #' macOS/Linux branch: forked parallelism
  #' ------------------------------------------------------------
  #' mclapply() is efficient on Unix-like systems because workers
  #' inherit the current R session through forking.
  #' ------------------------------------------------------------
  if (.Platform$OS.type != "windows") {

    #' Return the result.
    return(parallel::mclapply(X, 
                              FUN, ..., 
                              mc.cores = n_cores, 
                              mc.set.seed = TRUE))
  }

  #' ------------------------------------------------------------
  #' Windows branch: PSOCK parallelism
  #' ------------------------------------------------------------
  #' Windows cannot fork. A PSOCK cluster starts clean R sessions,
  #' so packages and functions must be explicitly made available.
  #' ------------------------------------------------------------
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(
    parallel::stopCluster(cl),
    add = TRUE
  )
  #' Reproducible independent RNG streams on PSOCK workers.
  parallel::clusterSetRNGStream(cl, iseed = .seed)
  #' Load the same package set on every worker.
  .parallel_packages <- unique(as.character(.packages))
  parallel::clusterExport(cl = cl, varlist = ".parallel_packages", envir   = environment()
  )

  parallel::clusterEvalQ(
    cl,
    {
      options(repos = c(CRAN = "https://cloud.r-project.org"))
      options(mc_doScale_quiet = TRUE)
      #' Calculate loaded.
      loaded <- vapply(
        .parallel_packages,
        function(pkg) {
          suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))
        },
        logical(1L)
      )
      missing_pkgs <- .parallel_packages[!loaded]
      #' Validate the condition.
      if (length(missing_pkgs) > 0L) {
        stop("The following packages could not be loaded on a parallel worker: ",
          paste(missing_pkgs, collapse = ", ")
        )
      }
      invisible(TRUE)
    }
  )

  #' Export framework functions and small constants to PSOCK workers.
  export_names <- get_parallel_export_names(.export_env)
  #' Validate the condition.
  if (length(export_names) > 0L) {
    parallel::clusterExport(cl = cl, varlist = export_names, envir = .export_env)
  }

  #' Load-balanced scheduling is safer when some features or sample sizes
  #' take longer than others.
  parallel::parLapplyLB(cl  = cl, X   = X,fun = FUN,...)
}

#' ---------------------------------------------------------------------------
#' Evaluate model-training code with a temporary caret parallel backend.
#'
#' Caret parallelism is used inside cross-validation and tuning. Sample sizes
#' should be processed sequentially while this backend is active to avoid
#' nested parallelism and excessive memory use.
#' ------------------------------------------------------------------------
with_caret_parallel <- function(
    n_cores,
    expr,
    seed = 12345,
    packages = getOption("normality_parallel_packages", character(0))
) {
  n_cores <- normalize_n_cores(n_cores)

  if (n_cores == 1L) {
    return(force(expr))
  }

  if (!requireNamespace("doParallel", quietly = TRUE) ||
      !requireNamespace("foreach", quietly = TRUE)) {
    stop("Packages 'doParallel' and 'foreach' are required for caret parallelism.")
  }

  cluster <- parallel::makeCluster(n_cores, type = "PSOCK")

  #' Continue the calculation.
  on.exit(
    {
      parallel::stopCluster(cluster)
      foreach::registerDoSEQ()
    },
    add = TRUE
  )

  parallel::clusterSetRNGStream(cluster, iseed = seed)

  worker_packages <- unique(c(as.character(packages), "caret"))
  worker_packages <- worker_packages[nzchar(worker_packages)]

  parallel::clusterExport(
    cluster,
    varlist = "worker_packages",
    envir = environment()
  )

  #' Continue the calculation.
  parallel::clusterEvalQ(
    cluster,
    {
      loaded <- vapply(
        worker_packages,
        function(package) {
          suppressPackageStartupMessages(
            require(package, character.only = TRUE, quietly = TRUE))
        },
        logical(1)
      )

      #' Validate the condition.
      if (any(!loaded)) {
        stop("Could not load worker package(s): ",
          paste(worker_packages[!loaded], collapse = ", ")
        )
      }

      invisible(TRUE)
    }
  )

  doParallel::registerDoParallel(cluster)
  force(expr)
}


#' ============================================================
#' SECTION 2: DEFAULT DISTRIBUTION SETS
#' ============================================================
#'
#' Each alternative distribution will be paired with its own
#' moment-matched Normal distribution.
#'
#' A user can add a new distribution in two ways:
#' 1. use a supported distribution name and parameters;
#' 2. provide normal_par = c(mean, sd) manually.
#' ============================================================


default_training_specs <- function() {
  list(
    #' alternatives.
    list(dist = "log_gamma", par = c(0.75, 1)),
    list(dist = "log_gamma", par = c(1.25, 1)),
    list(dist = "log_gamma", par = c(1.50, 1)),
    list(dist = "chi_square",  par = 7),
    list(dist = "exponential", par = 3),
    list(dist = "weibull",     par = c(2, 1)),
    list(dist = "laplace",     par = c(0, 4)),
    list(dist = "beta",        par = c(8, 2)),
    list(dist = "uniform",     par = c(0, 10)),
    list(dist = "t",           par = 3),
    list(dist = "f",           par = c(6, 15)),
    list(dist = "cauchy",      par = c(0, 1))
  )
}


#' ------------------------------------------------------------
#' Return the default distributions used for final model evaluation.
#' ------------------------------------------------------------
default_eval_specs <- function() {
  list(
    list(dist = "gumbel",       par = c(0, 1)),
    list(dist = "beta",         par = c(8, 2)),
    list(dist = "chi_square",   par = 3),
    list(dist = "laplace",      par = c(0, 4)),
    list(dist = "uniform",      par = c(0, 1)),
    list(dist = "cauchy",       par = c(0, 1)),
    list(dist = "contaminated", par = c(0.75, 0, 1, 5)),
    list(dist = "t",            par = 3)
  )
}


#' ============================================================
#' SECTION 3.1: MOMENT-MATCHED NORMAL SPECIFICATION
#' ============================================================
#'
#' moment_matched_normal()
#' Constructs the Normal distribution paired with a given alternative
#' distribution. The returned object has the form
#' list(dist = "normal", par = c(mean, sd)).
#'
#' Matching priority:
#' 1. If normal_par is supplied, use it directly.
#' 2. If the distribution is supported, use its exact theoretical
#' mean and standard deviation.
#' 3. If rfun is supplied, estimate the mean and standard deviation
#' by simulation.
#' 4. Otherwise, stop and require normal_par.
#'
#' Special case:
#' For Cauchy distributions, the mean and variance do not exist.
#' In that case, matching is based on the median and IQR rather than
#' the mean and standard deviation.
#'
#' ============================================================

moment_matched_normal <- function(spec, n_match = 1e6, seed = 12345) {

  #' Use supplied Normal parameters when available.
  if (!is.null(spec$normal_par)) {
    return(list(dist = "normal", par = spec$normal_par))
  }

  #' Estimate matching moments from a custom generator.
  if (!is.null(spec$rfun) && is.function(spec$rfun)) {
    #' Reproducibility.
    set.seed(seed)
    #' Simulated reference sample.
    x <- spec$rfun(n_match)
    #' Estimated mean and SD.
    return(list(dist = "normal", par = c(mean(x), sd(x))))
  }

  #' Standardized distribution name and parameters.
  dist <- tolower(trimws(spec$dist))
  par <- spec$par
  #' Calculate the matching Normal mean and SD.
  matched_par <- switch(
    dist,
    
    normal = {
      #' Existing Normal mean and SD.
      c(par[1], par[2])
    },

    gumbel = {
      #' Location and scale.
      loc <- par[1]
      scale <- par[2]
      #' Theoretical mean and SD.
      c(loc + scale * 0.5772156649, pi * scale / sqrt(6))
    },

    chi_square = {
      #' Degrees of freedom.
      df <- par[1]
      #' Theoretical mean and SD.
      c(df, sqrt(2 * df))
    },

    gamma = {
      #' Shape and rate.
      shape <- par[1]
      rate <- par[2]
      #' Theoretical mean and SD.
      c(shape / rate, sqrt(shape) / rate)
    },

    reflected_gamma = {
      #' Shape and rate.
      shape <- par[1]
      rate <- par[2]
      #' Reflected mean and unchanged SD.
      c(-shape / rate, sqrt(shape) / rate)
    },

    exponential = {
      #' Rate parameter.
      rate <- par[1]
      #' Theoretical mean and SD.
      c(1 / rate, 1 / rate)
    },

    reflected_exponential = {
      #' Rate parameter.
      rate <- par[1]
      #' Reflected mean and unchanged SD.
      c(-1 / rate, 1 / rate)
    },

    weibull = {
      #' Shape and scale.
      shape <- par[1]
      scale <- par[2]
      #' Theoretical mean.
      mu <- scale * gamma(1 + 1 / shape)
      #' Theoretical SD.
      sd <- scale * sqrt(gamma(1 + 2 / shape) - gamma(1 + 1 / shape)^2)
      c(mu, sd)
    },

    reflected_weibull = {
      #' Shape and scale.
      shape <- par[1]
      scale <- par[2]
      #' Original mean and SD.
      mu <- scale * gamma(1 + 1 / shape)
      sd <- scale * sqrt(gamma(1 + 2 / shape) - gamma(1 + 1 / shape)^2)
      #' Reflected mean and unchanged SD.
      c(-mu, sd)
    },

    laplace = {
      #' Location and scale.
      location <- par[1]
      scale <- par[2]
      #' Theoretical mean and SD.
      c(location, sqrt(2) * scale)
    },

    beta = {
      #' Shape parameters.
      a <- par[1]
      b <- par[2]
      #' Theoretical mean.
      mu <- a / (a + b)
      #' Theoretical SD.
      sd <- sqrt(a * b / ((a + b)^2 * (a + b + 1)))
      c(mu, sd)
    },

    uniform = {
      #' Lower and upper bounds.
      a <- par[1]
      b <- par[2]
      #' Theoretical mean and SD.
      c((a + b) / 2, (b - a) / sqrt(12))
    },

    logistic = {
      #' Location and scale.
      location <- par[1]
      scale <- par[2]
      #' Theoretical mean and SD.
      c(location, scale * pi / sqrt(3))
    },

    lognormal = {
      #' Log-scale mean and SD.
      meanlog <- par[1]
      sdlog <- par[2]
      #' Original-scale mean.
      mu <- exp(meanlog + sdlog^2 / 2)
      #' Original-scale SD.
      sd <- sqrt((exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2))
      c(mu, sd)
    },

    reflected_lognormal = {
      #' Log-scale mean and SD.
      meanlog <- par[1]
      sdlog <- par[2]
      #' Original mean and SD.
      mu <- exp(meanlog + sdlog^2 / 2)
      sd <- sqrt((exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2))
      #' Reflected mean and unchanged SD.
      c(-mu, sd)
    },

    t = {
      #' Degrees of freedom.
      df <- par[1]
      #' Finite-variance requirement.
      if (df <= 2) {
        stop("t distribution needs df > 2 for finite variance.")
      }
      #' Zero mean and theoretical SD.
      c(0, sqrt(df / (df - 2)))
    },

    f = {
      #' Numerator and denominator degrees of freedom.
      df1 <- par[1]
      df2 <- par[2]
      #' Finite-variance requirement.
      if (df2 <= 4) {
        stop("F distribution needs df2 > 4 for finite variance.")
      }

      #' Theoretical mean.
      mu <- df2 / (df2 - 2)
      #' Theoretical SD.
      sd <- sqrt(2 * df2^2 * (df1 + df2 - 2) / (df1 * (df2 - 2)^2 * (df2 - 4)))
      c(mu, sd)
    },

    pareto = {
      #' Shape and minimum value.
      alpha <- par[1]
      xm <- par[2]

      #' Finite-variance requirement.
      if (alpha <= 2) {
        stop("Pareto distribution needs alpha > 2 for finite variance.")
      }

      #' Theoretical mean and SD.
      mu <- xm * alpha / (alpha - 1)
      sd <- xm * sqrt(alpha / ((alpha - 1)^2 * (alpha - 2)))
      c(mu, sd)
    },

    cauchy = {
      #' Location and scale.
      location <- par[1]
      scale <- par[2]
      #' Match the median and IQR because moments do not exist.
      c(location, scale / qnorm(0.75))
    },

    contaminated = {
      #' Mixture weight, center, and component SDs.
      p <- par[1]
      mu <- par[2]
      sd1 <- par[3]
      sd2 <- par[4]

      #' Mixture SD around the common mean.
      sd_mix <- sqrt(p * sd1^2 + (1 - p) * sd2^2)
      c(mu, sd_mix)
    },

    log_gamma = {
      #' Shape and rate.
      shape <- par[1]
      rate <- par[2]
      #' Positive finite parameters required.
      if (length(par) != 2L || !is.finite(shape) || !is.finite(rate) || shape <= 0 || rate <= 0) {
        stop("For log_gamma, par must be c(shape, rate), with both parameters greater than zero.")
      }
      #' Mean and SD of X = -log(Y), where Y follows Gamma(shape, rate).
      c(log(rate) - digamma(shape), sqrt(trigamma(shape)))
    },

    #' Unsupported distribution.
    stop(
      "No automatic moment-matching rule for distribution: ", spec$dist, "\n",
      "Add normal_par = c(mean, sd) to the specification."
    )
  )

  #' Return the matched Normal specification.
  list(dist = "normal", par = matched_par)
}


#' ============================================================
#' SECTION 3.2:MOMENT-MATCHED DISTRIBUTION PAIRS
#' ============================================================
#'
#' Construct Normal-vs-alternative distribution pairs for training,
#' evaluation, or ROC comparison.
#'
#' For each non-normal alternative specification, this function creates
#' a corresponding Normal distribution with the same theoretical mean
#' and standard deviation.
#'
#' Args:
#' alt_specs : A list of alternative distribution specifications.
#'
#' Returns:
#' A list of paired specifications. Each element contains:
#' normal : the moment-matched Normal specification.
#' alt    : the original alternative distribution specification.
#'
#' Notes:
#' The matched Normal distribution is created by moment_matched_normal(),
#' which in turn uses get_distribution_moments() to obtain the
#' theoretical mean and standard deviation of the alternative.
#'
#' ============================================================
make_paired_specs <- function(alt_specs) {

  #' Continue the calculation.
  lapply(alt_specs, function(alt) {
    list(
      normal = moment_matched_normal(alt),
      alt    = alt
    )
  })
}

#' ============================================================
#' SECTION 3.3:THEORETICAL DISTRIBUTION MOMENTS
#' ============================================================
#'
#' Return the theoretical mean and standard deviation of a specified
#' distribution under the parameterization used in this framework.
#'
#' This function is used to support moment-matched comparisons between
#' each non-normal alternative distribution and a Normal distribution
#' with the same mean and standard deviation.
#'
#' Args:
#' dist : Character string giving the distribution name.
#' par  : Numeric vector of distribution parameters.
#'
#' Returns:
#' A list with two elements:
#' mean : theoretical mean of the distribution.
#' sd   : theoretical standard deviation of the distribution.
#'
#' Notes:
#' The returned second moment is the standard deviation, not the
#' variance. These values are used by moment_matched_normal() to
#' construct the corresponding Normal(mean, sd) distribution.
#'
#' ============================================================
get_distribution_moments <- function(dist, par) {

  #' Standardize the distribution name.
  dist <- tolower(trimws(dist))

  #' Apply the distribution-specific moment formulas.
  switch(
    dist,

    normal = {
      #' Location, scale, and Normal quartiles.
      list(
        mean = par[1], median = par[1], sd = par[2],
        iqr = qnorm(0.75, par[1], par[2]) - qnorm(0.25, par[1], par[2])
      )
    },

    gumbel = {
      #' Location and scale.
      loc <- par[1]
      scale <- par[2]
      #' Theoretical moments and quartile spread.
      list(
        mean = loc + scale * 0.5772156649,
        median = loc - scale * log(log(2)),
        sd = pi * scale / sqrt(6),
        iqr = qevd_gumbel(0.75, loc, scale) - qevd_gumbel(0.25, loc, scale)
      )
    },

    chi_square = {
      #' Degrees of freedom.
      df <- par[1]
      #' Theoretical moments and quantile-based spread.
      list(mean = df, median = qchisq(0.50, df), sd = sqrt(2 * df),
        iqr = qchisq(0.75, df) - qchisq(0.25, df))
    },

    gamma = {
      #' Shape and rate.
      shape <- par[1]
      rate <- par[2]

      #' Theoretical moments and Gamma quartiles.
      list(
        mean = shape / rate,
        median = qgamma(0.50, shape = shape, rate = rate),
        sd = sqrt(shape) / rate,
        iqr = qgamma(0.75, shape = shape, rate = rate) - qgamma(0.25, shape = shape, rate = rate))
    },

    reflected_gamma = {
      #' Shape and rate.
      shape <- par[1]
      rate <- par[2]

      #' Negate location measures; retain spread measures.
      list(
        mean = -shape / rate,
        median = -qgamma(0.50, shape = shape, rate = rate),
        sd = sqrt(shape) / rate,
        iqr = qgamma(0.75, shape = shape, rate = rate) - qgamma(0.25, shape = shape, rate = rate))
    },

    exponential = {
      #' Rate parameter.
      rate <- par[1]
      #' Theoretical moments and Exponential quartiles.
      list(
        mean = 1 / rate, median = log(2) / rate, sd = 1 / rate,
        iqr = qexp(0.75, rate = rate) - qexp(0.25, rate = rate))
    },

    reflected_exponential = {
      #' Rate parameter.
      rate <- par[1]
      #' Negate location measures; retain spread measures.
      list(
        mean = -1 / rate, median = -log(2) / rate, sd = 1 / rate,
        iqr = qexp(0.75, rate = rate) - qexp(0.25, rate = rate))
    },

    weibull = {
      #' Shape and scale.
      shape <- par[1]
      scale <- par[2]
      #' Theoretical moments and Weibull quartiles.
      list(
        mean = scale * gamma(1 + 1 / shape),
        median = scale * log(2)^(1 / shape),
        sd = scale * sqrt(gamma(1 + 2 / shape) - gamma(1 + 1 / shape)^2),
        iqr = qweibull(0.75, shape = shape, scale = scale) - qweibull(0.25, shape = shape, scale = scale))
    },

    reflected_weibull = {
      #' Shape and scale.
      shape <- par[1]
      scale <- par[2]

      #' Original Weibull mean and SD.
      mu <- scale * gamma(1 + 1 / shape)
      sd <- scale * sqrt(gamma(1 + 2 / shape) - gamma(1 + 1 / shape)^2)
      #' Negate location measures; retain spread measures.
      list(
        mean = -mu,
        median = -scale * log(2)^(1 / shape),
        sd = sd,
        iqr = qweibull(0.75, shape = shape, scale = scale) - qweibull(0.25, shape = shape, scale = scale))
    },

    laplace = {
      #' Location and scale.
      location <- par[1]
      scale <- par[2]
      #' Closed-form Laplace moments.
      list(mean = location, median = location, sd = sqrt(2) * scale, iqr = 2 * scale * log(2))
    },

    beta = {
      #' Two shape parameters.
      a <- par[1]
      b <- par[2]
      #' Theoretical moments and Beta quartiles.
      list(
        mean = a / (a + b),
        median = qbeta(0.50, a, b),
        sd = sqrt(a * b / ((a + b)^2 * (a + b + 1))),
        iqr = qbeta(0.75, a, b) - qbeta(0.25, a, b)
      )
    },

    uniform = {
      #' Lower and upper bounds.
      a <- par[1]
      b <- par[2]
      #' Closed-form Uniform moments.
      list(mean = (a + b) / 2, median = (a + b) / 2, sd = (b - a) / sqrt(12), iqr = (b - a) / 2)
    },

    logistic = {
      #' Location and scale.
      location <- par[1]
      scale <- par[2]
      #' Theoretical moments and Logistic quartiles.
      list(
        mean = location, median = location, sd = scale * pi / sqrt(3),
        iqr = qlogis(0.75, location, scale) - qlogis(0.25, location, scale))
    },

    lognormal = {
      #' Log-scale parameters.
      meanlog <- par[1]
      sdlog <- par[2]

      #' Theoretical moments and Lognormal quartiles.
      list(
        mean = exp(meanlog + sdlog^2 / 2),
        median = exp(meanlog),
        sd = sqrt((exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2)),
        iqr = qlnorm(0.75, meanlog, sdlog) - qlnorm(0.25, meanlog, sdlog))
    },

    reflected_lognormal = {
      #' Log-scale parameters.
      meanlog <- par[1]
      sdlog <- par[2]
      #' Original Lognormal mean and SD.
      mu <- exp(meanlog + sdlog^2 / 2)
      sd <- sqrt((exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2))
      #' Negate location measures; retain spread measures.
      list(
        mean = -mu, median = -exp(meanlog), sd = sd,
        iqr = qlnorm(0.75, meanlog, sdlog) - qlnorm(0.25, meanlog, sdlog))
    },

    t = {
      #' Degrees of freedom.
      df <- par[1]
      #' Mean and SD exist only above their respective thresholds.
      list(
        mean = if (df > 1) 0 else NA_real_,
        median = 0,
        sd = if (df > 2) sqrt(df / (df - 2)) else NA_real_,
        iqr = qt(0.75, df) - qt(0.25, df))
    },

    f = {
      #' Numerator and denominator degrees of freedom.
      df1 <- par[1]
      df2 <- par[2]
      #' Mean requires df2 > 2; SD requires df2 > 4.
      list(
        mean = if (df2 > 2) df2 / (df2 - 2) else NA_real_,
        median = qf(0.50, df1, df2),
        sd = if (df2 > 4) sqrt(2 * df2^2 * (df1 + df2 - 2) / (df1 * (df2 - 2)^2 * (df2 - 4))) else NA_real_,
        iqr = qf(0.75, df1, df2) - qf(0.25, df1, df2))
    },

    pareto = {
      #' Shape and minimum value.
      alpha <- par[1]
      xm <- par[2]
      #' Mean requires alpha > 1; SD requires alpha > 2.
      list(
        mean = if (alpha > 1) xm * alpha / (alpha - 1) else NA_real_,
        median = xm * 2^(1 / alpha),
        sd = if (alpha > 2) xm * sqrt(alpha / ((alpha - 1)^2 * (alpha - 2))) else NA_real_,
        iqr = xm * 0.25^(-1 / alpha) - xm * 0.75^(-1 / alpha))
    },

    cauchy = {
      #' Location and scale.
      location <- par[1]
      scale <- par[2]
      #' Cauchy mean and SD are undefined; median and IQR remain finite.
      list(mean = NA_real_, median = location, sd = NA_real_, iqr = 2 * scale)
    },

    contaminated = {
      #' Mixture weight, center, and component SDs.
      p <- par[1]
      mu <- par[2]
      sd1 <- par[3]
      sd2 <- par[4]

      #' Common-center mixture SD.
      sd_mix <- sqrt(p * sd1^2 + (1 - p) * sd2^2)
      #' Use the mixture center and Normal-scale IQR approximation.
      list(
        mean = mu, median = mu, sd = sd_mix,
        iqr = qnorm(0.75, mu, sd_mix) - qnorm(0.25, mu, sd_mix))
    },

    log_gamma = {
      #' Gamma shape and rate.
      shape <- par[1]
      rate <- par[2]
      #' Positive finite parameters required.
      if (length(par) != 2L || !is.finite(shape) || !is.finite(rate) || shape <= 0 || rate <= 0) {
        stop("For log_gamma, par must be c(shape, rate), with both parameters greater than zero.")
      }

      #' Reverse probabilities because X = -log(Y) is decreasing.
      q25 <- -log(qgamma(0.75, shape = shape, rate = rate))
      q50 <- -log(qgamma(0.50, shape = shape, rate = rate))
      q75 <- -log(qgamma(0.25, shape = shape, rate = rate))

      #' Transformation moments and IQR.
      list(mean = log(rate) - digamma(shape), median = q50, sd = sqrt(trigamma(shape)), iqr = q75 - q25)
    },

    #' Reject unsupported names.
    stop("Unsupported distribution in get_distribution_moments(): ", dist)
  )
}


#' ------------------------------------------------------------
#' Compute Gumbel quantiles using the framework parameterization.
#' ------------------------------------------------------------
qevd_gumbel <- function(p, loc, scale) {
  loc - scale * log(-log(p))
}


#' ============================================================
#' SECTION 4: DATA GENERATION
#' ============================================================
#'
#' generate_data()
#' Generates one sample from a named distribution.
#'
#' standardize = FALSE:
#' Returns the raw sample.
#'
#' standardize = TRUE:
#' Centers and scales using theoretical values where available.
#' If theoretical SD is unavailable, IQR is used.
#' ============================================================

generate_data <- function(n, dist, par = NULL, standardize = FALSE,
                          center_by = c("mean", "median")) {

  #' Standardize the requested distribution and centering rule.
  dist <- tolower(trimws(dist))
  center_by <- match.arg(center_by)

  #' Supply default parameters when none are provided.
  if (is.null(par)) {
    par <- switch(
      dist,
      normal       = c(0, 1),
      gumbel       = c(0, 1),
      log_gamma    = c(1.25, 1),
      chi_square   = 3,
      gamma        = c(3, 1),
      reflected_gamma = c(3, 1),
      exponential  = 1,
      reflected_exponential = 1,
      weibull      = c(2, 1),
      reflected_weibull = c(2, 1),
      laplace      = c(0, 1),
      beta         = c(2, 5),
      uniform      = c(0, 1),
      logistic     = c(0, 1),
      lognormal    = c(0, 1),
      reflected_lognormal = c(0, 1),
      t            = 3,
      f            = c(6, 15),
      pareto       = c(3, 1),
      cauchy       = c(0, 1),
      contaminated = c(0.75, 0, 1, 5),
      stop("Unsupported distribution: ", dist)
    )
  }

  #' Generate the requested random sample.
  x <- switch(
    dist,

    #' Normal sample.
    normal = {
      rnorm(n, mean = par[1], sd = par[2])
    },

    #' Gumbel sample.
    gumbel = {
      evd::rgumbel(n, loc = par[1], scale = par[2])
    },

    #' Chi-square sample.
    chi_square = {
      rchisq(n, df = par[1])
    },

    #' Gamma sample.
    gamma = {
      rgamma(n, shape = par[1], rate = par[2])
    },

    #' Reflected Gamma sample.
    reflected_gamma = {
      -rgamma(n, shape = par[1], rate = par[2])
    },

    #' Exponential sample.
    exponential = {
      rexp(n, rate = par[1])
    },

    #' Reflected Exponential sample.
    reflected_exponential = {
      -rexp(n, rate = par[1])
    },

    #' Weibull sample.
    weibull = {
      rweibull(n, shape = par[1], scale = par[2])
    },

    #' Reflected Weibull sample.
    reflected_weibull = {
      -rweibull(n, shape = par[1], scale = par[2])
    },

    #' Laplace sample.
    laplace = {
      LaplacesDemon::rlaplace( n, location = par[1], scale    = par[2])
    },

    #' Beta sample.
    beta = {
      rbeta(n, shape1 = par[1], shape2 = par[2])
    },

    #' Uniform sample.
    uniform = {
      runif(n, min = par[1], max = par[2])
    },

    #' Logistic sample.
    logistic = {
      rlogis(n, location = par[1], scale = par[2])
    },

    #' Lognormal sample.
    lognormal = {
      rlnorm(n, meanlog = par[1], sdlog = par[2])
    },

    #' Reflected Lognormal sample.
    reflected_lognormal = {
      -rlnorm(n, meanlog = par[1], sdlog = par[2])
    },

    #' Student t sample.
    t = {
      rt(n, df = par[1])
    },

    #' F sample.
    f = {
      rf(n, df1 = par[1], df2 = par[2])
    },

    #' Pareto sample.
    pareto = {
      VGAM::rpareto(n, shape = par[1], scale = par[2])
    },

    #' Cauchy sample.
    cauchy = {
      rcauchy(n, location = par[1], scale = par[2])
    },

    #' Contaminated Normal sample.
    contaminated = {
      p   <- par[1]
      mu  <- par[2]
      sd1 <- par[3]
      sd2 <- par[4]

      #' Assign mixture components.
      group <- rbinom(n, size = 1, prob = p)
      sds <- ifelse(group == 1, sd1, sd2)
      rnorm(n, mean = mu, sd = sds)
    },

    #' Log-Gamma sample.
    log_gamma = {
      #' Extract the Gamma shape and rate parameters.
      shape <- par[1]
      rate  <- par[2]

      #' Validate the distribution parameters.
      if (length(par) != 2L || !is.finite(shape) || !is.finite(rate) ||
          shape <= 0 || rate <= 0) {
        stop("For log_gamma, par must be c(shape, rate), ","with both parameters greater than zero."
        )
      }

      #' Generate Gamma values and apply the negative-log transformation.
      gamma_values <- rgamma(n = n, shape = shape, rate  = rate)

      #' Protect against taking the logarithm of numerical zero.
      gamma_values <- pmax(gamma_values, .Machine$double.xmin)
      -log(gamma_values)
    },

    stop("Unsupported distribution: ", dist)
  )

  #' Return raw values when standardization is disabled.
  if (!standardize) {
    return(x)
  }

  #' Obtain the theoretical center and spread.
  moments <- get_distribution_moments(dist, par)

  center_value <- switch(
    center_by,
    mean   = moments$mean,
    median = moments$median
  )

  if (!is.finite(center_value)) {
    center_value <- moments$median
  }

  scale_value <- moments$sd

  if (!is.finite(scale_value) || scale_value <= 0) {
    scale_value <- moments$iqr
  }
  if (!is.finite(scale_value) || scale_value <= 0) {
    scale_value <- sd(x)
  }
  if (!is.finite(scale_value) || scale_value <= 0) {
    scale_value <- 1
  }

  x_std <- (x - center_value) / scale_value

  #' Continue the calculation.
  attr(x_std, "distribution") <- dist
  attr(x_std, "parameters")   <- par
  attr(x_std, "center_by")    <- center_by
  attr(x_std, "center_value") <- center_value
  attr(x_std, "scale_value")  <- scale_value
  attr(x_std, "standardized") <- TRUE

  x_std
}

#' ============================================================
#' SECTION 5: FEATURE ROW GENERATION
#' ============================================================
#'
#' This function draws repeated samples from one distribution,
#' computes features, adds the class label, and returns a data
#' frame ready for model training or evaluation.
#'
#' calculate_features() will be defined in Part 2.
#' ============================================================

#' Calculate generate feature rows.
generate_feature_rows <- function(n, n_rep, spec, label, feature_set = NULL,
    standardize_sample = FALSE, center_by = NULL, progress_bar = NULL,
    progress_offset = 0L
) {
  rows <- vector("list", n_rep)
  progress_every <- max(1L, floor(n_rep / 100L))
  #' Generate each sample and calculate only the requested features.
  for (iteration in seq_len(n_rep)) {
    x <- generate_data( n = n, dist = spec$dist, par = spec$par, 
                        standardize = standardize_sample, center_by = center_by)
    features <- calculate_features( x = x, feature_set = feature_set)
    features$Label <- label
    rows[[iteration]] <- features
    #' Update the shared training progress bar at approximately 1% intervals.
    if (!is.null(progress_bar) &&
        (iteration %% progress_every == 0L || iteration == n_rep)) {
      update_progress_bar(progress_bar, progress_offset + iteration)
    }
  }
  do.call(rbind, rows)
}

#' ============================================================
#' PART 2: FEATURE ENGINEERING
#' ============================================================
#'
#' This section computes numerical summaries of a sample.
#' Each row returned by calculate_features() represents one
#' simulated sample.
#'
#' The goal is not to decide normality from one statistic.
#' The goal is to give the ML models many complementary
#' signals about shape, tails, spacing, entropy, and Q-Q fit.
#' ============================================================

#' ============================================================
#' SECTION 6: BASIC FEATURE HELPERS
#' ============================================================

#' Standardize a numeric sample
#'
standardize_sample <- function(x) {
  #' Center and scale the sample using the selected theoretical summaries.
  x <- as.numeric(na.omit(x))
  n <- length(x)
  if (n < 2L) {
    return(rep(NA_real_, n))
  }
  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(rep(NA_real_, n))
  }
  as.numeric((x - mean(x)) / s)
}


#' Safely compute a ratio.
safe_ratio <- function(num, den) {
  if (!is.finite(num) || !is.finite(den)) {
    return(NA_real_)
  }
  if (abs(den) < .Machine$double.eps) {
    return(NA_real_)
  }
  num / den
}


#' Mean log spacing between consecutive order statistics.
mean_log_spacing <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  #' Consecutive order-statistic spacings.
  d <- diff(x)
  #' Log spacing requires strictly positive spacings.
  d <- d[d > 0]
  if (length(d) == 0L) {
    return(NA_real_)
  }
  mean(log(d))
}


#' ============================================================
#' SECTION 7: NORMALITY TEST STATISTICS
#' ============================================================

#' ------------------------------------------------------------
#' Calculate the classical normality-test statistics used as features.
#' ------------------------------------------------------------
normality_test_stats <- function(x) {

  #' Build the result table.
  data.frame(
    SW_stat  = safe_calc(as.numeric(shapiro.test(x)$statistic)),
    SF_stat  = safe_calc(as.numeric(nortest::sf.test(x)$statistic)),
    AD_stat  = safe_calc(as.numeric(nortest::ad.test(x)$statistic)),
    LF_stat  = safe_calc(as.numeric(nortest::lillie.test(x)$statistic)),
    KS_stat  = safe_calc(as.numeric(ks.test(x, "pnorm", mean(x), sd(x))$statistic)),
    JB_stat  = safe_calc(as.numeric(tseries::jarque.bera.test(x)$statistic))
  )
}


#' ============================================================
#' SECTION 8: ENTROPY AND SPACING FEATURES
#' ============================================================
#'
#' 1. Lin-Mudholkar Zp statistic.
#'
#' The statistic is based on the correlation between the ordered
#' sample and leave-one-out variance terms. It is designed to
#' capture departures from Normality through the relationship
#' between order statistics and local dispersion.

calculate_zp_statistic <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)

  if (n < 4L || length(unique(x)) < 4L) {
    return(NA_real_)
  }
  #' Leave-one-out variance terms.
  sum_x  <- sum(x)
  sum_x2 <- sum(x^2)
  n1     <- n - 1L
  u <- ((sum_x2 - x^2) - ((sum_x - x)^2 / n1))^(1 / 3)
  safe_calc(atanh(cor(x, u)), default = NA_real_
  )
}

#' ---------------------------------------------
#' 2. Vasicek entropy-based statistic.
#'
#' This statistic estimates entropy using order-statistic spacings.
calculate_vasicek_kmn <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  #' Need enough samples to compute shape/entropy
  if (n < 4L) {
    return(NA_real_)
  }
  m <- floor(sqrt(n))
  s <- sd(x)
  #' Constant data has zero entropy variance
  if (!is.finite(s) || s <= 0.0) {
    return(NA_real_)
  }
  #' Fast vectorized index bounding
  i <- seq_len(n)
  upper_idx <- pmin(n, i + m)
  lower_idx <- pmax(1L, i - m)
  spacings <- x[upper_idx] - x[lower_idx]
  #' Handle heavy ties or duplicates gracefully
  if (any(spacings <= 0.0)) {
    #' Add a tiny epsilon to prevent log(0) 
    spacings[spacings <= 0.0] <- 1e-7
  }
  #' Compute entropy value
  entropy_val <- (n / (2 * m * s)) * exp(mean(log(spacings)))
  #' check against Inf or NaN values
  if (!is.finite(entropy_val)) {
    return(NA_real_)
  }
  return(entropy_val)
}


#' --------------------------------------------------------
#' 3.  Renyi entropy using a histogram approximation.
#'
#' This summarizes how concentrated or dispersed the empirical
#' distribution is across bins. It is included as a broad shape
#' feature rather than as a formal Normality test.
#' alpha = 0 Hartley Entropy
#' alpha = 1 Shannon Entropy
#' alpha = 2 Collision Entropy
#' alpha = infinity Min- Entropy
#' --------------------------------------------------------
renyi_entropy <- function(x, alpha = 2, bins = 10, global_range = NULL) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  #' Requires at least 3 data points to build a histogram
  if (n < 3L) {
    return(NA_real_)
  }
  #' Determine bin boundaries
  if (!is.null(global_range)) {
    breaks <- seq(global_range[1], global_range[2], length.out = bins + 1)
  } else {
    breaks <- bins
  }
  #' Compute histogram probabilities
  p <- table(cut(x, breaks = breaks, include.lowest = TRUE))
  p <- as.numeric(p) / n

  #' Filter out empty bins (0 * log(0) is 0 in entropy)
  p <- p[p > 0]
  if (length(p) == 0L) {
    return(NA_real_)
  }
  #' Handle special case alpha = 1 (Shannon Entropy limit)
  if (abs(alpha - 1) < 1e-9) {
    entropy_val <- -sum(p * log(p))
  } else {
    entropy_val <- (1 / (1 - alpha)) * log(sum(p^alpha))
  }
  #' check for numerical instability
  if (!is.finite(entropy_val)) {
    return(NA_real_)
  }
  return(entropy_val)
}


#' --------------------------------------------------------
#' 4. Approximate negentropy measures departure from Gaussianity
#' Pre-computed reference constant for negentropy_approx().
#' --------------------------------------------------------
negentropy_normal_ref <- local({
  set.seed(12345)
  #' Pre-compute with a large draw .
  mean(log(cosh(rnorm(1e6L))))
})

#' Approximate negentropy.
negentropy_approx <- function(x) {
  z <- standardize_sample(x)
  if (length(z) < 3L) {
    return(NA_real_)
  }
  #' log(cosh(z)) is mathematically equal to |z| - log(2) + log(1 + exp(-2*|z|))
  abs_z <- abs(z)
  log_cosh_z <- abs_z - log(2) + log(1 + exp(-2 * abs_z))
  #' Execute calculation
  raw_negentropy <- safe_calc(
    (mean(log_cosh_z) - negentropy_normal_ref)^2,
    default = NA_real_
  )
  #' Catch any boundary failures
  if (!is.finite(raw_negentropy)) {
    return(NA_real_)
  }
  return(raw_negentropy)
}

#' --------------------------------------------------------------
#' 5. Variance and skewness of log spacings.
#'
#' Consecutive order-statistic spacings describe how observations
#' are distributed across the support. The variance and skewness
#' of log spacings can detect clustering, gaps, and tail behavior.
#' --------------------------------------------------------------
calc_spacing_ratio_stats <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  #' Compute consecutive differences
  d <- diff(x)
  d <- d[d > 0]
  #' Need at least 3 points to compute skewness of differences safely
  if (length(d) < 3L) {
    return(list(
      var_log_spacings = NA_real_,
      skew_log_spacings = NA_real_
    ))
  }
  log_d <- log(d)
  #' Calculate variance
  v_val <- var(log_d)
  if (!is.finite(v_val)) {
    v_val <- NA_real_
  }
  #' Calculate skewness using safe_calc
  s_val <- safe_calc(
    as.numeric(moments::skewness(log_d)),
    default = NA_real_
  )
  #' Intercept NaN/Inf returns
  if (!is.finite(s_val)) {
    s_val <- NA_real_
  }

  list(
    var_log_spacings = v_val,
    skew_log_spacings = s_val
  )
}


#' --------------------------------------------------------------
#' 6. Greenwood Spacing Statistic
#'
#' Measures clustering and spacing unevenness.
#' -------------------------------------------------------------
greenwood_stat <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  d <- diff(x)
  d <- d[d > 0]
  #'  need at least 3 valid intervals (4 points)
  n_spacings <- length(d)
  if (n_spacings < 3L) {
    return(NA_real_)
  }
  total_d <- sum(d)
  if (!is.finite(total_d) || total_d <= 0) {
    return(NA_real_)
  }
  #' Raw Greenwood statistic
  g <- sum(d^2) / (total_d^2)

  #' EXACT uniform moments for m = n_spacings internal intervals sum to 1:
  #' E[G] = 2 / (m + 1)
  #' Var(G) = 4 * (m - 1) / ((m + 1)^2 * (m + 2) * (m + 3))
  m <- n_spacings
  expected_g <- 2 / (m + 1)
  var_g <- (4 * (m - 1)) / (((m + 1)^2) * (m + 2) * (m + 3))

  if (var_g <= 0) {
    return(NA_real_)
  }

  #' Standardize using safe_ratio helper
  standardized_g <- safe_ratio(g - expected_g, sqrt(var_g))
  #' Catch any numerical boundary failures
  if (!is.finite(standardized_g)) {
    return(NA_real_)
  }
  return(standardized_g)
}

#' --------------------------------------------------------------
#' 7. Rao Spacing Statistic
#'
#' Measures the total absolute deviation of spacings from equal intervals.
#' --------------------------------------------------------------
rao_spacing_stat <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  #' Need at least 4 observations
  if (n < 4L) {
    return(NA_real_)
  }
  #' 1. Capture ALL n spacings
  internal_gaps <- diff(x)
  boundary_gap <- x[n] - x[1L]

  #' Filter out duplicate ties
  internal_gaps <- internal_gaps[internal_gaps > 0]
  if (length(internal_gaps) == 0L || boundary_gap <= 0) {
    return(NA_real_)
  }
  #' Combine internal intervals with the wrap-around interval
  d <- c(internal_gaps, boundary_gap)
  m <- length(d) #' Total number of valid spacing intervals
  total_d <- sum(d)
  if (!is.finite(total_d) || total_d <= 0) {
    return(NA_real_)
  }
  #' Compute true normalized spacings (must sum to 1.0)
  s <- d / total_d
  #' Compute Rao's raw statistic
  raw_rao <- sum(abs(s - (1 / m)))
  #' 4. Apply exact asymptotic limits for the uniform distribution:
  #' Expected Mean -> 2 / e (~0.7357589)
  #' Expected Variance -> (2*e - 5) / (e^2 * m)
  e_const <- exp(1)
  mean0 <- 2 / e_const
  var0 <- (2 * e_const - 5) / ((e_const^2) * m)

  if (var0 <= 0) {
    return(NA_real_)
  }
  #' Standardize into a clean Z-score feature
  standardized_rao <- safe_ratio(raw_rao - mean0, sqrt(var0))
  #' Catch any numerical instabilities or boundary failures
  if (!is.finite(standardized_rao)) {
    return(NA_real_)
  }
  return(standardized_rao)
}


#' ============================================================
#' SECTION 9: TAIL, OUTLIER, AND ROBUSTNESS FEATURES
#' ============================================================
#' -----------------------------------------------------------
#' 8. Tail Asymmetry Ratio Feature
#'
#' Compares the absolute structural span of the upper 10% tail
#' to the lower 10% tail relative to the median.
#' -----------------------------------------------------------
calculate_tail_asymmetry <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  #' Need at least 10 observations to extract a valid 10% tail slice (1 element each)
  if (n < 10L) {
    return(NA_real_)
  }
  #' Calculate the Median
  md <- median(x)
  #' Calculate the exact number of elements that fit into a 10% tail slice
  k <- floor(n * 0.10)
  if (k < 1L) {
    return(NA_real_)
  }
  #' Extract structurally identical slice lengths directly using sorted index positions
  lower_tail <- x[1:k]
  upper_tail <- x[(n - k + 1L):n]
  #' Calculate absolute structural distance from the center
  upper_excess <- mean(upper_tail) - md
  lower_deficit <- md - mean(lower_tail)
  #' If the data is completely flat or has zero tail variance, return NA_real_
  if (!is.finite(upper_excess) || !is.finite(lower_deficit) ||
      upper_excess <= 0 || lower_deficit <= 0) {
    return(NA_real_)
  }
  #' Compute the ratio
  ratio_val <- safe_ratio(upper_excess, lower_deficit)
  #' Catch any numerical instabilities or boundary failures
  if (!is.finite(ratio_val)) {
    return(NA_real_)
  }
  return(ratio_val)
}

#' ---------------------------------------------------------------
#' 9. Moors Robust Kurtosis Feature
#'
#' Evaluates peak flatness and tail weight using octile points.
#' For a standard Normal distribution, this metric equals exactly 0.0.
#' ----------------------------------------------------------------
calculate_moors_kurtosis <- function(x) {

  x <- as.numeric(na.omit(x))
  if (length(x) < 8L) {
    return(NA_real_)
  }
  #' Extract the 7 precise octile boundaries (12.5% up to 87.5%)
  probs <- seq(0.125, 0.875, by = 0.125)
  q <- quantile(x, probs = probs, names = FALSE)
  #' Moors Formula calculation
  num <- (q[7L] - q[5L]) + (q[3L] - q[1L])
  den <- q[6L] - q[2L]
  if (!is.finite(num) || !is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  #' Standardize so that a normal baseline shape maps perfectly to 0.0
  #' Raw normal value for the Moors layout is ~1.2331
  moors_stat <- (num / den) - 1.233093
  if (!is.finite(moors_stat)) {
    #' Handles extreme tie breakdowns safely
    return(NA_real_)
  }
  return(moors_stat)
}


#' ----------------------------------------------------------------
#' 10. Excess Tail Proportion Feature
#'
#' Calculates the proportion of standardized observations beyond +/- 2 SD
#' relative to the exact mathematical standard Normal reference value.
#' -------------------------------------------------------------------
calculate_excess_tail_prop <- function(x) {
  #' Standardize the vector to mean 0, SD 1
  z <- standardize_sample(x)
  #' Strip out NAs to get the true count of valid standardized elements
  z <- as.numeric(na.omit(z))
  n <- length(z)
  #' Need an adequate sample size to calculate a meaningful proportion
  if (n < 5L) {
    return(NA_real_)
  }
  #' Calculate the exact theoretical normal reference: P(|Z| > 2) = 0.04550026...
  normal_reference <- 2 * pnorm(-2)
  #' Compute the empirical proportion of exceedances
  observed_proportion <- mean(abs(z) > 2.0)
  #' Execute the ratio
  ratio_val <- safe_ratio(observed_proportion, normal_reference)
  #' Catch any numerical instabilities or boundary failures
  if (!is.finite(ratio_val)) {
    return(NA_real_)
  }
  return(ratio_val)
}


#' -------------------------------------------------------------------------
#' 11. Tukey Outlier Proportion Feature
#'
#' Calculates the proportion of observations falling
#' outside the Tukey boxplot fences.
#' ------------------------------------------------------------------------
calculate_outlier_proportion <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  #' Need at least 4 elements to establish distinct quartile bounds
  if (n < 4L) {
    return(NA_real_)
  }
  #' Extract quartiles exactly once
  q <- quantile(x, probs = c(0.25, 0.75), names = FALSE)
  #' Derive the robust scale interval directly
  interquartile_range <- q[2L] - q[1L]
  #' If the data has heavy ties or is a flat line, IQR is 0.
  if (!is.finite(interquartile_range) || interquartile_range <= 0) {
    return(NA_real_)
  }
  #' Compute the standard 1.5 Tukey fence multiplier
  h <- 1.5 * interquartile_range
  lower_fence <- q[1L] - h
  upper_fence <- q[2L] + h
  #' Calculate the empirical proportion of points breaching the fences
  outlier_prop <- mean(x < lower_fence | x > upper_fence)
  if (!is.finite(outlier_prop)) {
    return(NA_real_)
  }
  return(outlier_prop)
}


#' ------------------------------------------------------------
#' 12. Qn robust scale ratio
#' ------------------------------------------------------------
#'
#' Computes the ratio of the Qn robust scale estimator to the
#' sample standard deviation. Under normality, this ratio is
#' approximately 1. Values below 1 suggest heavier tails relative
#' to the normal distribution, while values above 1 suggest
#' lighter tails or more bounded behavior.
#' ------------------------------------------------------------

calc_qn_robust_spread <- function(x) {
  #' Scale the robust Qn estimator relative to the sample SD.
  x <- as.numeric(na.omit(x))
  if (length(x) < 5L) {
    return(NA_real_)
  }
  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  qn_val <- safe_calc(
    as.numeric(robustbase::Qn(x, finite.corr = TRUE)),
    default = NA_real_
  )
  if (!is.finite(qn_val) || qn_val <= 0) {
    return(NA_real_)
  }
  safe_ratio(qn_val, s)
}


#' -----------------------------------------------------------------
#' 13. Medcouple Robust Skewness Feature
#'
#' Computes the medcouple estimator of skewness. Robust to outliers
calc_medcouple <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  #' Medcouple requires at least 3 distinct observations
  if (n < 5L) {
    return(NA_real_)
  }
  #' If the most common number takes up more than 50% of the row,
  #' the median collapses into a tie, causing robustbase::mc to freeze or crash.
  max_tie_count <- max(table(x))
  if (max_tie_count > (n / 2.0)) {
    return(NA_real_)
  }
  #' Compute the medcouple metric safely
  mc_val <- safe_calc(
    as.numeric(robustbase::mc(x)),
    default = NA_real_
  )
  #' Catch any unmanaged numerical instabilities or boundary failures
  if (!is.finite(mc_val)) {
    return(NA_real_)
  }
  return(mc_val)
}

#' -----------------------------------------------------------------------
#' 14. Robust Left and Right Tail Weights (Brys Methodology)
#'
#' Computes independent left and right tail heaviness metrics using octile ratios.
calc_tail_weights_brys <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  #' Brys octile tail metrics require at least 8 elements
  if (n < 8L) {
    return(list(
      LTW = NA_real_,
      RTW = NA_real_
    ))
  }

  #' Extract the specific quantiles needed for the Brys tail formulas:
  probs <- c(0.03125, 0.0625, 0.25, 0.50, 0.75, 0.9375, 0.96875)
  q <- quantile(x, probs = probs, names = FALSE)
  #' Denominators: Inner body spread from the median on each side
  left_body  <- q[4L] - q[3L]  #' Median - Q(0.25)
  right_body <- q[5L] - q[4L]  #' Q(0.75) - Median
  #' Guard against completely flat core bodies or heavy ties
  if (left_body <= 0 || right_body <= 0) {
    return(list(
      LTW = NA_real_,
      RTW = NA_real_
    ))
  }

  #' TRUE BRYS FORMULAS:
  #' Left Tail Weight (LTW): Ratio of extreme left spread to inner left spread
  ltw_val <- (q[2L] - q[1L]) / left_body
  #' Right Tail Weight (RTW): Ratio of extreme right spread to inner right spread
  rtw_val <- (q[7L] - q[6L]) / right_body
  #' Catch any floating-point or boundary anomalies
  if (!is.finite(ltw_val) || ltw_val < 0) ltw_val <- NA_real_
  if (!is.finite(rtw_val) || rtw_val < 0) rtw_val <- NA_real_

  list(
    LTW = ltw_val,
    RTW = rtw_val
  )
}



#' ============================================================
#' SECTION 10: Q-Q AND CORRELATION FEATURES
#' ============================================================
#'
#' ------------------------------------------------------------
#' 15. Q-Q cubic coefficient
#' ------------------------------------------------------------
#'
#' Fits a degree-3 polynomial of theoretical Normal quantiles
#' to the ordered sample values and returns the cubic coefficient.
#' Under Normality, the Q-Q relationship is approximately linear;
#' a nonzero cubic coefficient captures systematic curvature.
#' ------------------------------------------------------------

qq_cubic_coef <- function(x) {
  #' Fit a cubic QQ relationship and retain the nonlinear coefficient.
  x <- as.numeric(na.omit(x))
  n <- length(x)
  if (n < 6L || length(unique(x)) < 4L) {
    return(NA_real_)
  }

  theoretical <- qnorm(ppoints(n))
  empirical   <- sort(x)
  #' Fit the model.
  fit <- tryCatch(
    suppressWarnings(
      lm(empirical ~ poly(theoretical, 3, raw = TRUE))
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(NA_real_)
  }
  coefs <- coef(fit)
  if (length(coefs) < 4L || !is.finite(coefs[4L])) {
    return(NA_real_)
  }
  as.numeric(coefs[4L])
}


#' ----------------------------------------------------------------------
#' 16. de Wet-Venter Weighted Correlation Feature
#'
#' Computes the de Wet-Venter normality statistic, standardized into a
#' sample-size invariant Z-score.
dewet_venter_stat <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  #' The statistic requires enough order statistics to compute valid weights
  if (n < 8L) {
    return(NA_real_)
  }

  #' Blom plotting positions
  p <- (seq_len(n) - 0.375) / (n + 0.25)
  m <- qnorm(p)
  #' Density weights
  w <- dnorm(m)
  w_sum <- sum(w)

  if (!is.finite(w_sum) || w_sum <= 0) {
    return(NA_real_)
  }
  
  #' Weighted centering
  x_bar <- sum(w * x) / w_sum
  m_bar <- sum(w * m) / w_sum
  x_dev <- x - x_bar
  m_dev <- m - m_bar

  #' Calculate raw weighted correlation coefficient (r)
  numerator <- sum(w * x_dev * m_dev)
  denom_x <- sum(w * x_dev^2)
  denom_m <- sum(w * m_dev^2)

  if (denom_x <= 0 || denom_m <= 0) {
    return(NA_real_)
  }

  r <- safe_ratio(num = numerator, den = sqrt(denom_x * denom_m))
  if (!is.finite(r)) {
    return(NA_real_)
  }

  #' Bound correlation strictly to mathematical limits
  r <- max(-1, min(1, r))
  r2 <- r^2
  #' Convert raw R^2 into a standardized test statistic.
  #' The raw statistic is L = n * (1 - R^2).
  raw_L <- n * (1 - r2)
  #' De Wet & Venter's asymptotic normal moment corrections for sample size n:
  #' E[L] ~ ln(n) + 0.5772 (Euler-Mascheroni constant) - 1
  #' Var(L) ~ pi^2 / 6 - 1
  expected_L <- log(n) + 0.5772156649 - 1.0
  var_L <- (pi^2 / 6.0) - 1.0

  if (var_L <= 0) {
    return(NA_real_)
  }

  #' Convert to a standardized scale-invariant Z-score
  standardized_z <- safe_ratio(raw_L - expected_L, sqrt(var_L))
  if (!is.finite(standardized_z)) {
    return(NA_real_)
  }
  return(standardized_z)
}


#' ------------------------------------------------------------------------
#' 17. Standardized Ryan-Joiner Candidate Feature
#'
#' Computes the Ryan-Joiner normality statistic, transformed into a
#' sample-size invariant Z-score. Evaluates global Q-Q plot linearity evenly.
#' ------------------------------------------------------------------------
ryan_joiner_stat <- function(x) {
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  #' Requires adequate order statistics to resolve the correlation structure
  if (n < 8L) {
    return(NA_real_)
  }
  #' Blom plotting positions and expected Normal order statistics
  p <- (seq_len(n) - 0.375) / (n + 0.25)
  m <- qnorm(p)
  #' Calculate raw correlation coefficient (r)
  r <- safe_calc(cor(x, m), default = NA_real_)
  if (is.na(r) || !is.finite(r)) {
    return(NA_real_)
  }
  #' Bound numerically and compute R^2
  r <- max(-1, min(1, r))
  r2 <- r^2
  #' Transform raw R^2 into a size-invariant Z-score.
  #' Based on the Shapiro-Francia asymptotic normal properties for log(1 - R^2).
  raw_W <- 1 - r2
  #' If raw_W is exactly 0 (perfect correlation), bound it to prevent log(0)
  if (raw_W <= 0) raw_W <- 1e-10
  log_W <- log(raw_W)

  #' Asymptotic mean and SD formulas for log(1 - R^2) given sample size n
  expected_log_W <- -log(n) - 0.5772156649 + 1.4808
  sd_log_W       <- sqrt(1.5707963268 / n) #' sqrt(pi / 2n)
  #' Convert to a standardized Z-score
  standardized_z <- safe_ratio(log_W - expected_log_W, sd_log_W)
  if (!is.finite(standardized_z)) {
    return(NA_real_)
  }
  return(standardized_z)
}


#' ============================================================
#' SECTION 11: ECF AND DISTANCE FEATURES
#' ============================================================
#' 18. Epps-Pulley Geometric Normality Invariant
#'
#' Measures omnibus deviation from a normal distribution shape.
#' Returns a bounded metric: 0 means perfectly normal baseline shape.
epps_pulley_stat <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  if (n < 4L) return(0.0)

  s <- sd(x)
  if (!is.finite(s) || s <= 0.0) return(0.0)
  #' Center and scale data to standard normal space
  z <- (x - mean(x)) / s
  #' Double-sum loop vectorized using outer-product matrices
  z_diff_sq <- outer(z, z, "-")^2
  z_sum_sq  <- outer(z, z, "+")^2
  #' Epps-Pulley integration terms
  term_1 <- sum(exp(-z_diff_sq / 2.0)) / (n^2)
  term_2 <- (2.0 / n) * sum(exp(-z^2 / 4.0)) / sqrt(1.5)
  term_3 <- 1.0 / sqrt(2.0)
  statistic <- n * (term_1 - term_2 + term_3)
  if (!is.finite(statistic) || statistic < 0.0) return(0.0)
  return(statistic)
}

#' --------------------------------------------------------------
#' 19. Empirical Characteristic Function Discrete Deviations
#'
#' Extracts localized spatial shape deviations from the real and imaginary
#' coordinates of the empirical characteristic function.
#' Empirical Characteristic Function Discrete Deviations
ecf_deviations <- function(x, t_vals = c(1.0, 2.0, 3.0)) {
  x <- as.numeric(na.omit(x))
  #' Generate explicit named slots for automated matrix column binding
  out_names <- c(
    paste0("ecf_re_t", seq_along(t_vals)),
    paste0("ecf_im_t", seq_along(t_vals))
  )
  #' Initialize the target NA vector layout
  na_return <- setNames(rep(NA_real_, length(out_names)), out_names)

  if (length(x) < 3L) {
    return(na_return)
  }
  #' Swapped '||' for '|' to evaluate vectors
  if (any(t_vals > 3.0 | t_vals <= 0.0)) {
    return(na_return)
  }

  #' Standardize sample to mean 0, variance 1
  z <- standardize_sample(x)
  z <- as.numeric(na.omit(z))
  #' If the row has no variance, standardization returns NA
  if (length(z) < 3L) {
    return(na_return)
  }

  #' Real deviations: Evaluates localized variance and kurtosis shapes
  real_dev <- vapply(
    t_vals, function(t) {
      val <- mean(cos(t * z)) - exp(-t^2 / 2)
      if (!is.finite(val)) NA_real_ else val
    }, numeric(1L)
  )

  #' Imaginary deviations: robust asymmetry metric (0.0 for true normal)
  imag_dev <- vapply(
    t_vals, function(t) {
      val <- mean(sin(t * z))
      if (!is.finite(val)) NA_real_ else val
    }, numeric(1L)
  )

  final_vector <- c(real_dev, imag_dev)
  #' Intercept any rare floating-point
  if (any(!is.finite(final_vector))) {
    return(na_return)
  }
  setNames(final_vector, out_names)
}




#' ============================================================
#' SECTION 12: OTHER SHAPE FEATURES
#' ============================================================
#' 20. Excess Moors Kurtosis Feature
#'
#' Evaluates peak flatness and tail weight using octile intervals.
#' Calibrated so that a standard Normal distribution shape maps perfectly to 0.0.
calc_moors_kurtosis <- function(x) {
  x <- as.numeric(na.omit(x))
  if (length(x) < 8L) {
    return(NA_real_)
  }

  #' Extract the 7 precise octile points
  octiles <- quantile(
    x, probs = seq(1 / 8, 7 / 8, by = 1 / 8), names = FALSE
  )

  numerator <- (octiles[7L] - octiles[5L]) + (octiles[3L] - octiles[1L])
  denominator <- octiles[6L] - octiles[2L]

  #' Guard against zero-variance rows or dense identical ties
  if (!is.finite(numerator) || !is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }

  raw_moors <- safe_ratio(numerator, denominator)
  if (is.na(raw_moors) || !is.finite(raw_moors)) {
    return(NA_real_)
  }
  #' Subtract the theoretical normal constant
  #' This sets the baseline expectation to exactly 0.0 under normality
  excess_moors <- raw_moors - 1.233093
  return(excess_moors)
}


#' ============================================================
#' ADDITIONAL FEATURE FUNCTIONS
#' ============================================================
#'
#' These functions add classical, moment-based, robust-spread,
#' Q-Q, and tail features that complement the current ML
#' normality feature set.
#'
#' Important:
#' These functions do not standardize the sample before feature
#' extraction.  When z-scores are used, the standardization is
#' local to the feature definition.
#' ============================================================


#' ------------------------------------------------------------
#' Robust spread ratio: IQR / SD, normalized to equal 1 under
#' normality asymptotically.
#' ------------------------------------------------------------
calc_iqr_sd_ratio <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)

  if (n < 4L) {
    return(NA_real_)
  }

  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  #' For normal data, IQR / SD is approximately 1.34898.
  as.numeric(IQR(x) / s / 1.3489795)
}


#' ------------------------------------------------------------
#' Robust spread ratio: MAD / SD, normalized to equal 1 under
#' normality asymptotically.
#'
#' Note:
#' stats::mad() already multiplies by 1.4826 by default, so it
#' estimates SD under normality.  Therefore mad(x) / sd(x)
#' should be near 1 under normality.
#' ------------------------------------------------------------
calc_mad_sd_ratio <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)

  if (n < 4L) {
    return(NA_real_)
  }

  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  as.numeric(mad(x, constant = 1.4826) / s)
}


#' ------------------------------------------------------------
#' Mean absolute deviation ratio, normalized to equal 1 under
#' normality asymptotically.
#' ------------------------------------------------------------
calc_mean_abs_dev_ratio <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)

  if (n < 4L) {
    return(NA_real_)
  }

  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }

  mad_mean <- mean(abs(x - mean(x)))
  #' For normal data, E|X - mu| / sigma = sqrt(2 / pi).
  as.numeric((mad_mean / s) / sqrt(2 / pi))
}


#' ------------------------------------------------------------
#' Classical skewness and kurtosis features.
#' ------------------------------------------------------------
calc_classical_moments <- function(x) {

  #' Compute sample skewness and kurtosis with safe fallbacks.
  x <- as.numeric(na.omit(x))
  n <- length(x)

  out <- list(
    skewness = NA_real_,
    kurtosis = NA_real_
  )

  if (n < 4L) {
    return(out)
  }

  out$skewness <- safe_calc(
    as.numeric(e1071::skewness(x, type = 2)),
    default = NA_real_
  )

  out$kurtosis <- safe_calc(
    as.numeric(e1071::kurtosis(x, type = 2)),
    default = NA_real_
  )

  out
}


#' ------------------------------------------------------------
#' L-moment skewness and L-moment kurtosis.
#'
#' These are often more stable than ordinary moment skewness and
#' kurtosis under heavy-tailed or contaminated samples.
#' ------------------------------------------------------------
calc_l_moment_features <- function(x) {

  #' Compute L-moment shape summaries from the sample.
  x <- as.numeric(na.omit(x))
  n <- length(x)

  out <- list(
    l_skewness = NA_real_,
    l_kurtosis = NA_real_
  )

  if (n < 5L) {
    return(out)
  }

  lmom <- safe_calc(
    Lmoments::Lmoments(x),
    default = rep(NA_real_, 4L)
  )

  if (length(lmom) >= 4L) {
    out$l_skewness <- as.numeric(lmom[3L])
    out$l_kurtosis <- as.numeric(lmom[4L])
  }

  out
}


#' ------------------------------------------------------------
#' Q-Q correlation with standard normal scores.
#'
#' Values near 1 indicate approximate Q-Q linearity.
#' ------------------------------------------------------------
calc_qq_correlation <- function(x) {

  #' Measure linear agreement between sample and Normal quantiles.
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)

  if (n < 4L || length(unique(x)) < 3L) {
    return(NA_real_)
  }
  
  z_theory <- qnorm(ppoints(n))
  out <- safe_calc(
    cor(x, z_theory),
    default = NA_real_
  )
  as.numeric(out)
}


#' ------------------------------------------------------------
#' Standard deviation of Q-Q residuals from a linear Q-Q fit
#' ------------------------------------------------------------
#'
#' Fits a linear regression of ordered sample values on theoretical
#' Normal quantiles and returns the standard deviation of the
#' residuals. Larger values indicate greater departure from a
#' linear Normal Q-Q pattern.
#' ------------------------------------------------------------

calc_qq_resid_sd <- function(x) {
  #' Measure residual variation around the fitted Normal QQ line.
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)
  
  if (n < 5L || length(unique(x)) < 3L) {
    return(NA_real_)
  }

  z_theory <- qnorm(ppoints(n))
  fit <- safe_calc(
    lm(x ~ z_theory),
    default = NULL
  )

  if (is.null(fit)) {
    return(NA_real_)
  }

  out <- safe_calc(
    sd(residuals(fit)),
    default = NA_real_
  )

  if (!is.finite(out)) {
    return(NA_real_)
  }
  as.numeric(out)
}


#' ------------------------------------------------------------
#' Mean log spacing between consecutive order statistics.
#'
#' Useful as a simple companion to the variance and skewness of
#' log spacings already in the current feature set.
#' ------------------------------------------------------------
calc_mean_log_spacing <- function(x) {
  #' Summarize ordered-value spacings on the log scale.
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)

  if (n < 4L) {
    return(NA_real_)
  }

  d <- diff(x)
  d <- d[d > 0]
  if (length(d) < 2L) {
    return(NA_real_)
  }

  out <- mean(log(d))
  if (!is.finite(out)) {
    return(NA_real_)
  }
  as.numeric(out)
}


#' ------------------------------------------------------------
#' Pearson's second skewness coefficient:
#' 3 * (mean - median) / SD.
#' ------------------------------------------------------------
calc_pearson_skewness <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)

  if (n < 4L) {
    return(NA_real_)
  }

  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  as.numeric(3 * (mean(x) - median(x)) / s)
}


#' ------------------------------------------------------------
#' Normalized 5%-95% tail-weight ratio.
#'
#' The empirical 5%-95% spread of the locally standardized sample
#' is divided by the theoretical standard normal 5%-95% spread.
#' Values above 1 suggest heavier tails than normal.
#' ------------------------------------------------------------
calc_tail_weight_ratio <- function(x) {

  #' Compare outer-tail spread with central spread.
  x <- as.numeric(na.omit(x))
  n <- length(x)

  if (n < 8L) {
    return(NA_real_)
  }

  z <- standardize_sample(x)
  z <- as.numeric(na.omit(z))
  if (length(z) < 8L) {
    return(NA_real_)
  }

  observed <- safe_calc(
    as.numeric(quantile(z, 0.95, names = FALSE) - quantile(z, 0.05, names = FALSE)),
    default = NA_real_
  )

  theoretical <- qnorm(0.95) - qnorm(0.05)
  if (!is.finite(observed) || theoretical <= 0) {
    return(NA_real_)
  }
  as.numeric(observed / theoretical)
}


#' ------------------------------------------------------------
#' Classical normality test statistics.
#'
#' These are used as numerical features, not as p-values.
#' ------------------------------------------------------------
calc_classical_normality_stats <- function(x) {

  #' Compute the requested classical normality-test statistics.
  x <- as.numeric(na.omit(x))
  n <- length(x)

  #' Store the result.
  out <- list(
    sw_stat = NA_real_,
    ad_stat = NA_real_,
    jb_stat = NA_real_,
    ks_stat = NA_real_
  )

  if (n < 4L || length(unique(x)) < 3L) {
    return(out)
  }

  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(out)
  }

  out$sw_stat <- safe_calc(
    as.numeric(shapiro.test(x)$statistic),
    default = NA_real_
  )

  out$ad_stat <- safe_calc(
    as.numeric(nortest::ad.test(x)$statistic),
    default = NA_real_
  )

  out$jb_stat <- safe_calc(
    as.numeric(tseries::jarque.bera.test(x)$statistic),
    default = NA_real_
  )

  out$ks_stat <- safe_calc(
    as.numeric(ks.test(x, "pnorm", mean(x), sd(x))$statistic),
    default = NA_real_
  )
  out
}


#' ------------------------------------------------------------
#' Henze-Zirkler-type univariate ECF statistic.
#'
#' This is a weighted empirical characteristic-function distance
#' from fitted normality.  It complements the existing Epps-Pulley
#' and fixed-grid ECF deviation features.
#' ------------------------------------------------------------
henze_zirkler_stat <- function(x) {
  #' Compute the univariate Henze–Zirkler normality statistic.
  x <- as.numeric(na.omit(x))
  n <- length(x)

  if (n < 5L) {
    return(NA_real_)
  }

  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }

  z <- (x - mean(x)) / s
  beta <- (1 / sqrt(2)) * ((2 * n + 1) / 4)^(1 / 7)
  b    <- beta^2
  dmat <- outer(z, z, "-")^2
  term1 <- (1 / n) * sum(exp(-0.5 * b * dmat))
  term2 <- 2 * (1 + b)^(-0.5) * sum(exp(-b * z^2 / (2 * (1 + b))))
  term3 <- n * (1 + 2 * b)^(-0.5)

  hz <- term1 - term2 + term3
  if (!is.finite(hz)) {
    return(NA_real_)
  }
  as.numeric(hz)
}


#' ------------------------------------------------------------
#' Empirical moment-ratio distance from normal skewness/kurtosis.
#'
#' This combines skewness and kurtosis deviations into one
#' size-adjusted statistic.
#' ------------------------------------------------------------
calc_empirical_moment_ratio_distance <- function(x) {

  #' Compare empirical moment ratios with their Normal-reference values.
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  if (n < 10L) {
    return(NA_real_)
  }

  xc <- x - mean(x)
  m2 <- mean(xc^2)
  if (!is.finite(m2) || m2 <= 0) {
    return(NA_real_)
  }

  sqrt_b1 <- mean(xc^3) / m2^(3 / 2)
  b2      <- mean(xc^4) / m2^2
  var_sqrt_b1 <- 6 / n
  var_b2      <- 24 / n

  if (var_sqrt_b1 <= 0 || var_b2 <= 0) {
    return(NA_real_)
  }

  dist_val <- sqrt((sqrt_b1^2 / var_sqrt_b1) + ((b2 - 3)^2 / var_b2))
  if (!is.finite(dist_val)) {
    return(NA_real_)
  }
  as.numeric(dist_val)
}


#' ------------------------------------------------------------
#' Pickands tail-index estimator.
#'
#' This is useful for detecting heavy-tail behavior.  It can be
#' unstable at very small n, so it is guarded conservatively.
#' ------------------------------------------------------------
pickands_tail_index <- function(x, k = NULL) {

  #' Estimate tail heaviness from upper order statistics.
  x <- sort(as.numeric(na.omit(x)), decreasing = TRUE)
  n <- length(x)

  if (n < 12L) {
    return(NA_real_)
  }

  if (is.null(k)) {
    k <- max(1L, floor(n / 5L))
  }

  k <- min(k, floor(n / 4L))
  if (k < 1L || 4L * k > n) {
    return(NA_real_)
  }

  num <- x[k] - x[2L * k]
  den <- x[2L * k] - x[4L * k]
  if (!is.finite(num) || !is.finite(den) || den <= 0 || num <= 0) {
    return(NA_real_)
  }

  out <- log(num / den) / log(2)
  if (!is.finite(out)) {
    return(NA_real_)
  }
  as.numeric(out)
}


#' ------------------------------------------------------------
#' Mean excess slope.
#'
#' Positive slopes suggest heavier upper-tail behavior; negative
#' slopes suggest lighter upper-tail behavior.
#' ------------------------------------------------------------
mean_excess_slope <- function(x) {
  #' Estimate upper-tail growth from the empirical mean-excess function.
  x <- as.numeric(na.omit(x))
  n <- length(x)

  if (n < 12L) {
    return(NA_real_)
  }

  p_min <- 0.60
  p_max <- 1 - 2 / n

  if (p_max <= p_min) {
    return(NA_real_)
  }

  n_grid <- max(3L, min(8L, floor(n / 3L)))
  probs  <- seq(p_min, p_max, length.out = n_grid)
  thresholds <- quantile(x, probs = probs, names = FALSE)
  mean_excess <- vapply(
    thresholds,
    function(u) {
      exc <- x[x > u] - u

      if (length(exc) < 2L) {
        return(NA_real_)
      }
      mean(exc)
    },
    numeric(1L)
  )

  ok <- is.finite(thresholds) & is.finite(mean_excess)
  if (sum(ok) < 3L) {
    return(NA_real_)
  }

  fit <- safe_calc(
    lm(mean_excess[ok] ~ thresholds[ok]),
    default = NULL
  )

  if (is.null(fit)) {
    return(NA_real_)
  }

  out <- safe_calc(
    as.numeric(coef(fit)[2L]),
    default = NA_real_
  )

  if (!is.finite(out)) {
    return(NA_real_)
  }
  out
}

#' ------------------------------------------------
#' Gini coefficient of absolute deviations from center.
#' ----------------------------------------------------
calc_gini_abs_dev <- function(x, center = c("median", "mean")) {
  center <- match.arg(center)
  x <- x[is.finite(x)]
  if (length(x) < 3L) {
    return(NA_real_)
  }
  c0 <- switch(
    center,
    median = median(x),
    mean   = mean(x)
  )
  z <- abs(x - c0)
  if (all(z == 0)) {
    return(0)
  }
  safe_calc(as.numeric(ineq::ineq(z, type = "Gini")))
}

#' --------------------------------------------------
#' RMS-to-SD ratio.
#' --------------------------------------------------
calc_rms_sd_ratio <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) {
    return(NA_real_)
  }
  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  rms <- sqrt(mean(x^2))
  safe_ratio(rms, s)
}

#' -----------------------------------------------
#' Quantile-based upper-to-lower spread ratio.
#' ----------------------------------------------
calc_quantile_peak_trough_ratio <- function(x) {
  #' Compare central concentration with outer quantile spread.
  x <- x[is.finite(x)]
  if (length(x) < 5L) {
    return(NA_real_)
  }
  #' Calculate q.
  q <- quantile(
    x,
    probs = c(0.05, 0.50, 0.95),
    names = FALSE,
    type = 8
  )
  lower_spread <- q[2L] - q[1L]
  upper_spread <- q[3L] - q[2L]
  safe_ratio(upper_spread, lower_spread)
}

#' ------------------------------------------------------------
#'  tail_index_ratio()
#' ------------------------------------------------------------
#'
#' Purpose
#'   Estimate the ratio of the upper to lower Hill tail indices.
#'   The Hill estimator for the upper tail of a sample x is
#'     alpha_hat = 1 / mean(log(x[1:k]) - log(x[k+1]))
#'   where x is sorted in decreasing order and k = floor(sqrt(n)).
#'   The lower-tail index is estimated on -x.  A ratio near 1
#'   suggests symmetric tails; heavy-tailed distributions tend
#'   to produce ratios far from 1.
#' ------------------------------------------------------------
cal_tail_index_ratio <- function(samples, k = NULL) {
  x <- na.omit(samples)
  n <- length(x)
  if (n < 6L) return(NA_real_)
  if (is.null(k)) k <- max(2L, floor(sqrt(n)))

#' ------------------------------------------------------------
#' Compute one Hill tail-index estimate for a selected tail length.
#' ------------------------------------------------------------
  hill_one <- function(vals) {
    vals <- sort(vals, decreasing = TRUE)
    nv   <- length(vals)
    kv   <- min(k, nv - 1L)
    if (vals[nv] <= 0) vals <- vals - vals[nv] + 1
    anchor <- vals[kv + 1L]
    if (!is.finite(anchor) || anchor <= 0) return(NA_real_)
    h <- mean(log(vals[seq_len(kv)]) - log(anchor))
    if (!is.finite(h) || h <= 0) return(NA_real_)
    1 / h
  }
  au <- safe_calc(hill_one(x))
  al <- safe_calc(hill_one(-x))
  if (!is.finite(au) || !is.finite(al))  return(NA_real_)
  if (abs(al) < .Machine$double.eps)     return(NA_real_)

  au / al
}

#' -------------------------------------------------
#' Calculate Biweight_qq_corr
#' -------------------------------------------------
calc_biweight_qq_corr <- function(x) {
  #' Compute a robust correlation between sample and Normal quantiles.
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)

  if (n < 8L || length(unique(x)) < 4L) {
    return(NA_real_)
  }
  expected <- qnorm(ppoints(n))
  mad_x <- mad(x)
  mad_e <- mad(expected)
  if (!is.finite(mad_x) || !is.finite(mad_e) ||
      mad_x <= 0 || mad_e <= 0) {
    return(NA_real_)
  }

  c_val <- 9
  u <- (x - median(x)) / (c_val * mad_x)
  v <- (expected - median(expected)) / (c_val * mad_e)
  wx <- ifelse(abs(u) < 1, (1 - u^2)^2, 0)
  we <- ifelse(abs(v) < 1, (1 - v^2)^2, 0)
  w  <- wx * we

  if (sum(w) <= 0) {
    return(NA_real_)
  }

  x_c <- x - weighted.mean(x, w)
  e_c <- expected - weighted.mean(expected, w)
  num <- sum(w * x_c * e_c)
  den <- sqrt(sum(w * x_c^2) * sum(w * e_c^2))
  out <- safe_ratio(num, den)

  if (!is.finite(out)) {
    return(NA_real_)
  }
  as.numeric(out)
}

#' ------------------------------------------------------------
#' Geary statistic
#' ------------------------------------------------------------
#'
#' Geary's ratio compares the mean absolute deviation from the mean
#' to the sample standard deviation.
#'
#' For a Normal distribution:
#' E|X - mu| / sigma = sqrt(2 / pi).
#'
#' This feature is normalized so that the Normal reference value is 1.
#' Values below or above 1 indicate departures in tail/body structure.
#' ------------------------------------------------------------

calc_geary_stat <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  if (n < 4L) {
    return(NA_real_)
  }
  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  mean_abs_dev <- mean(abs(x - mean(x)))
  if (!is.finite(mean_abs_dev)) {
    return(NA_real_)
  }

  #' Normalized to equal 1 under normality.
  geary <- safe_ratio(mean_abs_dev / s, sqrt(2 / pi))
  if (!is.finite(geary)) {
    return(NA_real_)
  }
  as.numeric(geary)
}


#' ------------------------------------------------------------
#' Bonett-Seier-type robust kurtosis feature
#' ------------------------------------------------------------
#'
#' This is a robust kurtosis-style ratio based on the second central
#' moment divided by the squared mean absolute deviation from the mean.
#'
#' For a Normal distribution:
#' E[(X - mu)^2] / {E|X - mu|}^2 = pi / 2.
#'
#' The returned value is centered so that the Normal reference is 0:
#' bonett_seier = raw_ratio / (pi / 2) - 1.
#'
#' Positive values suggest heavier tails relative to the Normal reference.
#' Negative values suggest lighter or more bounded tails.
#' ------------------------------------------------------------

calc_bonett_seier <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)
  
  if (n < 4L) {
    return(NA_real_)
  }

  xc <- x - mean(x)
  m2 <- mean(xc^2)
  d1 <- mean(abs(xc))
  if (!is.finite(m2) || !is.finite(d1) || m2 <= 0 || d1 <= 0) {
    return(NA_real_)
  }
  raw_ratio <- safe_ratio(m2, d1^2)
  if (!is.finite(raw_ratio)) {
    return(NA_real_)
  }

  #' Center at the Normal reference value.
  out <- raw_ratio / (pi / 2) - 1
  if (!is.finite(out)) {
    return(NA_real_)
  }
  as.numeric(out)
}


#' ------------------------------------------------------------
#' Bowman-Shenton omnibus moment statistic
#' ------------------------------------------------------------
#'
#' This combines skewness and kurtosis departures from normality.
#' It is closely related to the Jarque-Bera form:
#'
#' n / 6 * [skewness^2 + (kurtosis - 3)^2 / 4]
#'
#' where kurtosis is Pearson kurtosis, not excess kurtosis.
#'
#' Larger values indicate stronger departure from Normal skewness
#' and/or Normal kurtosis.
#' ------------------------------------------------------------

calc_bowman_shenton <- function(x) {
  #' Combine sample skewness and kurtosis into the Bowman–Shenton statistic.
  x <- as.numeric(na.omit(x))
  n <- length(x)

  if (n < 8L) {
    return(NA_real_)
  }

  xc <- x - mean(x)
  m2 <- mean(xc^2)
  if (!is.finite(m2) || m2 <= 0) {
    return(NA_real_)
  }

  m3 <- mean(xc^3)
  m4 <- mean(xc^4)
  if (!is.finite(m3) || !is.finite(m4)) {
    return(NA_real_)
  }
  skew <- m3 / (m2^(3 / 2))
  kurt <- m4 / (m2^2)
  if (!is.finite(skew) || !is.finite(kurt)) {
    return(NA_real_)
  }

  bs <- (n / 6) * (skew^2 + ((kurt - 3)^2 / 4))
  if (!is.finite(bs)) {
    return(NA_real_)
  }
  as.numeric(bs)
}


#' ------------------------------------------------------------
#' Vasicek-type entropy estimator
#' ------------------------------------------------------------
#'
#' Estimates differential entropy using order-statistic spacings.
#' This helper is used by entropy_sd_gap().
#' ------------------------------------------------------------

calc_vasicek_entropy <- function(x, m = NULL) {
  x <- sort(as.numeric(na.omit(x)))
  n <- length(x)

  if (n < 8L) {
    return(NA_real_)
  }

  if (is.null(m)) {
    m <- floor(sqrt(n))
  }
  m <- max(1L, min(as.integer(m), floor((n - 1L) / 2L)))
  if (m < 1L) {
    return(NA_real_)
  }

  i <- seq_len(n)
  lower_idx <- pmax(1L, i - m)
  upper_idx <- pmin(n, i + m)
  spacing <- x[upper_idx] - x[lower_idx]
  if (any(!is.finite(spacing)) || any(spacing <= 0)) {
    spacing[!is.finite(spacing) | spacing <= 0] <- NA_real_
  }

  ok <- is.finite(spacing) & spacing > 0
  if (sum(ok) < max(3L, floor(n / 2L))) {
    return(NA_real_)
  }
  #' Boundary-corrected local window width.
  window_width <- upper_idx - lower_idx
  entropy_terms <- log((n / window_width[ok]) * spacing[ok])
  out <- mean(entropy_terms, na.rm = TRUE)
  if (!is.finite(out)) {
    return(NA_real_)
  }
  as.numeric(out)
}


#' ------------------------------------------------------------
#' Entropy-SD gap
#' ------------------------------------------------------------
#'
#' The Normal distribution has maximum differential entropy among
#' distributions with the same standard deviation.
#'
#' Normal entropy with sample SD s:
#' H_N = log(s * sqrt(2*pi*exp(1)))
#'
#' Feature:
#' entropy_sd_gap = H_N - H_hat
#'
#' Values near 0 are more Normal-like. Larger positive values indicate
#' lower estimated entropy than a Normal distribution with the same SD.
#' ------------------------------------------------------------

calc_entropy_sd_gap <- function(x) {
  #' Compare entropy-based and standard-deviation-based Normal scales.
  x <- as.numeric(na.omit(x))
  n <- length(x)

  if (n < 8L) {
    return(NA_real_)
  }
  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  h_hat <- calc_vasicek_entropy(x)
  if (!is.finite(h_hat)) {
    return(NA_real_)
  }
  h_normal <- log(s * sqrt(2 * pi * exp(1)))
  out <- h_normal - h_hat
  if (!is.finite(out)) {
    return(NA_real_)
  }
  as.numeric(out)
}


#' ------------------------------------------------------------
#' Energy-distance normality statistic
#' ------------------------------------------------------------
#'
#' This computes a one-sample energy goodness-of-fit statistic against
#' fitted standard normality.
#'
#' Steps:
#' 1. Standardize the sample using its sample mean and SD.
#' 2. Compare the standardized empirical distribution to N(0,1).
#'
#' For Z ~ N(0,1), the closed-form expectation is:
#' E|z - Z| = 2*phi(z) + z*(2*Phi(z) - 1).
#'
#' Also:
#' E|Z - Z'| = 2 / sqrt(pi).
#'
#' The statistic is:
#' n * [ 2/n sum_i E|z_i - Z|
#' - E|Z - Z'|
#' - 1/n^2 sum_{i,j} |z_i - z_j| ]
#'
#' Larger values indicate stronger departure from fitted Normality.
#' ------------------------------------------------------------

calc_energy_normal_stat <- function(x) {
  x <- as.numeric(na.omit(x))
  n <- length(x)

  if (n < 5L) {
    return(NA_real_)
  }

  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  z <- (x - mean(x)) / s
  if (any(!is.finite(z))) {
    return(NA_real_)
  }

  #' E|z_i - Z| where Z ~ N(0,1)
  ez <- 2 * dnorm(z) + z * (2 * pnorm(z) - 1)
  if (any(!is.finite(ez))) {
    return(NA_real_)
  }

  term_emp_norm <- (2 / n) * sum(ez)
  #' E|Z - Z'| for independent standard Normals.
  term_norm_norm <- 2 / sqrt(pi)
  #' Empirical pairwise distance term.
  pair_dist <- abs(outer(z, z, "-"))
  term_emp_emp <- mean(pair_dist)

  if (!is.finite(term_emp_emp)) {
    return(NA_real_)
  }

  stat <- n * (term_emp_norm - term_norm_norm - term_emp_emp)
  #' Small numerical negatives can occur from floating-point error.
  if (is.finite(stat) && stat < 0 && abs(stat) < 1e-10) {
    stat <- 0
  }

  if (!is.finite(stat)) {
    return(NA_real_)
  }
  as.numeric(stat)
}

#' ===================================================================
#' Master Feature Generation Engine
#'
#' Computes an expansive, highly rigorous pool of location, scale, spacing,
#' curvature, and robust distribution-shape features for a single sample row.
#' Optimized for high-velocity gating and designed for rigorous feature selection.
#'
#' @param x Numeric vector representing a single distribution sample row.
#' @param feature_set Character vector of allowed feature names. If NULL, all are generated.
#' @param global_range Numeric vector c(min, max) passed to histogram-based entropy functions.
calculate_features <- function(x, feature_set = NULL, global_range = NULL) {
  #' Clean the sample.
  x <- as.numeric(na.omit(x))
  n <- length(x)
  #' Minimum sample size.
  if (n < 8L) {
    stop("calculate_features() requires at least 8 observations for stable shape extraction.")
  }

  #' Available feature names.
  all_features <- c(
    "Mean", "Median", "Mean_Median_Diff", "Studentized_Range",
    "negentropy", "epps_pulley",
    "ecf_re_t1", "ecf_re_t2", "ecf_re_t3",
    "ecf_im_t1", "ecf_im_t2", "ecf_im_t3",
    "dewet_venter", "ryan_joiner",
    "zp_stat", "rao_spacing", "vasicek", "greenwood",
    "var_log_spacings", "skew_log_spacings", "qq_cubic_coef",
    "tail_asymmetry", "medcouple", "moors_kurtosis",
    "excess_tail_prop", "outlier_proportion", "qn_robust_spread",
    "ltw_brys", "rtw_brys", "biweight_qq_corr",
    "renyi_a0.5", "renyi_a1.0", "renyi_a2.0",
    "iqr_sd_ratio", "mad_sd_ratio", "mean_abs_dev_ratio",
    "skewness", "kurtosis", "l_skewness", "l_kurtosis",
    "geary_stat", "bonett_seier", "bowman_shenton",
    "entropy_sd_gap", "energy_normal_stat",
    "qq_corr", "qq_resid_sd", "mean_log_spacing",
    "sw_stat", "ad_stat", "jb_stat", "ks_stat",
    "henze_zirkler", "tail_weight_ratio", "moment_ratio_distance",
    "pickands_tail_index", "mean_excess_slope", "gini_absdev_med",
    "rms_sd_ratio", "q95_q05_asym_ratio",
    "inter_asym_vs_moors", "inter_negent_vs_epps"
  )

  #' Requested subset.
  requested_features <- if (is.null(feature_set)) all_features else unique(as.character(feature_set))

  #' Unsupported names.
  missing_features <- setdiff(requested_features, all_features)
  if (length(missing_features) > 0L) {
    stop("The following requested features are not supported: ", paste(missing_features, collapse = ", "))
  }

  #' Shared-value cache.
  cache <- new.env(parent = emptyenv())

  #' Retrieve or calculate once.
  get_cached <- function(name, calculation) {
    if (!exists(name, envir = cache, inherits = FALSE)) {
      assign(name, calculation(), envir = cache)
    }
    get(name, envir = cache, inherits = FALSE)
  }

  #' Force one numeric value.
  as_scalar <- function(value) {
    if (is.null(value) || length(value) == 0L) return(NA_real_)
    value <- suppressWarnings(as.numeric(value[1L]))
    if (length(value) == 0L) NA_real_ else value
  }

  #' Calculate one requested feature.
  calculate_one <- function(feature_name) {
    #' Reuse completed feature.
    if (exists(feature_name, envir = cache, inherits = FALSE)) {
      return(get(feature_name, envir = cache, inherits = FALSE))
    }

    #' Basic sample summaries.
    mean_x <- get_cached(".mean", function() mean(x))
    median_x <- get_cached(".median", function() median(x))
    sd_x <- get_cached(".sd", function() {
      value <- sd(x)
      if (is.finite(value) && value > 0) value else NA_real_
    })

    #' Select the requested calculation.
    value <- switch(
      feature_name,

      #' Location and scale.
      Mean = mean_x,
      Median = median_x,
      Mean_Median_Diff = mean_x - median_x,
      Studentized_Range = if (is.finite(sd_x)) diff(range(x)) / sd_x else NA_real_,

      #' Entropy and characteristic-function features.
      negentropy = negentropy_approx(x),
      epps_pulley = epps_pulley_stat(x),
      ecf_re_t1 = get_cached(".ecf", function() ecf_deviations(x, t_vals = c(1, 2, 3)))["ecf_re_t1"],
      ecf_re_t2 = get_cached(".ecf", function() ecf_deviations(x, t_vals = c(1, 2, 3)))["ecf_re_t2"],
      ecf_re_t3 = get_cached(".ecf", function() ecf_deviations(x, t_vals = c(1, 2, 3)))["ecf_re_t3"],
      ecf_im_t1 = get_cached(".ecf", function() ecf_deviations(x, t_vals = c(1, 2, 3)))["ecf_im_t1"],
      ecf_im_t2 = get_cached(".ecf", function() ecf_deviations(x, t_vals = c(1, 2, 3)))["ecf_im_t2"],
      ecf_im_t3 = get_cached(".ecf", function() ecf_deviations(x, t_vals = c(1, 2, 3)))["ecf_im_t3"],

      #' Quantile-correlation statistics.
      dewet_venter = dewet_venter_stat(x),
      ryan_joiner = ryan_joiner_stat(x),

      #' Spacing-based statistics.
      zp_stat = calculate_zp_statistic(x),
      rao_spacing = rao_spacing_stat(x),
      vasicek = calculate_vasicek_kmn(x),
      greenwood = greenwood_stat(x),
      var_log_spacings = get_cached(".spacing", function() calc_spacing_ratio_stats(x))$var_log_spacings,
      skew_log_spacings = get_cached(".spacing", function() calc_spacing_ratio_stats(x))$skew_log_spacings,

      #' Shape and tail features.
      qq_cubic_coef = qq_cubic_coef(x),
      tail_asymmetry = calculate_tail_asymmetry(x),
      medcouple = calc_medcouple(x),
      moors_kurtosis = calc_moors_kurtosis(x),
      excess_tail_prop = calculate_excess_tail_prop(x),
      outlier_proportion = calculate_outlier_proportion(x),
      qn_robust_spread = calc_qn_robust_spread(x),
      ltw_brys = get_cached(".tail_weights", function() calc_tail_weights_brys(x))$LTW,
      rtw_brys = get_cached(".tail_weights", function() calc_tail_weights_brys(x))$RTW,
      biweight_qq_corr = calc_biweight_qq_corr(x),

      #' Rényi entropy.
      renyi_a0.5 = renyi_entropy(x, alpha = 0.5, global_range = global_range),
      renyi_a1.0 = renyi_entropy(x, alpha = 1.0, global_range = global_range),
      renyi_a2.0 = renyi_entropy(x, alpha = 2.0, global_range = global_range),

      #' Robust spread ratios.
      iqr_sd_ratio = calc_iqr_sd_ratio(x),
      mad_sd_ratio = calc_mad_sd_ratio(x),

      #' Same Geary-type calculation, retained under two feature names.
      mean_abs_dev_ratio = get_cached(".geary_shared", function() calc_mean_abs_dev_ratio(x)),
      geary_stat = get_cached(".geary_shared", function() calc_mean_abs_dev_ratio(x)),

      #' Classical and L-moment shape summaries.
      skewness = get_cached(".classical_moments", function() calc_classical_moments(x))$skewness,
      kurtosis = get_cached(".classical_moments", function() calc_classical_moments(x))$kurtosis,
      l_skewness = get_cached(".l_moments", function() calc_l_moment_features(x))$l_skewness,
      l_kurtosis = get_cached(".l_moments", function() calc_l_moment_features(x))$l_kurtosis,

      #' Additional omnibus and QQ features.
      bonett_seier = calc_bonett_seier(x),
      bowman_shenton = calc_bowman_shenton(x),
      entropy_sd_gap = calc_entropy_sd_gap(x),
      energy_normal_stat = calc_energy_normal_stat(x),
      qq_corr = calc_qq_correlation(x),
      qq_resid_sd = calc_qq_resid_sd(x),
      mean_log_spacing = calc_mean_log_spacing(x),

      #' Classical normality-test statistics.
      sw_stat = safe_calc(as.numeric(shapiro.test(x)$statistic)),
      ad_stat = safe_calc(as.numeric(nortest::ad.test(x)$statistic)),
      jb_stat = safe_calc(as.numeric(tseries::jarque.bera.test(x)$statistic)),
      ks_stat = if (is.finite(sd_x)) {
        safe_calc(as.numeric(ks.test(x, "pnorm", mean_x, sd_x)$statistic))
      } else {
        NA_real_
      },

      #' Heavy-tail and asymmetry measures.
      henze_zirkler = henze_zirkler_stat(x),
      tail_weight_ratio = calc_tail_weight_ratio(x),
      moment_ratio_distance = calc_empirical_moment_ratio_distance(x),
      pickands_tail_index = pickands_tail_index(x),
      mean_excess_slope = mean_excess_slope(x),
      gini_absdev_med = calc_gini_abs_dev(x, center = "median"),
      rms_sd_ratio = calc_rms_sd_ratio(x),
      q95_q05_asym_ratio = calc_quantile_peak_trough_ratio(x),

      #' Interaction features.
      inter_asym_vs_moors = safe_calc(calculate_one("tail_asymmetry") * calculate_one("moors_kurtosis")),
      inter_negent_vs_epps = safe_calc(calculate_one("negentropy") * calculate_one("epps_pulley")),

      stop("Unsupported feature: ", feature_name)
    )

    #' Store the completed feature.
    value <- as_scalar(value)
    assign(feature_name, value, envir = cache)
    value
  }

  #' Calculate requested columns.
  feature_values <- lapply(requested_features, calculate_one)
  names(feature_values) <- requested_features

  #' One-row feature data frame.
  as.data.frame(feature_values, check.names = FALSE, stringsAsFactors = FALSE)
}

#' ============================================================
#' SECTION 14: FEATURE PREPROCESSING
#' ============================================================
#' Streamlined Preprocessing Engine
#'
#' Applies missing-value imputation followed by centering and scaling to
#' the numeric feature columns, while strictly protecting target labels from leakage.
preprocess_data <- function(train_data,
                            impute_method = c("median", "mean"),
                            scale_method  = c("center", "scale")) {

  impute_method <- match.arg(impute_method)

  if (!"Label" %in% names(train_data)) {
    stop("preprocess_data() requires a Label column.")
  }

  #' Keep the outcome separate to prevent leakage.
  y_train <- train_data$Label

  feature_df <- train_data[, setdiff(names(train_data), "Label"), drop = FALSE]
  numeric_cols <- vapply(feature_df, is.numeric, logical(1L))
  numeric_train <- feature_df[, numeric_cols, drop = FALSE]

  if (ncol(numeric_train) == 0L) {
    stop("No numeric feature columns found for preprocessing.")
  }

  #' Convert Inf and NaN to NA before imputation.
  numeric_train[] <- lapply(numeric_train, function(z) {
    z <- as.numeric(z)
    z[!is.finite(z)] <- NA_real_
    z
  })

  #' Drop columns that are entirely missing or constant.
  keep_cols <- vapply(numeric_train, function(z) {
    any(is.finite(z)) && length(unique(z[is.finite(z)])) > 1L
  }, logical(1L))

  numeric_train <- numeric_train[, keep_cols, drop = FALSE]

  if (ncol(numeric_train) == 0L) {
    stop("All feature columns were empty or constant after cleaning.")
  }

  #' Simple training-only imputation.
  impute_values <- vapply(numeric_train, function(z) {
    val <- switch(
      impute_method,
      median = median(z, na.rm = TRUE),
      mean   = mean(z, na.rm = TRUE)
    )

    if (!is.finite(val)) 0 else val
  }, numeric(1L))

  for (nm in names(numeric_train)) {
    numeric_train[[nm]][is.na(numeric_train[[nm]])] <- impute_values[[nm]]
  }

  #' Centering/scaling is fit on the training features only.
  preproc_obj <- caret::preProcess(
    numeric_train,
    method = scale_method
  )

  train_transformed <- predict(preproc_obj, numeric_train)
  train_transformed$Label <- y_train

  #' Assemble the result.
  list(
    train         = train_transformed,
    preProc      = preproc_obj,
    impute_values = impute_values,
    kept_features = names(numeric_train),
    feature_names = names(train_transformed)[names(train_transformed) != "Label"]
  )
}


#' Apply saved preprocessing to new feature rows.
apply_preprocess <- function(new_data, prep_obj) {

  #' Apply the saved centering and scaling transformations.
  if (!all(c("preProc", "impute_values", "kept_features") %in% names(prep_obj))) {
    stop("prep_obj must be the object returned by preprocess_data().")
  }

  x <- new_data[, prep_obj$kept_features, drop = FALSE]
  x[] <- lapply(x, function(z) {
    z <- as.numeric(z)
    z[!is.finite(z)] <- NA_real_
    z
  })

  for (nm in names(x)) {
    x[[nm]][is.na(x[[nm]])] <- prep_obj$impute_values[[nm]]
  }
  predict(prep_obj$preProc, x)
}


#' ============================================================
#' SECTION 15: TRAINING DATA GENERATION
#' ============================================================
#' Balanced Training Data Generation Engine
#'
#' Generates simulated feature matrices from paired moment-matched Normal and
#' Alternative distribution specifications to force shape-based classification learning.
#' -----------------------------------------------------------------------------------

make_training_data <- function(
    n,
    paired_specs,
    num_sim,
    feature_set = NULL,
    standardize_sample = FALSE,
    center_by = NULL,
    show_progress = TRUE
) {
  #' Each pair contributes one Normal block and one alternative block.
  total_blocks <- 2L * length(paired_specs)
  total_samples <- total_blocks * num_sim
  rows <- vector("list", total_blocks)
  row_id <- 0L

  progress_bar <- start_progress_bar(
    total = total_samples,
    label = paste0("Generating training samples for n = ", n),
    show_progress = show_progress
  )

  on.exit(close_progress_bar(progress_bar), add = TRUE)

  for (pair in paired_specs) {
    #' Generate the moment-matched Normal block.
    row_id <- row_id + 1L
    rows[[row_id]] <- generate_feature_rows(
      n = n,
      n_rep = num_sim,
      spec = pair$normal,
      label = "Normal",
      feature_set = feature_set,
      standardize_sample = standardize_sample,
      center_by = center_by,
      progress_bar = progress_bar,
      progress_offset = (row_id - 1L) * num_sim
    )

    #' Generate the corresponding nonnormal block.
    row_id <- row_id + 1L
    rows[[row_id]] <- generate_feature_rows(
      n = n,
      n_rep = num_sim,
      spec = pair$alt,
      label = "Non_Normal",
      feature_set = feature_set,
      standardize_sample = standardize_sample,
      center_by = center_by,
      progress_bar = progress_bar,
      progress_offset = (row_id - 1L) * num_sim
    )
  }

  train_data <- do.call(rbind, rows)

  #' Fix the class order required by caret's two-class ROC summary.
  train_data$Label <- factor(
    train_data$Label,
    levels = c("Non_Normal", "Normal")
  )

  rownames(train_data) <- NULL
  train_data
}

#' ============================================================
#' SECTION 16: MODEL TRAINING HELPERS
#' ============================================================
#' Caret Cross-Validation Settings
#'
#' Establishes balanced 10-fold cross-validation
make_cv_control <- function(
    k_folds = 5L,
    repeats = 1L,
    verbose = FALSE
) {
  #' Use ordinary cross-validation unless the user explicitly requests repeats.
  repeats <- max(1L, as.integer(repeats))
  cv_method <- if (repeats > 1L) "repeatedcv" else "cv"

  control_args <- list(
    method = cv_method,
    number = as.integer(k_folds),
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,

    #' Resampling predictions are not retained because final evaluation is separate.
    savePredictions = FALSE,

    #' Let caret use the temporary backend registered by train_models_all_n().
    allowParallel = TRUE,
    verboseIter = isTRUE(verbose)
  )

  if (repeats > 1L) {
    control_args$repeats <- repeats
  }

  do.call(caret::trainControl, control_args)
}

#' Train One Robust Classifier
train_one_model <- function(model_name, train_data, cv_ctrl) {
  #' Determine actual feature dimensions dynamically (subtracting the 'Label' column)
  num_features <- ncol(train_data) - 1L

  #' -------------------------------------------------------------------------
  #' Logistic Regression via Elastic Net (Glmnet)
  #' -------------------------------------------------------------------------
  if (model_name == "LR") {
    return(
      suppressWarnings(
        caret::train(
          Label ~ .,
          data      = train_data,
          method    = "glmnet",
          metric    = "ROC",
          trControl = cv_ctrl,
          tuneGrid  = expand.grid(
            alpha  = c(0.0, 0.5, 1.0),
            lambda = 10^seq(-4, 1, length.out = 20L)
          ),
          family    = "binomial"
        )
      )
    )
  }

  #' -------------------------------------------------------------------------
  #' Random Forest (Dynamic Tuning Bounds)
  #' -------------------------------------------------------------------------
  if (model_name == "RF") {
    #' Use a wider mtry grid so the forest is not locked into a
    #' small hand-picked set.  The grid is clipped to the available
    #' number of features.
    mtry_pool <- unique(round(c(
      1,
      sqrt(num_features),
      num_features / 4,
      num_features / 3,
      num_features / 2,
      num_features
    )))

    valid_mtry <- sort(unique(pmax(1L, pmin(num_features, mtry_pool))))

    #' Return the result.
    return(
      caret::train(
        Label ~ .,
        data       = train_data,
        method     = "rf",
        metric     = "ROC",
        trControl  = cv_ctrl,
        tuneGrid   = expand.grid(mtry = valid_mtry),
        ntree      = 1000L,
        importance = TRUE
      )
    )
  }

  #' -------------------------------------------------------------------------
  #' Artificial Neural Network (Nnet)
  #' -------------------------------------------------------------------------
  if (model_name == "ANN") {
    #' Slightly broader network-size and decay search.  Inputs are
    #' already centered and scaled by preprocess_data().
    return(
      caret::train(
        Label ~ .,
        data      = train_data,
        method    = "nnet",
        metric    = "ROC",
        trControl = cv_ctrl,
        tuneGrid  = expand.grid(
          size  = c(3, 5, 10),
          decay = c(0.001, 0.01, 0.10)
        ),
        MaxNWts = 10000L,
        maxit   = 500L,
        trace   = FALSE
      )
    )
  }

  #' -------------------------------------------------------------------------
  #' Gradient Boosting Machine (Gbm)
  #' -------------------------------------------------------------------------
  if (model_name == "GBM") {
    #' Broader boosting grid.  The extra depths and tree counts help
    #' capture nonlinear feature interactions without changing the
    #' preprocessing or outcome definition.
    min_node_size <- if (nrow(train_data) < 500L) c(2, 5) else c(5, 10)

    #' Return the result.
    return(
      caret::train(
        Label ~ .,
        data      = train_data,
        method    = "gbm",
        metric    = "ROC",
        trControl = cv_ctrl,
        tuneGrid  = expand.grid(
          n.trees           = c(100, 300, 600),
          interaction.depth = c(1, 3),
          shrinkage         = c(0.01, 0.05, 0.10),
          n.minobsinnode    = min_node_size
        ),
        verbose = FALSE
      )
    )
  }

  #' -------------------------------------------------------------------------
  #' Custom e1071 LIBSVM radial container (Platt Probability Calibrated)
  #' -------------------------------------------------------------------------
  if (model_name == "SVM") {
    svm_e1071 <- list(
      label      = "Support Vector Machine (RBF, e1071)",
      type       = "Classification",
      library    = "e1071",
      loop       = NULL,
      parameters = data.frame(
        parameter = c("cost", "gamma"),
        class     = c("numeric", "numeric"),
        label     = c("Cost", "Gamma")
      ),
      grid       = function(x, y, len = NULL, search = "grid") {
        expand.grid(cost = c(0.1, 1, 10), gamma = c(0.01, 0.1, 1))
      },
      fit        = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
        fit <- e1071::svm(
          x           = x,
          y           = y,
          kernel      = "radial",
          cost        = param$cost,
          gamma       = param$gamma,
          probability = TRUE,
          ...
        )

        #' Continue the calculation.
        fit$obsLevels <- lev
        fit
      },
      predict    = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
        predict(modelFit, newdata = newdata)
      },
      prob       = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
        prob_mat <- attr(
          predict(modelFit, newdata = newdata, probability = TRUE),
          "probabilities"
        )

        if (is.null(prob_mat)) {
          stop("SVM probability estimates were not returned.")
        }

        prob_mat[, modelFit$obsLevels, drop = FALSE]
      },
      sort       = function(x) x[order(x$cost, x$gamma), ]
    )

    #' Return the result.
    return(
      caret::train(
        Label ~ .,
        data      = train_data,
        method    = svm_e1071,
        metric    = "ROC",
        trControl = cv_ctrl,
        tuneGrid  = expand.grid(
          cost  = 2^seq(-3, 5, by = 2),
          gamma = 2^seq(-7, -1, by = 2)
        )
      )
    )
  }

  #' -------------------------------------------------------------------------
  #' K-Nearest Neighbors (Knn)
  #' -------------------------------------------------------------------------
  if (model_name == "KNN") {
    #' Wider neighbor grid.  The upper bound is tied to training size
    #' so KNN does not overfit to very small neighborhoods only.
    max_k <- max(5L, min(51L, floor(sqrt(nrow(train_data)))))
    k_sequence <- seq(3L, max_k, by = 2L)

    #' Return the result.
    return(
      caret::train(
        Label ~ .,
        data      = train_data,
        method    = "knn",
        metric    = "ROC",
        trControl = cv_ctrl,
        tuneGrid  = expand.grid(k = k_sequence)
      )
    )
  }
  stop("Unsupported model name: ", model_name)
}


#' ============================================================
#' SECTION 17: TRAIN MODELS FOR ONE SAMPLE SIZE
#' ============================================================

#' ------------------------------------------------------------
#' Train the requested machine-learning models for one sample size.
#' ------------------------------------------------------------
train_models_one_n <- function(
    n,
    paired_specs,
    num_sim = 100,
    models_to_train = c("LR", "RF", "ANN", "GBM", "SVM", "KNN"),
    feature_set = NULL,
    standardize_sample = FALSE,
    center_by = NULL,
    k_folds = 5L,
    cv_repeats = 1L,
    cv_verbose = FALSE,
    show_progress = TRUE
) {
  valid_models <- c("LR", "RF", "ANN", "GBM", "SVM", "KNN")
  bad_models <- setdiff(models_to_train, valid_models)

  if (length(bad_models) > 0L) {
    stop("Unsupported model(s): ", paste(bad_models, collapse = ", "))
  }

  center_by <- resolve_center_by(standardize_sample, center_by)
  start_time <- Sys.time()

  #' Continue the calculation.
  cat("\n", strrep("=", 65), "\n", sep = "")
  cat("Training models for n =", n, "\n")
  cat("Matched distribution pairs:", length(paired_specs), "\n")
  cat("Replicates per class per pair:", num_sim, "\n")
  cat("Selected features:", if (is.null(feature_set)) "all" else length(feature_set), "\n")
  cat("CV design:", k_folds, "folds x", cv_repeats, "repeat(s)\n")
  cat(strrep("=", 65), "\n")

  #' Generate the complete training matrix before fitting any classifier.
  raw_train <- make_training_data(
    n = n,
    paired_specs = paired_specs,
    num_sim = num_sim,
    feature_set = feature_set,
    standardize_sample = standardize_sample,
    center_by = center_by,
    show_progress = show_progress
  )

  #' Fit imputation, centering, and scaling on the training data only.
  prep <- preprocess_data(
    train_data = raw_train,
    impute_method = "median",
    scale_method = c("center", "scale")
  )

  #' Calculate train std.
  train_std <- prep$train
  cv_ctrl <- make_cv_control(
    k_folds = k_folds,
    repeats = cv_repeats,
    verbose = cv_verbose
  )

  models <- vector("list", length(models_to_train))
  names(models) <- models_to_train

  model_bar <- start_progress_bar(
    total = length(models_to_train),
    label = paste0("Fitting models for n = ", n),
    show_progress = show_progress
  )

  on.exit(close_progress_bar(model_bar), add = TRUE)

  #' Fit each requested model and report its elapsed time.
  for (model_index in seq_along(models_to_train)) {
    model_name <- models_to_train[model_index]
    model_start <- Sys.time()

    cat("\nFitting ", model_name, " for n = ", n, " ...\n", sep = "")

    models[[model_name]] <- train_one_model(
      model_name = model_name,
      train_data = train_std,
      cv_ctrl = cv_ctrl
    )

    cat(
      "Completed ", model_name, " in ",
      format_elapsed_time(model_start), ".\n",
      sep = ""
    )

    update_progress_bar(model_bar, model_index)
  }

  cat(
    "Completed all training for n = ", n, " in ",
    format_elapsed_time(start_time), ".\n",
    sep = ""
  )

  #' Assemble the result.
  list(
    n = n,
    models = models,
    prep_obj = prep,
    feature_names = prep$feature_names,
    paired_specs = paired_specs,
    standardize_sample = standardize_sample,
    center_by = center_by
  )
}

#' ============================================================
#' SECTION 18: TRAIN MODELS ACROSS SAMPLE SIZES
#' ============================================================
#'
#' Trains the requested models for each sample size.
#'
#' Returns:
#' A named list, where each element is the trained model bundle
#' for one sample size.
#' ============================================================

#' Calculate train models all n.
train_models_all_n <- function(
    sample_sizes,
    paired_specs,
    num_sim = 100,
    models_to_train = c("LR", "RF", "ANN", "GBM", "SVM", "KNN"),
    feature_set = NULL,
    standardize_sample = FALSE,
    center_by = NULL,
    k_folds = 5L,
    cv_repeats = 1L,
    cv_verbose = FALSE,
    n_cores = 1L,
    n_cores_cv = 1L,
    show_progress = TRUE
) {
  if (length(sample_sizes) == 0L) {
    stop("sample_sizes must contain at least one sample size.")
  }

  sample_sizes <- as.numeric(sample_sizes)
  n_cores <- normalize_n_cores(n_cores, max_tasks = length(sample_sizes))
  n_cores_cv <- normalize_n_cores(n_cores_cv)

  #' Continue the calculation.
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("TRAINING MODELS ACROSS SAMPLE SIZES\n")
  cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
  cat("Sample-size workers:", n_cores, "\n")
  cat("Caret CV workers:", n_cores_cv, "\n")
  cat(strrep("=", 70), "\n")

  #' Define the work for one sample size in a reusable local function.
  train_one_size <- function(n, local_progress) {

  #' Train the requested models for one sample size.
    train_models_one_n(
      n = n,
      paired_specs = paired_specs,
      num_sim = num_sim,
      models_to_train = models_to_train,
      feature_set = feature_set,
      standardize_sample = standardize_sample,
      center_by = center_by,
      k_folds = k_folds,
      cv_repeats = cv_repeats,
      cv_verbose = cv_verbose,
      show_progress = local_progress
    )
  }

  if (n_cores_cv > 1L) {
    #' Process sample sizes sequentially while caret parallelizes CV internally.
    if (n_cores > 1L) {
      message("n_cores_train is ignored because n_cores_cv > 1. ",
        "This prevents nested parallelism."
      )
    }

    trained_list <- with_caret_parallel(
      n_cores = n_cores_cv,
      expr = lapply(
        sample_sizes,
        function(n) train_one_size(n, show_progress)
      )
    )
  } else if (n_cores > 1L) {
    #' Parallel sample-size workers use messages instead of competing progress bars.
    trained_list <- safe_mclapply(
      X = sample_sizes,
      FUN = function(n) train_one_size(n, FALSE),
      n_cores = n_cores
    )
  } else {
    #' Sequential training displays the full progress bars.
    trained_list <- lapply(
      sample_sizes,
      function(n) train_one_size(n, show_progress)
    )
  }
  names(trained_list) <- as.character(sample_sizes)
  trained_list
}

#' ============================================================
#' PART 4: PREDICTION AND EVALUATION
#' ============================================================

#' ============================================================
#' SECTION 19: PREDICTION HELPERS
#' ============================================================

#' Resolve centering choice for optional sample standardization.
#'
#' If standardize_sample = FALSE, no centering is applied to the raw sample.
#' If standardize_sample = TRUE and center_by is NULL, mean-centering is used.
resolve_center_by <- function(standardize_sample = FALSE, center_by = NULL) {
  if (!isTRUE(standardize_sample)) {
    return("mean")  #' ignored by generate_data() when standardize = FALSE
  }
  if (is.null(center_by)) {
    return("mean")
  }
  match.arg(center_by, choices = c("mean", "median"))
}

#' Prepare one sample for prediction.
#'
#' Computes only the features used during training and applies the
#' preprocessing object fitted on the training data.
prepare_prediction_features <- function(x, model_bundle) {

  #' Calculate and clean the feature values required for prediction.
  if (is.null(model_bundle$feature_names)) {
    stop("model_bundle$feature_names is missing.")
  }
  feats <- calculate_features(
    x           = x,
    feature_set = model_bundle$feature_names
  )
  missing_features <- setdiff(model_bundle$feature_names, names(feats))
  #' Validate the condition.
  if (length(missing_features) > 0L) {
    stop("Missing feature(s) from calculate_features(): ",
      paste(missing_features, collapse = ", ")
    )
  }
  feats <- feats[, model_bundle$feature_names, drop = FALSE]
  if (!is.null(model_bundle$prep_obj)) {
    feats <- apply_preprocess(feats, model_bundle$prep_obj)
  }
  feats
}


#' Prepare a batch of feature rows for prediction.
prepare_prediction_matrix <- function(feature_df, model_bundle) {
  #' Align prediction columns with the features used during training.
  if (is.null(model_bundle$feature_names)) {
    stop("model_bundle$feature_names is missing.")
  }
  missing_features <- setdiff(model_bundle$feature_names, names(feature_df))
  #' Validate the condition.
  if (length(missing_features) > 0L) {
    stop("Evaluation feature matrix is missing feature(s): ",
      paste(missing_features, collapse = ", ")
    )
  }

  feature_df <- feature_df[, model_bundle$feature_names, drop = FALSE]
  feature_df[] <- lapply(feature_df, function(z) {
    z <- as.numeric(z)
    z[!is.finite(z)] <- NA_real_
    z
  })

  if (!is.null(model_bundle$prep_obj)) {
    feature_df <- apply_preprocess(feature_df, model_bundle$prep_obj)
  }
  feature_df
}

#' Predict P(Non_Normal) from all trained models.
#'
#' Returns a data frame with one probability column per model.
predict_model_probs <- function(model_bundle, feats_std) {
  #' Dispatch probability prediction to the selected fitted model.
  model_names <- names(model_bundle$models)
  prob_df <- lapply(model_names, function(model_name) {
    pred <- predict( model_bundle$models[[model_name]], newdata = feats_std, type = "prob")

    if (!"Non_Normal" %in% names(pred)) {
      stop("Model ", model_name, " did not return a 'Non_Normal' probability column.")
    }
    as.numeric(pred[, "Non_Normal"])
  })
  prob_df <- as.data.frame(prob_df)
  names(prob_df) <- model_names
  prob_df
}


#' ------------------------------------------------------------
#' Create ML normality-score function for the utility framework
#' Returns P(Normal) for one sample and one trained model.
#' ------------------------------------------------------------
predict_prob_normal_fast <- function(x, model_bundle, model_name = "RF") {
  #' Return the fitted model probability assigned to the Normal class.
  if (!model_name %in% names(model_bundle$models)) {
    stop("Model not found in model_bundle: ", model_name)
  }

  feats_std <- prepare_prediction_features(x = x, model_bundle = model_bundle)

  pred_probs <- predict(
    model_bundle$models[[model_name]],
    newdata = feats_std,
    type    = "prob"
  )

  if (!"Normal" %in% names(pred_probs)) {
    stop("Model ", model_name, " did not return a 'Normal' probability column.")
  }
  as.numeric(pred_probs[, "Normal"])
}

#' Convert P(Non_Normal) to a class label.
class_from_prob <- function(prob, threshold = 0.50) {
  ifelse(prob >= threshold, "Non_Normal", "Normal")
}

#' Create a readable distribution label.
make_distribution_label <- function(spec) {
  #' Continue the calculation.
  paste0( spec$dist, "(", paste(round(spec$par, 4), collapse = ", "),")" )
}


#' ============================================================
#' SECTION 20: EVALUATION DATA GENERATION
#' ============================================================

#' Create flattened evaluation specifications.
#'
#' For each alternative distribution, this returns:
#' 1. the moment-matched Normal distribution;
#' 2. the corresponding alternative distribution.
make_eval_specs <- function(alt_specs) {
  pairs <- make_paired_specs(alt_specs)
  eval_specs <- vector("list", 2L * length(pairs))
  row_id <- 0L

  #' Retain a pair identifier so ROC curves can reuse evaluation predictions.
  for (pair_id in seq_along(pairs)) {
    pair <- pairs[[pair_id]]

    #' Calculate row id.
    row_id <- row_id + 1L
    eval_specs[[row_id]] <- list(
      spec = pair$normal,
      label = "Normal",
      pair_id = pair_id,
      pair_role = "Normal"
    )

    #' Calculate row id.
    row_id <- row_id + 1L
    eval_specs[[row_id]] <- list(
      spec = pair$alt,
      label = "Non_Normal",
      pair_id = pair_id,
      pair_role = "Alternative"
    )
  }
  eval_specs
}


#' Generate evaluation features for one distribution specification.
#'
#' This returns a feature matrix plus the corresponding true labels
#' and distribution labels.
generate_eval_block <- function(n, n_iter, eval_item, feature_names, 
                                standardize_sample = FALSE, center_by = NULL, 
                                progress_bar = NULL, progress_offset = 0L) {
  spec <- eval_item$spec
  true_class <- eval_item$label
  dist_label <- make_distribution_label(spec)
  center_by_use <- resolve_center_by(standardize_sample = standardize_sample,center_by = center_by)

  rows <- vector("list", n_iter)
  progress_every <- max(1L, floor(n_iter / 100L))

  #' Generate one feature row per independent evaluation sample.
  for (iteration in seq_len(n_iter)) {
    x <- generate_data(n = n,dist = spec$dist,par = spec$par,standardize = standardize_sample,center_by = center_by_use)

    rows[[iteration]] <- calculate_features(x = x, feature_set = feature_names)

    #' Update the shared evaluation progress bar at approximately 1% intervals.
    if (!is.null(progress_bar) &&
        (iteration %% progress_every == 0L || iteration == n_iter)) {
      update_progress_bar(progress_bar, progress_offset + iteration)
    }
  }

  feature_df <- do.call(rbind, rows)
  rownames(feature_df) <- NULL

  #' Assemble the result.
  list(
    features = feature_df,
    true_class = rep(true_class, n_iter),
    distribution = rep(dist_label, n_iter),
    pair_id = rep(eval_item$pair_id, n_iter),
    pair_role = rep(eval_item$pair_role, n_iter)
  )
}


#' ============================================================
#' SECTION 21: MODEL EVALUATION FOR ONE SAMPLE SIZE
#' ============================================================

#' Evaluate trained models for one sample size.
#'
#' For each evaluation distribution:
#' 1. generate n_iter samples;
#' 2. compute the trained feature set;
#' 3. apply the training preprocessing object;
#' 4. predict P(Non_Normal);
#' 5. return prediction data frames for each model and MajorityVote.
evaluate_one_n <- function(
    n,
    model_bundle,
    eval_specs,
    n_iter = 1000,
    threshold = 0.50,
    majority_threshold = 0.40,
    standardize_sample = FALSE,
    center_by = NULL,
    store_holdout_features = FALSE,
    show_progress = TRUE
) {
  if (is.null(model_bundle$models) || length(model_bundle$models) == 0L) {
    stop("model_bundle$models is missing or empty.")
  }
  if (is.null(model_bundle$feature_names)) {
    stop("model_bundle$feature_names is missing.")
  }

  model_names <- names(model_bundle$models)
  all_names <- c(model_names, "MajorityVote")
  start_time <- Sys.time()

  cat("\nEvaluating n =", n, "\n")
  cat("Evaluation blocks:", length(eval_specs), "\n")
  cat("Replicates per block:", n_iter, "\n")

  block_results <- vector("list", length(eval_specs))
  total_samples <- length(eval_specs) * n_iter

  progress_bar <- start_progress_bar(
    total = total_samples,
    label = paste0("Generating evaluation samples for n = ", n),
    show_progress = show_progress
  )

  on.exit(close_progress_bar(progress_bar), add = TRUE)
  
  #' Generate features and probabilities one distribution block at a time.
  for (block_index in seq_along(eval_specs)) {
    eval_item <- eval_specs[[block_index]]

    #' Calculate eval block.
    eval_block <- generate_eval_block(
      n = n,
      n_iter = n_iter,
      eval_item = eval_item,
      feature_names = model_bundle$feature_names,
      standardize_sample = standardize_sample,
      center_by = center_by,
      progress_bar = progress_bar,
      progress_offset = (block_index - 1L) * n_iter
    )

    feats_std <- prepare_prediction_matrix(feature_df = eval_block$features, model_bundle = model_bundle)

    probs <- predict_model_probs(model_bundle = model_bundle,feats_std = feats_std)

    #' Continue the calculation.
    block_results[[block_index]] <- list(
      probs = probs,
      true_class = eval_block$true_class,
      distribution = eval_block$distribution,
      pair_id = eval_block$pair_id,
      pair_role = eval_block$pair_role,
      features = if (store_holdout_features) feats_std else NULL
    )
  }

  true_class_all <- unlist(lapply(block_results, `[[`, "true_class"), use.names = FALSE)

  distribution_all <- unlist(lapply(block_results, `[[`, "distribution"),use.names = FALSE)

  pair_id_all <- unlist(lapply(block_results, `[[`, "pair_id"),use.names = FALSE)

  pair_role_all <- unlist(lapply(block_results, `[[`, "pair_role"),use.names = FALSE)

  prob_all <- do.call(rbind, lapply(block_results, `[[`, "probs"))
  
  rownames(prob_all) <- NULL

  pred_store <- vector("list", length(all_names))
  names(pred_store) <- all_names

  #' Store one prediction table for each fitted model.
  for (model_name in model_names) {
    prob_nn <- prob_all[[model_name]]

    #' Build the result table.
    pred_store[[model_name]] <- data.frame(
      True_Class = true_class_all,
      Predicted_Class = class_from_prob(prob_nn, threshold),
      Prob_Non_Normal = prob_nn,
      Distribution = distribution_all,
      Pair_ID = pair_id_all,
      Pair_Role = pair_role_all,
      stringsAsFactors = FALSE
    )
  }

  #' Average model probabilities for the majority-vote classifier.
  mean_prob <- rowMeans(prob_all[, model_names, drop = FALSE], na.rm = TRUE)

  #' Build the result table.
  pred_store[["MajorityVote"]] <- data.frame(
    True_Class = true_class_all,
    Predicted_Class = class_from_prob(mean_prob, majority_threshold),
    Prob_Non_Normal = mean_prob,
    Distribution = distribution_all,
    Pair_ID = pair_id_all,
    Pair_Role = pair_role_all,
    stringsAsFactors = FALSE
  )

  holdout_features <- NULL

  #' Validate the condition.
  if (isTRUE(store_holdout_features)) {
    holdout_features <- do.call( rbind, lapply(block_results, `[[`, "features"))
    rownames(holdout_features) <- NULL
  }

  cat("Completed evaluation for n = ", n, " in ", format_elapsed_time(start_time), ".\n", sep = "")
  list(
    predictions = pred_store,
    holdout_features = holdout_features
  )
}

#' The ML classifier returns an estimated probability that a sample belongs
#' to the non-normal class. For the primary analysis, a fixed threshold of 0.50
#' was used: samples with estimated probability at least 0.50 were classified as
#' non-normal, and samples below 0.50 were classified as normal. This choice
#' avoids tuning the classification rule on the evaluation distributions and
#' gives a transparent decision rule. Alternative thresholds may be examined as
#' sensitivity analyses.

#' ============================================================
#' SECTION 22: PERFORMANCE METRICS
#' ============================================================

#' Compute classification metrics from prediction data.
#'
#' The positive class is "Non_Normal".
compute_metrics <- function(pred_df, positive_class = "Non_Normal", digits = NULL) {

  required_cols <- c("True_Class", "Predicted_Class", "Prob_Non_Normal")
  missing_cols  <- setdiff(required_cols, names(pred_df))

  #' Validate the condition.
  if (length(missing_cols) > 0L) {
    stop("pred_df is missing required column(s): ",paste(missing_cols, collapse = ", ")
    )
  }

  class_levels <- c("Non_Normal", "Normal")

  pred_df$True_Class <- factor( pred_df$True_Class, levels = class_levels)
  pred_df$Predicted_Class <- factor( pred_df$Predicted_Class,levels = class_levels)

  cm <- caret::confusionMatrix(
    data      = pred_df$Predicted_Class,
    reference = pred_df$True_Class,
    positive  = positive_class
  )

  accuracy    <- as.numeric(cm$overall["Accuracy"])
  sensitivity <- as.numeric(cm$byClass["Sensitivity"])
  specificity <- as.numeric(cm$byClass["Specificity"])
  precision   <- as.numeric(cm$byClass["Precision"])
  f1          <- as.numeric(cm$byClass["F1"])

  auc <- NA_real_

  if (length(unique(pred_df$True_Class)) == 2L && all(is.finite(pred_df$Prob_Non_Normal))) {

    #' Calculate auc.
    auc <- safe_calc(
      as.numeric(
        pROC::auc(
          response  = pred_df$True_Class,
          predictor = pred_df$Prob_Non_Normal,
          levels    = c("Normal", "Non_Normal"),
          direction = "<",
          quiet     = TRUE
        )
      ),
      default = NA_real_
    )
  }

  #' Calculate metric values.
  metric_values <- c(
    Accuracy    = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision   = precision,
    F1          = f1,
    AUC         = auc
  )

  if (!is.null(digits)) {
    metric_values <- round(metric_values, digits)
  }

  #' Assemble the result.
  list(
    ConfusionMatrix = cm,
    Metrics         = metric_values,
    Predictions     = pred_df
  )
}

#' -----------------------------------------------------
#' Summarize metrics for all models at one sample size.
#' ----------------------------------------------------
get_metrics_summary <- function(model_metrics, digits = 4) {
  model_names <- names(model_metrics)
  summary_df <- do.call(
    rbind,
    lapply(model_names, function(model_name) {
      metrics <- model_metrics[[model_name]]$Metrics
      #' Build the result table.
      data.frame(
        Model       = model_name,
        Accuracy    = metrics["Accuracy"],
        Sensitivity = metrics["Sensitivity"],
        Specificity = metrics["Specificity"],
        Precision   = metrics["Precision"],
        F1          = metrics["F1"],
        AUC         = metrics["AUC"],
        stringsAsFactors = FALSE
      )
    })
  )

  rownames(summary_df) <- NULL
  numeric_cols <- setdiff(names(summary_df), "Model")
  summary_df[numeric_cols] <- lapply(summary_df[numeric_cols],
    function(z) round(as.numeric(z), digits)
  )
  summary_df
}


#' ============================================================
#' SECTION 23: EVALUATE MODELS ACROSS SAMPLE SIZES
#' ============================================================

#' Evaluate trained models across all sample sizes.
#'
#' Returns one list element per sample size. Each element contains:
#' - metrics: per-model metric objects;
#' - summary: compact metric table;
#' - predictions: raw prediction data frames;
#' - holdout_features: optional processed holdout feature matrix.
evaluate_models_all_n <- function(
    trained_models,
    eval_specs,
    n_iter = 1000,
    threshold = 0.50,
    majority_threshold = 0.50,
    standardize_sample = FALSE,
    center_by = NULL,
    store_holdout_features = FALSE,
    digits = 4,
    n_cores = 1L,
    show_progress = TRUE
) {
  if (length(trained_models) == 0L) {
    stop("trained_models is empty.")
  }

  model_keys <- names(trained_models)
  n_cores <- normalize_n_cores(n_cores, max_tasks = length(model_keys))

  cat("\n", strrep("=", 65), "\n", sep = "")
  cat("EVALUATING MODELS ACROSS SAMPLE SIZES\n")
  cat("Sample sizes:", paste(model_keys, collapse = ", "), "\n")
  cat("Parallel workers:", n_cores, "\n")
  cat(strrep("=", 65), "\n")

#' ------------------------------------------------------------
#' Evaluate trained models for one sample size.
#' ------------------------------------------------------------
  evaluate_one_size <- function(n_key, local_progress) {
    model_bundle <- trained_models[[n_key]]

    if (is.null(model_bundle$n)) {
      stop("trained_models[['", n_key, "']] is missing element 'n'.")
    }

    #' Calculate eval raw.
    eval_raw <- evaluate_one_n(
      n = model_bundle$n,
      model_bundle = model_bundle,
      eval_specs = eval_specs,
      n_iter = n_iter,
      threshold = threshold,
      majority_threshold = majority_threshold,
      standardize_sample = standardize_sample,
      center_by = center_by,
      store_holdout_features = store_holdout_features,
      show_progress = local_progress
    )

    #' Calculate held-out performance metrics.
    model_metrics <- lapply(eval_raw$predictions, compute_metrics)
    summary_df <- get_metrics_summary(model_metrics, digits = digits)

    cat("\nPerformance summary for n =", model_bundle$n, "\n")
    print(summary_df)

    #' Assemble the result.
    list(
      n = model_bundle$n,
      metrics = model_metrics,
      summary = summary_df,
      predictions = eval_raw$predictions,
      holdout_features = eval_raw$holdout_features
    )
  }
  if (n_cores > 1L) {
    message("Evaluation progress bars are disabled inside parallel workers. ",
      "Completion messages will still be printed."
    )
    #' Store the results.
    results <- safe_mclapply(X = model_keys,
      FUN = function(n_key) evaluate_one_size(n_key, FALSE),
      n_cores = n_cores
    )
  } else {
    results <- lapply(model_keys,
      function(n_key) evaluate_one_size(n_key, show_progress)
    )
  }
  names(results) <- model_keys
  results
}

#' ============================================================
#' PART 5: ROC PLOTTING AND CLASSICAL-VS-ML COMPARISON
#' ============================================================

#' ============================================================
#' SECTION 24: ROC HELPER FUNCTIONS
#' ============================================================

#' Compute AUC using the trapezoidal rule.
#'
#' This is mainly used for manually constructed ROC curves, such as
#' classical-test ROC curves computed over a grid of alpha cutoffs.
compute_auc <- function(fpr,
                        tpr,
                        ensure_endpoints = TRUE) {

  ok <- is.finite(fpr) & is.finite(tpr)
  fpr <- fpr[ok]
  tpr <- tpr[ok]

  if (length(fpr) < 2L) {
    return(NA_real_)
  }

  ord <- order(fpr, tpr)
  fpr <- fpr[ord]
  tpr <- tpr[ord]

  if (isTRUE(ensure_endpoints)) {

    if (fpr[1L] > 0 || tpr[1L] > 0) {
      fpr <- c(0, fpr)
      tpr <- c(0, tpr)
    }

    last <- length(fpr)

    if (fpr[last] < 1 || tpr[last] < 1) {
      fpr <- c(fpr, 1)
      tpr <- c(tpr, 1)
    }
  }

  sum(diff(fpr) * (head(tpr, -1L) + tail(tpr, -1L)) / 2)
}


#' Compute ROC curve and AUC from ML prediction data.
#'
#' The positive class is "Non_Normal". Larger values of
#' Prob_Non_Normal indicate stronger evidence against Normality.
get_ml_roc <- function(pred_df) {
  #' Extract labels and scores before constructing the ML ROC curve.
  required_cols <- c("True_Class", "Prob_Non_Normal")
  missing_cols  <- setdiff(required_cols, names(pred_df))
  if (length(missing_cols) > 0L) {
    stop("pred_df is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  #' Calculate the ROC curve.
  roc_obj <- pROC::roc(
    response = factor(pred_df$True_Class, levels = c("Normal", "Non_Normal")),
    predictor = pred_df$Prob_Non_Normal,
    levels    = c("Normal", "Non_Normal"),
    direction = "<",
    quiet     = TRUE
  )

  #' Assemble the result.
  list(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    auc = as.numeric(pROC::auc(roc_obj))
  )
}


#' Select plotting color.
get_plot_color <- function(name, color_map = NULL, default = "gray40") {
  if (!is.null(color_map) && name %in% names(color_map)) {
    return(color_map[[name]])
  }
  default
}


#' ============================================================
#' SECTION 25: ML-ONLY ROC PLOTS
#' ============================================================

#' Plot ML ROC curves for one sample size.
#'
#' This function uses stored prediction results. It does not rerun
#' simulations.
plot_ml_roc_one_n <- function(eval_results_n, n, ml_colors = NULL, main_title = NULL) {

  if (is.null(eval_results_n$predictions)) {
    stop("eval_results_n must contain a 'predictions' element.")
  }

  prediction_list <- eval_results_n$predictions
  model_names     <- names(prediction_list)
  if (length(model_names) == 0L) {
    stop("No model prediction data found.")
  }
  if (is.null(main_title)) {
    main_title <- paste0("ML ROC Curves (n = ", n, ")")
  }

  #' Initialize the ROC panel.
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "False Positive Rate",
    ylab = "True Positive Rate",
    main = main_title,
    las  = 1
  )

  #' Chance-performance reference.
  abline(0, 1, lty = 2, col = "gray60")

  legend_labels <- character(0)
  legend_cols   <- character(0)

  #' Add one ROC curve per model.
  for (model_name in model_names) {
    
    #' Model-specific ROC coordinates.
    roc_vals <- get_ml_roc(prediction_list[[model_name]])
    col_m    <- get_plot_color(model_name, ml_colors)

    lines(roc_vals$fpr, roc_vals$tpr, col = col_m, lwd = 2)

    #' Store the AUC label.
    legend_labels <- c(legend_labels, sprintf("%s (AUC = %.3f)", model_name, roc_vals$auc))
    legend_cols <- c(legend_cols, col_m)
  }

  #' Display model labels and AUC values.
  legend("bottomright", legend = legend_labels, col = legend_cols, lwd = 2, bty = "o", cex = 0.85)
  invisible(NULL)
}


#' Save ML ROC plots for all sample sizes.
save_ml_roc_plots <- function(eval_results_by_n, sample_sizes = names(eval_results_by_n),
                              output_dir = "results",ml_colors = NULL) {

  make_output_dir(output_dir)
  
  for (n in sample_sizes) {
    n_key <- as.character(n)
    if (!n_key %in% names(eval_results_by_n)) {
      warning("Skipping n = ", n, ": no evaluation result found.")
      next
    }

    pdf_file <- file.path( output_dir,paste0("ml_roc_n", n_key, ".pdf"))
    pdf(pdf_file, width = 6, height = 6)

    plot_ml_roc_one_n(
      eval_results_n = eval_results_by_n[[n_key]],
      n              = n_key,
      ml_colors      = ml_colors
    )

    dev.off()
    cat("Saved ML ROC plot:", pdf_file, "\n")
  }
  invisible(NULL)
}


#' ============================================================
#' SECTION 26: CLASSICAL NORMALITY-TEST ROC CURVES
#' ============================================================

#' Compute one classical normality-test p-value.
#'
#' Small p-values indicate evidence against Normality.
run_classical_test <- function(x, test, normal_par = NULL) {
  test <- toupper(trimws(test))
  if (test == "SW") {
    return(safe_calc(shapiro.test(x)$p.value))
  }
  if (test == "SF") {
    return(safe_calc(nortest::sf.test(x)$p.value))
  }
  if (test == "AD") {
    return(safe_calc(nortest::ad.test(x)$p.value))
  }
  if (test == "LF") {
    return(safe_calc(nortest::lillie.test(x)$p.value))
  }
  if (test == "JB") {
    return(safe_calc(tseries::jarque.bera.test(x)$p.value))
  }
  if (test == "KS") {
    #' Validate the condition.
    if (is.null(normal_par)) {
      mu <- mean(x)
      s  <- sd(x)
    } else {
      mu <- normal_par[1L]
      s  <- normal_par[2L]
    }
    if (!is.finite(s) || s <= 0) {
      return(NA_real_)
    }
    #' Return the result.
    return(safe_calc(ks.test(x, "pnorm", mean = mu, sd = s)$p.value))
  }
  stop("Unsupported classical test: ", test)
}


#' Compute p-values for one classical test under Normal and alternative samples.
simulate_all_classical_pvalues <- function(n, normal_spec, alt_spec, tests, 
                                           n_sim = 1000,show_progress = FALSE) {
  tests <- toupper(trimws(tests))
  #' Calculate p norm.
  p_norm <- matrix( NA_real_, nrow = n_sim, ncol = length(tests), dimnames = list(NULL, tests))
  p_alt <- p_norm
  progress_bar <- start_progress_bar(
    total = n_sim,
    label = "Classical ROC simulation",
    show_progress = show_progress
  )

  on.exit(close_progress_bar(progress_bar), add = TRUE)

  #' Generate each Normal-alternative sample pair once, then apply every test.
  for (iteration in seq_len(n_sim)) {
    x_norm <- generate_data( n = n, dist = normal_spec$dist, par = normal_spec$par)
    x_alt <- generate_data( n = n, dist = alt_spec$dist,par = alt_spec$par)

    #' Process each item.
    for (test_name in tests) {
      p_norm[iteration, test_name] <- run_classical_test(x = x_norm, test = test_name, normal_par = normal_spec$par)
      #' Continue the calculation.
      p_alt[iteration, test_name] <- run_classical_test(x = x_alt, test = test_name,normal_par = normal_spec$par)
    }
    update_progress_bar(progress_bar, iteration)
  }
  list(p_norm = p_norm, p_alt = p_alt)
}


#' Compute one test's p-values while retaining backward compatibility.
simulate_classical_pvalues <- function(n, normal_spec, alt_spec,test, n_sim = 1000) {
  pvalues <- simulate_all_classical_pvalues(n = n, normal_spec = normal_spec,
    alt_spec = alt_spec, tests = test, n_sim = n_sim,show_progress = FALSE)
  list(
    p_norm = pvalues$p_norm[, 1L],
    p_alt = pvalues$p_alt[, 1L]
  )
}


#' Convert Normal and alternative p-values into ROC coordinates.
pvalues_to_roc <- function(p_norm, p_alt, alpha_grid) {
  fpr <- vapply(alpha_grid, function(alpha) mean(p_norm < alpha, na.rm = TRUE),
    numeric(1L)
  )
  tpr <- vapply(alpha_grid, function(alpha) mean(p_alt < alpha, na.rm = TRUE),
    numeric(1L)
  )

  list(
    fpr = fpr,
    tpr = tpr
  )
}


#' Compute classical ROC curves for one Normal/alternative pair.
#'
#' For each classical normality test:
#' FPR = P(reject Normality | matched Normal)
#' TPR = P(reject Normality | alternative)
#'
#' Rejection rule:
#' p-value < alpha.
compute_classical_roc <- function(n, normal_spec, alt_spec, tests = c("SW", "AD", "JB"),
    alpha_grid = seq(0, 1, by = 0.05), n_sim = 1000, show_progress = FALSE
) {
  
  tests <- toupper(trimws(tests))
  
  #' Use common Monte Carlo samples across all classical tests.
  pvalues <- simulate_all_classical_pvalues(n = n, normal_spec = normal_spec, alt_spec = alt_spec,
                                          tests = tests, n_sim = n_sim, show_progress = show_progress)

  #' Calculate fpr.
  fpr <- matrix( NA_real_, nrow = length(tests), ncol = length(alpha_grid), dimnames = list(tests, NULL))
  tpr <- fpr

  #' Convert each test's stored p-values to ROC coordinates.
  for (test_name in tests) {
    roc_values <- pvalues_to_roc(
      p_norm = pvalues$p_norm[, test_name],
      p_alt = pvalues$p_alt[, test_name],
      alpha_grid = alpha_grid
    )

    fpr[test_name, ] <- roc_values$fpr
    tpr[test_name, ] <- roc_values$tpr
  }

  #' Assemble the result.
  list(
    FPR = fpr,
    TPR = tpr,
    alpha = alpha_grid,
    tests = tests
  )
}


#' ============================================================
#' SECTION 27: ML ROC CURVES FOR ONE DISTRIBUTION PAIR
#' ============================================================

#' Generate feature rows for one distribution pair.
generate_pair_prediction_features <- function(n, normal_spec, alt_spec, model_bundle, n_sim = 1000) {

  specs <- list(
    list(spec = normal_spec, label = "Normal"),
    list(spec = alt_spec,    label = "Non_Normal")
  )

  rows        <- vector("list", 2L * n_sim)
  true_class  <- character(2L * n_sim)
  row_id      <- 0L

  for (item in specs) {
    for (i in seq_len(n_sim)) {
      row_id <- row_id + 1L
      x <- generate_data(n = n, dist = item$spec$dist, par  = item$spec$par)
      rows[[row_id]] <- calculate_features(x = x,feature_set = model_bundle$feature_names)
      true_class[row_id] <- item$label
    }
  }

  feature_df <- do.call(rbind, rows)
  rownames(feature_df) <- NULL

  list(
    features   = feature_df,
    true_class = true_class
  )
}


#' Compute ML probabilities for one Normal/alternative pair.
compute_ml_pair_probabilities <- function(n, normal_spec, alt_spec, model_bundle,
                                          ml_methods = names(model_bundle$models),
                                          n_sim = 1000) {

  ml_methods <- intersect(ml_methods, names(model_bundle$models))

  if (length(ml_methods) == 0L) {
    stop("None of the requested ML methods are available in model_bundle.")
  }
  #' Calculate pair data.
  pair_data <- generate_pair_prediction_features(
    n            = n,
    normal_spec  = normal_spec,
    alt_spec     = alt_spec,
    model_bundle = model_bundle,
    n_sim        = n_sim
  )

  feats_std <- prepare_prediction_matrix(
    feature_df   = pair_data$features,
    model_bundle = model_bundle
  )

  all_probs <- predict_model_probs(
    model_bundle = model_bundle,
    feats_std    = feats_std
  )

  all_probs <- all_probs[, ml_methods, drop = FALSE]
  pred_store <- lapply(ml_methods, function(model_name) {
    
    #' Build the result table.
    data.frame(
      True_Class      = pair_data$true_class,
      Prob_Non_Normal = all_probs[[model_name]],
      stringsAsFactors = FALSE
    )
  })

  names(pred_store) <- ml_methods

  pred_store$MajorityVote <- data.frame(
    True_Class      = pair_data$true_class,
    Prob_Non_Normal = rowMeans(all_probs, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  pred_store
}


#' Compute ML ROC curves for one Normal/alternative pair.
#'
#' This function reruns pair-specific simulations and obtains classifier
#' probabilities. It is intended for direct classical-vs-ML ROC comparison
#' on the same distribution pair.
compute_ml_roc_pair <- function(n, normal_spec, alt_spec, model_bundle,
                                ml_methods = names(model_bundle$models),
                                n_sim = 1000) {

  #' Calculate pred store.
  pred_store <- compute_ml_pair_probabilities(
    n            = n,
    normal_spec  = normal_spec,
    alt_spec     = alt_spec,
    model_bundle = model_bundle,
    ml_methods   = ml_methods,
    n_sim        = n_sim
  )
  lapply(pred_store, get_ml_roc)
}


#' ============================================================
#' SECTION 28: COMBINED CLASSICAL-VS-ML ROC PLOT
#' ============================================================

#' Plot classical and ML ROC curves in one panel.
#'
#' Classical tests are dashed.
#' ML classifiers are solid.
#' Legend is placed inside each panel.
plot_classical_ml_roc <- function(classical_roc, ml_roc, main_title, test_colors = NULL,
                                  ml_colors = NULL,legend_cex = 0.70) {

  #' Initialize the ROC panel.
  plot(
    NA,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "False Positive Rate (1 - Specificity)",
    ylab = "True Positive Rate (Sensitivity)",
    main = main_title,
    las  = 1,
    cex.main = 0.85,
    font.main = 2
  )

  abline(0, 1, lty = 3, col = "gray60")

  legend_labels <- character(0)
  legend_cols   <- character(0)
  legend_lty    <- numeric(0)
  legend_lwd    <- numeric(0)

  #' ------------------------------------------------------------
  #' Classical ROC curves: dashed lines
  #' ------------------------------------------------------------

  for (test_name in rownames(classical_roc$FPR)) {
    auc_val <- compute_auc(
      fpr = classical_roc$FPR[test_name, ],
      tpr = classical_roc$TPR[test_name, ]
    )

    col_t <- get_plot_color(
      name      = test_name,
      color_map = test_colors,
      default   = "gray40"
    )

    #' Continue the calculation.
    lines(
      classical_roc$FPR[test_name, ],
      classical_roc$TPR[test_name, ],
      col = col_t,
      lwd = 1.5,
      lty = 2
    )

    #' Store the AUC label.
    legend_labels <- c( legend_labels,
      sprintf("%s (AUC = %.3f)", test_name, auc_val)
    )

    legend_cols <- c(legend_cols, col_t)
    legend_lty  <- c(legend_lty, 2)
    legend_lwd  <- c(legend_lwd, 1.5)
  }

  #' ------------------------------------------------------------
  #' ML ROC curves: solid lines
  #' ------------------------------------------------------------

  for (model_name in names(ml_roc)) {
    col_m <- get_plot_color(
      name      = model_name,
      color_map = ml_colors,
      default   = "gray40"
    )

    lines(ml_roc[[model_name]]$fpr, ml_roc[[model_name]]$tpr, col = col_m, lwd = 1.8, lty = 1)

    #' Store the AUC label.
    legend_labels <- c(legend_labels,
      sprintf("%s (AUC = %.3f)", model_name, ml_roc[[model_name]]$auc)
    )

    legend_cols <- c(legend_cols, col_m)
    legend_lty  <- c(legend_lty, 1)
    legend_lwd  <- c(legend_lwd, 1.8)
  }

  #' Display model labels and AUC values.
  legend("bottomright", legend = legend_labels,
    col    = legend_cols,
    lty    = legend_lty,
    lwd    = legend_lwd,
    bty    = "o",
    cex    = legend_cex,
    bg     = "white"
  )

  invisible(NULL)
}

#' ============================================================
#' SECTION 29: RUN ROC COMPARISON ACROSS PAIRS AND SAMPLE SIZES
#' ============================================================

#' Create a safe label for file names.
safe_dist_label <- function(spec) {
  label <- paste0( spec$dist, "_", paste(spec$par, collapse = "_"))
  gsub("[^A-Za-z0-9_]+", "_", label)
}

#' Create a readable distribution label.
roc_dist_label <- function(spec, digits = 3) {
  dist_name <- tools::toTitleCase( gsub("_", " ", spec$dist))
  par_text <- paste(round(spec$par, digits),collapse = ", ")
  paste0(dist_name, "(", par_text, ")")
}

#' ============================================================
#' SECTION 29: RUN ROC COMPARISON ACROSS PAIRS AND SAMPLE SIZES
#' ============================================================

#' Build pair-specific ML ROC curves from saved evaluation predictions.
#'
#' Reusing held-out predictions avoids regenerating the same evaluation samples
#' and recalculating their features during the classical-versus-ML ROC stage.
get_ml_roc_from_evaluation <- function(eval_results_n, pair_id, normal_spec,
    alt_spec, ml_methods) {
  if (is.null(eval_results_n$predictions)) {
    stop("Evaluation results do not contain prediction tables.")
  }
  prediction_list <- eval_results_n$predictions
  requested_methods <- unique(c(ml_methods, "MajorityVote"))
  requested_methods <- intersect(requested_methods, names(prediction_list))

  if (length(requested_methods) == 0L) {
    stop("None of the requested ML methods are available in evaluation results.")
  }

  normal_label <- make_distribution_label(normal_spec)
  alt_label <- make_distribution_label(alt_spec)

  roc_results <- lapply(
    requested_methods,
    function(model_name) {
      prediction_df <- prediction_list[[model_name]]

      #' Prefer the exact pair identifier stored by the optimized evaluation.
      if ("Pair_ID" %in% names(prediction_df)) {
        pair_df <- prediction_df[prediction_df$Pair_ID == pair_id, , drop = FALSE]
      } else {
        #' Fall back to distribution labels for evaluation files made earlier.
        pair_df <- prediction_df[prediction_df$Distribution %in% c(normal_label, alt_label),
          ,drop = FALSE]
      }

      #' Validate the condition.
      if (nrow(pair_df) == 0L || length(unique(pair_df$True_Class)) < 2L) {
        stop("Could not recover both classes for ROC pair ", pair_id," and model ", model_name, "."
        )
      }
      get_ml_roc(pair_df)
    }
  )

  names(roc_results) <- requested_methods
  roc_results
}

#' ------------------------------------------------------------
#' Compare classical normality tests and ML classifiers using ROC curves.
#' ------------------------------------------------------------
run_roc_comparison <- function(
    trained_models,
    roc_alt_specs,
    sample_sizes = names(trained_models),
    classical_tests = c("SW", "AD", "JB"),
    ml_methods = c("RF", "GBM", "ANN"),
    alpha_grid = seq(0, 1, by = 0.05),
    n_sim = 1000,
    output_dir = "results",
    test_colors = NULL,
    ml_colors = NULL,
    n_cores = 1L,
    eval_results_by_n = NULL,
    reuse_evaluation_predictions = TRUE,
    show_progress = TRUE
) {
  if (length(trained_models) == 0L) {
    stop("trained_models is empty.")
  }

  if (length(roc_alt_specs) == 0L) {
    stop("roc_alt_specs is empty.")
  }

  make_output_dir(output_dir)
  
  roc_pairs <- make_paired_specs(roc_alt_specs)
  roc_classical <- list()
  roc_ml <- list()
  sample_sizes <- as.numeric(sample_sizes)
  n_cores <- normalize_n_cores(n_cores, max_tasks = length(sample_sizes))

  pair_bar <- start_progress_bar(
    total = length(roc_pairs),
    label = "Classical-versus-ML ROC pairs",
    show_progress = show_progress
  )

  on.exit(close_progress_bar(pair_bar), add = TRUE)

  #' Process each alternative pair and save one PDF per pair.
  for (pair_id in seq_along(roc_pairs)) {
    pair <- roc_pairs[[pair_id]]
    pair_start <- Sys.time()

    alt_name <- safe_dist_label(pair$alt)
    alt_label <- roc_dist_label(pair$alt)
    norm_label <- roc_dist_label(pair$normal)

    cat("\nRunning ROC comparison for ", alt_label, " vs ", norm_label, " ...\n",sep = "")

    #' Compute the sample-size ROC results in parallel when requested.
    roc_by_n <- safe_mclapply(
      X = sample_sizes,
      FUN = function(n) {
        n_key <- as.character(n)

        if (!n_key %in% names(trained_models)) {
          stop("No trained model found for n = ", n_key)
        }

        #' Calculate classical result.
        classical_result <- compute_classical_roc(
          n = n,
          normal_spec = pair$normal,
          alt_spec = pair$alt,
          tests = classical_tests,
          alpha_grid = alpha_grid,
          n_sim = n_sim,
          show_progress = show_progress && n_cores == 1L
        )

        can_reuse <- isTRUE(reuse_evaluation_predictions) &&
          !is.null(eval_results_by_n) &&
          n_key %in% names(eval_results_by_n)

        if (can_reuse) {
          ml_result <- get_ml_roc_from_evaluation(
            eval_results_n = eval_results_by_n[[n_key]],
            pair_id = pair_id,
            normal_spec = pair$normal,
            alt_spec = pair$alt,
            ml_methods = ml_methods
          )
        } else {
          #' Retain the original simulation path as a compatibility fallback.
          ml_result <- compute_ml_roc_pair(
            n = n,
            normal_spec = pair$normal,
            alt_spec = pair$alt,
            model_bundle = trained_models[[n_key]],
            ml_methods = ml_methods,
            n_sim = n_sim
          )
        }

        #' Assemble the result.
        list(
          n = n,
          n_key = n_key,
          classical = classical_result,
          ml = ml_result
        )
      },
      n_cores = n_cores
    )

    names(roc_by_n) <- as.character(sample_sizes)

    #' Store the numerical ROC results before plotting.
    for (result in roc_by_n) {
      result_key <- paste0(alt_name, "_n", result$n_key)
      roc_classical[[result_key]] <- result$classical
      roc_ml[[result_key]] <- result$ml
    }

    pdf_file <- file.path(output_dir, paste0("roc_comparison_", alt_name, ".pdf"))

    grDevices::pdf( file = pdf_file,
      width = 5.8 * length(sample_sizes),
      height = 5.8
    )

    #' Continue the calculation.
    graphics::par(
      mfrow = c(1, length(sample_sizes)),
      mar = c(4.4, 4.5, 3.0, 1.0),
      oma = c(1.0, 1.0, 3.0, 1.0),
      mgp = c(2.6, 0.8, 0)
    )

    #' Draw each sample-size panel on the open PDF device.
    for (result in roc_by_n) {
      panel_title <- paste0(alt_label, "  vs  ", norm_label,"  |  n = ", result$n)

      #' Continue the calculation.
      plot_classical_ml_roc(
        classical_roc = result$classical,
        ml_roc = result$ml,
        main_title = panel_title,
        test_colors = test_colors,
        ml_colors = ml_colors
      )
    }

    #' Continue the calculation.
    graphics::mtext(
      paste0("ROC Comparison: Classical vs ML  |  ",alt_label, "  vs  ", norm_label),
      outer = TRUE,
      cex = 1.0,
      font = 1,
      line = 1.0
    )

    grDevices::dev.off()
    cat("Saved ROC comparison: ", pdf_file," (", format_elapsed_time(pair_start), ")\n",sep = "")
    update_progress_bar(pair_bar, pair_id)
  }
  invisible(list(classical = roc_classical, ml = roc_ml))
}


#' ============================================================
#' PART 6: PERMUTATION AUC VARIABLE IMPORTANCE
#' ============================================================
#' 
#' ============================================================
#' SECTION 30: PERMUTATION VIP HELPERS
#' ============================================================
#' 
#' Compute AUC from true labels and predicted probabilities.
#'
#' The positive class is "Non_Normal".
compute_prob_auc <- function(true_class,prob_non_normal) {

  if (length(unique(true_class)) < 2L) {
    return(NA_real_)
  }
  if (!all(is.finite(prob_non_normal))) {
    return(NA_real_)
  }
  #' Continue the calculation.
  safe_calc(
    as.numeric(
      pROC::auc(
        pROC::roc(
          response = factor(
            true_class,
            levels = c("Normal", "Non_Normal")
          ),
          predictor = prob_non_normal,
          levels    = c("Normal", "Non_Normal"),
          direction = "<",
          quiet     = TRUE
        )
      )
    ),
    default = NA_real_
  )
}


#' Extract held-out labels from one sample-size evaluation object.
#'
#' All classifiers are evaluated on the same held-out samples, so
#' the labels can be taken from the first prediction data frame.
get_holdout_labels <- function(eval_results_n) {
  if (is.null(eval_results_n$predictions)) {
    stop("eval_results_n must contain a 'predictions' element.")
  }
  prediction_list <- eval_results_n$predictions
  if (length(prediction_list) == 0L) {
    stop("No prediction data found in eval_results_n$predictions.")
  }
  first_model <- names(prediction_list)[1L]
  prediction_list[[first_model]]$True_Class
}


#' Predict P(Non_Normal) for one trained model.
predict_one_model_prob <- function(model,new_data) {
  pred <- predict(
    model,
    newdata = new_data,
    type    = "prob"
  )
  if (!"Non_Normal" %in% names(pred)) {
    stop("Model prediction does not contain a 'Non_Normal' probability column.")
  }
  as.numeric(pred[, "Non_Normal"])
}



#' ============================================================
#' SECTION 31: PERMUTATION VIP FOR ONE MODEL
#' ============================================================
#' Compute permutation AUC-drop importance for one model.
#'
#' For each feature:
#' 1. permute the feature in the held-out data(from training set of dist.);
#' 2. recompute predicted probabilities;
#' 3. recompute AUC;
#' 4. record baseline AUC - permuted AUC.
compute_vip_one_model <- function(model,
                                  model_name,
                                  holdout_x,
                                  holdout_y,
                                  n_permutations = 20,
                                  digits = 4,
                                  n_cores = 1) {

  if (ncol(holdout_x) == 0L) {
    stop("holdout_x has no feature columns.")
  }

  feature_names <- names(holdout_x)

  base_prob <- predict_one_model_prob(
    model    = model,
    new_data = holdout_x
  )

  baseline_auc <- compute_prob_auc(
    true_class      = holdout_y,
    prob_non_normal = base_prob
  )

  #' Continue the calculation.
  cat(
    sprintf("  %s baseline AUC = %.4f | %d features x %d permutations | cores = %d\n",
      model_name,
      baseline_auc,
      length(feature_names),
      n_permutations,
      n_cores
    )
  )

  feature_results <- safe_mclapply(
    X = feature_names,
    FUN = function(feature) {
      perm_auc <- numeric(n_permutations)
      for (b in seq_len(n_permutations)) {

        perm_x <- holdout_x
        perm_x[[feature]] <- sample(perm_x[[feature]])
        #' Calculate perm prob.
        perm_prob <- safe_calc(
          predict_one_model_prob(
            model    = model,
            new_data = perm_x
          ),
          default = rep(NA_real_, nrow(perm_x))
        )

        perm_auc[b] <- compute_prob_auc(
          true_class      = holdout_y,
          prob_non_normal = perm_prob
        )
      }

      #' Build the result table.
      data.frame(
        Feature  = feature,
        AUC_Drop = baseline_auc - mean(perm_auc, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    },
    n_cores = n_cores
  )

  imp_df <- do.call(rbind, feature_results)
  imp_df$AUC_Drop <- round(imp_df$AUC_Drop, digits)
  imp_df <- imp_df[order(imp_df$AUC_Drop, decreasing = TRUE), ,drop = FALSE]

  imp_df$Rank <- seq_len(nrow(imp_df))
  imp_df <- imp_df[, c("Rank", "Feature", "AUC_Drop")]
  rownames(imp_df) <- NULL

  #' Assemble the result.
  list(
    model_name   = model_name,
    baseline_auc = baseline_auc,
    imp_df       = imp_df
  )
}


#' ============================================================
#' SECTION 32: PERMUTATION VIP FOR ALL MODELS AND SAMPLE SIZES
#' ============================================================

#' Compute permutation AUC-drop importance for selected models.
#'
#' This requires evaluate_models_all_n(..., store_holdout_features = TRUE).
#' The function is optional because it can be time-consuming.
compute_permutation_vip <- function(trained_models,
                                    eval_results_by_n,
                                    n_permutations = 20,
                                    models_for_vip = NULL,
                                    digits = 4,
                                    n_cores = 1) {

  if (length(trained_models) == 0L) {
    stop("trained_models is empty.")
  }

  vip_results <- vector("list", length(trained_models))
  names(vip_results) <- names(trained_models)

  for (n_key in names(trained_models)) {
    cat("\nPermutation VIP for n =", n_key, "\n")
    if (!n_key %in% names(eval_results_by_n)) {
      stop("No evaluation results found for n = ", n_key)
    }

    model_bundle <- trained_models[[n_key]]
    eval_n       <- eval_results_by_n[[n_key]]

    #' Validate the condition.
    if (is.null(eval_n$holdout_features)) {
      stop("Missing holdout_features for n = ",n_key,
           ". Rerun evaluate_models_all_n(..., store_holdout_features = TRUE).")
    }

    holdout_x <- eval_n$holdout_features
    holdout_y <- get_holdout_labels(eval_n)
    model_names <- names(model_bundle$models)

    if (!is.null(models_for_vip)) {
      model_names <- intersect(model_names, models_for_vip)
    }

    if (length(model_names) == 0L) {
      stop("No requested VIP models found for n = ", n_key)
    }

    vip_results[[n_key]] <- vector("list", length(model_names))
    names(vip_results[[n_key]]) <- model_names

    for (model_name in model_names) {

      #' Continue the calculation.
      vip_results[[n_key]][[model_name]] <- compute_vip_one_model(
        model          = model_bundle$models[[model_name]],
        model_name     = model_name,
        holdout_x      = holdout_x,
        holdout_y      = holdout_y,
        n_permutations = n_permutations,
        digits         = digits,
        n_cores        = n_cores
      )
      vip_results[[n_key]][[model_name]]$n <- model_bundle$n
    }
  }

  vip_results
}


#' ============================================================
#' SECTION 33: SUMMARIZE AND SELECT PERMUTATION VIP FEATURES
#' ============================================================

#' Convert permutation VIP results to long format.
vip_results_to_long <- function(vip_results) {
  #' Reshape model-specific importance results for plotting and export.
  rows <- list()
  for (n_key in names(vip_results)) {
    vip_n <- vip_results[[n_key]]
    if (is.null(vip_n)) {
      next
    }
    for (model_name in names(vip_n)) {
      entry <- vip_n[[model_name]]
      if (is.null(entry)) {
        next
      }

      #' Build the result table.
      rows[[length(rows) + 1L]] <- data.frame(
        n        = entry$n,
        Model    = model_name,
        Feature  = entry$imp_df$Feature,
        AUC_Drop = entry$imp_df$AUC_Drop,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0L) {
    return(NULL)
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Summarize permutation VIP across models and sample sizes.
summarize_permutation_vip <- function(vip_results,digits = 4) {
  vip_long <- vip_results_to_long(vip_results)
  if (is.null(vip_long) || nrow(vip_long) == 0L) {
    return(NULL)
  }

  split_rows <- split(vip_long, vip_long$Feature)
  
  summary_df <- do.call(
    rbind,
    lapply(names(split_rows), function(feature) {
      z <- split_rows[[feature]]$AUC_Drop
      #' Build the result table.
      data.frame(
        Feature         = feature,
        Mean_AUC_Drop   = round(mean(z, na.rm = TRUE), digits),
        Median_AUC_Drop = round(median(z, na.rm = TRUE), digits),
        Min_AUC_Drop    = round(min(z, na.rm = TRUE), digits),
        Max_AUC_Drop    = round(max(z, na.rm = TRUE), digits),
        Pct_Positive    = round(mean(z > 0, na.rm = TRUE), digits),
        stringsAsFactors = FALSE
      )
    })
  )

  summary_df <- summary_df[order(summary_df$Mean_AUC_Drop, decreasing = TRUE),,drop = FALSE]
  rownames(summary_df) <- NULL
  summary_df
}


#' Select top features from a permutation VIP summary table.
select_top_vip_features <- function(vip_summary,
                                    top_n = 30,
                                    rank_by = c("mean", "median", "min"),
                                    min_pct_positive = 0,
                                    require_positive_mean = TRUE) {

  rank_by <- match.arg(rank_by)

  if (is.null(vip_summary) || nrow(vip_summary) == 0L) {
    stop("vip_summary is empty.")
  }

  #' Calculate needed cols.
  needed_cols <- c("Feature", "Mean_AUC_Drop", "Median_AUC_Drop", "Min_AUC_Drop", "Pct_Positive")
  missing_cols <- setdiff(needed_cols, names(vip_summary))
  
  #' Validate the condition.
  if (length(missing_cols) > 0L) {
    stop("vip_summary is missing required column(s): ", paste(missing_cols, collapse = ", ")
    )
  }

  vip_keep <- vip_summary[vip_summary$Pct_Positive >= min_pct_positive, ,drop = FALSE]

  #' Validate the condition.
  if (isTRUE(require_positive_mean)) {
    vip_keep <- vip_keep[vip_keep$Mean_AUC_Drop > 0, , drop = FALSE]
  }

  if (nrow(vip_keep) == 0L) {
    warning("No features met the selection criteria.")
    return(character(0))
  }

  #' Calculate rank col.
  rank_col <- switch(
    rank_by,
    mean   = "Mean_AUC_Drop",
    median = "Median_AUC_Drop",
    min    = "Min_AUC_Drop"
  )

  vip_keep <- vip_keep[order(vip_keep[[rank_col]], decreasing = TRUE), ,drop = FALSE]
  head(vip_keep$Feature, top_n)
}



#' ------------------------------------------------------------
#' Remove highly correlated features from a feature matrix
#' ------------------------------------------------------------
#'
#' This is used for the pre-VIP correlation screening step.
#' It is intentionally unsupervised and should be used with a
#' conservative cutoff such as 0.97 or 0.98.
#' ------------------------------------------------------------

remove_highly_correlated_features <- function(feature_df,
                                              candidate_features = NULL,
                                              cutoff = 0.98,
                                              impute_method = c("median", "mean"),
                                              cor_method = c("spearman", "pearson")) {

  impute_method <- match.arg(impute_method)
  cor_method    <- match.arg(cor_method)

  if (is.null(feature_df) || nrow(feature_df) == 0L) {
    stop("feature_df is empty or NULL.")
  }

  if (!is.null(candidate_features)) {
    missing_features <- setdiff(candidate_features, names(feature_df))
    #' Validate the condition.
    if (length(missing_features) > 0L) {
      stop("The following candidate features are missing from feature_df: ",
        paste(missing_features, collapse = ", ")
      )
    }
    feature_df <- feature_df[, candidate_features, drop = FALSE]
  }

  numeric_cols <- vapply(feature_df, is.numeric, logical(1L))
  feature_df   <- feature_df[, numeric_cols, drop = FALSE]

  if (ncol(feature_df) == 0L) {
    stop("No numeric feature columns available for correlation filtering.")
  }

  feature_df[] <- lapply(feature_df, function(z) {
    z <- as.numeric(z)
    z[!is.finite(z)] <- NA_real_
    z
  })

  keep_cols <- vapply(feature_df, function(z) {
    finite_z <- z[is.finite(z)]
    length(finite_z) > 1L && length(unique(finite_z)) > 1L
  }, logical(1L))

  feature_df <- feature_df[, keep_cols, drop = FALSE]

  if (ncol(feature_df) <= 1L) {
    return(names(feature_df))
  }

  feature_df <- impute_na(feature_df, method = impute_method)

  #' Calculate cor mat.
  cor_mat <- suppressWarnings(cor(feature_df, use    = "pairwise.complete.obs", method = cor_method))
  cor_mat[!is.finite(cor_mat)] <- 0
  diag(cor_mat) <- 1

  #' Calculate drop idx.
  drop_idx <- caret::findCorrelation(cor_mat, cutoff = cutoff, names  = FALSE, exact  = TRUE)

  if (length(drop_idx) == 0L) {
    return(names(feature_df))
  }

  names(feature_df)[-drop_idx]
}


#' ------------------------------------------------------------
#' VIP-priority correlation filter
#' ------------------------------------------------------------
#'
#' Keeps features in the order provided by vip_ordered_features.
#' A feature is retained if it is not highly correlated with any
#' feature already retained.
#'
#' This avoids the problem where caret::findCorrelation() removes
#' a highly ranked VIP feature and keeps a lower-ranked correlated one.
#' ------------------------------------------------------------

remove_correlated_features_by_vip_priority <- function(feature_df,
                                                       vip_ordered_features,
                                                       cutoff = 0.95,
                                                       impute_method = c("median", "mean"),
                                                       cor_method = c("spearman", "pearson")) {

  impute_method <- match.arg(impute_method)
  cor_method    <- match.arg(cor_method)

  if (is.null(feature_df) || nrow(feature_df) == 0L) {
    stop("feature_df is empty or NULL.")
  }

  if (length(vip_ordered_features) == 0L) {
    stop("vip_ordered_features is empty.")
  }

  missing_features <- setdiff(vip_ordered_features, names(feature_df))
  #' Validate the condition.
  if (length(missing_features) > 0L) {
    stop("The following VIP features are missing from feature_df: ",
      paste(missing_features, collapse = ", ")
    )
  }

  x <- feature_df[, vip_ordered_features, drop = FALSE]
  x[] <- lapply(x, function(z) {
    z <- as.numeric(z)
    z[!is.finite(z)] <- NA_real_
    z
  })

  #' Remove all-missing or constant columns.
  keep_cols <- vapply(x, function(z) {
    finite_z <- z[is.finite(z)]
    length(finite_z) > 1L && length(unique(finite_z)) > 1L
  }, logical(1L))

  x <- x[, keep_cols, drop = FALSE]
  vip_ordered_features <- vip_ordered_features[vip_ordered_features %in% names(x)]
  if (length(vip_ordered_features) <= 1L) {
    return(vip_ordered_features)
  }

  x <- impute_na(x, method = impute_method)
  #' Calculate cor mat.
  cor_mat <- suppressWarnings(cor(x, use = "pairwise.complete.obs", method = cor_method ))
  cor_mat[!is.finite(cor_mat)] <- 0
  diag(cor_mat) <- 1

  retained <- character(0)
  removed  <- character(0)

  for (f in vip_ordered_features) {
    if (length(retained) == 0L) {
      retained <- c(retained, f)
      next
    }
    max_abs_cor <- max(abs(cor_mat[f, retained]), na.rm = TRUE)
    #' Validate the condition.
    if (is.finite(max_abs_cor) && max_abs_cor >= cutoff) {
      removed <- c(removed, f)
    } else {
      retained <- c(retained, f)
    }
  }

  attr(retained, "removed_by_correlation") <- removed
  retained
}


#' ------------------------------------------------------------
#' Screen highly correlated features across sample sizes
#' ------------------------------------------------------------
#'
#' This function generates a pilot feature matrix for each sample size,
#' removes near-duplicate features within each sample size, and returns
#' the union of retained features.
#'
#' The union is used because a feature may be redundant at n = 50 but
#' useful at n = 10, or vice versa.
#' ------------------------------------------------------------

#' Calculate screen correlated features across n.
screen_correlated_features_across_n <- function(sample_sizes,
                                                paired_specs,
                                                candidate_features,
                                                num_sim = 300,
                                                cutoff = 0.98,
                                                standardize_sample = FALSE,
                                                center_by = NULL,
                                                impute_method = "median",
                                                cor_method = "spearman") {

  if (length(sample_sizes) == 0L) {
    stop("sample_sizes must contain at least one value.")
  }
  if (length(candidate_features) == 0L) {
    stop("candidate_features must contain at least one feature.")
  }

  retained_by_n <- vector("list", length(sample_sizes))
  names(retained_by_n) <- as.character(sample_sizes)

  for (n in sample_sizes) {
    n_key <- as.character(n)
    cat("\n")
    cat("Correlation screening for n =", n, "\n")
    #' Calculate pilot data.
    pilot_data <- make_training_data(
      n                  = n,
      paired_specs       = paired_specs,
      num_sim            = num_sim,
      feature_set        = candidate_features,
      standardize_sample = standardize_sample,
      center_by          = center_by
    )

    feature_df <- pilot_data[, setdiff(names(pilot_data), "Label"), drop = FALSE]

    #' Calculate retained.
    retained <- remove_highly_correlated_features(
      feature_df          = feature_df,
      candidate_features  = candidate_features,
      cutoff              = cutoff,
      impute_method       = impute_method,
      cor_method          = cor_method
    )

    retained_by_n[[n_key]] <- retained
    
    cat("  Candidate features:", length(candidate_features), "\n")
    cat("  Retained features :", length(retained), "\n")
    cat("  Removed features  :", length(setdiff(candidate_features, retained)), "\n")
  }

  #' Use the union so features useful at one sample size are not lost.
  retained_union <- Reduce(union, retained_by_n)
  retained_union <- candidate_features[candidate_features %in% retained_union]
  cat("\n")
  cat("Final retained feature count across all n:", length(retained_union), "\n")
  attr(retained_union, "retained_by_n") <- retained_by_n
  retained_union
}

#' Return selected features together with their VIP metrics.
get_selected_feature_metrics <- function(vip_summary, selected_features) {
  
  out <- vip_summary[match(selected_features, vip_summary$Feature),, drop = FALSE]
  out <- data.frame(Rank = seq_len(nrow(out)), out, row.names = NULL)
  out
}


#' ============================================================
#' FINAL FEATURE FILTER
#' ============================================================

#' Select a stable final feature set from the permutation VIP summary.
#'
#' This is intended for the final feature-selection stage, where both
#' magnitude and stability of the AUC drop are considered.
final_feature_filter <- function(vip_summary,
                                 min_pct_positive = 0.75,
                                 min_median_drop = 0,
                                 allow_small_negative_min = TRUE,
                                 min_mean_for_negative_min = 0.01,
                                 top_n = 30) {

  #' Calculate required cols.
  required_cols <- c("Feature", "Mean_AUC_Drop", "Median_AUC_Drop", "Min_AUC_Drop", "Pct_Positive")
  missing_cols <- setdiff(required_cols, names(vip_summary))

  #' Validate the condition.
  if (length(missing_cols) > 0L) {
    stop("vip_summary is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  out <- vip_summary

  #' Store the result.
  out <- out[out$Pct_Positive >= min_pct_positive & out$Median_AUC_Drop > min_median_drop & out$Mean_AUC_Drop > 0, , drop = FALSE]
  if (isTRUE(allow_small_negative_min)) {
    #' Store the result.
    out <- out[out$Min_AUC_Drop >= 0 | out$Mean_AUC_Drop >= min_mean_for_negative_min, ,drop = FALSE]
  } else {
    #' Store the result.
    out <- out[out$Min_AUC_Drop >= 0, ,drop = FALSE ]
  }
  out <- out[order(out$Mean_AUC_Drop, decreasing = TRUE), , drop = FALSE]
  out <- head(out, top_n)
  rownames(out) <- NULL
  out
}

#' ============================================================
#' SECTION 34: PERMUTATION VIP PLOTS
#' ============================================================

#' Plot permutation VIP for one model and one sample size.
plot_permutation_vip <- function(vip_entry, top_n = 30, selected_features = NULL,main_title = NULL) {

  if (is.null(vip_entry)) {
    stop("vip_entry is NULL.")
  }

  imp_df <- vip_entry$imp_df
  if (is.null(imp_df) || nrow(imp_df) == 0L) {
    stop("vip_entry$imp_df is empty.")
  }

  #' Plot the final selected features when supplied.
  if (!is.null(selected_features)) {
    selected_features <- unique(as.character(selected_features))

    missing_features <- setdiff(selected_features,imp_df$Feature)

    #' Validate the condition.
    if (length(missing_features) > 0L) {
      warning("Selected features unavailable in this VIP table: ", paste(missing_features, collapse = ", "))
    }

    #' Calculate imp df.
    imp_df <- imp_df[imp_df$Feature %in% selected_features, , drop = FALSE]
  } else {
    imp_df <- head(imp_df, top_n)
  }
  
  if (nrow(imp_df) == 0L) {
    stop("No features are available for the VIP plot.")
  }
  #' Sort so the largest AUC drop appears at the top.
  imp_df <- imp_df[order(imp_df$AUC_Drop, decreasing = FALSE), ,drop = FALSE]
  #' Validate the condition.
  if (is.null(main_title)) {
    title_prefix <- if (is.null(selected_features)) {
      "Permutation VIP: "
    } else {
      "Permutation VIP for final selected features: "
    }
    #' Calculate main title.
    main_title <- paste0(title_prefix, vip_entry$model_name," (n = ", vip_entry$n,")")
  }

  old_mar <- par("mar")
  on.exit(par(mar = old_mar), add = TRUE)
  par(mar = c(5, 12, 4, 2))

  #' Draw the plot.
  barplot(
    height = imp_df$AUC_Drop,
    names.arg = imp_df$Feature,
    horiz = TRUE,
    las = 1,
    border = NA,
    col = "steelblue",
    xlab = "AUC Drop After Permutation",
    main = main_title,
    cex.names = 0.82,
    cex.main = 0.95
  )

  abline(v = 0, lty = 2, col = "gray50")
  invisible(imp_df)
}


#' ============================================================
#' COMBINED PERMUTATION VIP PLOTS
#' ============================================================

#' Save permutation AUC VIP plots for two sample sizes in a 1 x 2 layout.
save_permutation_vip_1x2_plot <- function(
    vip_results,
    sample_sizes = c(10, 50),
    model = "RF",
    top_n = 25,
    selected_features = NULL,
    output_dir = "results",
    file_name = NULL,
    width = 15,
    height = NULL,
    minimum_height = 8,
    height_per_feature = 0.33
) {

  if (length(sample_sizes) != 2L) {
    stop("sample_sizes must contain exactly two sample sizes.")
  }

  make_output_dir(output_dir)

  #' Validate the condition.
  if (is.null(file_name)) {
    file_name <- paste0("vip_perm_1x2_", model, "_n", paste(sample_sizes, collapse = "_n"), ".pdf")
  }

  pdf_file <- file.path(output_dir, file_name)

  #' Increase height as the number of displayed features increases.
  displayed_feature_count <- if (is.null(selected_features)) {
    top_n
  } else {
    length(unique(selected_features))
  }

  #' Validate the condition.
  if (is.null(height)) {
    height <- max( minimum_height, 2.5 + height_per_feature * displayed_feature_count)
  }

  pdf(file = pdf_file, width = width, height = height)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  #' Continue the calculation.
  par(
    mfrow = c(1, 2),
    mar = c(5.0, 9.5, 3.5, 1.0),
    oma = c(1.0, 1.0, 2.5, 1.0),
    mgp = c(2.6, 0.8, 0)
  )

  for (n in sample_sizes) {
    n_key <- as.character(n)

    #' Validate the condition.
    if (!n_key %in% names(vip_results)) {
      warning("Skipping n = ", n, ": not found in vip_results.")
      plot.new()
      title(main = paste0("n = ", n, " not available"))
      next
    }

    #' Validate the condition.
    if (!model %in% names(vip_results[[n_key]])) {
      warning("Skipping model ", model, " for n = ", n, ".")
      plot.new()
      title(main = paste0(model, ", n = ", n, " not available"))
      next
    }

    imp_df <- vip_results[[n_key]][[model]]$imp_df

    #' Validate the condition.
    if (is.null(imp_df) || nrow(imp_df) == 0L) {
      warning("No VIP table found for model ", model, " and n = ", n, ".")
      plot.new()
      title(main = paste0(model, ", n = ", n, " empty"))
      next
    }

    #' Use the exact final selected feature set when supplied.
    if (!is.null(selected_features)) {
      imp_df <- imp_df[imp_df$Feature %in% selected_features, , drop = FALSE]
    } else {
      imp_df <- imp_df[seq_len(min(top_n, nrow(imp_df))), , drop = FALSE]
    }

    #' Validate the condition.
    if (nrow(imp_df) == 0L) {
      warning("No selected features were found for n = ", n, ".")
      plot.new()
      title(main = paste0(model, ", n = ", n, " empty"))
      next
    }

    imp_df <- imp_df[order(imp_df$AUC_Drop, decreasing = FALSE), , drop = FALSE]

    #' Draw the plot.
    barplot(
      height = imp_df$AUC_Drop,
      names.arg = imp_df$Feature,
      horiz = TRUE,
      las = 1,
      cex.names = 0.80,
      xlab = "AUC Drop",
      main = paste0(model, " permutation VIP, n = ", n),
      border = NA,
      col = "steelblue"
    )

    abline(v = 0, lty = 3, col = "gray50")
  }

  #' Continue the calculation.
  mtext(paste0("Permutation AUC Variable Importance: ", model),
    outer = TRUE,
    cex = 1.05,
    font = 2,
    line = 1.0
  )

  dev.off()
  cat("\nSaved combined permutation VIP plot to:", pdf_file, "\n")
  invisible(pdf_file)
}


#' ---------------------------------------------------
#' Save VIP plots for selected models and sample sizes.
save_permutation_vip_plots <- function(
    vip_results,
    sample_sizes = names(vip_results),
    models = c("RF"),
    top_n = 30,
    selected_features = NULL,
    output_dir = "results",
    width = 9,
    height = NULL,
    minimum_height = 8,
    height_per_feature = 0.33
) {

  make_output_dir(output_dir)

  for (n in sample_sizes) {
    n_key <- as.character(n)
    if (!n_key %in% names(vip_results)) {
      warning("Skipping n = ", n_key, ": no VIP result found.")
      next
    }

    for (model_name in models) {
      vip_entry <- vip_results[[n_key]][[model_name]]
      if (is.null(vip_entry)) {
        next
      }

      pdf_file <- file.path(output_dir, paste0("vip_perm_", model_name, "_n", n_key, ".pdf"))

      displayed_feature_count <- if (is.null(selected_features)) {
        top_n
      } else {
        length(unique(selected_features))
      }

      #' Calculate plot height.
      plot_height <- if (is.null(height)) {
        max(minimum_height, 2.5 + height_per_feature * displayed_feature_count)
      } else {
        height
      }

      pdf(file = pdf_file, width = width, height = plot_height)

      plot_permutation_vip(
        vip_entry = vip_entry,
        top_n = top_n,
        selected_features = selected_features
      )

      dev.off()

      cat("Saved permutation VIP plot:", pdf_file, "\n")
    }
  }

  invisible(NULL)
}


#' ============================================================
#' PART 7: DENSITY PLOTS FOR MOMENT-MATCHED PAIRS
#' ============================================================


#' ============================================================
#' SECTION 35: DENSITY HELPER FUNCTIONS
#' ============================================================

#' Evaluate the theoretical density for a distribution specification.
#'
#' Used only for plotting.
#'
#' Required structure:
#' spec$dist
#' spec$par
density_value <- function(x_grid, spec) {

  #' Standardize the distribution name and parameters.
  dist <- tolower(trimws(spec$dist))
  par <- spec$par

  #' Evaluate the distribution-specific density.
  out <- switch(
    dist,
    #' Normal density.
    normal = dnorm(x_grid, mean = par[1L], sd = par[2L]),
    #' Gumbel density.
    gumbel = evd::dgumbel(x_grid, loc = par[1L], scale = par[2L]),
    #' Chi-square density on positive support.
    chi_square = ifelse(x_grid >= 0, dchisq(x_grid, df = par[1L]), 0),
    #' Gamma density on positive support.
    gamma = ifelse(x_grid >= 0, dgamma(x_grid, shape = par[1L], rate = par[2L]), 0),
    #' Reflected Gamma density on negative support.
    reflected_gamma = ifelse(x_grid <= 0, dgamma(-x_grid, shape = par[1L], rate = par[2L]), 0),
    #' Exponential density on positive support.
    exponential = ifelse(x_grid >= 0, dexp(x_grid, rate = par[1L]), 0),
    #' Reflected Exponential density on negative support.
    reflected_exponential = ifelse(x_grid <= 0, dexp(-x_grid, rate = par[1L]), 0),
    #' Weibull density on positive support.
    weibull = ifelse(x_grid >= 0, dweibull(x_grid, shape = par[1L], scale = par[2L]), 0),
    #' Reflected Weibull density on negative support.
    reflected_weibull = ifelse(x_grid <= 0, dweibull(-x_grid, shape = par[1L], scale = par[2L]), 0),
    #' Laplace density.
    laplace = LaplacesDemon::dlaplace(x_grid, location = par[1L], scale = par[2L]),
    #' Beta density on the unit interval.
    beta = ifelse(x_grid >= 0 & x_grid <= 1, dbeta(x_grid, shape1 = par[1L], shape2 = par[2L]), 0),
    #' Uniform density.
    uniform = dunif(x_grid, min = par[1L], max = par[2L]),
    #' Logistic density.
    logistic = dlogis(x_grid, location = par[1L], scale = par[2L]),
    #' Lognormal density on positive support.
    lognormal = ifelse(x_grid > 0, dlnorm(x_grid, meanlog = par[1L], sdlog = par[2L]), 0),
    #' Reflected Lognormal density on negative support.
    reflected_lognormal = ifelse(x_grid < 0, dlnorm(-x_grid, meanlog = par[1L], sdlog = par[2L]), 0),
    #' Student t density.
    t = dt(x_grid, df = par[1L]),
    #' F density on positive support.
    f = ifelse(x_grid >= 0, df(x_grid, df1 = par[1L], df2 = par[2L]), 0),
    #' Pareto density.
    pareto = VGAM::dpareto(x_grid, shape = par[1L], scale = par[2L]),
    #' Cauchy density.
    cauchy = dcauchy(x_grid, location = par[1L], scale = par[2L]),
    #' Two-component contaminated Normal density.
    contaminated = {
      p <- par[1L]
      mu <- par[2L]
      sd1 <- par[3L]
      sd2 <- par[4L]
      #' Weighted component densities.
      p * dnorm(x_grid, mean = mu, sd = sd1) + (1 - p) * dnorm(x_grid, mean = mu, sd = sd2)
    },
    stop("Unsupported distribution in density_value(): ", dist)
  )
  #' Remove undefined or infinite density values.
  out[!is.finite(out)] <- NA_real_
  out
}

#' Choose an x-axis range for a moment-matched density comparison.
#'
#' The range combines:
#' 1. a distribution-specific plotting range for the alternative;
#' 2. mean +/- 4 SD for the matched Normal.
density_xlim <- function(alt_spec, normal_spec) {

  #' Standardize the alternative name and parameters.
  dist <- tolower(trimws(alt_spec$dist))
  par <- alt_spec$par
  #' Matched Normal center and spread.
  normal_mu <- normal_spec$par[1L]
  normal_sd <- normal_spec$par[2L]
  #' Distribution-specific alternative range.
  alt_range <- switch(
    dist,
    #' Gumbel: asymmetric location-scale range.
    gumbel = c(par[1L] - 3 * par[2L], par[1L] + 7 * par[2L]),
    #' Chi-square: positive support through the 99.9th percentile.
    chi_square = c(0, qchisq(0.999, df = par[1L])),
    #' Gamma: positive support through the 99.9th percentile.
    gamma = c(0, qgamma(0.999, shape = par[1L], rate = par[2L])),
    #' Reflected Gamma: mirror the Gamma range around zero.
    reflected_gamma = c(-qgamma(0.999, shape = par[1L], rate = par[2L]), 0),
    #' Exponential: positive support through the 99.9th percentile.
    exponential = c(0, qexp(0.999, rate = par[1L])),
    #' Reflected Exponential: mirror the upper quantile around zero.
    reflected_exponential = c(-qexp(0.999, rate = par[1L]), 0),
    #' Weibull: positive support through the 99.9th percentile.
    weibull = c(0, qweibull(0.999, shape = par[1L], scale = par[2L])),
    #' Reflected Weibull: mirror the Weibull range around zero.
    reflected_weibull = c(-qweibull(0.999, shape = par[1L], scale = par[2L]), 0),
    #' Laplace: symmetric range around the location.
    laplace = c(par[1L] - 8 * par[2L], par[1L] + 8 * par[2L]),
    #' Beta: unit interval with a small plotting margin.
    beta = c(-0.05, 1.05),
    #' Uniform: support limits with a small plotting margin.
    uniform = c(par[1L] - 0.05, par[2L] + 0.05),
    #' Logistic: central 99.8% probability range.
    logistic = qlogis(c(0.001, 0.999), location = par[1L], scale = par[2L]),
    #' Lognormal: positive support through the 99.9th percentile.
    lognormal = c(0, qlnorm(0.999, meanlog = par[1L], sdlog = par[2L])),
    #' Reflected Lognormal: mirror the upper quantile around zero.
    reflected_lognormal = c(-qlnorm(0.999, meanlog = par[1L], sdlog = par[2L]), 0),
    #' Student t: central 99.8% probability range.
    t = qt(c(0.001, 0.999), df = par[1L]),
    #' F: positive support through the 99.9th percentile.
    f = c(0, qf(0.999, df1 = par[1L], df2 = par[2L])),
    #' Pareto: lower support limit through the 99.9th percentile.
    pareto = c(par[2L], VGAM::qpareto(0.999, shape = par[1L], scale = par[2L])),
    #' Cauchy: broad symmetric range around the location.
    cauchy = c(par[1L] - 8 * par[2L], par[1L] + 8 * par[2L]),
    #' Contaminated Normal: range based on the wider component.
    contaminated = c(par[2L] - 4 * par[4L], par[2L] + 4 * par[4L]),
    stop("Unsupported distribution in density_xlim(): ", dist)
  )

  #' Matched Normal range.
  normal_range <- c(normal_mu - 4 * normal_sd, normal_mu + 4 * normal_sd)
  #' Combine both ranges.
  xlim <- range(c(alt_range, normal_range), finite = TRUE)
  #' Require finite, increasing limits.
  if (length(xlim) != 2L || !all(is.finite(xlim)) || xlim[1L] >= xlim[2L]) {
    stop("Could not construct a valid x-axis range for ", dist)
  }
  xlim
}

#' Create a readable distribution label.
dist_label <- function(spec, digits = 3) {
  dist_name <- tools::toTitleCase(gsub("_", " ", spec$dist))
  par_text  <- paste(round(spec$par, digits), collapse = ", ")
  paste0(dist_name, "(", par_text, ")")
}



#' ============================================================
#' SECTION 36: DENSITY PLOT FUNCTION
#' ============================================================

#' Plot alternative distributions against their moment-matched Normals.
#'
#' Each panel compares one alternative distribution to its own
#' moment-matched Normal distribution.
#'
#' Input:
#' paired_specs from make_paired_specs().
#'
#' Output:
#' A saved PDF file.
plot_density_pairs <- function(paired_specs, output_dir = "results",
                               file_name = "density_pairs.pdf", n_grid = 2000,
                               n_col = 4, alt_col = "brown", normal_col = "steelblue") {

  if (length(paired_specs) == 0L) {
    stop("paired_specs is empty.")
  }

  if (n_grid < 100L) {
    stop("n_grid should be at least 100.")
  }

  make_output_dir(output_dir)

  n_pairs <- length(paired_specs)
  n_col   <- min(n_col, n_pairs)
  n_row   <- ceiling(n_pairs / n_col)

  pdf_file <- file.path(output_dir, file_name)
  pdf(file   = pdf_file, width  = 4.2 * n_col, height = 4.0 * n_row)
  old_par <- par(no.readonly = TRUE)

  on.exit({par(old_par)
    dev.off()
  }, add = TRUE)

  #' Continue the calculation.
  par(
    mfrow = c(n_row, n_col),
    mar   = c(4.2, 4.2, 3.0, 1.2),
    oma   = c(1.0, 1.0, 3.0, 0.5),
    mgp   = c(2.5, 0.7, 0)
  )

  for (pair in paired_specs) {
    normal_spec <- pair$normal
    alt_spec    <- pair$alt
    xlim <- density_xlim( alt_spec = alt_spec, normal_spec = normal_spec)
    x_grid <- seq(from = xlim[1L], to = xlim[2L], length.out = n_grid)

    d_alt <- density_value(x_grid = x_grid,spec   = alt_spec)
    d_norm <- density_value(x_grid = x_grid, spec   = normal_spec)
    ymax <- max(c(d_alt, d_norm), na.rm = TRUE)

    if (!is.finite(ymax) || ymax <= 0) {
      ymax <- 1
    }

    #' Initialize the density panel.
    plot(
      NA,
      xlim = xlim,
      ylim = c(0, 1.10 * ymax),
      xlab = "x",
      ylab = "Density",
      main = paste0(tools::toTitleCase(gsub("_", " ", alt_spec$dist))," vs Matched Normal"),
      las      = 1,
      cex.main = 1.2
    )

    #' Shade the alternative density.
    polygon(
      x      = c(x_grid, rev(x_grid)),
      y      = c(d_alt, rep(0, length(x_grid))),
      col    = adjustcolor(alt_col, alpha.f = 0.18),
      border = NA
    )

    #' Shade the matched Normal density.
    polygon(
      x      = c(x_grid, rev(x_grid)),
      y      = c(d_norm, rep(0, length(x_grid))),
      col    = adjustcolor(normal_col, alpha.f = 0.18),
      border = NA
    )
    
    #' Alternative density curve.
    lines(x_grid, d_alt, col = alt_col, lwd = 3)
    #' Matched Normal density curve.
    lines(x_grid, d_norm, col = normal_col, lwd = 3, lty = 2)
    #' Matched Normal mean.
    abline(v = normal_spec$par[1L], col = "gray50", lty = 3, lwd = 1.5)
    #' Distribution legend.
    legend("topright", legend = c(dist_label(alt_spec), dist_label(normal_spec)),
           col = c(alt_col, normal_col), lwd = 3, lty = c(1, 2), bty = "o", cex = 0.95)
  }
  #' Overall plot title.
  mtext("Alternative Distributions vs Moment-Matched Normal Distributions",
        outer = TRUE, cex = 1.2, font = 3, line = 1.0)
  cat("Saved density plots:", pdf_file, "\n")
  invisible(pdf_file)
}


#' ============================================================
#' PART 8: MAIN RUN FUNCTION
#' ============================================================


#' ============================================================
#' SECTION 37: SMALL PIPELINE HELPERS
#' ============================================================

print_pipeline_step <- function(title) {
  cat("\n")
  cat(strrep("=", 60), "\n")
  cat(title, "\n")
  cat(strrep("=", 60), "\n")
}


#' ------------------------------------------------------------
#' Save an available result table to a CSV file.
#' ------------------------------------------------------------
save_csv_if_not_null <- function(x, file) {
  if (!is.null(x)) {
    write.csv(x, file = file, row.names = FALSE)
    cat("Saved:", file, "\n")
  }
  invisible(file)
}


#' ------------------------------------------------------------
#' Check whether a model bundle contains held-out feature data.
#' ------------------------------------------------------------
has_holdout_features <- function(eval_results_by_n,
                                 sample_sizes) {
  if (is.null(eval_results_by_n)) {
    return(FALSE)
  }
  for (n in sample_sizes) {
    n_key <- as.character(n)
    if (!n_key %in% names(eval_results_by_n)) {
      return(FALSE)
    }
    if (is.null(eval_results_by_n[[n_key]]$holdout_features)) {
      return(FALSE)
    }
  }

  TRUE
}


#' ============================================================
#' COMBINE EVALUATION SUMMARIES
#' ============================================================
combine_eval_summaries <- function(eval_results_by_n) {
  #' Combine sample-size evaluation summaries into one table.
  if (is.null(eval_results_by_n) || length(eval_results_by_n) == 0L) {
    return(NULL)
  }
  out <- lapply(names(eval_results_by_n), function(n_key) {
    summary_n <- eval_results_by_n[[n_key]]$summary
    if (is.null(summary_n)) {
      return(NULL)
    }
    summary_n$n <- as.numeric(n_key)
    #' Continue the calculation.
    summary_n[, c("n", setdiff(names(summary_n), "n")), drop = FALSE]
  })
  out <- Filter(Negate(is.null), out)
  if (length(out) == 0L) {
    return(NULL)
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}


#' ============================================================
#' SECTION 38: MAIN PIPELINE WRAPPER
#' ============================================================

#' Run the ML-based normality testing framework.
#'
#' This wrapper allows the user to:
#' 1. train or load models;
#' 2. generate or load held-out validation/evaluation results;
#' 3. save ML ROC plots;
#' 4. run permutation AUC VIP;
#' 5. select top VIP features;
#' 6. plot VIP for selected models only;
#' 7. optionally create density plots;
#' 8. optionally run classical-vs-ML ROC comparisons.
run_normality_framework <- function(
  #' Core simulation settings
  sample_sizes = c(10, 50),
  num_sim      = 500,
  n_iter_eval  = 1000,

  #' Parallel processing
  #' These defaults work on macOS, Linux, and Windows.
  #' They reserve one core for the operating system and also avoid
  #' requesting more workers than useful for sample-size-level tasks.
  n_cores_train = default_parallel_cores(max_tasks = length(sample_sizes)),
  n_cores_eval  = default_parallel_cores(max_tasks = length(sample_sizes)),
  n_cores_vip   = default_parallel_cores(),
  n_cores_roc   = default_parallel_cores(max_tasks = length(sample_sizes)),
  n_cores_cv    = default_parallel_cores(),

  #' ------------------------------------------------------------
  #' Model settings
  #' ------------------------------------------------------------
  models_to_train = c("LR", "RF", "ANN", "GBM", "SVM", "KNN"),
  #' Fixed probability threshold
  threshold          = 0.50,
  majority_threshold = 0.40,

  #' ------------------------------------------------------------
  #' Distribution settings
  #' ------------------------------------------------------------
  #' eval_alt_specs supplies held-out data for validation or final evaluation.
  train_alt_specs = default_training_specs(),
  eval_alt_specs  = default_eval_specs(),
  roc_alt_specs   = default_eval_specs(),
  #' ------------------------------------------------------------
  #' Feature settings
  #' ------------------------------------------------------------
  feature_set        = NULL,
  standardize_sample = FALSE,
  center_by          = NULL,
  #' ------------------------------------------------------------
  #' Model training/loading options
  #' ------------------------------------------------------------
  train_models    = TRUE,
  model_save_path = NULL,
  k_folds         = 5L,
  cv_repeats     = 1L,
  cv_verbose     = FALSE,
  #' ------------------------------------------------------------
  #' Held-out validation/evaluation options
  #' ------------------------------------------------------------
  run_evaluation       = TRUE,
  eval_save_path       = NULL,
  load_eval_if_exists  = FALSE,
  save_eval_results    = TRUE,

  #' Store holdout feature matrix during evaluation.
  #' This is required for permutation AUC VIP.
  store_holdout_features = NULL,
  #' ------------------------------------------------------------
  #' Optional analysis steps
  #' ------------------------------------------------------------
  make_density_plots  = FALSE,
  run_ml_roc_plots    = TRUE,
  run_permutation_vip = FALSE,
  run_roc_analysis    = FALSE,
  reuse_roc_evaluation_predictions = TRUE,
  show_progress = TRUE,
  #' ------------------------------------------------------------
  #' Permutation VIP settings
  #' ------------------------------------------------------------
  n_permutations = 20,
  #' These models are used for feature selection.
  models_for_vip = c("RF", "GBM", "ANN"),
  vip_top_n = 30,
  vip_rank_by = c("mean", "median", "min"),
  vip_min_pct_positive = 0.50,
  min_mean_for_negative_min = 0.005,
  vip_require_positive_mean = TRUE,

  #' These models are used only for plotting VIP.
  vip_plot_models = c("RF"),
  vip_plot_sample_sizes = sample_sizes,
  #' ------------------------------------------------------------
  #' ROC comparison settings
  #' ------------------------------------------------------------
  classical_tests = c("SW", "AD", "JB"),
  ml_methods_roc  = c("RF", "GBM", "ANN"),
  alpha_grid      = seq(0, 1, by = 0.05),
  n_sim_roc       = 1000,

  #' ------------------------------------------------------------
  #' Plot colors
  #' ------------------------------------------------------------
  test_colors = c(SW = "#CC79A7", AD = "#FF7F00", JB = "#56B4E9"),
  ml_colors = c(LR = "black", RF = "red", GBM = "blue", ANN = "forestgreen", SVM = "green", KNN = "purple", MajorityVote = "brown"),
  #' ------------------------------------------------------------
  #' Output
  #' ------------------------------------------------------------
  output_dir = "results"
) {
  #' ------------------------------------------------------------
  #' Initial setup
  #' ------------------------------------------------------------
  make_output_dir(output_dir)
  vip_rank_by <- match.arg(vip_rank_by)
  if (is.null(model_save_path)) {
    model_save_path <- file.path(output_dir, "trained_models.RData")
  }
  if (is.null(eval_save_path)) {
    eval_save_path <- file.path(output_dir, "eval_results.RData")
  }
  if (is.null(store_holdout_features)) {
    store_holdout_features <- isTRUE(run_permutation_vip)
  }
  sample_sizes <- as.numeric(sample_sizes)
  train_pairs <- make_paired_specs(train_alt_specs)
  eval_specs  <- make_eval_specs(eval_alt_specs)

  #' Calculate run settings.
  run_settings <- list(
    sample_sizes              = sample_sizes,
    num_sim                   = num_sim,
    n_iter_eval               = n_iter_eval,
    models_to_train           = models_to_train,
    threshold                 = threshold,
    majority_threshold        = majority_threshold,
    train_alt_specs           = train_alt_specs,
    eval_alt_specs            = eval_alt_specs,
    roc_alt_specs             = roc_alt_specs,
    feature_set               = feature_set,
    standardize_sample        = standardize_sample,
    center_by                 = center_by,
    k_folds                   = k_folds,
    cv_repeats                 = cv_repeats,
    cv_verbose                 = cv_verbose,
    n_permutations            = n_permutations,
    models_for_vip            = models_for_vip,
    vip_top_n                 = vip_top_n,
    vip_rank_by               = vip_rank_by,
    vip_min_pct_positive      = vip_min_pct_positive,
    vip_require_positive_mean = vip_require_positive_mean,
    vip_plot_models           = vip_plot_models,
    vip_plot_sample_sizes     = vip_plot_sample_sizes,
    classical_tests           = classical_tests,
    ml_methods_roc            = ml_methods_roc,
    alpha_grid                = alpha_grid,
    n_sim_roc                 = n_sim_roc,
    n_cores_train             = n_cores_train,
    n_cores_eval              = n_cores_eval,
    n_cores_vip               = n_cores_vip,
    n_cores_roc               = n_cores_roc,
    n_cores_cv                = n_cores_cv,
    parallel_backend_train    = get_parallel_backend(n_cores_train),
    parallel_backend_eval     = get_parallel_backend(n_cores_eval),
    parallel_backend_vip      = get_parallel_backend(n_cores_vip),
    parallel_backend_roc      = get_parallel_backend(n_cores_roc),
    parallel_backend_cv       = if (n_cores_cv > 1L) "PSOCK" else "sequential",
    reuse_roc_evaluation_predictions = reuse_roc_evaluation_predictions,
    show_progress             = show_progress,
    output_dir                = output_dir
  )

  #' Calculate density plot file.
  density_plot_file         <- NULL
  trained_models            <- NULL
  eval_results_by_n         <- NULL
  eval_summary_all          <- NULL
  vip_results               <- NULL
  vip_summary               <- NULL
  top_features              <- NULL
  selected_feature_metrics  <- NULL
  final_feature_table       <- NULL
  final_selected_features   <- NULL
  roc_results               <- NULL

  #' ------------------------------------------------------------
  #' Optional density plots
  #' ------------------------------------------------------------
  if (isTRUE(make_density_plots)) {
    print_pipeline_step("DENSITY PLOTS")
    #' Calculate density plot file.
    density_plot_file <- plot_density_pairs(
      paired_specs = make_paired_specs(eval_alt_specs),
      output_dir   = output_dir,
      file_name    = "density_pairs.pdf"
    )
  }

  #' ------------------------------------------------------------
  #' Step 1: Train or load models
  #' ------------------------------------------------------------
  print_pipeline_step("STEP 1: MODEL TRAINING / MODEL LOADING")

  if (isTRUE(train_models)) {
    #' Calculate trained models.
    trained_models <- train_models_all_n(
      sample_sizes       = sample_sizes,
      paired_specs       = train_pairs,
      num_sim            = num_sim,
      models_to_train    = models_to_train,
      feature_set        = feature_set,
      standardize_sample = standardize_sample,
      center_by          = center_by,
      k_folds            = k_folds,
      cv_repeats        = cv_repeats,
      cv_verbose        = cv_verbose,
      n_cores            = n_cores_train,
      n_cores_cv         = n_cores_cv,
      show_progress      = show_progress
    )

    save(
      trained_models,
      run_settings,
      file = model_save_path
    )
    
    cat("\nSaved trained models to:", model_save_path, "\n")
  } else {

    #' Validate the condition.
    if (!file.exists(model_save_path)) {
      stop("Model file not found: ", model_save_path,
        "\nSet train_models = TRUE to train models first."
      )
    }
    
    load(model_save_path)
    if (!exists("trained_models")) {
      stop("Loaded model file does not contain an object named 'trained_models'.")
    }
    cat("Loaded trained models from:", model_save_path, "\n")
  }

  #' Restrict to requested sample sizes if the loaded model file has more.
  model_keys <- as.character(sample_sizes)
  missing_models <- setdiff(model_keys, names(trained_models))
  #' Validate the condition.
  if (length(missing_models) > 0L) {
    stop("The following sample sizes are missing from trained_models: ",
      paste(missing_models, collapse = ", ")
    )
  }

  trained_models <- trained_models[model_keys]
  #' ------------------------------------------------------------
  #' Step 2: Evaluate or load evaluation results
  #' ------------------------------------------------------------
  print_pipeline_step("STEP 2: HELD-OUT EVALUATION")
  if (isTRUE(load_eval_if_exists) && file.exists(eval_save_path)) {
    load(eval_save_path)
    if (!exists("eval_results_by_n")) {
      stop("Loaded evaluation file does not contain 'eval_results_by_n'.")
    }
    cat("Loaded evaluation results from:", eval_save_path, "\n")
  } else if (isTRUE(run_evaluation)) {

    #' Calculate eval results by n.
    eval_results_by_n <- evaluate_models_all_n(
      trained_models         = trained_models,
      eval_specs             = eval_specs,
      n_iter                 = n_iter_eval,
      threshold              = threshold,
      majority_threshold     = majority_threshold,
      standardize_sample     = standardize_sample,
      center_by              = center_by,
      store_holdout_features = store_holdout_features,
      n_cores                = n_cores_eval,
      show_progress          = show_progress
    )
    if (isTRUE(save_eval_results)) {
      save(
        eval_results_by_n,
        run_settings,
        file = eval_save_path
      )

      cat("\nSaved evaluation results to:", eval_save_path, "\n")
    }

  } else {

    cat("Skipping evaluation because run_evaluation = FALSE.\n")
  }

  #' ------------------------------------------------------------
  #' Evaluation metrics summary: print and save
  #' ------------------------------------------------------------

  if (!is.null(eval_results_by_n)) {

    eval_summary_all <- combine_eval_summaries(eval_results_by_n)

    if (!is.null(eval_summary_all)) {

      cat("\n")
      cat(strrep("=", 70), "\n")
      cat("COMBINED EVALUATION METRICS SUMMARY\n")
      cat(strrep("=", 70), "\n")
      print(eval_summary_all)

      eval_summary_path <- file.path(
        output_dir,
        "evaluation_metrics_summary.csv"
      )

      write.csv(
        eval_summary_all,
        eval_summary_path,
        row.names = FALSE
      )

      cat("\nSaved evaluation metrics summary to:", eval_summary_path, "\n")

      #' Save the evaluation object again, now including the combined summary.
      if (isTRUE(save_eval_results)) {

        #' Save the results.
        save(
          eval_results_by_n,
          eval_summary_all,
          run_settings,
          file = eval_save_path
        )

        cat("Updated evaluation RData with eval_summary_all:", eval_save_path, "\n")
      }
    }
  }

  #' ------------------------------------------------------------
  #' Step 3: ML-only ROC plots
  #' ------------------------------------------------------------

  if (isTRUE(run_ml_roc_plots)) {

    print_pipeline_step("STEP 3: ML ROC PLOTS")

    if (is.null(eval_results_by_n)) {
      warning("Skipping ML ROC plots because eval_results_by_n is NULL.")
    } else {

      #' Continue the calculation.
      save_ml_roc_plots(
        eval_results_by_n = eval_results_by_n,
        sample_sizes      = sample_sizes,
        output_dir        = output_dir,
        ml_colors         = ml_colors
      )
    }

  } else {

    cat("\nSkipping ML ROC plots because run_ml_roc_plots = FALSE.\n")
  }

  #' ------------------------------------------------------------
  #' Step 4: Permutation AUC VIP
  #' ------------------------------------------------------------

  if (isTRUE(run_permutation_vip)) {

    print_pipeline_step("STEP 4: PERMUTATION AUC VARIABLE IMPORTANCE")

    #' VIP requires holdout_features. If the current evaluation object
    #' does not have them, rerun evaluation only to store the feature matrix.
    if (!has_holdout_features(eval_results_by_n, sample_sizes)) {

      cat(
        "Evaluation results do not contain holdout_features.\n",
        "Rerunning evaluation only to store holdout features for VIP.\n",
        sep = ""
      )

      #' Calculate eval results by n.
      eval_results_by_n <- evaluate_models_all_n(
        trained_models         = trained_models,
        eval_specs             = eval_specs,
        n_iter                 = n_iter_eval,
        threshold              = threshold,
        majority_threshold     = majority_threshold,
        standardize_sample     = standardize_sample,
        center_by              = center_by,
        store_holdout_features = TRUE,
        n_cores                = n_cores_eval,
        show_progress          = show_progress
      )

      save(
        eval_results_by_n,
        run_settings,
        file = eval_save_path
      )

      cat("\nSaved VIP-ready evaluation results to:", eval_save_path, "\n")
    }

    #' Calculate vip results.
    vip_results <- compute_permutation_vip(
      trained_models     = trained_models,
      eval_results_by_n  = eval_results_by_n,
      n_permutations     = n_permutations,
      models_for_vip     = models_for_vip,
      n_cores            = n_cores_vip
    )

    vip_summary <- summarize_permutation_vip(vip_results)

    #' Calculate top features.
    top_features <- select_top_vip_features(
      vip_summary           = vip_summary,
      top_n                 = vip_top_n,
      rank_by               = vip_rank_by,
      min_pct_positive      = vip_min_pct_positive,
      require_positive_mean = vip_require_positive_mean
    )

    #' Calculate final feature table.
    final_feature_table <- final_feature_filter(
      vip_summary,
      min_pct_positive           = vip_min_pct_positive,
      min_median_drop            = 0,
      allow_small_negative_min   = TRUE,
      min_mean_for_negative_min  = min_mean_for_negative_min,
      top_n                      = vip_top_n
    )

    #' ------------------------------------------------------------
    #' Final selected features
    #' ------------------------------------------------------------
    #' Post-VIP correlation filtering is intentionally skipped.
    #' The final selected feature set is exactly the set retained by
    #' final_feature_filter().

    final_selected_features <- final_feature_table$Feature


    #' Keep exported metrics consistent with the final feature set.
    selected_feature_metrics <- get_selected_feature_metrics(
      vip_summary = vip_summary,
      selected_features = final_selected_features
    )

    rownames(final_feature_table) <- NULL

    excluded_after_final_filter <- setdiff(top_features, final_selected_features)

    if (length(excluded_after_final_filter) > 0L) {
      cat("\nTop VIP features excluded by final_feature_filter():\n")
      print(excluded_after_final_filter)
    }

    write.csv(
      final_feature_table,
      file.path(output_dir, "final_feature_table.csv"),
      row.names = FALSE
    )

    write.csv(
      data.frame(Feature = final_selected_features),
      file.path(output_dir, "final_selected_features.csv"),
      row.names = FALSE
    )

    cat("\nSelected features from permutation AUC VIP:\n")
    print(selected_feature_metrics)

    cat("\nFinal selected features after final_feature_filter():\n")
    print(final_feature_table)

    selected_metrics_path <- file.path(
      output_dir,
      "selected_feature_metrics.csv"
    )

    save_csv_if_not_null(
      selected_feature_metrics,
      selected_metrics_path
    )

    vip_summary_path <- file.path(
      output_dir,
      "vip_summary.csv"
    )

    save_csv_if_not_null(
      vip_summary,
      vip_summary_path
    )

    vip_save_path <- file.path(
      output_dir,
      "permutation_vip_results.RData"
    )

    #' Save the results.
    save(
      eval_results_by_n,
      vip_results,
      vip_summary,
      top_features,
      selected_feature_metrics,
      final_feature_table,
      final_selected_features,
      run_settings,
      file = vip_save_path
    )

    cat("\nSaved permutation VIP results to:", vip_save_path, "\n")

    #' Continue the calculation.
    save_permutation_vip_plots(
      vip_results = vip_results,
      sample_sizes = vip_plot_sample_sizes,
      models = vip_plot_models,
      top_n = vip_top_n,
      selected_features = final_selected_features,
      output_dir = output_dir
    )

    #' Draw the plot.
    save_permutation_vip_1x2_plot(
      vip_results = vip_results,
      sample_sizes = c(10, 50),
      model = "RF",
      top_n = vip_top_n,
      selected_features = final_selected_features,
      output_dir = output_dir
    )

  } else {

    cat("\nSkipping permutation VIP because run_permutation_vip = FALSE.\n")
  }

  #' ------------------------------------------------------------
  #' Step 5: Classical-vs-ML ROC comparison
  #' ------------------------------------------------------------

  if (isTRUE(run_roc_analysis)) {

    print_pipeline_step("STEP 5: CLASSICAL VS ML ROC COMPARISON")

    #' Calculate roc results.
    roc_results <- run_roc_comparison(
      trained_models  = trained_models,
      roc_alt_specs   = roc_alt_specs,
      sample_sizes    = sample_sizes,
      classical_tests = classical_tests,
      ml_methods      = ml_methods_roc,
      alpha_grid      = alpha_grid,
      n_sim           = n_sim_roc,
      output_dir      = output_dir,
      test_colors     = test_colors,
      ml_colors       = ml_colors,
      n_cores         = n_cores_roc,
      eval_results_by_n = eval_results_by_n,
      reuse_evaluation_predictions = reuse_roc_evaluation_predictions,
      show_progress   = show_progress
    )

    roc_save_path <- file.path(output_dir, "roc_comparison_results.RData")

    save(
      roc_results,
      run_settings,
      file = roc_save_path
    )

    cat("\nSaved ROC comparison results to:", roc_save_path, "\n")
  } else {
    cat("\nSkipping classical-vs-ML ROC comparison because run_roc_analysis = FALSE.\n")
  }

  #' ------------------------------------------------------------
  #' Completion
  #' ------------------------------------------------------------

  print_pipeline_step("PIPELINE COMPLETE")
  cat("Output directory:", output_dir, "\n")

  #' Continue the calculation.
  invisible(
    list(
      density_plot_file         = density_plot_file,
      trained_models            = trained_models,
      eval_results_by_n         = eval_results_by_n,
      eval_summary_all          = eval_summary_all,
      vip_results               = vip_results,
      vip_summary               = vip_summary,
      top_features              = top_features,
      selected_feature_metrics  = selected_feature_metrics,
      final_feature_table       = final_feature_table,
      final_selected_features   = final_selected_features,
      roc_results               = roc_results,
      run_settings              = run_settings
    )
  )
}
