#' =============================================================================
#' RUN SELECTED USER-FRAMEWORK ANALYSES
#' =============================================================================
#'
#' This script provides one clearly marked run-control panel for selecting:
#' 1. the downstream comparison;
#' 2. the normality-pretest approach;
#' 3. the analytical phases;
#' 4. the threshold mode; and
#' 5. the simulation scale.
#'
#' Statistical functions remain in the separate user-framework function file.
#' Only Step 1 normally needs to be edited before a run.
#' =============================================================================


#' =============================================================================
#' SUPPORTING FILE AND EXECUTION HELPERS
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Return the project directory supplied through --wd= or use the current folder.
#' -----------------------------------------------------------------------------
resolve_project_directory <- function() {
  command_arguments <- commandArgs(trailingOnly = TRUE)
  working_directory_argument <- grep(
    pattern = "^--wd=",
    x = command_arguments,
    value = TRUE
  )

  project_directory <- if (length(working_directory_argument) == 1L) {
    sub("^--wd=", "", working_directory_argument)
  } else {
    getwd()
  }

  #' Prepare the output location.
  project_directory <- normalizePath(
    project_directory,
    winslash = "/",
    mustWork = TRUE
  )

  project_directory
}


#' -----------------------------------------------------------------------------
#' Set the random-number generator and seed used by one analysis run.
#' -----------------------------------------------------------------------------
set_reproducible_rng <- function(seed,
                                 rng_version = "4.2.0") {
  RNGversion(rng_version)

  set.seed(
    seed,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )

  invisible(seed)
}


#' -----------------------------------------------------------------------------
#' Create a directory and return its normalized path.
#' -----------------------------------------------------------------------------
make_run_directory <- function(directory) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!dir.exists(directory)) {
    stop("Could not create directory: ", directory, call. = FALSE)
  }

  normalizePath(
    directory,
    winslash = "/",
    mustWork = TRUE
  )
}


#' -----------------------------------------------------------------------------
#' Open a file browser and return one validated file path.
#' -----------------------------------------------------------------------------
browse_for_required_file <- function(file_description,
                                     allowed_extensions) {
  if (!interactive()) {
    stop(
      "Selecting required files needs an interactive R session.",
      call. = FALSE
    )
  }

  #' Tell the user which file is required.
  message(
    "\nBrowse for ",
    file_description,
    "."
  )

  #' Open the system file browser.
  selected_file <- tryCatch(
    file.choose(new = FALSE),
    error = function(error_condition) {
      stop(
        "No file was selected for ",
        file_description,
        ".",
        call. = FALSE
      )
    }
  )

  #' Confirm that the selected file exists.
  if (!file.exists(selected_file)) {
    stop(
      "The selected file does not exist: ",
      selected_file,
      call. = FALSE
    )
  }

  #' Standardize the allowed extensions.
  allowed_extensions <- tolower(
    sub("^\\.", "", as.character(allowed_extensions))
  )

  #' Read the selected extension.
  selected_extension <- tolower(
    tools::file_ext(selected_file)
  )

  #' Reject an incorrect file type.
  if (!selected_extension %in% allowed_extensions) {
    stop(
      "Invalid file type for ",
      file_description,
      ". Expected: ",
      paste0(".", allowed_extensions, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  #' Return a normalized path.
  normalizePath(
    selected_file,
    winslash = "/",
    mustWork = TRUE
  )
}


#' -----------------------------------------------------------------------------
#' Confirm that a selected R file defines the required functions.
#' -----------------------------------------------------------------------------
validate_selected_r_file <- function(file,
                                     required_functions,
                                     file_description) {
  #' Read the source code without executing it.
  source_lines <- readLines(
    file,
    warn = FALSE
  )

  #' Check for each required function definition.
  function_is_defined <- vapply(
    required_functions,
    function(function_name) {
      definition_pattern <- paste0(
        "^[[:space:]]*",
        function_name,
        "[[:space:]]*<-[[:space:]]*function[[:space:]]*\\("
      )

      any(grepl(definition_pattern, source_lines))
    },
    logical(1)
  )

  #' Identify missing definitions.
  missing_functions <- required_functions[!function_is_defined]

  if (length(missing_functions) > 0L) {
    stop(
      "The selected ",
      file_description,
      " is missing function definition(s): ",
      paste(missing_functions, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}



#' =============================================================================
#' STEP 1: USER RUN-CONTROL PANEL
#'
#' EDIT THIS STEP TO SELECT WHAT THE SCRIPT SHOULD RUN.
#' The remaining steps build and execute the requested analyses automatically.
#' =============================================================================


#' -----------------------------------------------------------------------------
#' 1A. WORKFLOW AND OUTPUT CONTROLS
#' -----------------------------------------------------------------------------

#' TRUE executes the run plan. FALSE creates and prints the plan only.
run_workflow <- TRUE

#' TRUE previews the run plan without starting any simulation.
dry_run <- FALSE

#' TRUE stops after the first failed analysis.
#' FALSE records the failure and continues with the remaining analyses.
stop_on_error <- TRUE

#' TRUE retains completed result objects in memory.
#' Set to FALSE for large runs when only the saved RData files are needed.
keep_results_in_memory <- TRUE

#' TRUE saves a lightweight RDS and CSV manifest of all requested analyses.
write_run_manifest <- TRUE

#' Install missing packages only when explicitly enabled.
install_missing_packages <- FALSE

#' Save the structured result object produced by each analysis.
save_results <- TRUE

#' Available compression values: FALSE, TRUE, "gzip", "bzip2", or "xz".
save_compression <- "gzip"


#' -----------------------------------------------------------------------------
#' 1B. SIMULATION SCALE AND REPRODUCIBILITY
#' -----------------------------------------------------------------------------

#' Select "pilot", "full", or "custom".
run_mode <- "custom"

#' A different reproducible seed is derived for each requested analysis.
base_seed <- 12345L
rng_version <- "4.2.0"

#' Set simulation and resampling counts for each run mode.
run_settings <- switch(
  run_mode,

  pilot = list(
    Nsim = 200L,
    N_tradeoff = 200L,
    nresample = 499L
  ),

  full = list(
    Nsim = 10000L,
    N_tradeoff = 10000L,
    nresample = 1000L
  ),

  custom = list(
    Nsim = 10000L,
    N_tradeoff = 10000L,
    nresample = 1000L
  ),

  stop("run_mode must be 'pilot', 'full', or 'custom'.", call. = FALSE)
)

#' Sample sizes used in the across-sample-size analysis.
sample_sizes <- c(10, 20, 30, 40, 50)

#' Sample size used by phases that require one fixed n.
single_n <- 10

#' Effect-size grid used in Phase 5.
effect_sizes_plot <- c(0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0)

#' Threshold grids for classical and ML normality procedures.
threshold_grid_classical <- seq(0.00, 1.00, by = 0.005)
threshold_grid_ml <- seq(0.00, 1.00, by = 0.005)

#' ROC grids for classical and ML normality procedures.
normality_roc_grid_classical <- seq(0.00, 1.00, by = 0.025)
normality_roc_grid_ml <- seq(0.00, 1.00, by = 0.025)

#' Downstream significance-level grid used in Phase 3.
significance_levels <- seq(0.00, 1.00, by = 0.005)

#' Nominal downstream test level.
test_alpha <- 0.05

#' Threshold-selection tolerances used in Phase 2.
type1_inflation_tolerance <- 0.005
power_loss_tolerance <- 0.01

#' Internal downstream method codes retained throughout the framework.
downstream_methods <- c("test_1", "test_2", "adaptive")

#' Optional plot controls. Leave NULL for automatic behavior.
type1_ylim <- NULL
reference_effect <- NULL


#' -----------------------------------------------------------------------------
#' 1C. SELECT THE DOWNSTREAM COMPARISON
#' -----------------------------------------------------------------------------

#' Available downstream-analysis identifiers:
#'
#' ONE-SAMPLE
#' "onesample_t_vs_sign"
#' "onesample_t_vs_permutation"
#' "onesample_t_vs_boxcox"
#'
#' TWO-SAMPLE
#' "twosample_t_vs_mann_whitney"
#' "twosample_t_vs_permutation"
#'
#' ONE-WAY ANOVA
#' "anova_vs_kruskal_wallis"
#' "anova_vs_permutation"
#'
#' SIMPLE LINEAR REGRESSION
#' "regression_ols_vs_rank"
#' "regression_ols_vs_permutation"
#' "regression_ols_vs_mm"
#' "regression_ols_vs_m"

#' Enter one or more identifiers from the list above.
selected_downstream_analyses <- c(
  "onesample_t_vs_sign"
  #' "anova_vs_kruskal_wallis"
)


#' -----------------------------------------------------------------------------
#' 1D. SELECT THE ANALYTICAL PHASES
#' -----------------------------------------------------------------------------

#' Phase guide:
#' 1 = normality-test ROC analysis
#' 2 = threshold trade-off analysis
#' 3 = power-versus-Type-I-error ROC analysis
#' 4 = power and Type I error across sample sizes
#' 5 = power across effect sizes
#'
#' Phases 3, 4, and 5 require a threshold.
#' Include Phase 2 to estimate it, or supply pretest_threshold below.

#' Run classical normality-pretest analyses.
run_classical_analyses <- TRUE

#' Classical normality tests used for adaptive routing.
classical_methods <- c("SW")

#' Analytical phases run for each selected classical method.
classical_phases <- 1:5

#' TRUE calculates the common multi-test classical ROC plot once per design.
run_classical_roc_once <- TRUE

#' Run the Fisher-combined Shapiro-Wilk and Anderson-Darling procedure.
run_fisher_analyses <- FALSE

#' Analytical phases run for the Fisher-combined procedure.
fisher_phases <- 1:5


#' -----------------------------------------------------------------------------
#' 1E. SELECT THE THRESHOLD MODE
#' -----------------------------------------------------------------------------

#' Allowed values:
#' "single" = estimate one threshold at single_n;
#' "per_n"  = estimate a separate threshold for each value in sample_sizes.
#'
#' Use c("single", "per_n") to run both modes.
threshold_modes_to_run <- c("per_n")

#' Supply this only when Phase 2 is skipped but Phase 3, 4, or 5 is requested.
#'
#' One common threshold:
#' pretest_threshold <- 0.05
#'
#' Sample-size-specific thresholds:
#' pretest_threshold <- c(
#'   "10" = 0.05,
#'   "20" = 0.05,
#'   "30" = 0.05,
#'   "40" = 0.05,
#'   "50" = 0.05
#' )
pretest_threshold <- NULL


#' -----------------------------------------------------------------------------
#' 1F. NORMALITY ROC DISPLAY SETTINGS
#' -----------------------------------------------------------------------------

#' Classical tests calculated in the Phase 1 ROC simulation.
normality_tests_for_roc <- c(
  "SW", "SF", "LF", "KS", "JB",
  "SKEW", "KURT", "DAP", "AD", "CVM"
)

#' Classical tests displayed in the final Phase 1 ROC plot.
selected_normality_tests_for_roc <- c(
  "SW", "SF", "LF", "JB", "SKEW", "DAP", "AD", "CVM"
)


#' -----------------------------------------------------------------------------
#' 1G. OPTIONAL ML NORMALITY APPLICATION
#' -----------------------------------------------------------------------------

#' TRUE compares the selected fitted ML models in Phase 1.
run_ml_roc <- TRUE
ml_roc_phases <- 1L

#' TRUE uses one selected ML model for adaptive routing.
run_ml_adaptive <- TRUE
ml_adaptive_phases <- 2:5

#' Select the downstream comparisons that should use the ML normality procedure.
#'
#' These identifiers must exist in the downstream-analysis catalog. The ML
#' classifier is applied to the numeric object returned by each design's
#' fn_to_get_norm_obj function, including raw samples or model residuals.
ml_downstream_analyses <- c(
  "onesample_t_vs_sign"
  #'"anova_vs_kruskal_wallis"
)

#' To apply ML to every downstream comparison selected in Section 1C, use:
#' ml_downstream_analyses <- selected_downstream_analyses

#' Model used for adaptive routing.
ml_chosen_model <- "SVM"

#' Models shown in each ML ROC comparison.
ml_models_for_roc <- c("RF", "GBM", "ANN", "SVM")

#' Include the mean model probability in the ML ROC comparison.
ml_use_majority_vote <- TRUE

#' These paths remain NULL because the script opens file browsers when ML is used.
ml_framework_file <- NULL
ml_model_file <- NULL

#' Reusable flag for validation, package checks, and run-plan construction.
ml_requested <- isTRUE(run_ml_roc) ||
  isTRUE(run_ml_adaptive)

ml_resources <- NULL
ml_configs <- NULL


#' =============================================================================
#' STEP 2: SELECT, SOURCE, AND VALIDATE INPUTS
#' =============================================================================

#' Use --wd= when supplied; otherwise keep the directory from which R was started.
project_dir <- resolve_project_directory()
setwd(project_dir)

#' Ask the user to select the cleaned user-framework function file.
framework_file <- browse_for_required_file(
  file_description = paste0(
    "the cleaned user-framework function file ",
    "(for example, user_framework_rpkg_revised.R)"
  ),
  allowed_extensions = "R"
)

#' Check the selected framework before sourcing it.
validate_selected_r_file(
  file = framework_file,
  required_functions = "run_simulation",
  file_description = "user-framework file"
)

#' Display the selected framework path.
cat(
  "\nSelected user-framework file:\n  ",
  framework_file,
  "\n",
  sep = ""
)

#' Store all analysis-specific result folders under this parent directory.
output_dir <- file.path(
  project_dir,
  "results"
)

#' Stop before running when the cleaned framework file cannot be found.
if (!file.exists(framework_file)) {
  stop(
    "Framework file not found: ",
    framework_file,
    call. = FALSE
  )
}

#' Source all reusable statistical and simulation functions.
source(
  framework_file,
  local = .GlobalEnv,
  chdir = TRUE
)

#' Always require the central simulation function.
required_framework_functions <- "run_simulation"

#' Require ML resource helpers only when an ML analysis is requested.
if (isTRUE(run_ml_roc) ||
    isTRUE(run_ml_adaptive)) {
  required_framework_functions <- c(
    required_framework_functions,
    "load_ml_resources",
    "make_ml_normality_configs"
  )
}

#' Identify missing framework functions.
missing_framework_functions <- required_framework_functions[
  !vapply(
    required_framework_functions,
    exists,
    logical(1),
    mode = "function",
    inherits = TRUE
  )
]

if (length(missing_framework_functions) > 0L) {
  stop(
    "The sourced framework is missing function(s): ",
    paste(missing_framework_functions, collapse = ", "),
    call. = FALSE
  )
}

#' Create the common output parent directory before any analyses begin.
output_dir <- make_run_directory(output_dir)

#' =============================================================================
#' SUPPORT FOR STEP 2: PACKAGE VALIDATION HELPERS
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Install explicitly requested missing packages or stop with a clear message.
#' -----------------------------------------------------------------------------
ensure_required_packages <- function(packages,
                                     install_missing = FALSE) {
  packages <- sort(unique(as.character(packages)))
  packages <- packages[nzchar(packages)]

  if (length(packages) == 0L) {
    return(invisible(character(0)))
  }

  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) == 0L) {
    return(invisible(packages))
  }

  if (isTRUE(install_missing)) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    utils::install.packages(missing_packages)

    still_missing <- missing_packages[
      !vapply(missing_packages, requireNamespace, logical(1), quietly = TRUE)
    ]

    if (length(still_missing) > 0L) {
      stop(
        "Package installation failed for: ",
        paste(still_missing, collapse = ", "),
        call. = FALSE
      )
    }

    return(invisible(packages))
  }

  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them with install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "),
    ")) or set install_missing_packages <- TRUE.",
    call. = FALSE
  )
}


#' -----------------------------------------------------------------------------
#' Return packages required by the selected classical normality-test codes.
#' -----------------------------------------------------------------------------
normality_test_packages <- function(test_codes) {
  package_map <- list(
    SW = character(0),
    SF = "nortest",
    LF = "nortest",
    KS = character(0),
    AD = "nortest",
    AD2 = "DescTools",
    CVM = "nortest",
    JB = "tseries",
    DAG = "moments",
    DAP = "moments",
    ANS = "moments",
    SKEW = "moments",
    KURT = "moments"
  )

  #' Compute unknown codes.
  unknown_codes <- setdiff(toupper(test_codes), names(package_map))

  if (length(unknown_codes) > 0L) {
    stop(
      "Unknown classical normality-test code(s): ",
      paste(unknown_codes, collapse = ", "),
      call. = FALSE
    )
  }

  unique(unlist(package_map[toupper(test_codes)], use.names = FALSE))
}


#' -----------------------------------------------------------------------------
#' Return packages required by named data-generating distributions.
#' -----------------------------------------------------------------------------
distribution_packages <- function(distributions) {
  package_map <- c(laplace = "LaplacesDemon", pareto = "VGAM", gumbel = "evd")

  matched_packages <- unname(
    package_map[intersect(tolower(distributions), names(package_map))]
  )

  unique(matched_packages[!is.na(matched_packages)])
}


#' -----------------------------------------------------------------------------
#' Return packages required when the external ML framework is requested.
#' -----------------------------------------------------------------------------
ml_required_packages <- function() {
  c(
    "pacman",
    "LaplacesDemon",
    "VGAM",
    "evd",
    "nortest",
    "DescTools",
    "moments",
    "tseries",
    "Lmoments",
    "robustbase",
    "pracma",
    "ineq",
    "caret",
    "glmnet",
    "randomForest",
    "gbm",
    "nnet",
    "kernlab",
    "e1071",
    "foreach",
    "doParallel",
    "pROC"
  )
}


#' =============================================================================
#' STEP 3: BUILD ANALYSIS CATALOG AND CONFIGURATIONS
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Build all supported downstream comparisons with their functions and settings.
#' -----------------------------------------------------------------------------
build_downstream_catalog <- function(nresample) {

  #' Use the selected resampling count for the one-sample sign-flip test.
  one_sample_permutation <- function(data) {
    one_sample_perm_test(
      data = data,
      nresample = nresample
    )
  }

  #' Use the selected resampling count for the two-sample permutation test.
  two_sample_permutation <- function(data) {
    two_sample_perm_test(
      data = data,
      nresample = nresample
    )
  }

  #' Use the selected resampling count for the permutation-based ANOVA.
  anova_permutation <- function(data) {
    permutation_anova(
      data = data,
      nresample = nresample
    )
  }

  #' Use the selected resampling count for the permutation-regression test.
  regression_permutation <- function(data) {
    perm_regression(
      data = data,
      nresample = nresample
    )
  }

  list(
    onesample_t_vs_sign = list(
      label = "One-sample t-test versus sign test",
      test_type = "onesample_ttest_vs_sign",
      distributions = c("exponential", "normal"),
      effect_size = 0.5,
      center_by = "median",
      gen_data = onesample_data,
      get_parameters = onesample_parameters,
      fn_to_get_norm_obj = raw_data,
      fn_for_ds_test_1 = one_sample_t_test,
      fn_for_ds_test_2 = sign_test,
      method_labels = c(test_1 = "One-sample t-test", test_2 = "Sign test", adaptive = "Adaptive procedure"),
      packages = character(0),
      extra_args = list()
    ),

    onesample_t_vs_permutation = list(
      label = "One-sample t-test versus sign-flip permutation test",
      test_type = "onesample_ttest_vs_permutation",
      distributions = c("laplace", "normal"),
      effect_size = 0.5,
      center_by = "mean",
      gen_data = onesample_data,
      get_parameters = onesample_parameters,
      fn_to_get_norm_obj = raw_data,
      fn_for_ds_test_1 = one_sample_t_test,
      fn_for_ds_test_2 = one_sample_permutation,
      method_labels = c(
        test_1 = "One-sample t-test",
        test_2 = "Permutation test",
        adaptive = "Adaptive procedure"
      ),
      packages = "LaplacesDemon",
      extra_args = list()
    ),

    onesample_t_vs_boxcox = list(
      label = "One-sample t-test versus Box-Cox t-test",
      test_type = "onesample_ttest_vs_boxcox",
      distributions = c("exponential", "normal"),
      effect_size = 0.5,
      center_by = "median",
      gen_data = onesample_data,
      get_parameters = onesample_parameters,
      fn_to_get_norm_obj = raw_data,
      fn_for_ds_test_1 = one_sample_t_test,
      fn_for_ds_test_2 = boxcox_t_test,
      method_labels = c(
        test_1 = "One-sample t-test",
        test_2 = "Box-Cox t-test",
        adaptive = "Adaptive procedure"
      ),
      packages = character(0),
      extra_args = list()
    ),

    twosample_t_vs_mann_whitney = list(
      label = "Two-sample t-test versus Mann-Whitney U test",
      test_type = "twosample_ttest_vs_mann_whitney",
      distributions = c("exponential", "normal"),
      effect_size = 0.5,
      center_by = "median",
      gen_data = two_sample_data,
      get_parameters = twosample_parameters,
      fn_to_get_norm_obj = raw_data,
      fn_for_ds_test_1 = twosample_t_test,
      fn_for_ds_test_2 = Mann_whitney_U_test,
      method_labels = c(
        test_1 = "Two-sample t-test",
        test_2 = "Mann-Whitney U test",
        adaptive = "Adaptive procedure"
      ),
      packages = character(0),
      extra_args = list()
    ),

    twosample_t_vs_permutation = list(
      label = "Two-sample t-test versus permutation test",
      test_type = "twosample_ttest_vs_permutation",
      distributions = c("exponential", "normal"),
      effect_size = 0.5,
      center_by = "mean",
      gen_data = two_sample_data,
      get_parameters = twosample_parameters,
      fn_to_get_norm_obj = raw_data,
      fn_for_ds_test_1 = twosample_t_test,
      fn_for_ds_test_2 = two_sample_permutation,
      method_labels = c(
        test_1 = "Two-sample t-test",
        test_2 = "Permutation test",
        adaptive = "Adaptive procedure"
      ),
      packages = character(0),
      extra_args = list()
    ),

    anova_vs_kruskal_wallis = list(
      label = "One-way ANOVA versus Kruskal-Wallis test",
      test_type = "anova_vs_kruskal_wallis",
      distributions = c("exponential", "normal"),
      effect_size = list(
        H0 = c(0, 0, 0),
        H1 = c(0, 0, 0.5)
      ),
      center_by = "median",
      gen_data = anova_gen_data,
      get_parameters = anova_parameters,
      fn_to_get_norm_obj = anova_residuals,
      fn_for_ds_test_1 = one_way_anova,
      fn_for_ds_test_2 = kruskal_wallis_test,
      method_labels = c(
        test_1 = "One-way ANOVA",
        test_2 = "Kruskal-Wallis test",
        adaptive = "Adaptive procedure"
      ),
      packages = character(0),
      extra_args = list()
    ),

    anova_vs_permutation = list(
      label = "One-way ANOVA versus permutation ANOVA",
      test_type = "anova_vs_permutation",
      distributions = c("exponential", "normal"),
      effect_size = list(
        H0 = c(0, 0, 0),
        H1 = c(0, 0, 0.5)
      ),
      center_by = "median",
      gen_data = anova_gen_data,
      get_parameters = anova_parameters,
      fn_to_get_norm_obj = anova_residuals,
      fn_for_ds_test_1 = one_way_anova,
      fn_for_ds_test_2 = anova_permutation,
      method_labels = c(
        test_1 = "One-way ANOVA",
        test_2 = "Permutation ANOVA",
        adaptive = "Adaptive procedure"
      ),
      packages = "coin",
      extra_args = list()
    ),

    regression_ols_vs_rank = list(
      label = "OLS regression versus rank-based regression",
      test_type = "regression_ols_vs_rank",
      distributions = c("exponential", "normal"),
      effect_size = 0.5,
      center_by = "mean",
      gen_data = reg_data,
      get_parameters = reg_parameters,
      fn_to_get_norm_obj = reg_residuals,
      fn_for_ds_test_1 = simple_linear_reg,
      fn_for_ds_test_2 = rank_regression,
      method_labels = c(
        test_1 = "OLS slope test",
        test_2 = "Rank-based regression",
        adaptive = "Adaptive procedure"
      ),
      packages = "Rfit",
      extra_args = list(
        x_dist = "exponential",
        x_par = 1,
        error_par = NULL
      )
    ),

    regression_ols_vs_permutation = list(
      label = "OLS regression versus permutation regression",
      test_type = "regression_ols_vs_permutation",
      distributions = c("exponential", "normal"),
      effect_size = 0.5,
      center_by = "mean",
      gen_data = reg_data,
      get_parameters = reg_parameters,
      fn_to_get_norm_obj = reg_residuals,
      fn_for_ds_test_1 = simple_linear_reg,
      fn_for_ds_test_2 = regression_permutation,
      method_labels = c(
        test_1 = "OLS slope test",
        test_2 = "Permutation regression",
        adaptive = "Adaptive procedure"
      ),
      packages = character(0),
      extra_args = list(
        x_dist = "exponential",
        x_par = 1,
        error_par = NULL
      )
    ),

    regression_ols_vs_mm = list(
      label = "OLS regression versus MM-estimation regression",
      test_type = "regression_ols_vs_mm",
      distributions = c("exponential", "normal"),
      effect_size = 0.5,
      center_by = "mean",
      gen_data = reg_data,
      get_parameters = reg_parameters,
      fn_to_get_norm_obj = reg_residuals,
      fn_for_ds_test_1 = simple_linear_reg,
      fn_for_ds_test_2 = mm_regression,
      method_labels = c(
        test_1 = "OLS slope test",
        test_2 = "MM-estimation regression",
        adaptive = "Adaptive procedure"
      ),
      packages = "robustbase",
      extra_args = list(
        x_dist = "exponential",
        x_par = 1,
        error_par = NULL
      )
    ),

    regression_ols_vs_m = list(
      label = "OLS regression versus M-estimation regression",
      test_type = "regression_ols_vs_m",
      distributions = c("exponential", "normal"),
      effect_size = 0.5,
      center_by = "mean",
      gen_data = reg_data,
      get_parameters = reg_parameters,
      fn_to_get_norm_obj = reg_residuals,
      fn_for_ds_test_1 = simple_linear_reg,
      fn_for_ds_test_2 = m_regression,
      method_labels = c(
        test_1 = "OLS slope test",
        test_2 = "M-estimation regression",
        adaptive = "Adaptive procedure"
      ),
      packages = "MASS",
      extra_args = list(
        x_dist = "exponential",
        x_par = 1,
        error_par = NULL
      )
    )
  )
}


#' Build the reusable catalog after the resampling count is resolved.
downstream_catalog <- build_downstream_catalog(
  nresample = run_settings$nresample
)

#' Stop when the user selected an analysis identifier absent from the catalog.
unknown_downstream_analyses <- setdiff(
  selected_downstream_analyses,
  names(downstream_catalog)
)

if (length(unknown_downstream_analyses) > 0L) {
  stop(
    "Unknown downstream analysis identifier(s): ",
    paste(unknown_downstream_analyses, collapse = ", "),
    ". Available identifiers: ",
    paste(names(downstream_catalog), collapse = ", "),
    call. = FALSE
  )
}


#' Validate ML downstream identifiers independently from classical selections.
unknown_ml_downstream_analyses <- setdiff(
  ml_downstream_analyses,
  names(downstream_catalog)
)

if (length(unknown_ml_downstream_analyses) > 0L) {
  stop(
    "Unknown ML downstream analysis identifier(s): ",
    paste(unknown_ml_downstream_analyses, collapse = ", "),
    ". Available identifiers: ",
    paste(names(downstream_catalog), collapse = ", "),
    call. = FALSE
  )
}

if (isTRUE(ml_requested) && length(ml_downstream_analyses) == 0L) {
  stop(
    "Select at least one value in ml_downstream_analyses when ML is requested.",
    call. = FALSE
  )
}


#' =============================================================================
#' SUPPORT FOR STEP 3: NORMALITY CONFIGURATIONS
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Combine Shapiro-Wilk and Anderson-Darling p-values using Fisher's method.
#' -----------------------------------------------------------------------------
fisher_combined_test <- function(x) {
  p_sw <- stats::shapiro.test(x)$p.value
  p_ad <- nortest::ad.test(x)$p.value

  #' Protect the logarithm from numerical zero.
  p_values <- pmax(
    c(p_sw, p_ad),
    .Machine$double.xmin
  )

  fisher_statistic <- -2 * sum(log(p_values))

  c(p.value = stats::pchisq(fisher_statistic, df = 4, lower.tail = FALSE))
}


#' -----------------------------------------------------------------------------
#' Create the standard configuration for one classical normality test.
#' -----------------------------------------------------------------------------
make_classical_normality_config <- function(test_code) {
  list(
    method = "classical",
    config = list(
      norm_test = toupper(test_code)
    )
  )
}


#' -----------------------------------------------------------------------------
#' Create the custom configuration for the Fisher-combined procedure.
#' -----------------------------------------------------------------------------
make_fisher_normality_config <- function() {
  list(
    method = "custom",
    config = list(
      fn = fisher_combined_test,
      label = "Fisher SW + AD"
    )
  )
}


#' =============================================================================
#' SUPPORT FOR STEP 3: BUILD THE ANALYSIS RUN PLAN
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Create a valid file-safe suffix from an analysis method label.
#' -----------------------------------------------------------------------------
clean_run_suffix <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}


#' -----------------------------------------------------------------------------
#' Create one fully specified request for run_simulation().
#' -----------------------------------------------------------------------------
make_analysis_request <- function(design_id,
                                  method_id,
                                  norm_config,
                                  phases,
                                  per_n_thresholds,
                                  threshold_grid,
                                  normality_roc_grid) {
  design <- downstream_catalog[[design_id]]
  method_suffix <- clean_run_suffix(method_id)

  test_type <- paste(
    design$test_type,
    method_suffix,
    sep = "_"
  )

  if (isTRUE(per_n_thresholds) && 2L %in% phases) {
    test_type <- paste0(test_type, "_per_n")
  }

  list(
    id = test_type,
    design_id = design_id,
    design_label = design$label,
    method_id = method_id,
    norm_config = norm_config,
    phases = sort(unique(as.integer(phases))),
    per_n_thresholds = isTRUE(per_n_thresholds),
    threshold_grid = threshold_grid,
    normality_roc_grid = normality_roc_grid
  )
}


#' -----------------------------------------------------------------------------
#' Return validated threshold-selection modes for analyses that use a threshold.
#' -----------------------------------------------------------------------------
resolve_threshold_modes <- function() {
  valid_modes <- c("single", "per_n")
  selected_modes <- unique(tolower(trimws(threshold_modes_to_run)))
  invalid_modes <- setdiff(selected_modes, valid_modes)

  if (length(invalid_modes) > 0L) {
    stop(
      "Invalid threshold mode(s): ",
      paste(invalid_modes, collapse = ", "),
      ". Allowed values are 'single' and 'per_n'.",
      call. = FALSE
    )
  }

  #' Validate the current condition.
  if (length(selected_modes) == 0L) {
    stop("Select at least one threshold mode.", call. = FALSE)
  }

  valid_modes[valid_modes %in% selected_modes]
}


#' -----------------------------------------------------------------------------
#' Build all requests implied by the selected analyses and approach switches.
#' -----------------------------------------------------------------------------
build_run_plan <- function() {
  requests <- list()
  request_index <- 1L
  threshold_modes <- resolve_threshold_modes()

  for (design_id in selected_downstream_analyses) {
    classical_roc_added <- FALSE

    #' Add one request for every selected classical normality method and mode.
    if (isTRUE(run_classical_analyses)) {
      classical_methods <- toupper(classical_methods)

      for (method_code in classical_methods) {
        for (threshold_mode in threshold_modes) {
          phases <- classical_phases

          #' Calculate the common classical ROC comparison once per design.
          if (isTRUE(run_classical_roc_once) && classical_roc_added) {
            phases <- setdiff(phases, 1L)
          }

          requests[[request_index]] <- make_analysis_request(
            design_id = design_id,
            method_id = method_code,
            norm_config = make_classical_normality_config(method_code),
            phases = phases,
            per_n_thresholds = identical(threshold_mode, "per_n"),
            threshold_grid = threshold_grid_classical,
            normality_roc_grid = normality_roc_grid_classical
          )

          #' Validate the current condition.
          if (1L %in% phases) {
            classical_roc_added <- TRUE
          }

          request_index <- request_index + 1L
        }
      }
    }

    #' Add the Fisher-combined procedure for every selected threshold mode.
    if (isTRUE(run_fisher_analyses)) {
      fisher_roc_added <- FALSE

      for (threshold_mode in threshold_modes) {
        phases <- fisher_phases

        #' Avoid repeating the same Fisher ROC curve across threshold modes.
        if (fisher_roc_added) {
          phases <- setdiff(phases, 1L)
        }

        requests[[request_index]] <- make_analysis_request(
          design_id = design_id,
          method_id = "Fisher",
          norm_config = make_fisher_normality_config(),
          phases = phases,
          per_n_thresholds = identical(threshold_mode, "per_n"),
          threshold_grid = threshold_grid_classical,
          normality_roc_grid = normality_roc_grid_classical
        )

        #' Validate the current condition.
        if (1L %in% phases) {
          fisher_roc_added <- TRUE
        }

        request_index <- request_index + 1L
      }
    }

  }

  #' Add one multi-model ML ROC request for each selected ML design.
  if (isTRUE(run_ml_roc)) {
    for (design_id in ml_downstream_analyses) {
      requests[[request_index]] <- make_analysis_request(
        design_id = design_id,
        method_id = "ML_ROC",
        norm_config = ml_configs$roc,
        phases = ml_roc_phases,
        per_n_thresholds = FALSE,
        threshold_grid = threshold_grid_ml,
        normality_roc_grid = normality_roc_grid_ml
      )

      request_index <- request_index + 1L
    }
  }

  #' Add adaptive ML requests for every selected design and threshold mode.
  if (isTRUE(run_ml_adaptive)) {
    for (design_id in ml_downstream_analyses) {
      for (threshold_mode in threshold_modes) {
        requests[[request_index]] <- make_analysis_request(
          design_id = design_id,
          method_id = paste0("ML_", ml_chosen_model),
          norm_config = ml_configs$decision,
          phases = ml_adaptive_phases,
          per_n_thresholds = identical(threshold_mode, "per_n"),
          threshold_grid = threshold_grid_ml,
          normality_roc_grid = normality_roc_grid_ml
        )

        request_index <- request_index + 1L
      }
    }
  }

  if (length(requests) == 0L) {
    stop(
      "No analyses were selected. Enable at least one normality approach.",
      call. = FALSE
    )
  }

  #' Compute request ids.
  request_ids <- vapply(requests, `[[`, character(1), "id")

  if (anyDuplicated(request_ids)) {
    duplicated_ids <- unique(request_ids[duplicated(request_ids)])

    stop(
      "The run plan contains duplicate analysis identifiers: ",
      paste(duplicated_ids, collapse = ", "),
      call. = FALSE
    )
  }

  names(requests) <- request_ids
  requests
}

#' -----------------------------------------------------------------------------
#' Return the package set required by the current selections.
#' -----------------------------------------------------------------------------
collect_required_packages <- function() {
  required_packages <- character(0)

  for (design_id in selected_downstream_analyses) {
    design <- downstream_catalog[[design_id]]

    required_packages <- c(required_packages, design$packages, distribution_packages(design$distributions))
  }

  #' Validate the current condition.
  if (isTRUE(run_classical_analyses)) {
    required_packages <- c(required_packages, normality_test_packages(classical_methods))

    if (1L %in% classical_phases) {
      required_packages <- c(required_packages, normality_test_packages(normality_tests_for_roc))
    }
  }

  #' Validate the current condition.
  if (isTRUE(run_fisher_analyses)) {
    required_packages <- c(required_packages, "nortest")
  }

  if (isTRUE(ml_requested)) {
    #' Include design-specific requirements for every selected ML procedure.
    for (design_id in ml_downstream_analyses) {
      ml_design <- downstream_catalog[[design_id]]

      required_packages <- c(
        required_packages,
        ml_design$packages,
        distribution_packages(ml_design$distributions)
      )
    }

    required_packages <- c(
      required_packages,
      ml_required_packages()
    )
  }

  sort(unique(required_packages))
}


#' Validate the packages needed by the selected analyses before loading ML files.
required_packages <- collect_required_packages()

ensure_required_packages(
  packages = required_packages,
  install_missing = install_missing_packages
)

#' Load trained ML models only after package validation and only when requested.
if (isTRUE(ml_requested)) {
  #' Ask for the exact ML function file used during training.
  ml_framework_file <- browse_for_required_file(
    file_description = paste0(
      "the ML function file used to train the models ",
      "(select ML_framework_func.R, not the ML training runner)"
    ),
    allowed_extensions = "R"
  )

  #' Reject an ML runner or incomplete function file before sourcing it.
  validate_selected_r_file(
    file = ml_framework_file,
    required_functions = c("calculate_features", "prepare_prediction_features", "predict_prob_normal_fast"),
    file_description = "ML function file"
  )

  #' Ask for the saved trained-model bundles.
  ml_model_file <- browse_for_required_file(
    file_description = paste0(
      "the trained-model RData file containing sample-size-specific bundles"
    ),
    allowed_extensions = c("RData", "rda")
  )

  #' Display the selected ML resources.
  cat(
    "\nSelected ML function file:\n  ",
    ml_framework_file,
    "\nSelected trained-model file:\n  ",
    ml_model_file,
    "\n",
    sep = ""
  )

  #' Load and validate the selected ML resources.
  ml_resources <- load_ml_resources(
    framework_file = ml_framework_file,
    model_file = ml_model_file,
    project_dir = project_dir,
    required_sample_sizes = NULL
  )

  #' Display the trained sample sizes available for nearest-size selection.
  cat(
    "\nAvailable trained ML sample sizes: ",
    paste(ml_resources$available_sample_sizes, collapse = ", "),
    "\n",
    sep = ""
  )

  #' Create the ROC and adaptive-decision configurations.
  ml_configs <- make_ml_normality_configs(
    ml_resources = ml_resources,
    chosen_model = ml_chosen_model,
    models_for_roc = ml_models_for_roc,
    use_majority_vote = ml_use_majority_vote
  )

  #' Reset the analysis RNG because some external ML files set their own seed.
  set_reproducible_rng(
    seed = base_seed,
    rng_version = rng_version
  )
}

#' Build the final list of analysis requests after optional ML setup is complete.
run_plan <- build_run_plan()


#' =============================================================================
#' SUPPORT FOR STEP 4: RUN-PLAN DISPLAY AND EXECUTION HELPERS
#' =============================================================================

#' -----------------------------------------------------------------------------
#' Convert the run plan into a concise table for review and record keeping.
#' -----------------------------------------------------------------------------
run_plan_table <- function(run_plan) {
  do.call(
    rbind,
    lapply(seq_along(run_plan), function(i) {
      request <- run_plan[[i]]
      design <- downstream_catalog[[request$design_id]]

      data.frame(
        Order = i,
        Analysis_ID = request$id,
        Downstream_Comparison = request$design_label,
        Normality_Method = request$method_id,
        Phases = paste(request$phases, collapse = ","),
        Per_N_Thresholds = request$per_n_thresholds,
        Distributions = paste(design$distributions, collapse = " vs "),
        Seed = base_seed + i - 1L,
        stringsAsFactors = FALSE
      )
    })
  )
}


#' -----------------------------------------------------------------------------
#' Print the complete analysis plan before any simulations begin.
#' -----------------------------------------------------------------------------
print_run_plan <- function(run_plan) {
  plan_table <- run_plan_table(run_plan)

  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("USER-FRAMEWORK ANALYSIS PLAN\n")
  cat(strrep("=", 78), "\n", sep = "")
  print(plan_table, row.names = FALSE)
  cat("\nOutput parent directory:", output_dir, "\n")
  cat("Run mode:", run_mode, "\n")
  cat("Nsim:", run_settings$Nsim, "\n")
  cat("N_tradeoff:", run_settings$N_tradeoff, "\n")
  cat("Resamples:", run_settings$nresample, "\n")
  cat(strrep("=", 78), "\n\n", sep = "")

  invisible(plan_table)
}


#' -----------------------------------------------------------------------------
#' Run one request by combining shared settings with its downstream design.
#' -----------------------------------------------------------------------------
run_one_analysis <- function(request,
                             seed) {
  design <- downstream_catalog[[request$design_id]]

  #' Reset the reproducible seed immediately before this analysis.
  set_reproducible_rng(
    seed = seed,
    rng_version = rng_version
  )

  shared_arguments <- list(
    Nsim = run_settings$Nsim,
    N_tradeoff = run_settings$N_tradeoff,
    test_type = request$id,
    distributions = design$distributions,
    norm_config = request$norm_config,
    threshold_grid = request$threshold_grid,
    normality_roc_grid = request$normality_roc_grid,
    tol_pos = type1_inflation_tolerance,
    loss_tol = power_loss_tolerance,
    test_alpha = test_alpha,
    center_by = design$center_by,
    effect_size = design$effect_size,
    sample_sizes = sample_sizes,
    single_n = single_n,
    effect_sizes_plot = effect_sizes_plot,
    sig_levels = significance_levels,
    ds_test_methods = downstream_methods,
    method_labels = design$method_labels,
    type1_ylim = type1_ylim,
    reference_effect = reference_effect,
    norm_test = normality_tests_for_roc,
    selected_tests = selected_normality_tests_for_roc,
    phases = request$phases,
    pretest_threshold = pretest_threshold,
    per_n_thresholds = request$per_n_thresholds,
    output_dir = output_dir,
    save_results = save_results,
    save_compression = save_compression,
    gen_data = design$gen_data,
    get_parameters = design$get_parameters,
    fn_to_get_norm_obj = design$fn_to_get_norm_obj,
    fn_for_ds_test_1 = design$fn_for_ds_test_1,
    fn_for_ds_test_2 = design$fn_for_ds_test_2
  )

  #' Pass only design-specific generator arguments through the ellipsis.
  do.call(
    run_simulation,
    c(shared_arguments, design$extra_args)
  )
}


#' -----------------------------------------------------------------------------
#' Execute every request, record status information, and optionally retain results.
#' -----------------------------------------------------------------------------
run_selected_analyses <- function(run_plan,
                                  keep_results = TRUE,
                                  stop_on_error = TRUE) {
  stored_results <- list()
  manifest_rows <- vector("list", length(run_plan))

  for (i in seq_along(run_plan)) {
    request <- run_plan[[i]]
    analysis_seed <- base_seed + i - 1L
    analysis_start <- Sys.time()

    cat("\n", strrep("#", 78), "\n", sep = "")
    cat("RUNNING ANALYSIS ", i, " OF ", length(run_plan), "\n", sep = "")
    cat("Analysis ID: ", request$id, "\n", sep = "")
    cat(strrep("#", 78), "\n", sep = "")

    analysis_result <- tryCatch(
      run_one_analysis(
        request = request,
        seed = analysis_seed
      ),
      error = function(error_condition) {
        structure(
          list(
            message = conditionMessage(error_condition),
            call = conditionCall(error_condition)
          ),
          class = "framework_run_error"
        )
      }
    )

    analysis_end <- Sys.time()
    elapsed_minutes <- as.numeric(
      difftime(
        analysis_end,
        analysis_start,
        units = "mins"
      )
    )

    failed <- inherits(analysis_result, "framework_run_error")

    if (failed) {
      status <- "failed"
      error_message <- analysis_result$message
      result_file <- NA_character_
      output_directory <- file.path(output_dir, request$id)

      message(
        "Analysis failed: ",
        request$id,
        "\n",
        error_message
      )
    } else {
      status <- "completed"
      error_message <- NA_character_
      result_file <- analysis_result$files$results_rdata %||% NA_character_
      output_directory <- analysis_result$root_dir %||%
        file.path(output_dir, request$id)

      if (isTRUE(keep_results)) {
        stored_results[[request$id]] <- analysis_result
      }
    }

    manifest_rows[[i]] <- data.frame(
      Order = i,
      Analysis_ID = request$id,
      Downstream_Comparison = request$design_label,
      Normality_Method = request$method_id,
      Phases = paste(request$phases, collapse = ","),
      Per_N_Thresholds = request$per_n_thresholds,
      Seed = analysis_seed,
      Status = status,
      Elapsed_Minutes = elapsed_minutes,
      Output_Directory = output_directory,
      Results_File = result_file,
      Error_Message = error_message,
      stringsAsFactors = FALSE
    )

    if (failed && isTRUE(stop_on_error)) {
      manifest <- do.call(rbind, manifest_rows[seq_len(i)])

      stop(
        "Stopped after failed analysis '",
        request$id,
        "'.\n",
        error_message,
        call. = FALSE
      )
    }

    if (!isTRUE(keep_results) && !failed) {
      rm(analysis_result)
      invisible(gc())
    }
  }

  manifest <- do.call(rbind, manifest_rows)

  list(
    manifest = manifest,
    results = stored_results,
    run_plan = run_plan,
    settings = list(
      project_dir = project_dir,
      framework_file = framework_file,
      ml_requested = ml_requested,
      ml_downstream_analyses = ml_downstream_analyses,
      ml_chosen_model = ml_chosen_model,
      ml_framework_file = ml_framework_file,
      ml_model_file = ml_model_file,
      output_dir = output_dir,
      run_mode = run_mode,
      run_settings = run_settings,
      sample_sizes = sample_sizes,
      single_n = single_n,
      test_alpha = test_alpha,
      threshold_grid_classical = threshold_grid_classical,
      threshold_grid_ml = threshold_grid_ml,
      normality_roc_grid_classical = normality_roc_grid_classical,
      normality_roc_grid_ml = normality_roc_grid_ml,
      type1_inflation_tolerance = type1_inflation_tolerance,
      power_loss_tolerance = power_loss_tolerance,
      completed_at = Sys.time()
    )
  )
}


#' -----------------------------------------------------------------------------
#' Save the lightweight run manifest and settings without duplicating full results.
#' -----------------------------------------------------------------------------
save_run_manifest <- function(run_output,
                              output_dir) {
  manifest_file <- file.path(
    output_dir,
    "framework_analysis_run_manifest.rds"
  )

  manifest_csv <- file.path(
    output_dir,
    "framework_analysis_run_manifest.csv"
  )

  lightweight_manifest <- list(
    manifest = run_output$manifest,
    settings = run_output$settings,
    run_plan = lapply(run_output$run_plan, function(request) {
      request$norm_config <- list(
        method = request$norm_config$method,
        config_names = names(request$norm_config$config)
      )
      request
    })
  )

  saveRDS(
    lightweight_manifest,
    file = manifest_file,
    compress = "gzip"
  )

  utils::write.csv(
    run_output$manifest,
    file = manifest_csv,
    row.names = FALSE
  )

  c(rds = manifest_file, csv = manifest_csv)
}


#' =============================================================================
#' STEP 4: EXECUTE THE SELECTED WORKFLOW
#' =============================================================================

#' Display the final plan before the first simulation starts.
run_plan_preview <- print_run_plan(run_plan)

if (!isTRUE(run_workflow) || isTRUE(dry_run)) {
  message(
    "The run plan was created but not executed because run_workflow = FALSE ",
    "or dry_run = TRUE."
  )

  framework_run <- list(
    manifest = run_plan_preview,
    results = list(),
    run_plan = run_plan
  )

} else {
  framework_run <- run_selected_analyses(
    run_plan = run_plan,
    keep_results = keep_results_in_memory,
    stop_on_error = stop_on_error
  )

  if (isTRUE(write_run_manifest)) {
    manifest_files <- save_run_manifest(
      run_output = framework_run,
      output_dir = output_dir
    )

    cat("\nRun manifest files:\n")
    print(manifest_files)
  }

  cat("\nFinal run summary:\n")
  print(framework_run$manifest, row.names = FALSE)
}


#' =============================================================================
#' STEP 5: READY-TO-USE RUN EXAMPLES
#'
#' Copy the settings from ONE example into Step 1.
#' Keep these examples commented so they do not overwrite the active controls.
#' =============================================================================

#' -----------------------------------------------------------------------------
#' EXAMPLE 1: Complete classical one-sample t-test versus sign-test analysis
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- "onesample_t_vs_sign"
#' run_classical_analyses <- TRUE
#' classical_methods <- "SW"
#' classical_phases <- 1:5
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- FALSE
#' run_ml_adaptive <- FALSE
#' threshold_modes_to_run <- "per_n"
#' pretest_threshold <- NULL


#' -----------------------------------------------------------------------------
#' EXAMPLE 2: Classical two-sample t-test versus Mann-Whitney U test
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- "twosample_t_vs_mann_whitney"
#' run_classical_analyses <- TRUE
#' classical_methods <- "SW"
#' classical_phases <- 1:5
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- FALSE
#' run_ml_adaptive <- FALSE
#' threshold_modes_to_run <- "per_n"
#' pretest_threshold <- NULL


#' -----------------------------------------------------------------------------
#' EXAMPLE 3: Classical ANOVA versus Kruskal-Wallis analysis
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- "anova_vs_kruskal_wallis"
#' run_classical_analyses <- TRUE
#' classical_methods <- "SW"
#' classical_phases <- 1:5
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- FALSE
#' run_ml_adaptive <- FALSE
#' threshold_modes_to_run <- "per_n"
#' pretest_threshold <- NULL


#' -----------------------------------------------------------------------------
#' EXAMPLE 4: Classical OLS versus rank-based regression analysis
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- "regression_ols_vs_rank"
#' run_classical_analyses <- TRUE
#' classical_methods <- "SW"
#' classical_phases <- 1:5
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- FALSE
#' run_ml_adaptive <- FALSE
#' threshold_modes_to_run <- "per_n"
#' pretest_threshold <- NULL


#' -----------------------------------------------------------------------------
#' EXAMPLE 5: Phase 1 classical normality ROC analysis only
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- "onesample_t_vs_sign"
#' run_classical_analyses <- TRUE
#' classical_methods <- "SW"
#' classical_phases <- 1L
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- FALSE
#' run_ml_adaptive <- FALSE
#' threshold_modes_to_run <- "single"
#' pretest_threshold <- NULL


#' -----------------------------------------------------------------------------
#' EXAMPLE 6: Phase 2 threshold trade-off analysis only
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- "anova_vs_kruskal_wallis"
#' run_classical_analyses <- TRUE
#' classical_methods <- "SW"
#' classical_phases <- 2L
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- FALSE
#' run_ml_adaptive <- FALSE
#' threshold_modes_to_run <- "per_n"
#' pretest_threshold <- NULL


#' -----------------------------------------------------------------------------
#' EXAMPLE 7: Phase 4 using supplied sample-size-specific thresholds
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- "twosample_t_vs_mann_whitney"
#' run_classical_analyses <- TRUE
#' classical_methods <- "SW"
#' classical_phases <- 4L
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- FALSE
#' run_ml_adaptive <- FALSE
#' threshold_modes_to_run <- "per_n"
#' pretest_threshold <- c(
#'   "10" = 0.05,
#'   "20" = 0.05,
#'   "30" = 0.05,
#'   "40" = 0.05,
#'   "50" = 0.05
#' )


#' -----------------------------------------------------------------------------
#' EXAMPLE 8: Phase 5 using one supplied threshold
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- "onesample_t_vs_permutation"
#' run_classical_analyses <- TRUE
#' classical_methods <- "SW"
#' classical_phases <- 5L
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- FALSE
#' run_ml_adaptive <- FALSE
#' threshold_modes_to_run <- "single"
#' pretest_threshold <- 0.05


#' -----------------------------------------------------------------------------
#' EXAMPLE 9: Run several classical downstream comparisons
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- c(
#'   "onesample_t_vs_sign",
#'   "twosample_t_vs_mann_whitney",
#'   "anova_vs_kruskal_wallis",
#'   "regression_ols_vs_rank"
#' )
#' run_classical_analyses <- TRUE
#' classical_methods <- "SW"
#' classical_phases <- 1:5
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- FALSE
#' run_ml_adaptive <- FALSE
#' threshold_modes_to_run <- "per_n"
#' pretest_threshold <- NULL


#' -----------------------------------------------------------------------------
#' EXAMPLE 10: Fisher-combined procedure only
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- "onesample_t_vs_sign"
#' run_classical_analyses <- FALSE
#' run_fisher_analyses <- TRUE
#' fisher_phases <- 1:5
#' run_ml_roc <- FALSE
#' run_ml_adaptive <- FALSE
#' threshold_modes_to_run <- "per_n"
#' pretest_threshold <- NULL


#' -----------------------------------------------------------------------------
#' EXAMPLE 11: ML ROC comparison for one-sample data only
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- character(0)
#' run_classical_analyses <- FALSE
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- TRUE
#' ml_roc_phases <- 1L
#' run_ml_adaptive <- FALSE
#' ml_downstream_analyses <- "onesample_t_vs_sign"
#' threshold_modes_to_run <- "single"
#' pretest_threshold <- NULL


#' -----------------------------------------------------------------------------
#' EXAMPLE 12: Complete adaptive ML analysis for a two-sample procedure
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- character(0)
#' run_classical_analyses <- FALSE
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- TRUE
#' ml_roc_phases <- 1L
#' run_ml_adaptive <- TRUE
#' ml_adaptive_phases <- 2:5
#' ml_downstream_analyses <- "twosample_t_vs_mann_whitney"
#' ml_chosen_model <- "SVM"
#' threshold_modes_to_run <- "per_n"
#' pretest_threshold <- NULL


#' -----------------------------------------------------------------------------
#' EXAMPLE 13: Complete adaptive ML analysis using ANOVA residuals
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- character(0)
#' run_classical_analyses <- FALSE
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- TRUE
#' ml_roc_phases <- 1L
#' run_ml_adaptive <- TRUE
#' ml_adaptive_phases <- 2:5
#' ml_downstream_analyses <- "anova_vs_kruskal_wallis"
#' ml_chosen_model <- "SVM"
#' threshold_modes_to_run <- "per_n"
#' pretest_threshold <- NULL


#' -----------------------------------------------------------------------------
#' EXAMPLE 14: Complete adaptive ML analysis using regression residuals
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- character(0)
#' run_classical_analyses <- FALSE
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- TRUE
#' ml_roc_phases <- 1L
#' run_ml_adaptive <- TRUE
#' ml_adaptive_phases <- 2:5
#' ml_downstream_analyses <- "regression_ols_vs_rank"
#' ml_chosen_model <- "SVM"
#' threshold_modes_to_run <- "per_n"
#' pretest_threshold <- NULL


#' -----------------------------------------------------------------------------
#' EXAMPLE 15: Apply ML to several raw-sample and residual-based procedures
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- character(0)
#' run_classical_analyses <- FALSE
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- TRUE
#' run_ml_adaptive <- TRUE
#' ml_downstream_analyses <- c(
#'   "onesample_t_vs_sign",
#'   "twosample_t_vs_mann_whitney",
#'   "anova_vs_kruskal_wallis",
#'   "regression_ols_vs_rank"
#' )
#' ml_roc_phases <- 1L
#' ml_adaptive_phases <- 2:5
#' ml_chosen_model <- "SVM"
#' threshold_modes_to_run <- "per_n"
#' pretest_threshold <- NULL


#' -----------------------------------------------------------------------------
#' EXAMPLE 16: Compare SW and ML on the same ANOVA procedure
#' -----------------------------------------------------------------------------
#' selected_downstream_analyses <- "anova_vs_kruskal_wallis"
#' run_classical_analyses <- TRUE
#' classical_methods <- "SW"
#' classical_phases <- 1:5
#' run_fisher_analyses <- FALSE
#' run_ml_roc <- TRUE
#' run_ml_adaptive <- TRUE
#' ml_downstream_analyses <- "anova_vs_kruskal_wallis"
#' ml_chosen_model <- "SVM"
#' threshold_modes_to_run <- "per_n"
#' pretest_threshold <- NULL
