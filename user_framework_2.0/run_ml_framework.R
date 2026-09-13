#' =============================================================================
#' FLEXIBLE ML NORMALITY WORKFLOW RUNNER
#' =============================================================================
#'
#' This runner organizes the existing functions in ML_framework_func.R.
#' It separates VIP validation distributions from final evaluation distributions
#' without changing model fitting, permutation VIP, evaluation, or ROC formulas.
#'
#' Workflows
#' - "full"          : feature selection, training, evaluation, and ROC analyses
#' - "training_only" : train selected models and save them
#' - "roc_only"      : load trained models and run only classical-versus-ML ROC
#' - "custom"        : choose individual stages with ml_stage_options()
#'
#' File modes
#' - "paths"  : use paths supplied through ml_file_paths()
#' - "browse" : prompt for each required external input file
#'
#' The examples at the end do not run automatically.
#' =============================================================================


#' =============================================================================
#' SECTION 1: USER-CONFIGURATION HELPERS
#' =============================================================================

#' Store optional external file paths.
ml_file_paths <- function(
    framework = NULL,
    trained_models = NULL,
    evaluation_results = NULL,
    final_features = NULL
) {
  list(
    framework = framework,
    trained_models = trained_models,
    evaluation_results = evaluation_results,
    final_features = final_features
  )
}


#' Override one or more workflow stages; NULL keeps the workflow preset.
ml_stage_options <- function(
    feature_selection = NULL,
    training = NULL,
    evaluation = NULL,
    ml_roc_plots = NULL,
    permutation_vip = NULL,
    classical_ml_roc = NULL,
    density_plots = NULL,
    load_evaluation = NULL,
    reuse_evaluation_for_roc = NULL
) {
  list(
    feature_selection = feature_selection,
    training = training,
    evaluation = evaluation,
    ml_roc_plots = ml_roc_plots,
    permutation_vip = permutation_vip,
    classical_ml_roc = classical_ml_roc,
    density_plots = density_plots,
    load_evaluation = load_evaluation,
    reuse_evaluation_for_roc = reuse_evaluation_for_roc
  )
}


#' Return a stage plan with every stage disabled.
empty_stage_options <- function() {
  list(
    feature_selection = FALSE,
    training = FALSE,
    evaluation = FALSE,
    ml_roc_plots = FALSE,
    permutation_vip = FALSE,
    classical_ml_roc = FALSE,
    density_plots = FALSE,
    load_evaluation = FALSE,
    reuse_evaluation_for_roc = FALSE
  )
}


#' Resolve a workflow preset and optional stage overrides.
resolve_stage_options <- function(workflow, run_mode, stage_options = NULL) {
  stages <- empty_stage_options()

  #' Validate the condition.
  if (workflow == "full") {
    stages$feature_selection <- TRUE
    stages$training <- TRUE
    stages$evaluation <- TRUE
    stages$ml_roc_plots <- TRUE
    stages$classical_ml_roc <- TRUE
    stages$density_plots <- run_mode == "pilot"
    stages$reuse_evaluation_for_roc <- TRUE
  }

  if (workflow == "training_only") {
    stages$training <- TRUE
  }
  if (workflow == "roc_only") {
    stages$classical_ml_roc <- TRUE
  }
  if (is.null(stage_options)) {
    return(stages)
  }
  if (!is.list(stage_options) || is.null(names(stage_options))) {
    stop("stage_options must be a named list from ml_stage_options().")
  }
  unknown <- setdiff(names(stage_options), names(stages))
  if (length(unknown) > 0L) {
    stop("Unknown stage option(s): ", paste(unknown, collapse = ", "))
  }
  for (name in names(stage_options)) {
    value <- stage_options[[name]]
    if (is.null(value)) {
      next
    }
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("Stage option '", name, "' must be TRUE, FALSE, or NULL.")
    }
    stages[[name]] <- value
  }
  stages
}


#' Return the original pilot or full simulation settings.
make_run_config <- function(run_mode) {
  switch(
    run_mode,
    pilot = list(
      screen_num_sim = 100L,
      vip_num_sim = 50L,
      vip_n_iter_eval = 100L,
      final_num_sim = 100L,
      final_n_iter_eval = 200L,
      n_permutations = 10L,
      n_sim_roc = 200L,
      k_folds = 3L,
      cv_repeats = 1L
    ),
    full = list(
      screen_num_sim = 500L,
      vip_num_sim = 500L,
      vip_n_iter_eval = 1000L,
      final_num_sim = 500L,
      final_n_iter_eval = 10000L,
      n_permutations = 20L,
      n_sim_roc = 1000L,
      k_folds = 10L,
      cv_repeats = 3L
    )
  )
}


#' =============================================================================
#' SECTION 2: FILE AND INPUT HELPERS
#' =============================================================================

#' Browse for one required input file.
browse_for_file <- function(prompt, extensions) {
  if (!interactive()) {
    stop(
      "Browsing requires an interactive R session. ",
      "Use file_mode = 'paths' for Rscript or a cluster."
    )
  }

  message("\n", prompt)
  file <- file.choose(new = FALSE)

  if (!file.exists(file)) {
    stop("The selected file does not exist.")
  }

  extension <- tolower(tools::file_ext(file))
  extensions <- tolower(sub("^\\.", "", extensions))
  
  if (!extension %in% extensions) {
    stop("Expected file type: ", paste0(".", extensions, collapse = ", "))
  }
  normalizePath(file, winslash = "/", mustWork = TRUE)
}


#' Resolve a required input from a path or browser.
resolve_input_file <- function(file_mode, path, prompt, extensions) {
  if (file_mode == "browse") {
    return(browse_for_file(prompt, extensions))
  }
  if (is.null(path) || length(path) != 1L || !nzchar(path)) {
    stop("A path is required when file_mode = 'paths'.\n", prompt)
  }

  path <- path.expand(path)

  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  extension <- tolower(tools::file_ext(path))
  extensions <- tolower(sub("^\\.", "", extensions))

  if (!extension %in% extensions) {
    stop("Expected file type: ", paste0(".", extensions, collapse = ", "))
  }

  normalizePath(path, winslash = "/", mustWork = TRUE)
}


#' Resolve an output file, creating its parent directory.
resolve_output_file <- function(path, default_name, output_dir) {
  file <- if (is.null(path) || !nzchar(path)) {
    file.path(output_dir, default_name)
  } else {
    path.expand(path)
  }

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(dirname(file))) {
    stop("Could not create output directory: ", dirname(file))
  }

  normalizePath(file, winslash = "/", mustWork = FALSE)
}


#' Create and normalize an output directory.
make_output_directory <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(path)) {
    stop("Could not create output directory: ", path)
  }

  normalizePath(path, winslash = "/", mustWork = TRUE)
}


#' Validate sample sizes.
validate_sample_sizes <- function(sample_sizes) {
  sample_sizes <- suppressWarnings(as.integer(sample_sizes))
  sample_sizes <- sort(unique(sample_sizes[is.finite(sample_sizes)]))

  if (length(sample_sizes) == 0L || any(sample_sizes < 8L)) {
    stop("sample_sizes must contain integer values of at least 8.")
  }

  sample_sizes
}


#' Validate model codes.
validate_model_names <- function(models, argument_name) {
  allowed <- c("LR", "RF", "ANN", "GBM", "SVM", "KNN")
  models <- unique(toupper(trimws(as.character(models))))
  models <- models[nzchar(models)]
  invalid <- setdiff(models, allowed)

  if (length(models) == 0L) {
    stop(argument_name, " must contain at least one model.")
  }

  #' Validate the condition.
  if (length(invalid) > 0L) {
    stop("Unsupported model(s) in ", argument_name, ": ",
      paste(invalid, collapse = ", ")
    )
  }
  models
}


#' Source the ML framework into its own environment.
load_ml_framework <- function(framework_file) {
  ml_env <- new.env(parent = globalenv())
  sys.source(framework_file, envir = ml_env, chdir = TRUE)

  #' Calculate required.
  required <- c("run_normality_framework", "screen_correlated_features_across_n",
                  "make_paired_specs", "plot_density_pairs" )

  #' Calculate missing.
  missing <- required[!vapply(required, exists, logical(1), envir = ml_env, mode = "function", inherits = FALSE)]

  if (length(missing) > 0L) {
    stop("The ML framework is missing: ", paste(missing, collapse = ", "))
  }
  ml_env
}


#' Load one named object from an RData file.
load_rdata_object <- function(file, object_name) {
  env <- new.env(parent = emptyenv())
  loaded <- load(file, envir = env)

  if (!object_name %in% loaded || !exists(object_name, envir = env, inherits = FALSE)) {
    stop("The file does not contain an object named '", object_name, "'.")
  }
  get(object_name, envir = env, inherits = FALSE)
}


#' Validate model bundles for selected sample sizes and methods.
validate_model_file <- function( model_file, sample_sizes, required_models = character(0)) {
  trained_models <- load_rdata_object(model_file, "trained_models")

  if (!is.list(trained_models) || length(trained_models) == 0L) {
    stop("trained_models must be a nonempty list.")
  }

  #' Validate the condition.
  if (is.null(names(trained_models)) || any(!nzchar(names(trained_models)))) {
    bundle_sizes <- vapply( trained_models, function(bundle) as.numeric(bundle$n %||% NA_real_),
      numeric(1)
    )

    if (any(!is.finite(bundle_sizes))) {
      stop("Could not identify sample sizes in trained_models.")
    }
    names(trained_models) <- as.character(bundle_sizes)
  }

  keys <- as.character(sample_sizes)
  missing_sizes <- setdiff(keys, names(trained_models))
  if (length(missing_sizes) > 0L) {
    stop("Missing trained models for n = ", paste(missing_sizes, collapse = ", "))
  }

  bundles <- trained_models[keys]
  #' Calculate invalid sizes.
  invalid_sizes <- keys[!vapply( bundles, function(bundle) {
        is.list(bundle) && !is.null(bundle$models) && length(bundle$models) > 0L && 
        !is.null(bundle$prep_obj) && !is.null(bundle$feature_names)
      },
      logical(1)
    )
  ]

  if (length(invalid_sizes) > 0L) {
    stop("Invalid model bundle for n = ", paste(invalid_sizes, collapse = ", "))
  }

  common_models <- Reduce(
    intersect,
    lapply(bundles, function(bundle) names(bundle$models))
  )

  missing_models <- setdiff(required_models, common_models)
  #' Validate the condition.
  if (length(missing_models) > 0L) {
    stop(
      "Unavailable requested model(s): ", paste(missing_models, collapse = ", "),
      "\nAvailable at every selected n: ", paste(common_models, collapse = ", ")
    )
  }
  invisible(common_models)
}


#' Validate saved evaluation results.
validate_evaluation_file <- function(evaluation_file, sample_sizes) {
  results <- load_rdata_object(evaluation_file, "eval_results_by_n")
  missing_sizes <- setdiff(as.character(sample_sizes), names(results))
  
  if (length(missing_sizes) > 0L) {
    stop("Missing evaluation results for n = ", paste(missing_sizes, collapse = ", "))
  }
  invisible(TRUE)
}


#' Read a final-feature CSV.
read_feature_file <- function(feature_file, candidate_features) {
  feature_table <- utils::read.csv(feature_file, stringsAsFactors = FALSE)

  if (!"Feature" %in% names(feature_table)) {
    stop("The final-feature CSV must contain a column named Feature.")
  }

  features <- unique(trimws(as.character(feature_table$Feature)))
  features <- features[nzchar(features)]

  if (length(features) == 0L) {
    stop("The final-feature CSV contains no feature names.")
  }
  unsupported <- setdiff(features, candidate_features)
  if (length(unsupported) > 0L) {
    stop("Unsupported feature(s): ", paste(unsupported, collapse = ", "))
  }
  features
}


#' Set reproducible random-number controls.
set_workflow_seed <- function(seed, rng_version = "4.2.0") {
  RNGversion(rng_version)

  #' Continue the calculation.
  set.seed(
    as.integer(seed),
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  invisible(seed)
}


#' Return conservative cross-platform core settings.
make_core_settings <- function(sample_sizes) {
  available <- parallel::detectCores(logical = TRUE)

  if (!is.finite(available) || available < 1L) {
    available <- 1L
  }

  available <- as.integer(available)
  #' Assemble the result.
  list(
    available = available,
    cv = min(4L, max(1L, available - 2L)),
    training = 1L,
    evaluation = min(2L, max(1L, available - 1L)),
    vip = min(2L, max(1L, available - 1L)),
    roc = min(length(sample_sizes), min(2L, max(1L, available - 1L)))
  )
}


#' Create a compact sample-size label.
sample_size_tag <- function(sample_sizes) {
  paste0("n", paste(sample_sizes, collapse = "_"))
}


#' Return the right value unless it is NULL.
`%||%` <- function(left, right) {
  if (is.null(left)) right else left
}


#' =============================================================================
#' SECTION 3: STUDY SPECIFICATIONS
#' =============================================================================
#' 
#' Return the fixed evaluation distributions.
evaluation_specs <- function() {
  list(
    list(dist = "gumbel", par = c(0, 1)),
    list(dist = "beta", par = c(8, 2)),
    list(dist = "chi_square", par = 3),
    list(dist = "laplace", par = c(0, 4)),
    list(dist = "uniform", par = c(0, 1)),
    list(dist = "cauchy", par = c(0, 1)),
    list(dist = "contaminated", par = c(0.75, 0, 1, 5)),
    list(dist = "t", par = 3)
  )
}

#' Return the out-of-bag-safe training alternatives.
training_specs <- function() {
  list(
    #' Right-skewed alternatives.
    list(dist = "exponential", par = 3),
    list(dist = "gamma", par = c(3, 1)),
    list(dist = "weibull", par = c(1.5, 1)),
    list(dist = "lognormal", par = c(0, 0.5)),
    list(dist = "lognormal", par = c(0, 1)),
    list(dist = "f", par = c(6, 15)),
    list(dist = "pareto", par = c(3, 1)),

    #' Gumbel-like alternatives.
    list(dist = "log_gamma", par = c(0.75, 1)),
    list(dist = "log_gamma", par = c(1.25, 1)),
    list(dist = "log_gamma", par = c(1.50, 1)),

    #' Left-skewed alternatives.
    list(dist = "reflected_exponential", par = 3),
    list(dist = "reflected_gamma", par = c(3, 1)),
    list(dist = "reflected_weibull", par = c(1.5, 1)),
    list(dist = "reflected_lognormal", par = c(0, 0.5)),
    list(dist = "reflected_lognormal", par = c(0, 1)),

    #' Bounded alternatives.
    list(dist = "beta", par = c(2, 8)),
    list(dist = "beta", par = c(8, 2)),
    list(dist = "beta", par = c(2, 2)),
    list(dist = "uniform", par = c(0, 1)),

    #' Symmetric alternatives.
    list(dist = "laplace", par = c(0, 4)),
    list(dist = "logistic", par = c(0, 1)),
    list(dist = "t", par = 3),
    list(dist = "t", par = 5)
  )
}


#' Return distributions for held-out permutation VIP.
#'
#' These specifications match the training design. The framework generates
#' new samples, so the VIP data remain independent of model-fitting samples.
vip_validation_specs <- function() {
  training_specs()
}

#' Return the original candidate-feature set.
candidate_features <- function() {
  c(
    "zp_stat", "qq_corr", "sw_stat", "vasicek", "ad_stat",
    "qq_cubic_coef", "negentropy", "l_kurtosis", "mean_log_spacing",
    "moment_ratio_distance", "ecf_re_t3", "mean_abs_dev_ratio",
    "ecf_im_t1", "dewet_venter", "Studentized_Range", "skewness",
    "henze_zirkler", "jb_stat", "qq_resid_sd", "l_skewness",
    "ryan_joiner", "ecf_im_t3", "epps_pulley", "iqr_sd_ratio",
    "mad_sd_ratio", "ecf_re_t1", "ecf_re_t2", "ecf_im_t2",
    "renyi_a2.0", "tail_asymmetry", "gini_absdev_med", "rms_sd_ratio",
    "greenwood", "skew_log_spacings", "biweight_qq_corr", "geary_stat",
    "bonett_seier", "bowman_shenton", "entropy_sd_gap", "energy_normal_stat"
  )
}


#' Verify that out-of-bag distributions are reserved for final evaluation.
validate_oob_split <- function(train_specs, vip_specs, eval_specs) {
  oob <- c("gumbel", "chi_square", "contaminated", "cauchy")

#' ------------------------------------------------------------
#' Extract normalized distribution names from specification lists.
#' ------------------------------------------------------------
  get_names <- function(specs) {
    tolower(
      vapply(specs, `[[`, character(1),"dist")
    )
  }

  #' Calculate development names.
  development_names <- unique(
    c(get_names(train_specs), get_names(vip_specs))
  )

  eval_names <- get_names(eval_specs)
  leakage <- intersect(development_names, oob)
  missing <- setdiff(oob, eval_names)

  #' Validate the condition.
  if (length(leakage) > 0L) {
    stop("Out-of-bag distribution(s) found in model development: ",
      paste(leakage, collapse = ", ")
    )
  }

  #' Validate the condition.
  if (length(missing) > 0L) {
    stop("Out-of-bag distribution(s) missing from final evaluation: ",
      paste(missing, collapse = ", ")
    )
  }

  invisible(TRUE)
}


#' =============================================================================
#' SECTION 4: PLAN VALIDATION AND DISPLAY
#' =============================================================================

#' Validate dependencies among selected stages.
validate_stage_plan <- function(stages) {
  selected <- c(
    stages$feature_selection,
    stages$training,
    stages$evaluation,
    stages$ml_roc_plots,
    stages$permutation_vip,
    stages$classical_ml_roc,
    stages$density_plots
  )

  if (!any(selected)) {
    stop("No workflow stages were selected.")
  }

  if (stages$ml_roc_plots && !stages$evaluation && !stages$load_evaluation) {
    stop("ML ROC plots require evaluation = TRUE or load_evaluation = TRUE.")
  }

  if (stages$permutation_vip && !stages$evaluation && !stages$load_evaluation) {
    stop("Permutation VIP requires evaluation = TRUE or load_evaluation = TRUE.")
  }

  #' Validate the condition.
  if (stages$classical_ml_roc && stages$reuse_evaluation_for_roc && 
      !stages$evaluation && !stages$load_evaluation) {
    stop("Reusing evaluation predictions requires evaluation = TRUE ", "or load_evaluation = TRUE.")
  }
  invisible(TRUE)
}


#' Print the resolved workflow before execution.
print_workflow_plan <- function(
    workflow,
    run_mode,
    file_mode,
    stages,
    sample_sizes,
    models_to_train,
    models_for_vip,
    models_for_roc,
    train_specs,
    vip_specs,
    eval_specs,
    framework_file,
    model_file,
    evaluation_file,
    feature_file,
    output_dir
) {
  cat("\n", strrep("=", 76), "\n", sep = "")
  cat("RESOLVED ML WORKFLOW\n")
  cat(strrep("=", 76), "\n", sep = "")
  cat("Workflow           :", workflow, "\n")
  cat("Run mode           :", run_mode, "\n")
  cat("File mode          :", file_mode, "\n")
  cat("Sample sizes       :", paste(sample_sizes, collapse = ", "), "\n")
  cat("Training specs     :", length(train_specs), "\n")
  cat("VIP validation     :", length(vip_specs), "training-set specs\n")
  cat("Final evaluation   :", length(eval_specs), "reserved specs\n")
  cat("Models to train    :",
    if (stages$training || stages$feature_selection) {
      paste(models_to_train, collapse = ", ")
    } else {
      "not used"
    },
    "\n"
  )
  cat("Models for VIP     :",
    if (stages$feature_selection || stages$permutation_vip) {
      paste(models_for_vip, collapse = ", ")
    } else {
      "not used"
    },
    "\n"
  )
  cat(
    "Models for ROC     :",
    if (stages$classical_ml_roc) {
      paste(models_for_roc, collapse = ", ")
    } else {
      "not used"
    },
    "\n"
  )
  cat("Framework file     :", framework_file, "\n")

  if (!is.null(model_file)) {
    cat("Model file         :", model_file, "\n")
  }

  if (!is.null(evaluation_file)) {
    cat("Evaluation file    :", evaluation_file, "\n")
  }

  if (!is.null(feature_file)) {
    cat("Feature file       :", feature_file, "\n")
  }

  cat("Output directory   :", output_dir, "\n")
  cat("\nStages\n")

  for (name in names(stages)) {
    cat("  ", format(name, width = 27), ": ", stages[[name]], "\n", sep = "")
  }

  cat(strrep("=", 76), "\n\n", sep = "")
}


#' =============================================================================
#' SECTION 5: MAIN RUN FUNCTION
#' =============================================================================

#' Run the selected ML workflow.
run_ml_framework <- function(
    workflow = c("full", "training_only", "roc_only", "custom"),
    run_mode = c("pilot", "full"),
    file_mode = c("browse", "paths"),
    files = ml_file_paths(),
    sample_sizes = c(10, 50),
    models_to_train = c("LR", "RF", "ANN", "GBM", "SVM", "KNN"),
    models_for_vip = c("RF", "GBM", "ANN", "SVM"),
    models_for_roc = c("RF", "GBM", "ANN", "SVM"),
    classical_tests = c("SW", "AD", "JB"),
    stage_options = NULL,
    feature_source = c("auto", "file", "vector"),
    final_features = character(0),
    alpha_grid = seq(0, 1, by = 0.05),
    n_sim_roc = NULL,
    classification_threshold = 0.50,
    majority_threshold = 0.50,
    test_colors = c(SW = "#CC79A7", AD = "#FF7F00", JB = "#56B4E9"),
    ml_colors = c(LR = "black", RF = "red", GBM = "blue", ANN = "forestgreen", SVM = "green", KNN = "purple", MajorityVote = "brown"),
    output_parent = ".",
    output_name = NULL,
    seed = 12345L,
    rng_version = "4.2.0",
    show_progress = TRUE,
    dry_run = FALSE
) {
  workflow <- match.arg(workflow)
  run_mode <- match.arg(run_mode)
  file_mode <- match.arg(file_mode)
  feature_source <- match.arg(feature_source)

  if (!is.list(files)) {
    stop("files must be created with ml_file_paths().")
  }

  #' Calculate sample sizes.
  sample_sizes <- validate_sample_sizes(sample_sizes)
  models_to_train <- validate_model_names(models_to_train, "models_to_train")
  models_for_vip <- validate_model_names(models_for_vip, "models_for_vip")
  models_for_roc <- validate_model_names(models_for_roc, "models_for_roc")
  classical_tests <- unique(toupper(trimws(as.character(classical_tests))))
  classical_tests <- classical_tests[nzchar(classical_tests)]

  if (length(classical_tests) == 0L) {
    stop("classical_tests must contain at least one test.")
  }

  stages <- resolve_stage_options(workflow, run_mode, stage_options)
  validate_stage_plan(stages)
  run_config <- make_run_config(run_mode)

  if (is.null(n_sim_roc)) {
    n_sim_roc <- run_config$n_sim_roc
  }

  if (!is.numeric(n_sim_roc) || length(n_sim_roc) != 1L || !is.finite(n_sim_roc) || n_sim_roc < 1L) {
    stop("n_sim_roc must be a positive integer.")
  }

  n_sim_roc <- as.integer(n_sim_roc)
  output_parent <- make_output_directory(output_parent)

  if (is.null(output_name) || !nzchar(output_name)) {
    output_name <- paste("ml", workflow, run_mode, sample_size_tag(sample_sizes), sep = "_")
  }

  output_dir <- make_output_directory(file.path(output_parent, output_name))

  #' Calculate framework file.
  framework_file <- resolve_input_file(
    file_mode = file_mode,
    path = files$framework,
    prompt = paste0("Select the ML framework function file.\n", "Choose ML_framework_func.R."),
    extensions = "R"
  )

  ml_env <- load_ml_framework(framework_file)

  #' The framework sets its own seed when sourced; restore this run's seed.
  set_workflow_seed(seed, rng_version)

  cores <- make_core_settings(sample_sizes)
  #' Define separate roles for development and final evaluation.
  train_specs <- training_specs()
  vip_specs <- vip_validation_specs()
  eval_specs <- evaluation_specs()

  feature_candidates <- candidate_features()

  validate_oob_split(
    train_specs = train_specs,
    vip_specs = vip_specs,
    eval_specs = eval_specs
  )

  #' Restrict stage-specific methods to newly trained models.
  if (stages$feature_selection || stages$training) {
    models_for_vip <- intersect(models_for_vip, models_to_train)
  }

  if (stages$training) {
    models_for_roc <- intersect(models_for_roc, models_to_train)
  }

  if (stages$feature_selection && length(models_for_vip) == 0L) {
    stop("Feature selection requires at least one model in models_for_vip.")
  }

  #' Validate the condition.
  if ((stages$feature_selection || stages$permutation_vip) &&
      !"RF" %in% models_for_vip) {
    stop("Permutation VIP currently requires RF because the framework ",
      "creates an RF-specific VIP plot."
    )
  }

  if (stages$classical_ml_roc && length(models_for_roc) == 0L) {
    stop("Classical-versus-ML ROC requires at least one model.")
  }

  #' Calculate model stages.
  model_stages <- any(
    c(
      stages$training,
      stages$evaluation,
      stages$ml_roc_plots,
      stages$permutation_vip,
      stages$classical_ml_roc
    )
  )

  #' Calculate existing model required.
  existing_model_required <- !stages$training &&
    any(
      c(
        stages$evaluation,
        stages$ml_roc_plots,
        stages$permutation_vip,
        stages$classical_ml_roc
      )
    )

  model_file <- NULL

  #' Validate the condition.
  if (existing_model_required) {
    model_file <- resolve_input_file(
      file_mode = file_mode,
      path = files$trained_models,
      prompt = paste0("Select the trained-model RData file.\n", "It must contain trained_models."),
      extensions = c("RData", "rda")
    )

    #' Calculate required models.
    required_models <- unique(
      c(
        if (stages$classical_ml_roc) models_for_roc else character(0),
        if (stages$permutation_vip) models_for_vip else character(0)
      )
    )

    validate_model_file(model_file, sample_sizes, required_models)

  #' Continue the calculation.
  } else if (stages$training) {
    model_file <- resolve_output_file(
      path = files$trained_models,
      default_name = paste0("trained_models_", sample_size_tag(sample_sizes), ".RData"
      ),
      output_dir = output_dir
    )
  }

  evaluation_file <- NULL

  #' Validate the condition.
  if (stages$load_evaluation) {
    evaluation_file <- resolve_input_file(
      file_mode = file_mode,
      path = files$evaluation_results,
      prompt = paste0("Select the saved evaluation-results RData file.\n",
        "It must contain eval_results_by_n."
      ),
      extensions = c("RData", "rda")
    )

    validate_evaluation_file(evaluation_file, sample_sizes)

  #' Continue the calculation.
  } else if (stages$evaluation || stages$permutation_vip) {
    evaluation_file <- resolve_output_file(
      path = files$evaluation_results,
      default_name = "eval_results.RData",
      output_dir = output_dir
    )
  }

  feature_file <- NULL
  selected_features <- character(0)

  if (stages$training && !stages$feature_selection) {
    source <- feature_source

    if (source == "auto") {
      source <- if (length(final_features) > 0L) "vector" else "file"
    }

    #' Validate the condition.
    if (source == "file") {
      feature_file <- resolve_input_file(
        file_mode = file_mode,
        path = files$final_features,
        prompt = paste0("Select the final-feature CSV file.\n",
          "It must contain a column named Feature."
        ),
        extensions = "csv"
      )

      selected_features <- read_feature_file(feature_file, feature_candidates)

    } else {
      selected_features <- unique(trimws(as.character(final_features)))
      selected_features <- selected_features[nzchar(selected_features)]

      if (length(selected_features) == 0L) {
        stop("feature_source = 'vector' requires final_features.")
      }

      unsupported <- setdiff(selected_features, feature_candidates)
      
      if (length(unsupported) > 0L) {
        stop("Unsupported feature(s): ", paste(unsupported, collapse = ", "))
      }
    }
  }

  #' Continue the calculation.
  print_workflow_plan(
    workflow = workflow,
    run_mode = run_mode,
    file_mode = file_mode,
    stages = stages,
    sample_sizes = sample_sizes,
    models_to_train = models_to_train,
    models_for_vip = models_for_vip,
    models_for_roc = models_for_roc,
    train_specs = train_specs,
    vip_specs = vip_specs,
    eval_specs = eval_specs,
    framework_file = framework_file,
    model_file = model_file,
    evaluation_file = evaluation_file,
    feature_file = feature_file,
    output_dir = output_dir
  )

  #' Validate the condition.
  if (isTRUE(dry_run)) {
    return(
      invisible(
        list(
          stages = stages,
          sample_sizes = sample_sizes,
          distribution_roles = list(
            training = train_specs,
            vip_validation = vip_specs,
            final_evaluation = eval_specs
          ),
          files = list(
            framework = framework_file,
            trained_models = model_file,
            evaluation_results = evaluation_file,
            final_features = feature_file
          ),
          output_dir = output_dir
        )
      )
    )
  }

  start_time <- Sys.time()
  discovery_result <- NULL
  discovery_dir <- NULL
  final_result <- NULL
  density_file <- NULL

  #' Run density plots independently of model stages.
  if (stages$density_plots) {
    density_file <- ml_env$plot_density_pairs(
      paired_specs = ml_env$make_paired_specs(eval_specs),
      output_dir = output_dir,
      file_name = "density_pairs.pdf"
    )
  }

  #' Run screening and permutation VIP feature selection.
  if (stages$feature_selection) {
    discovery_dir <- make_output_directory(
      file.path(output_dir, "stage1_vip_discovery")
    )

    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("FEATURE SELECTION: CORRELATION SCREENING\n")
    cat(strrep("=", 70), "\n", sep = "")

    #' Calculate screened features.
    screened_features <- ml_env$screen_correlated_features_across_n(
      sample_sizes = sample_sizes,
      paired_specs = ml_env$make_paired_specs(train_specs),
      candidate_features = feature_candidates,
      num_sim = run_config$screen_num_sim,
      cutoff = 0.98,
      standardize_sample = FALSE,
      center_by = NULL,
      impute_method = "median",
      cor_method = "spearman"
    )

    utils::write.csv(
      data.frame(Feature = screened_features),
      file.path(discovery_dir, "screened_features_pre_vip.csv"),
      row.names = FALSE
    )

    cat("\n", strrep("=", 70), "\n", sep = "")
    cat("FEATURE SELECTION: PERMUTATION VIP\n")
    cat(strrep("=", 70), "\n", sep = "")

    discovery_result <- ml_env$run_normality_framework(
      sample_sizes = sample_sizes,

      #' Use the configured pilot or full iteration counts without reduction.
      num_sim = run_config$vip_num_sim,
      n_iter_eval = run_config$vip_n_iter_eval,

      #' Fit discovery models on the training distributions.
      train_alt_specs = train_specs,

      #' Generate independent VIP holdout samples from the same specifications.
      eval_alt_specs = vip_specs,

      #' ROC analysis is disabled during feature selection.
      roc_alt_specs = eval_specs,

      n_cores_train = cores$training,
      n_cores_eval = cores$evaluation,
      n_cores_vip = cores$vip,
      n_cores_roc = cores$roc,
      n_cores_cv = cores$cv,

      train_models = TRUE,
      models_to_train = models_to_train,
      k_folds = run_config$k_folds,
      cv_repeats = run_config$cv_repeats,
      cv_verbose = FALSE,

      feature_set = screened_features,
      threshold = classification_threshold,
      majority_threshold = majority_threshold,
      standardize_sample = FALSE,
      center_by = NULL,

      #' Store independent holdout features for permutation AUC importance.
      run_evaluation = TRUE,
      save_eval_results = TRUE,
      store_holdout_features = TRUE,
      run_permutation_vip = TRUE,

      #' Final evaluation and ROC outputs are produced only in the final run.
      run_ml_roc_plots = FALSE,
      run_roc_analysis = FALSE,
      reuse_roc_evaluation_predictions = FALSE,
      make_density_plots = FALSE,
      show_progress = show_progress,

      #' Continue the calculation.
      n_permutations = run_config$n_permutations,
      models_for_vip = models_for_vip,
      vip_top_n = 30L,
      vip_rank_by = "mean",
      vip_min_pct_positive = 0.60,
      min_mean_for_negative_min = 0.005,
      vip_require_positive_mean = TRUE,
      vip_plot_models = "RF",
      vip_plot_sample_sizes = sample_sizes,

      #' Continue the calculation.
      classical_tests = classical_tests,
      ml_methods_roc = models_for_roc,
      alpha_grid = alpha_grid,
      n_sim_roc = n_sim_roc,
      test_colors = test_colors,
      ml_colors = ml_colors,
      output_dir = discovery_dir
    )

    selected_features <- discovery_result$final_selected_features

    if (is.null(selected_features) || length(selected_features) == 0L) {
      selected_file <- file.path(discovery_dir, "final_selected_features.csv")

      if (!file.exists(selected_file)) {
        stop("Could not recover the selected feature set.")
      }

      selected_features <- utils::read.csv( selected_file, stringsAsFactors = FALSE)$Feature
    }

    selected_features <- unique(as.character(selected_features))

    #' Continue the calculation.
    utils::write.csv(
      data.frame(Feature = selected_features),
      file.path(output_dir, "final_selected_features_for_refit.csv"),
      row.names = FALSE
    )
  }

  #' Release the large discovery object before the final run.
  if (model_stages && !is.null(discovery_result)) {
    discovery_result <- NULL
    invisible(gc())
  }

  #' Refit selected models and evaluate only on the reserved distributions.
  if (model_stages) {
    final_dir <- make_output_directory(file.path(output_dir, "final_run"))

    #' Validate the condition.
    if (stages$training && !stages$feature_selection) {
      utils::write.csv(data.frame(Feature = selected_features),
        file.path(output_dir, "final_selected_features_for_refit.csv"),
        row.names = FALSE
      )
    }

    final_result <- ml_env$run_normality_framework(
      sample_sizes = sample_sizes,
      num_sim = run_config$final_num_sim,
      n_iter_eval = run_config$final_n_iter_eval,
      #' Refit on training distributions and reserve eval_specs for testing.
      train_alt_specs = train_specs,
      eval_alt_specs = eval_specs,
      roc_alt_specs = eval_specs,
      n_cores_train = cores$training,
      n_cores_eval = cores$evaluation,
      n_cores_vip = cores$vip,
      n_cores_roc = cores$roc,
      n_cores_cv = cores$cv,
      train_models = stages$training,
      model_save_path = model_file,
      models_to_train = models_to_train,
      k_folds = run_config$k_folds,
      cv_repeats = run_config$cv_repeats,
      cv_verbose = FALSE,
      feature_set = if (stages$training) selected_features else NULL,
      threshold = classification_threshold,
      majority_threshold = majority_threshold,
      standardize_sample = FALSE,
      center_by = NULL,
      run_evaluation = stages$evaluation,
      eval_save_path = evaluation_file,
      load_eval_if_exists = stages$load_evaluation,
      save_eval_results = stages$evaluation || stages$permutation_vip,
      store_holdout_features = stages$permutation_vip,
      make_density_plots = FALSE,
      run_ml_roc_plots = stages$ml_roc_plots,
      run_permutation_vip = stages$permutation_vip,
      run_roc_analysis = stages$classical_ml_roc,
      reuse_roc_evaluation_predictions = stages$reuse_evaluation_for_roc,
      show_progress = show_progress,
      n_permutations = run_config$n_permutations,
      models_for_vip = models_for_vip,
      vip_top_n = 30L,
      vip_rank_by = "mean",
      vip_min_pct_positive = 0.60,
      min_mean_for_negative_min = 0.005,
      vip_require_positive_mean = TRUE,
      vip_plot_models = "RF",
      vip_plot_sample_sizes = sample_sizes,
      classical_tests = classical_tests,
      ml_methods_roc = models_for_roc,
      alpha_grid = alpha_grid,
      n_sim_roc = n_sim_roc,
      test_colors = test_colors,
      ml_colors = ml_colors,
      output_dir = final_dir
    )
  }

  #' Calculate elapsed time.
  elapsed_time <- if (
    exists("format_elapsed_time", envir = ml_env, mode = "function",inherits = FALSE)
  ) {
    ml_env$format_elapsed_time(start_time)
  } else {
    paste0(round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2),
      " minutes"
    )
  }

  #' Summarize the results.
  summary <- list(
    workflow = workflow,
    run_mode = run_mode,
    file_mode = file_mode,
    stages = stages,
    sample_sizes = sample_sizes,
    distribution_roles = list(
      training = train_specs,
      vip_validation = vip_specs,
      final_evaluation = eval_specs
    ),
    models_to_train = models_to_train,
    models_for_vip = models_for_vip,
    models_for_roc = models_for_roc,
    classical_tests = classical_tests,
    alpha_grid = alpha_grid,
    n_sim_roc = n_sim_roc,
    classification_threshold = classification_threshold,
    majority_threshold = majority_threshold,
    test_colors = test_colors,
    ml_colors = ml_colors,
    files = list(
      framework = framework_file,
      trained_models = model_file,
      evaluation_results = evaluation_file,
      final_features = feature_file
    ),
    selected_features = selected_features,
    discovery_dir = discovery_dir,
    cores = cores,
    output_dir = output_dir,
    elapsed_time = elapsed_time,
    completed_at = Sys.time()
  )

  summary_file <- file.path(output_dir, "ml_workflow_summary.rds")
  saveRDS(summary, summary_file, compress = "gzip")

  cat("\nWorkflow completed successfully.\n")
  cat("Elapsed time:", elapsed_time, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Workflow summary:", summary_file, "\n")

  #' Continue the calculation.
  invisible(
    list(
      discovery_result = discovery_result,
      discovery_dir = discovery_dir,
      final_result = final_result,
      density_file = density_file,
      summary = summary,
      files = list(
        summary = summary_file,
        trained_models = model_file,
        evaluation_results = evaluation_file
      )
    )
  )
}


#' =============================================================================
#' SECTION 6: EXAMPLE RUNS
#' =============================================================================
#'
#' Copy one example to the console and edit it. Nothing below runs automatically.
#' =============================================================================

if (FALSE) {

  #' Example 1: Browse for files and run a quick full pilot at n = 10 and 50.
  pilot_full_browse <- run_ml_framework(
    workflow = "full",
    run_mode = "pilot",
    file_mode = "browse",
    sample_sizes = c(10, 50),
    models_to_train = c("LR", "RF", "ANN", "GBM", "SVM", "KNN"),
    models_for_vip = c("RF", "GBM", "ANN", "SVM"),
    models_for_roc = c("RF", "GBM", "ANN", "SVM"),
    output_parent = "."
  )


  #' Example 2: Hardcoded paths and a full run for all models at n = 10 and 50.
  full_all_models_paths <- run_ml_framework(
    workflow = "full",
    run_mode = "full",
    file_mode = "paths",
    files = ml_file_paths(
      framework = "/replace/path/ML_framework_func.R"
    ),
    sample_sizes = c(10, 50),
    models_to_train = c("LR", "RF", "ANN", "GBM", "SVM", "KNN"),
    models_for_vip = c("RF", "GBM", "ANN", "SVM"),
    models_for_roc = c("RF", "GBM", "ANN", "SVM"),
    output_parent = "/replace/path/results"
  )


  #' Example 3: Browse and train selected models for all sample sizes.
  training_all_sizes_browse <- run_ml_framework(
    workflow = "training_only",
    run_mode = "full",
    file_mode = "browse",
    sample_sizes = c(10, 20, 30, 40, 50),
    models_to_train = c("RF", "ANN", "GBM", "SVM"),
    feature_source = "file",
    output_parent = "."
  )


  #' Example 4: Hardcoded paths, selected models, and selected sample sizes.
  training_selected_paths <- run_ml_framework(
    workflow = "training_only",
    run_mode = "full",
    file_mode = "paths",
    files = ml_file_paths(
      framework = "/replace/path/ML_framework_func.R",
      trained_models = "/replace/path/trained_models_n10_30_50.RData",
      final_features = "/replace/path/final_selected_features.csv"
    ),
    sample_sizes = c(10, 30, 50),
    models_to_train = c("RF", "SVM"),
    feature_source = "file",
    output_parent = "/replace/path/results"
  )


  #' Example 5: Browse and run only classical-versus-ML ROC comparisons.
  #' This skips training, evaluation, ML-only ROC plots, and VIP.
  roc_only_browse <- run_ml_framework(
    workflow = "roc_only",
    run_mode = "full",
    file_mode = "browse",
    sample_sizes = c(10, 20, 30, 40, 50),
    models_for_roc = c("RF", "GBM", "ANN"),
    classical_tests = c("SW", "AD", "JB"),
    output_parent = "."
  )


  #' Example 6: Hardcoded ROC-only run using saved evaluation predictions.
  roc_only_reuse_paths <- run_ml_framework(
    workflow = "roc_only",
    run_mode = "full",
    file_mode = "paths",
    files = ml_file_paths(
      framework = "/replace/path/ML_framework_func.R",
      trained_models = "/replace/path/trained_models.RData",
      evaluation_results = "/replace/path/eval_results.RData"
    ),
    sample_sizes = c(10, 50),
    models_for_roc = c("RF", "GBM", "ANN"),
    stage_options = ml_stage_options(
      load_evaluation = TRUE,
      reuse_evaluation_for_roc = TRUE
    ),
    output_parent = "/replace/path/results"
  )


  #' Example 7: Skip training; run evaluation and both ROC outputs.
  evaluation_and_roc_browse <- run_ml_framework(
    workflow = "custom",
    run_mode = "pilot",
    file_mode = "browse",
    sample_sizes = c(10, 50),
    models_for_roc = c("RF", "GBM", "ANN"),
    stage_options = ml_stage_options(
      evaluation = TRUE,
      ml_roc_plots = TRUE,
      classical_ml_roc = TRUE,
      reuse_evaluation_for_roc = TRUE
    ),
    output_parent = "."
  )


  #' Example 8: Run only screening and permutation VIP feature selection.
  feature_selection_only <- run_ml_framework(
    workflow = "custom",
    run_mode = "pilot",
    file_mode = "browse",
    sample_sizes = c(10, 50),
    models_to_train = c("RF", "GBM"),
    models_for_vip = c("RF", "GBM"),
    stage_options = ml_stage_options(
      feature_selection = TRUE
    ),
    output_parent = "."
  )


  #' Example 9: Train with a feature vector instead of a CSV.
  training_with_vector <- run_ml_framework(
    workflow = "training_only",
    run_mode = "pilot",
    file_mode = "paths",
    files = ml_file_paths(
      framework = "/replace/path/ML_framework_func.R"
    ),
    sample_sizes = c(10, 50),
    models_to_train = c("RF", "ANN"),
    feature_source = "vector",
    final_features = c(
      "sw_stat",
      "qq_corr",
      "skewness",
      "l_kurtosis",
      "tail_asymmetry"
    ),
    output_parent = "/replace/path/results"
  )


  #' Example 10: Override the full preset and skip selected stages.
  full_without_ml_roc_or_density <- run_ml_framework(
    workflow = "full",
    run_mode = "full",
    file_mode = "browse",
    sample_sizes = c(10, 50),
    models_to_train = c("RF", "GBM", "ANN"),
    models_for_vip = c("RF", "GBM", "ANN"),
    models_for_roc = c("RF", "GBM", "ANN"),
    stage_options = ml_stage_options(
      ml_roc_plots = FALSE,
      density_plots = FALSE
    ),
    output_parent = "."
  )
}
