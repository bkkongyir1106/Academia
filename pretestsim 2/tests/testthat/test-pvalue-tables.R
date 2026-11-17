# Remove elapsed times and timestamps when comparing numerical result objects.
without_runtime <- function(value) {
  if (!is.null(attr(value, "settings"))) {
    attr(value, "settings") <- without_runtime(attr(value, "settings"))
  }
  if (!is.list(value) || is.data.frame(value)) return(value)
  value[c("timing", "timestamp")] <- NULL
  lapply(value, without_runtime)
}

# Run small studies covering each generator and both normality interfaces.
test_that("all phases save four-column tables and reuse without testing again", {
  skip_if_not_installed("nortest")
  skip_if_not_installed("Rfit")

  cases <- list(
    one_sample = c("onesample_data", "onesample_parameters", "raw_data",
      "one_sample_t_test", "sign_test"),
    two_sample = c("two_sample_data", "twosample_parameters", "raw_data",
      "twosample_t_test", "Mann_whitney_U_test"),
    anova = c("anova_gen_data", "anova_parameters", "anova_residuals",
      "one_way_anova", "kruskal_wallis_test"),
    regression = c("reg_data", "reg_parameters", "reg_residuals",
      "simple_linear_reg", "rank_regression"),
    fisher = c("onesample_data", "onesample_parameters", "raw_data",
      "one_sample_t_test", "sign_test")
  )
  function_arguments <- c("gen_data", "get_parameters", "fn_to_get_norm_obj",
    "fn_for_ds_test_1", "fn_for_ds_test_2")

  for (case in names(cases)) {
    fw <- load_framework()
    functions <- mget(cases[[case]], envir = fw)
    names(functions) <- function_arguments
    custom <- case == "fisher"
    config <- if (custom) {
      list(method = "custom", config = list(fn = fw$fisher_combined, label = "Fisher"))
    } else {
      list(method = "classical", config = list(norm_test = "SW"))
    }
    effect <- if (case == "anova") list(H0 = c(0, 0, 0), H1 = c(0, 0, .5)) else .5
    directory <- tempfile(paste0("tables_", case))
    args <- c(functions, list(
      Nsim              = 12L,
      N_tradeoff        = 12L,
      test_type         = case,
      norm_config       = config,
      norm_test         = "SW",
      selected_tests    = if (custom) "Fisher" else "SW",
      effect_size       = effect,
      sample_sizes      = c(10, 20),
      single_n          = 10,
      effect_sizes_plot = c(0, .5),
      threshold_grid    = c(0, .5, 1),
      normality_roc_grid = c(0, .5, 1),
      sig_levels        = c(0, .5, 1),
      per_n_thresholds  = case %in% c("one_sample", "two_sample"),
      phases            = 1:5,
      output_dir        = directory
    ))

    set.seed(901)
    invisible(capture.output(result <- do.call(fw$run_simulation, args)))
    tables <- readRDS(result$files$pvalue_tables)
    expect_length(tables, 5L)
    expect_false("results_rdata" %in% names(result$files))
    expect_length(list.files(directory, pattern = "[.]RData$", recursive = TRUE), 0L)
    expect_length(list.files(directory, pattern = "[.]rds$", recursive = TRUE), 1L)

    for (phase in names(tables)) {
      for (batch in tables[[phase]]) {
        for (hypothesis in c("H0", "H1")) {
          table <- batch[[hypothesis]]
          expect_named(table, c("Normality_Score", "Test1_P_Value",
            "Test2_P_Value", "Adaptive_P_Value"))
          if (phase == "phase1_normality_roc") {
            expect_true(all(is.na(table$Adaptive_P_Value)))
          } else {
            threshold <- attr(table, "pretest_threshold")
            if (length(threshold) == 1L && is.finite(threshold)) {
              expected <- fw$route_saved_pvalues(table, threshold)
              expect_identical(table$Adaptive_P_Value, expected$Adaptive_P_Value)
            }
          }
        }
      }
    }
    if (case == "two_sample") {
      expect_length(tables$phase4_sample_sizes$batch_0001$H0$Normality_Score[[1]], 2L)
    }

    # Any recomputation of samples, normality tests, or fixed tests must fail.
    fw$generate_data <- function(...) stop("Unexpected sample generation")
    fw$generate_tests <- function(...) stop("Unexpected normality testing")
    fw$apply_ds_test <- function(...) stop("Unexpected downstream testing")
    args$reuse_pvalues <- TRUE
    set.seed(901)
    invisible(capture.output(replay <- do.call(fw$run_simulation, args)))
    expect_identical(tables, readRDS(replay$files$pvalue_tables))
    expect_identical(without_runtime(result$objects), without_runtime(replay$objects))

    # Changed seed and simulation counts must not be accepted as the same run.
    set.seed(902)
    expect_error(capture.output(do.call(fw$run_simulation, args)), "Random-number state differs")
    set.seed(901)
    args$Nsim <- 13L
    expect_error(capture.output(do.call(fw$run_simulation, args)), "Saved batch settings differ")
    unlink(directory, recursive = TRUE)
  }
})

# The adaptive column must use the two-sample routing rule at boundary values.
test_that("routing preserves both group scores and the strict threshold boundary", {
  fw <- load_framework()
  table <- data.frame(
    Normality_Score   = I(list(c(.4, .6), c(.4, .05), c(NA_real_, NA_real_))),
    Test1_P_Value    = c(.01, .02, .03),
    Test2_P_Value    = c(.10, .20, .30),
    Adaptive_P_Value = NA_real_
  )
  routed <- fw$route_saved_pvalues(table, .05)
  expect_identical(routed$Adaptive_P_Value, c(.01, .20, .30))
  expect_identical(routed[1:3], table[1:3])
})
