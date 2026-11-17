test_that("startup and input changes do not invoke simulation", {
  calls <- 0L
  fw <- new.env(parent = baseenv())
  for (name in c("onesample_data", "onesample_parameters", "raw_data", "one_sample_t_test", "sign_test")) fw[[name]] <- function(...) NULL
  fw$run_simulation <- function(...) {
    calls <<- calls + 1L
    list(params = list(phases_run = integer()), objects = list())
  }
  local_mocked_bindings(load_framework = function(...) fw, .package = "pretestsim")
  shiny::testServer(pretestsim:::app_server, {
    messages <- logical()
    session$sendCustomMessage <- function(type, message) {
      if (type == "pretestsim-run-state") messages <<- c(messages, message)
    }
    session$setInputs(run = 0L, family = "one_sample", downstream_alt = "sign",
      approach = "classical", norm_test = "SW", phases = "1", Nsim = 10,
      N_tradeoff = 10, grid_by = .5, sig_by = .5, tol_pos = .005,
      loss_tol = .01, test_alpha = .05, center_by = "median",
      effect_size = ".5", sample_sizes = "10,20", single_n = 10,
      effect_sizes_plot = "0,.5", dist = "exponential", battery = "SW", per_n = FALSE)
    session$flushReact()
    expect_equal(calls, 0L)
    expect_length(messages, 0L)
    session$setInputs(Nsim = 20)
    expect_equal(calls, 0L)
    expect_length(messages, 0L)
    before <- getwd()
    session$setInputs(run = 1L)
    expect_equal(calls, 1L)
    expect_identical(messages, c(TRUE, FALSE))
    expect_identical(getwd(), before)
    session$setInputs(phases = character(), run = 2L)
    expect_equal(calls, 1L)
    expect_identical(tail(messages, 2), c(TRUE, FALSE))
    fw$run_simulation <- function(...) stop("Test failure")
    session$setInputs(phases = "1", run = 3L)
    expect_null(run_state$results)
    expect_match(run_state$log, "Test failure")
    expect_identical(tail(messages, 2), c(TRUE, FALSE))
  })
})

test_that("UI exposes Fisher alone and never binds running to generic busy events", {
  old <- options(sass.cache = FALSE)
  on.exit(options(old), add = TRUE)
  html <- as.character(pretestsim:::app_ui())
  expect_match(html, "Fisher SW \\+ AD")
  expect_false(grepl("Machine learning|ml_model|shiny:busy", html))
  expect_match(html, "pretestsim-run-state")
})

# Exercise the real server and framework, including the reuse path.
test_that("dashboard saves and reuses the previous run's four-column tables", {
  shiny::testServer(pretestsim:::app_server, {
    session$setInputs(
      run = 0L, family = "one_sample", downstream_alt = "sign",
      approach = "classical", norm_test = "SW", phases = c("2", "4"),
      Nsim = 12, N_tradeoff = 12, grid_by = .5, sig_by = .5,
      tol_pos = .005, loss_tol = .01, test_alpha = .05,
      center_by = "median", effect_size = ".5", sample_sizes = "10,20",
      single_n = 10, effect_sizes_plot = "0,.5", dist = "exponential",
      battery = "SW", per_n = TRUE, seed = 12345, reuse_pvalues = FALSE
    )
    session$flushReact()
    expect_null(run_state$results)
    session$setInputs(run = 1L)
    expect_false(is.null(run_state$results))
    first_directory <- run_state$wd
    file <- run_state$results$files$pvalue_tables
    expect_true(file.exists(file))
    saved <- readRDS(file)
    expect_false(run_state$args$save_results)

    framework <- fw_rv()
    framework$generate_data <- function(...) stop("Dashboard unexpectedly generated data")
    session$setInputs(reuse_pvalues = TRUE)
    session$setInputs(run = 2L)
    expect_false(is.null(run_state$results))
    expect_true(run_state$args$reuse_pvalues)
    expect_identical(run_state$wd, first_directory)
    expect_identical(readRDS(file), saved)
    unlink(first_directory, recursive = TRUE)
  })
})
