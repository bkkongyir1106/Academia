DIST_CHOICES <- c(
  "exponential", "lognormal", "chi_square", "gamma", "t", "f", "uniform",
  "laplace", "cauchy", "gumbel", "weibull", "beta", "logistic", "pareto",
  "contaminated", "normal"
)

TEST_CODES <- c("SW", "SF", "LF", "KS", "JB", "SKEW", "KURT", "DAP", "AD", "CVM")

FAMILIES <- list(
  onesample_t_vs_sign = list(
    label = "One-sample t-test vs sign test",
    test_type = "onesample_ttest_vs_sign",
    functions = c(
      gen_data = "onesample_data",
      get_parameters = "onesample_parameters",
      fn_to_get_norm_obj = "raw_data",
      fn_for_ds_test_1 = "one_sample_t_test",
      fn_for_ds_test_2 = "sign_test"
    ),
    effect_size = 0.5,
    center_by = "median",
    distributions = c("exponential", "normal"),
    extra = list()
  ),
  twosample_t_vs_mann_whitney = list(
    label = "Two-sample t-test vs Mann-Whitney U",
    test_type = "twosample_ttest_vs_mann_whitney",
    functions = c(
      gen_data = "two_sample_data",
      get_parameters = "twosample_parameters",
      fn_to_get_norm_obj = "raw_data",
      fn_for_ds_test_1 = "twosample_t_test",
      fn_for_ds_test_2 = "Mann_whitney_U_test"
    ),
    effect_size = 0.5,
    center_by = "median",
    distributions = c("exponential", "normal"),
    extra = list()
  ),
  anova_vs_kruskal_wallis = list(
    label = "One-way ANOVA vs Kruskal-Wallis",
    test_type = "anova_vs_kruskal_wallis",
    functions = c(
      gen_data = "anova_gen_data",
      get_parameters = "anova_parameters",
      fn_to_get_norm_obj = "anova_residuals",
      fn_for_ds_test_1 = "one_way_anova",
      fn_for_ds_test_2 = "kruskal_wallis_test"
    ),
    effect_size = "0, 0, 0.5",
    center_by = "median",
    distributions = c("exponential", "normal"),
    extra = list()
  ),
  regression_ols_vs_rank = list(
    label = "OLS regression vs rank-based regression",
    test_type = "regression_ols_vs_rank",
    functions = c(
      gen_data = "reg_data",
      get_parameters = "reg_parameters",
      fn_to_get_norm_obj = "reg_residuals",
      fn_for_ds_test_1 = "simple_linear_reg",
      fn_for_ds_test_2 = "rank_regression"
    ),
    effect_size = 0.5,
    center_by = "mean",
    distributions = c("exponential", "normal"),
    extra = list(x_dist = "exponential", x_par = 1, error_par = NULL)
  )
)

parse_nums <- function(txt, default = numeric(0)) {
  if (is.null(txt) || !nzchar(trimws(txt))) return(default)
  values <- suppressWarnings(as.numeric(strsplit(txt, "[,;[:space:]]+")[[1]]))
  values[is.finite(values)]
}

parse_threshold <- function(txt) {
  txt <- trimws(txt %||% "")
  if (!nzchar(txt)) return(NULL)
  if (grepl("=", txt, fixed = TRUE)) {
    parts <- strsplit(txt, ",")[[1]]
    names <- trimws(sub("=.*", "", parts))
    values <- suppressWarnings(as.numeric(sub(".*=", "", parts)))
    return(stats::setNames(values, names))
  }
  suppressWarnings(as.numeric(txt))
}

build_theme <- function() {
  bslib::bs_add_rules(
    bslib::bs_theme(
      version = 5,
      bg = "#0f1419",
      fg = "#e6edf3",
      primary = "#4c8bf5",
      success = "#3fb950",
      "card-bg" = "#161b22",
      "card-border-color" = "#30363d"
    ),
    "
    .navbar { background-color: #21262d !important; border-bottom: 2px solid #4c8bf5 !important; }
    .navbar .nav-link { color: #fff !important; font-weight: 600; border-radius: 7px; margin: 4px 3px; }
    .navbar .nav-link.active { background-color: #4c8bf5 !important; }
    pre, .shiny-text-output { color: #e6edf3; }
    "
  )
}

app_ui <- function(raw_base = .default_raw_base()) {
  bslib::page_navbar(
    title = "Pretest Simulation Studio",
    theme = build_theme(),
    fillable = TRUE,
    sidebar = bslib::sidebar(
      width = 390,
      title = "Setup & Controls",
      bslib::accordion(
        open = c("assets", "approach", "design", "run"),
        bslib::accordion_panel(
          "GitHub assets",
          value = "assets",
          textInput("raw_base", "Raw GitHub asset folder", value = raw_base),
          checkboxInput("refresh_assets", "Refresh cached assets", FALSE),
          actionButton("load_fw", "Load framework", class = "btn-secondary w-100"),
          uiOutput("asset_status")
        ),
        bslib::accordion_panel(
          "Normality approach",
          value = "approach",
          radioButtons(
            "approach",
            NULL,
            c("Classical" = "classical", "Fisher SW + AD" = "fisher", "Machine learning" = "ml"),
            selected = "classical"
          ),
          conditionalPanel(
            "input.approach == 'classical'",
            selectInput("norm_test", "Adaptive-routing pretest", TEST_CODES, selected = "SW"),
            checkboxGroupInput("battery", "Phase 1 ROC battery", TEST_CODES,
                               selected = c("SW", "SF", "LF", "JB", "SKEW", "DAP", "AD", "CVM"),
                               inline = TRUE)
          ),
          conditionalPanel(
            "input.approach == 'ml'",
            selectInput("ml_model", "Adaptive ML model", c("SVM", "RF", "GBM", "ANN", "LR", "KNN"), selected = "SVM"),
            checkboxGroupInput("ml_roc_models", "Phase 1 ML ROC models",
                               c("RF", "GBM", "ANN", "SVM", "LR", "KNN"),
                               selected = c("RF", "GBM", "ANN", "SVM"),
                               inline = TRUE),
            checkboxInput("ml_vote", "Include majority-vote ROC score", TRUE)
          )
        ),
        bslib::accordion_panel(
          "Design & data",
          value = "design",
          selectInput("family", "Downstream comparison",
                      choices = stats::setNames(names(FAMILIES), vapply(FAMILIES, `[[`, character(1), "label"))),
          selectInput("dist", "Non-normal distribution", setdiff(DIST_CHOICES, "normal"), selected = "exponential"),
          radioButtons("center_by", "Center by", c("median", "mean"), selected = "median", inline = TRUE),
          textInput("effect_size", "Effect size", value = "0.5")
        ),
        bslib::accordion_panel(
          "Simulation parameters",
          value = "params",
          numericInput("Nsim", "Nsim", 200, min = 10, step = 50),
          numericInput("N_tradeoff", "N_tradeoff", 200, min = 10, step = 50),
          numericInput("test_alpha", "test_alpha", 0.05, min = 0.001, max = 0.5, step = 0.01),
          numericInput("single_n", "single_n", 10, min = 3, step = 1),
          textInput("sample_sizes", "sample_sizes", value = "10, 20, 30, 40, 50"),
          textInput("effect_sizes_plot", "effect_sizes_plot", value = "0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0"),
          numericInput("tol_pos", "Type I inflation tolerance", 0.005, min = 0, step = 0.001),
          numericInput("loss_tol", "Power-loss tolerance", 0.01, min = 0, step = 0.001),
          numericInput("grid_by", "Threshold grid step", 0.025, min = 0.001, max = 0.5, step = 0.005),
          numericInput("sig_by", "Significance-level grid step", 0.005, min = 0.001, max = 0.2, step = 0.001),
          checkboxInput("per_n", "Estimate threshold separately for each sample size", TRUE)
        ),
        bslib::accordion_panel(
          "Phases & run",
          value = "run",
          checkboxGroupInput(
            "phases",
            "Phases",
            c(
              "1 - Normality ROC" = 1,
              "2 - Threshold trade-off" = 2,
              "3 - Power vs Type I ROC" = 3,
              "4 - Power and Type I vs n" = 4,
              "5 - Power vs effect" = 5
            ),
            selected = c(1, 2, 4)
          ),
          uiOutput("threshold_ui"),
          actionButton("run", "Run simulation", class = "btn-primary w-100")
        )
      )
    ),
    bslib::nav_panel("Run Summary", bslib::layout_columns(
      bslib::card(bslib::card_header("Status"), verbatimTextOutput("run_status")),
      bslib::card(bslib::card_header("Resolved arguments"), verbatimTextOutput("arg_dump"))
    )),
    bslib::nav_panel("Phase 1", bslib::layout_columns(
      bslib::card(bslib::card_header("Normality ROC"), plotOutput("p1_plot", height = 430)),
      bslib::card(bslib::card_header("AUC table"), DT::DTOutput("p1_auc"))
    )),
    bslib::nav_panel("Phase 2", bslib::layout_columns(
      bslib::value_box("Threshold", textOutput("p2_star"), theme = "primary"),
      bslib::value_box("Power gain", textOutput("p2_gain"), theme = "success"),
      bslib::value_box("Power loss", textOutput("p2_loss"), theme = "secondary")
    ), bslib::card(bslib::card_header("Trade-off"), plotOutput("p2_plot", height = 460))),
    bslib::nav_panel("Phase 3", bslib::card(bslib::card_header("Power vs Type I error"), plotOutput("p3_plot", height = 520))),
    bslib::nav_panel("Phase 4", bslib::card(bslib::card_header("Power and Type I vs n"), plotOutput("p4_plot", height = 460)),
                     bslib::card(bslib::card_header("AUC tables"), uiOutput("p4_tables"))),
    bslib::nav_panel("Phase 5", bslib::card(bslib::card_header("Power vs effect size"), plotOutput("p5_plot", height = 480))),
    bslib::nav_panel("Log & Downloads",
                     bslib::card(bslib::card_header("Console output"), verbatimTextOutput("console_log")),
                     bslib::card(bslib::card_header("Generated files"), uiOutput("downloads")))
  )
}

app_server <- function(input, output, session, raw_base = .default_raw_base()) {
  fw_rv <- reactiveVal(NULL)
  ml_rv <- reactiveVal(NULL)
  status <- reactiveVal("Framework not loaded yet.")
  run_state <- reactiveValues(results = NULL, log = "", args = NULL, wd = NULL)

  load_fw <- function() {
    status("Downloading/loading framework...")
    fw <- tryCatch(
      load_framework(refresh = isTRUE(input$refresh_assets), raw_base = trimws(input$raw_base)),
      error = function(e) {
        status(paste("Framework load failed:", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(fw)) {
      fw_rv(fw)
      ml_rv(NULL)
      status(paste("Loaded framework from", attr(fw, "framework_file")))
    }
  }

  observeEvent(TRUE, load_fw(), once = TRUE)
  observeEvent(input$load_fw, load_fw())

  output$asset_status <- renderUI({
    cls <- if (grepl("^Loaded", status())) "text-success" else "text-warning"
    tags$div(class = paste("small", cls), status())
  })

  needs_threshold <- reactive(any(c("3", "4", "5") %in% input$phases) && !("2" %in% input$phases))
  output$threshold_ui <- renderUI({
    if (!needs_threshold()) return(NULL)
    tagList(
      textInput("pretest_threshold", "pretest_threshold", value = "0.05"),
      helpText("Use a scalar, or per-n pairs like 10=0.05, 20=0.04.")
    )
  })

  family_functions <- function(fw, family_id) {
    spec <- FAMILIES[[family_id]]
    out <- lapply(spec$functions, function(name) {
      if (!exists(name, envir = fw, mode = "function", inherits = FALSE)) {
        stop("Framework function not found: ", name, call. = FALSE)
      }
      get(name, envir = fw, mode = "function", inherits = FALSE)
    })
    c(out, spec$extra)
  }

  fisher_config <- function() {
    list(
      method = "custom",
      config = list(
        fn = function(x) {
          p_sw <- stats::shapiro.test(x)$p.value
          p_ad <- nortest::ad.test(x)$p.value
          p_values <- pmax(c(p_sw, p_ad), .Machine$double.xmin)
          stat <- -2 * sum(log(p_values))
          c(p.value = stats::pchisq(stat, df = 4, lower.tail = FALSE))
        },
        label = "Fisher SW + AD"
      )
    )
  }

  norm_config <- function(fw) {
    if (identical(input$approach, "classical")) {
      return(list(method = "classical", config = list(norm_test = input$norm_test)))
    }
    if (identical(input$approach, "fisher")) {
      return(fisher_config())
    }

    resources <- ml_rv()
    if (is.null(resources)) {
      resources <- load_ml_resources_from_github(
        framework_env = fw,
        refresh = isTRUE(input$refresh_assets),
        raw_base = trimws(input$raw_base),
        required_sample_sizes = parse_nums(input$sample_sizes, c(10, 20, 30, 40, 50))
      )
      ml_rv(resources)
    }
    fw$make_ml_normality_configs(
      ml_resources = resources,
      chosen_model = input$ml_model,
      models_for_roc = input$ml_roc_models,
      use_majority_vote = isTRUE(input$ml_vote)
    )$decision
  }

  roc_config <- function(fw) {
    if (!identical(input$approach, "ml")) return(norm_config(fw))
    resources <- ml_rv()
    if (is.null(resources)) {
      resources <- load_ml_resources_from_github(
        framework_env = fw,
        refresh = isTRUE(input$refresh_assets),
        raw_base = trimws(input$raw_base),
        required_sample_sizes = parse_nums(input$sample_sizes, c(10, 20, 30, 40, 50))
      )
      ml_rv(resources)
    }
    fw$make_ml_normality_configs(
      ml_resources = resources,
      chosen_model = input$ml_model,
      models_for_roc = input$ml_roc_models,
      use_majority_vote = isTRUE(input$ml_vote)
    )$roc
  }

  build_args <- function(fw) {
    spec <- FAMILIES[[input$family]]
    phases <- as.integer(input$phases)
    selected_norm_config <- if (identical(input$approach, "ml") && identical(phases, 1L)) roc_config(fw) else norm_config(fw)
    args <- c(
      family_functions(fw, input$family),
      list(
        Nsim = input$Nsim,
        N_tradeoff = input$N_tradeoff,
        test_type = paste(spec$test_type, input$approach, sep = "_"),
        distributions = c(input$dist, "normal"),
        norm_config = selected_norm_config,
        threshold_grid = seq(0, 1, by = input$grid_by),
        normality_roc_grid = seq(0, 1, by = max(input$grid_by, 0.025)),
        tol_pos = input$tol_pos,
        loss_tol = input$loss_tol,
        test_alpha = input$test_alpha,
        center_by = input$center_by,
        effect_size = parse_nums(input$effect_size, spec$effect_size),
        sample_sizes = parse_nums(input$sample_sizes, c(10, 20, 30, 40, 50)),
        single_n = input$single_n,
        effect_sizes_plot = parse_nums(input$effect_sizes_plot, c(0, 0.2, 0.5, 0.8)),
        sig_levels = seq(0, 1, by = input$sig_by),
        norm_test = TEST_CODES,
        selected_tests = if (identical(input$approach, "classical")) input$battery else input$norm_test,
        phases = phases,
        pretest_threshold = if (needs_threshold()) parse_threshold(input$pretest_threshold) else NULL,
        per_n_thresholds = isTRUE(input$per_n),
        save_results = TRUE,
        save_compression = "gzip"
      )
    )
    args
  }

  observeEvent(input$run, {
    fw <- fw_rv()
    if (is.null(fw)) {
      showNotification("Load the framework first.", type = "error")
      return()
    }
    if (length(input$phases) == 0L) {
      showNotification("Select at least one phase.", type = "warning")
      return()
    }

    args <- tryCatch(build_args(fw), error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = NULL)
      NULL
    })
    if (is.null(args)) return()
    run_state$args <- args

    wd <- file.path(tempdir(), paste0("pretestsim_", as.integer(Sys.time())))
    dir.create(wd, recursive = TRUE, showWarnings = FALSE)
    old <- setwd(wd)
    on.exit(setwd(old), add = TRUE)

    log <- capture.output({
      res <- tryCatch(
        do.call(fw$run_simulation, args),
        error = function(e) {
          cat("ERROR:", conditionMessage(e), "\n")
          NULL
        }
      )
    }, type = "output")

    run_state$results <- if (exists("res")) res else NULL
    run_state$log <- paste(log, collapse = "\n")
    run_state$wd <- wd

    if (is.null(run_state$results)) {
      showNotification("Run failed. Check the log.", type = "error", duration = NULL)
    } else {
      showNotification("Run complete.", type = "message")
    }
  })

  ran <- function(phase) {
    res <- run_state$results
    !is.null(res) && phase %in% res$params$phases_run
  }
  obj <- function() run_state$results$objects
  par <- function() run_state$results$params
  dark <- function(expr) {
    op <- graphics::par(bg = "#0f1419", fg = "#e6edf3", col.axis = "#c9d1d9", col.lab = "#e6edf3", col.main = "#e6edf3")
    on.exit(graphics::par(op), add = TRUE)
    expr
  }

  output$run_status <- renderText({
    paste(status(), if (!is.null(run_state$results)) "\nLast run complete." else "\nNo completed run yet.")
  })
  output$arg_dump <- renderText({
    args <- run_state$args
    if (is.null(args)) return("Run a simulation to see resolved arguments.")
    keep <- setdiff(names(args), c("gen_data", "get_parameters", "fn_to_get_norm_obj", "fn_for_ds_test_1", "fn_for_ds_test_2", "norm_config"))
    paste(vapply(keep, function(name) {
      value <- args[[name]]
      paste0(name, " = ", paste(value, collapse = ", "))
    }, character(1)), collapse = "\n")
  })
  output$console_log <- renderText(run_state$log %||% "Run a simulation to see console output.")

  output$p1_plot <- renderPlot({
    validate(need(ran(1), "Run Phase 1 to see ROC curves."))
    roc <- obj()$roc_pval_ds_test
    validate(need(!is.null(roc), "Phase 1 did not return ROC data."))
    tests <- rownames(roc$FPR)
    curves <- lapply(tests, function(test) list(fpr = as.numeric(roc$FPR[test, ]), tpr = as.numeric(roc$TPR[test, ]), label = test, lty = 1))
    dark(fw_rv()$plot_norm_roc_curve(curves = curves, dist_name = par()$distributions[1]))
  })
  output$p1_auc <- DT::renderDT({
    validate(need(ran(1), ""))
    roc <- obj()$roc_pval_ds_test
    validate(need(!is.null(roc), ""))
    tests <- rownames(roc$FPR)
    auc <- vapply(tests, function(test) fw_rv()$compute_auc(as.numeric(roc$FPR[test, ]), as.numeric(roc$TPR[test, ])), numeric(1))
    DT::datatable(data.frame(Test = tests, AUC = round(auc, 4)), rownames = FALSE, options = list(dom = "t"))
  })

  output$p2_star <- renderText({ validate(need(ran(2), "-")); sprintf("%.4f", obj()$param_star) })
  output$p2_gain <- renderText({ validate(need(ran(2), "-")); sprintf("%.4f", obj()$optimal_result$power_gain %||% NA_real_) })
  output$p2_loss <- renderText({ validate(need(ran(2), "-")); sprintf("%.4f", obj()$optimal_result$power_loss %||% NA_real_) })
  output$p2_plot <- renderPlot({
    validate(need(ran(2), "Run Phase 2 to see trade-off results."))
    dark(fw_rv()$plot_tradeoff_results_generic(obj()$optimal_result))
  })

  output$p3_plot <- renderPlot({
    validate(need(ran(3), "Run Phase 3 to see power vs Type I error."))
    dark(fw_rv()$power_vs_error_roc_plot(obj()$roc_data, nominal_alpha = par()$test_alpha, optimal_result = obj()$optimal_result))
  })

  output$p4_plot <- renderPlot({
    validate(need(ran(4), "Run Phase 4 to see power and Type I error by sample size."))
    dark(fw_rv()$plot_power_type1_results(
      combined_power = obj()$combined_power_n,
      combined_type1 = obj()$combined_type1_n,
      ds_test_methods = par()$ds_test_methods,
      distributions = par()$distributions,
      sample_sizes = par()$sample_sizes,
      test_alpha = par()$test_alpha,
      optimal_result = obj()$optimal_result,
      norm_config = par()$norm_config
    ))
  })
  output$p4_tables <- renderUI({
    validate(need(ran(4), ""))
    tables <- obj()$auc_tables
    do.call(tagList, lapply(names(tables), function(name) {
      tagList(h5(name), tableOutput(paste0("auc_", name)))
    }))
  })
  observe({
    tables <- if (!is.null(run_state$results)) obj()$auc_tables else NULL
    if (is.null(tables)) return()
    for (name in names(tables)) {
      local({
        id <- paste0("auc_", name)
        table <- tables[[name]]
        output[[id]] <- renderTable(table, striped = TRUE, bordered = TRUE)
      })
    }
  })

  output$p5_plot <- renderPlot({
    validate(need(ran(5), "Run Phase 5 to see power by effect size."))
    eff_len <- length(par()$effect_size_input)
    dark(fw_rv()$plot_power_by_effect_size(
      combined_power = obj()$combined_power_effect,
      fixed_n = par()$single_n,
      optimal_result = obj()$optimal_result,
      norm_config = par()$norm_config,
      eff_len = eff_len
    ))
  })

  output$downloads <- renderUI({
    if (is.null(run_state$results) || is.null(run_state$wd)) return(helpText("Run a simulation first."))
    files <- list.files(file.path(run_state$wd, "results"), recursive = TRUE, full.names = TRUE)
    if (!length(files)) return(helpText("No output files found."))
    do.call(tagList, lapply(seq_along(files), function(i) {
      id <- paste0("dl_", i)
      output[[id]] <- downloadHandler(
        filename = function() basename(files[i]),
        content = function(file) file.copy(files[i], file, overwrite = TRUE)
      )
      tags$div(downloadButton(id, basename(files[i]), class = "btn-sm mb-1"))
    }))
  })
}
