# =============================================================================
# app.R  —  Pretest Simulation Studio (drives the ORIGINAL framework)
# -----------------------------------------------------------------------------
# This dashboard runs the user's UNMODIFIED framework function `run_simulation()`
# (from user_framework_rpk_v3_0.R) through its six phases:
#
#   Phase 1  Normality-pretest ROC / AUROC comparison of tests
#   Phase 2  Trade-off analysis  ->  optimal pretest threshold (param_star)
#   Phase 3  Power vs Type I error ROC-like curves (at the optimal threshold)
#   Phase 4  Power & Type I error vs sample size  + AUC-Power / AUC-Type I tables
#   Phase 5  Power vs effect size (fixed n)
#   Phase 6  Save all results (RData) and the phase PDFs
#
# The user chooses Classical or Custom/ML normality assessment, supplies custom
# functions and optional trained-model RData, sets all simulation parameters
# (Nsim, N_tradeoff, tol_pos, loss_tol, effect_size, ...), and selects which
# phases to run. If Phase 2 is skipped while Phases 3/4/5 are requested, a
# pretest threshold must be supplied — exactly as the framework requires.
#
# Run:  shiny::runApp("this-folder")
#       (keep app.R, framework_loader.R and user_framework_rpk_v3_0.R together)
# =============================================================================

required <- c("shiny", "bslib", "DT")
miss <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")

library(shiny)
library(bslib)
library(DT)

# Allow large uploads (e.g. trained-model RData) up to ~2 GB; for very large
# files prefer loading by path on disk (see the ML panel).
options(shiny.maxRequestSize = 2 * 1024^3)

source("framework_loader.R", local = TRUE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# distributions supported by generate_data() in the framework
DIST_CHOICES <- c("exponential", "lognormal", "chi_square", "gamma", "t", "f",
                  "uniform", "laplace", "cauchy", "gumbel", "weibull", "beta",
                  "logistic", "pareto", "contaminated", "normal")

# classical normality-test codes recognised by generate_tests()
TEST_CODES <- c("SW", "SF", "LF", "KS", "JB", "SKEW", "KURT", "DAP", "AD", "CVM")

# test families -> (fns-bundle base name, label, non-parametric alternatives)
FAMILIES <- list(
  one_sample = list(label = "One-sample (t-test vs …)",
                    base = "onesample", perm = "onesample_perm",
                    alts = c("Sign test" = "onesample", "Permutation test" = "onesample_perm")),
  two_sample = list(label = "Two-sample (t-test vs …)",
                    base = "twosample", perm = "twosample_perm",
                    alts = c("Mann-Whitney U" = "twosample", "Permutation test" = "twosample_perm")),
  anova      = list(label = "One-way ANOVA (F-test vs …)",
                    base = "anova", perm = "anova_perm",
                    alts = c("Kruskal-Wallis" = "anova", "Permutation ANOVA" = "anova_perm")),
  regression = list(label = "Regression (OLS vs …)",
                    base = "regression", perm = "regression_perm",
                    alts = c("Rank-based" = "regression", "Permutation" = "regression_perm"))
)

parse_nums <- function(txt, default = numeric(0)) {
  if (is.null(txt) || !nzchar(trimws(txt))) return(default)
  v <- suppressWarnings(as.numeric(strsplit(txt, "[,;[:space:]]+")[[1]]))
  v[!is.na(v)]
}

PLOT_BG <- "#0f1419"

# Theme + explicit navbar rules so the top tabs are clearly visible (light text,
# highlighted active/hover) on the dark bar.
app_theme <- bs_add_rules(
  bs_theme(version = 5, bg = "#0f1419", fg = "#e6edf3",
           primary = "#4c8bf5", success = "#3fb950",
           base_font = font_google("Inter", local = FALSE),
           "card-bg" = "#161b22", "card-border-color" = "#30363d"),
  "
  /* ---- Top navigation: a clearly visible strip with button-style tabs ---- */
  .navbar {
    background-color: #21262d !important;
    border-bottom: 2px solid #4c8bf5 !important;
    min-height: 58px;
    --bs-nav-link-color: #ffffff;
    --bs-nav-link-hover-color: #ffffff;
  }
  .navbar-brand { color: #ffffff !important; font-weight: 700; }

  .navbar .navbar-nav .nav-link,
  .navbar .nav-link {
    color: #ffffff !important;
    font-weight: 600 !important;
    background-color: rgba(255,255,255,0.07);
    border: 1px solid #3a4250;
    border-radius: 7px;
    padding: 0.4rem 0.85rem !important;
    margin: 4px 3px;
  }
  .navbar .nav-link:hover,
  .navbar .nav-link:focus {
    color: #ffffff !important;
    background-color: rgba(76,139,245,0.40) !important;
    border-color: #4c8bf5;
  }
  .navbar .nav-link.active,
  .navbar .show > .nav-link {
    color: #ffffff !important;
    background-color: #4c8bf5 !important;
    border-color: #4c8bf5 !important;
    box-shadow: 0 0 0 2px rgba(76,139,245,0.35);
  }
  .navbar .navbar-toggler {
    border-color: #4c8bf5;
    background-color: rgba(76,139,245,0.15);
  }
  "
)

# =============================================================================
# UI
# =============================================================================
ui <- page_navbar(
  title = span(icon("flask-vial"), "Pretest Simulation Studio"),
  theme = app_theme,
  fillable = TRUE,

  # ----- shared, always-visible control sidebar ----------------------------
  sidebar = sidebar(
    width = 360, title = "Setup & Controls", class = "fw-controls",
    accordion(
      open = c("fw", "ap", "fam", "ph"),

      accordion_panel(
        "Framework file", value = "fw", icon = icon("file-code"),
        textInput("fw_path", "Path to your framework .R file", value = ""),
        actionButton("load_fw", "Load / Reload", icon = icon("rotate"),
                     class = "btn-secondary w-100 mb-2"),
        fileInput("fw_upload", "…or upload the framework .R", accept = ".R"),
        uiOutput("fw_status")
      ),

      accordion_panel(
        "Normality approach", value = "ap", icon = icon("scale-balanced"),
        radioButtons("approach", NULL,
          c("Classical" = "classical", "Custom / ML" = "custom"), inline = TRUE),
        conditionalPanel("input.approach == 'classical'",
          selectInput("norm_test_single", "Pretest used for the adaptive rule",
                      TEST_CODES, selected = "SW"),
          checkboxGroupInput("battery", "Phase-1 ROC battery (tests to compare)",
                      TEST_CODES, selected = c("SW", "SF", "LF", "JB", "AD"),
                      inline = TRUE)
        ),
        conditionalPanel("input.approach == 'custom'",
          checkboxInput("ml_mode", "Machine-Learning mode (load trained models)", FALSE),
          conditionalPanel("!input.ml_mode",
            textAreaInput("custom_code",
              "Custom normality function (returns a single named numeric scalar)",
              value = paste(
                "function(x) {",
                "  p_sw <- shapiro.test(x)$p.value",
                "  p_ad <- nortest::ad.test(x)$p.value",
                "  stat <- -2 * sum(log(c(p_sw, p_ad)))",
                "  c(p.value = pchisq(stat, df = 4, lower.tail = FALSE))",
                "}", sep = "\n"),
              height = "150px"),
            textInput("custom_label", "Curve label", value = "Fisher (SW+AD)")
          ),
          conditionalPanel("input.ml_mode",
            helpText(class = "small",
              "Point to files on disk (recommended — avoids the browser upload limit)."),
            textInput("ml_rfile",
              "ML R file to source into ml_env (path)", value = ""),
            textInput("ml_rdata_path",
              "trained_models .RData (path)", value = ""),
            tags$details(tags$summary(class = "small text-muted",
              "…or upload a small .RData instead"),
              fileInput("ml_rdata", NULL, accept = c(".RData", ".rdata", ".rda"))),
            radioButtons("ml_fn", "ML scoring function",
              c("ml_fn_roc (compare models — use for Phase 1)" = "roc",
                "ml_fn_decision (single model — use for Phases 2-6)" = "decision")),
            textInput("ml_models", "model_name(s), comma-separated",
                      value = "RF, GBM, ANN, SVM, LR"),
            checkboxInput("ml_vote", "use_majority_vote", TRUE),
            numericInput("ml_dthr", "decision_threshold", 0.50, 0, 1, 0.05)
          )
        )
      ),

      accordion_panel(
        "Test family & data", value = "fam", icon = icon("code-branch"),
        selectInput("family", "Downstream test family",
                    setNames(names(FAMILIES), vapply(FAMILIES, `[[`, "", "label"))),
        uiOutput("alt_ui"),
        selectInput("dist", "Non-normal distribution (H1)",
                    setdiff(DIST_CHOICES, "normal"), selected = "exponential"),
        helpText("The reference distribution (H0) is always 'normal'."),
        radioButtons("center_by", "center_by", c("median", "mean"), inline = TRUE),
        textInput("effect_size", "effect_size (scalar; ANOVA: comma-separated means)",
                  value = "0.5")
      ),

      accordion_panel(
        "Simulation parameters", value = "par", icon = icon("sliders"),
        div(class = "row g-2",
          div(class = "col", numericInput("Nsim", "Nsim", 200, 20, 1e5, 50)),
          div(class = "col", numericInput("N_tradeoff", "N_tradeoff", 300, 20, 1e5, 50))),
        div(class = "row g-2",
          div(class = "col", numericInput("test_alpha", "test_alpha", 0.05, 0.001, 0.5, 0.01)),
          div(class = "col", numericInput("single_n", "single_n", 20, 3, 1e4, 1))),
        textInput("sample_sizes", "sample_sizes", value = "10, 20, 30, 40, 50"),
        textInput("effect_sizes_plot", "effect_sizes_plot", value = "0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0"),
        div(class = "row g-2",
          div(class = "col", numericInput("tol_pos", "tol_pos", 0.005, 0, 1, 0.001)),
          div(class = "col", numericInput("loss_tol", "loss_tol", 0.01, 0, 1, 0.001))),
        h6("threshold_grid  (seq from / to / by)", class = "mt-2"),
        div(class = "row g-2",
          div(class = "col", numericInput("tg_from", "from", 0, 0, 1, 0.005)),
          div(class = "col", numericInput("tg_to", "to", 1, 0, 1, 0.005)),
          div(class = "col", numericInput("tg_by", "by", 0.025, 0.001, 0.5, 0.005))),
        numericInput("sig_by", "sig_levels step (seq 0..1 by)", 0.005, 0.001, 0.1, 0.005),
        checkboxInput("per_n", "per_n_thresholds (a threshold for each sample size)", FALSE)
      ),

      accordion_panel(
        "Phases & run", value = "ph", icon = icon("list-check"),
        checkboxGroupInput("phases", "Phases to run",
          c("1 · Normality ROC"        = 1,
            "2 · Trade-off (threshold)" = 2,
            "3 · Power–Error ROC"       = 3,
            "4 · Power & Type I vs n"   = 4,
            "5 · Power vs effect size"  = 5,
            "6 · Save outputs"          = 6),
          selected = c(1, 2, 4)),
        uiOutput("thr_ui"),
        uiOutput("phase_warn"),
        actionButton("run", "Run simulation", icon = icon("play"),
                     class = "btn-primary w-100"),
        helpText("Tip: start with small Nsim / N_tradeoff; raise for final runs.")
      )
    )
  ),

  # ----- phase result panels ------------------------------------------------
  nav_panel("Phase 1 · Normality ROC", icon = icon("chart-line"),
    layout_columns(col_widths = c(7, 5),
      card(card_header("AUROC curves — normality tests"),
           plotOutput("p1_plot", height = 420)),
      card(card_header("Area under ROC (per test)"), DTOutput("p1_auc")))
  ),
  nav_panel("Phase 2 · Trade-off", icon = icon("scale-balanced"),
    layout_columns(col_widths = c(3, 3, 3, 3),
      value_box("Optimal threshold", textOutput("p2_star"),
                showcase = icon("bullseye"), theme = "primary"),
      value_box("Power gain", textOutput("p2_gain"), theme = "success"),
      value_box("Power loss", textOutput("p2_loss"), theme = "secondary"),
      value_box("Type I inflation", textOutput("p2_infl"), theme = "secondary")),
    card(card_header("Trade-off analysis"), plotOutput("p2_plot", height = 430)),
    card(card_header("Selection note"), verbatimTextOutput("p2_note"))
  ),
  nav_panel("Phase 3 · Power–Error ROC", icon = icon("bezier-curve"),
    card(card_header("Power vs Type I error (ROC-like) at the optimal threshold"),
         plotOutput("p3_plot", height = 520))
  ),
  nav_panel("Phase 4 · Power & Type I vs n", icon = icon("bolt"),
    card(card_header("Power and Type I error across sample sizes"),
         plotOutput("p4_plot", height = 440)),
    card(card_header("AUC-Power / AUC-Type I tables"), uiOutput("p4_tables"))
  ),
  nav_panel("Phase 5 · Power vs effect", icon = icon("arrow-trend-up"),
    card(card_header("Power vs effect size (fixed n)"),
         plotOutput("p5_plot", height = 480))
  ),
  nav_panel("Log & Downloads", icon = icon("download"),
    layout_columns(col_widths = c(6, 6),
      card(card_header("Resolved run arguments"), verbatimTextOutput("arg_dump")),
      card(card_header("Framework load report"), verbatimTextOutput("load_report"))),
    card(card_header("Console output"), verbatimTextOutput("console_log")),
    card(card_header("Output files (Phase 6)"), uiOutput("downloads"))
  ),
  nav_panel("About", icon = icon("circle-info"),
    card(card_header("How this app uses your framework"), markdown(
"This dashboard loads your **original** `user_framework_rpk_v3_0.R` *unedited* and
calls its `run_simulation()` engine. The loader skips only the machine-specific
side-effects (`setwd`, the ML `source()`, the `pacman` install block, the
`load()` of trained models) and the example `run_simulation()` calls at the end
of the file. Every function, every `norm_config_*` preset, and every
`*_fns` / `*_params` bundle is used exactly as written.

**Phases** (select any subset in the sidebar):

1. **Normality ROC** — compares normality tests (classical battery, a custom
   combined test, or ML scores) via AUROC curves.
2. **Trade-off** — sweeps the pretest threshold and picks the optimal
   `param_star` via `select_optimal_parameter()` (governed by `tol_pos`,
   `loss_tol`).
3. **Power–Error ROC** — ROC-like power vs Type I error at the chosen threshold.
4. **Power & Type I vs n** — across `sample_sizes`, with AUC-Power and
   AUC-Type I tables.
5. **Power vs effect size** — at the fixed `single_n` over `effect_sizes_plot`.
6. **Save** — writes the phase PDFs and `*_results.RData`.

**Phase dependency.** Phases 3, 4 and 5 need a pretest threshold. It normally
comes from Phase 2; if Phase 2 is skipped you must supply `pretest_threshold`
(a scalar, or `n=value` pairs for per-n thresholds) — the same rule the
framework enforces.

**Custom / ML.** Choose *Custom / ML*, then either paste an R function that
returns a single named numeric scalar (e.g. the Fisher combined test), or enable
**ML mode** and upload a `.RData` containing `trained_models`; the app wires it
into the framework's `ml_fn_roc` / `ml_fn_decision` scorers."))
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  # ---- framework environment (auto-detect, manual path, or upload) ----
  fw_rv     <- reactiveVal(NULL)
  fw_status <- reactiveVal("Looking for the framework file…")

  detect_fw <- function() {
    names <- c("user_framework_rpk_v3_0.R", "user_framework_rpk_v3.0.R",
               "user_framework_rpkg_v3.0.R", "user_framework_rpkg_v3_0.R")
    cand <- c(names,
              file.path("..", names),
              file.path("../one_sample", names),
              file.path("/mnt/project", names),
              Sys.glob("*framework*v3*.R"),
              Sys.glob("../*framework*v3*.R"),
              Sys.glob("../*/*framework*v3*.R"))
    cand <- cand[file.exists(cand)]
    if (length(cand)) normalizePath(cand[1]) else NA_character_
  }

  try_load <- function(path) {
    if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
      fw_rv(NULL); fw_status(paste0("File not found: ", path %||% "(empty)"))
      return(invisible())
    }
    fw <- tryCatch(load_framework(path),
                   error = function(e) { fw_status(paste("Load error:", conditionMessage(e))); NULL })
    if (!is.null(fw)) {
      r <- attr(fw, "load_report")
      fw_rv(fw)
      fw_status(sprintf("Loaded \u2713  (%d definitions kept, %d skipped)\n%s",
                        r$kept, r$skipped, path))
      showNotification("Framework loaded.", type = "message")

      # Try to pre-fill the ML file paths from the framework's folder tree.
      d <- dirname(path)
      roots <- unique(c(d, dirname(d)))
      ml_r <- character(0); ml_d <- character(0)
      for (root in roots) {
        if (!length(ml_r)) ml_r <- list.files(root, "ML_Approach.*\\.R$",
                                               full.names = TRUE, recursive = TRUE)
        if (!length(ml_d)) ml_d <- list.files(root, "trained_models.*\\.RData$",
                                               full.names = TRUE, recursive = TRUE)
      }
      if (length(ml_r)) updateTextInput(session, "ml_rfile", value = normalizePath(ml_r[1]))
      if (length(ml_d)) updateTextInput(session, "ml_rdata_path", value = normalizePath(ml_d[1]))
    } else fw_rv(NULL)
  }

  observeEvent(TRUE, {                       # startup auto-detect
    p <- detect_fw()
    updateTextInput(session, "fw_path", value = if (is.na(p)) "" else p)
    if (!is.na(p)) try_load(p)
    else fw_status("Framework not found automatically. Enter its path or upload the .R file.")
  }, once = TRUE)

  observeEvent(input$load_fw, try_load(trimws(input$fw_path)))

  observeEvent(input$fw_upload, {
    f <- input$fw_upload; if (is.null(f)) return()
    dest <- file.path(tempdir(), f$name)
    file.copy(f$datapath, dest, overwrite = TRUE)
    updateTextInput(session, "fw_path", value = dest)
    try_load(dest)
  })

  output$fw_status <- renderUI({
    s <- fw_status()
    cls <- if (grepl("Loaded", s)) "text-success"
           else if (grepl("error|not found|not loaded", s, ignore.case = TRUE)) "text-warning"
           else "text-muted"
    div(class = paste("small", cls), HTML(gsub("\n", "<br>", s)))
  })

  # non-parametric alternative dropdown depends on family
  output$alt_ui <- renderUI({
    fam <- FAMILIES[[input$family]]
    selectInput("alt", "Non-parametric alternative (test_2)", choices = fam$alts)
  })

  # phase-dependency: threshold input + warning
  needs_thr <- reactive(any(c("3","4","5") %in% input$phases) && !("2" %in% input$phases))
  output$thr_ui <- renderUI({
    if (!needs_thr()) return(NULL)
    tagList(
      textInput("pretest_threshold",
                "pretest_threshold (Phase 2 skipped — required)",
                value = "0.05"),
      helpText("Scalar (e.g. 0.05) or per-n pairs (e.g. 10=0.05, 20=0.04).")
    )
  })
  output$phase_warn <- renderUI({
    if (needs_thr() && !nzchar(trimws(input$pretest_threshold %||% "")))
      div(class = "text-warning small mb-2",
          icon("triangle-exclamation"),
          " Supply a pretest threshold or include Phase 2.")
  })

  # ---- assemble the norm_config from the chosen approach ----
  build_norm_config <- function(fw) {
    if (input$approach == "classical")
      return(list(method = "classical", config = list(norm_test = input$norm_test_single)))

    # custom / ML
    if (isTRUE(input$ml_mode)) {
      # 1) Source the ML R file into the framework's ml_env so the ML scorers
      #    (ml_fn_roc / ml_fn_decision / classify_sample_prob) have their helpers.
      rfile <- trimws(input$ml_rfile %||% "")
      if (nzchar(rfile)) {
        if (!file.exists(rfile)) stop("ML R file not found: ", rfile)
        if (!is.environment(fw$ml_env))
          assign("ml_env", new.env(parent = globalenv()), envir = fw)
        tryCatch(sys.source(rfile, envir = fw$ml_env),
                 error = function(e) stop("Error sourcing ML R file: ", conditionMessage(e)))
      }
      # 2) Load trained_models — from a disk path (no upload limit) or an upload.
      rpath <- trimws(input$ml_rdata_path %||% "")
      if (nzchar(rpath)) {
        if (!file.exists(rpath)) stop("trained_models .RData not found: ", rpath)
        load(rpath, envir = fw)
      } else if (!is.null(input$ml_rdata)) {
        load(input$ml_rdata$datapath, envir = fw)
      } else {
        stop("Provide the trained_models .RData path (recommended) or upload it.")
      }
      if (!exists("trained_models", envir = fw, inherits = FALSE))
        stop("The .RData must contain an object named 'trained_models'.")
      models <- trimws(strsplit(input$ml_models, ",")[[1]])
      fn <- if (input$ml_fn == "roc") fw$ml_fn_roc else fw$ml_fn_decision
      if (!is.function(fn))
        stop("ml_fn_roc / ml_fn_decision not found — check the framework / ML R file.")
      return(list(method = "custom", config = list(
        fn = fn,
        trained_models = get("trained_models", envir = fw),
        model_name = if (input$ml_fn == "decision") models[1] else models,
        use_majority_vote = isTRUE(input$ml_vote),
        decision_threshold = input$ml_dthr)))
    }
    fn <- tryCatch(eval(parse(text = input$custom_code), envir = fw),
                   error = function(e) stop("Custom code error: ", conditionMessage(e)))
    if (!is.function(fn)) stop("Custom code must evaluate to a function.")
    list(method = "custom", config = list(fn = fn, label = input$custom_label))
  }

  # ---- assemble full argument list for run_simulation ----
  build_args <- function(fw) {
    fam   <- FAMILIES[[input$family]]
    fns   <- get(paste0(input$alt, "_fns"), envir = fw)   # the chosen bundle
    es    <- parse_nums(input$effect_size, 0.5)
    es    <- if (length(es) == 1) es[1] else es           # scalar or vector (anova)
    phases <- as.integer(input$phases)
    nc    <- build_norm_config(fw)

    test_type <- paste0(input$family, "_", input$alt, "_",
                        if (input$approach == "classical") input$norm_test_single else "custom")

    args <- c(fns, list(
      Nsim           = input$Nsim,
      N_tradeoff     = input$N_tradeoff,
      test_type      = test_type,
      distributions  = c(input$dist, "normal"),
      norm_config    = nc,
      threshold_grid = seq(input$tg_from, input$tg_to, by = input$tg_by),
      tol_pos        = input$tol_pos,
      loss_tol       = input$loss_tol,
      test_alpha     = input$test_alpha,
      center_by      = input$center_by,
      effect_size    = es,
      sample_sizes   = parse_nums(input$sample_sizes, c(10,20,30,40,50)),
      single_n       = input$single_n,
      effect_sizes_plot = parse_nums(input$effect_sizes_plot, c(0,0.2,0.5,0.8)),
      sig_levels     = seq(0, 1, by = input$sig_by),
      norm_test      = if (input$approach == "classical") input$battery else input$norm_test_single,
      selected_tests = if (input$approach == "classical") input$battery else input$norm_test_single,
      phases         = phases,
      per_n_thresholds = isTRUE(input$per_n)
    ))

    # pretest_threshold when Phase 2 is skipped
    if (needs_thr()) {
      raw <- trimws(input$pretest_threshold %||% "")
      if (grepl("=", raw)) {                       # per-n pairs "10=0.05, 20=0.04"
        pr <- strsplit(raw, ",")[[1]]
        nm <- trimws(sub("=.*", "", pr)); vl <- as.numeric(sub(".*=", "", pr))
        thr <- setNames(vl, nm)
      } else thr <- as.numeric(raw)
      args$pretest_threshold <- thr
    }
    args
  }

  # ---- run ----
  run_state <- reactiveValues(results = NULL, log = "", args = NULL)

  observeEvent(input$run, {
    fw <- fw_rv()
    if (is.null(fw)) { showNotification("Framework not loaded.", type = "error"); return() }
    if (length(input$phases) == 0) { showNotification("Select at least one phase.", type = "warning"); return() }
    if (needs_thr() && !nzchar(trimws(input$pretest_threshold %||% ""))) {
      showNotification("Phase 2 is skipped — supply a pretest threshold.", type = "error"); return()
    }

    args <- tryCatch(build_args(fw), error = function(e) {
      showNotification(conditionMessage(e), type = "error"); NULL })
    if (is.null(args)) return()
    run_state$args <- args

    withProgress(message = "Running run_simulation()", value = 0.1, {
      wd <- file.path(tempdir(), paste0("fwrun_", as.integer(Sys.time())))
      dir.create(wd, showWarnings = FALSE, recursive = TRUE)
      old <- setwd(wd); on.exit(setwd(old), add = TRUE)

      log <- capture.output({
        res <- tryCatch(suppressWarnings(do.call(fw$run_simulation, args)),
                        error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
      }, type = "output")
      incProgress(0.8)
      run_state$results <- if (exists("res")) res else NULL
      run_state$log <- paste(log, collapse = "\n")
      attr(run_state$results, "wd") <- wd
    })

    if (is.null(run_state$results))
      showNotification("Run failed — see Log tab.", type = "error", duration = NULL)
    else
      showNotification("Done. See the phase tabs.", type = "message")
  })

  ran <- function(p) {
    r <- run_state$results
    !is.null(r) && p %in% r$params$phases_run
  }
  obj <- function() run_state$results$objects
  par <- function() run_state$results$params

  # base-R plot styling for the dark theme
  with_dark <- function(expr) {
    op <- par_save <- graphics::par(bg = PLOT_BG, fg = "#e6edf3",
        col.axis = "#c9d1d9", col.lab = "#e6edf3", col.main = "#e6edf3")
    on.exit(graphics::par(par_save))
    expr
  }

  # ---------------- Phase 1 ----------------
  output$p1_plot <- renderPlot({
    validate(need(ran(1), "Run Phase 1 to see normality-test ROC curves."))
    fw <- fw_rv(); roc <- obj()$roc_pval_ds_test
    validate(need(!is.null(roc), "Phase 1 produced no ROC object (check Log tab)."))
    tests <- rownames(roc$FPR)
    curves <- lapply(tests, function(t) list(fpr = as.numeric(roc$FPR[t, ]),
                  tpr = as.numeric(roc$TPR[t, ]), label = t, lty = 1))
    with_dark(fw$plot_norm_roc_curve(curves = curves, dist_name = par()$distributions[1]))
  })
  output$p1_auc <- renderDT({
    validate(need(ran(1), ""))
    fw <- fw_rv(); roc <- obj()$roc_pval_ds_test
    validate(need(!is.null(roc), ""))
    tests <- rownames(roc$FPR)
    auc <- vapply(tests, function(t)
      fw$compute_auc(as.numeric(roc$FPR[t, ]), as.numeric(roc$TPR[t, ])), numeric(1))
    df <- data.frame(Test = tests, AUC = round(auc, 4))
    df <- df[order(-df$AUC), ]
    datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  # ---------------- Phase 2 ----------------
  output$p2_star <- renderText({ validate(need(ran(2), "—")); sprintf("%.4f", obj()$param_star) })
  output$p2_gain <- renderText({ validate(need(ran(2), "—"))
    v <- obj()$optimal_result$power_gain; if (is.null(v)||is.na(v)) "n/a" else sprintf("%.4f", v) })
  output$p2_loss <- renderText({ validate(need(ran(2), "—"))
    v <- obj()$optimal_result$power_loss; if (is.null(v)||is.na(v)) "n/a" else sprintf("%.4f", v) })
  output$p2_infl <- renderText({ validate(need(ran(2), "—"))
    v <- obj()$optimal_result$inflation_non_normal %||% obj()$optimal_result$inflation_normal
    if (is.null(v)||is.na(v)) "n/a" else sprintf("%.4f", v) })
  output$p2_plot <- renderPlot({
    validate(need(ran(2), "Run Phase 2 for the trade-off analysis."))
    fw <- fw_rv()
    with_dark(fw$plot_tradeoff_results_generic(
      optimal_result = obj()$optimal_result,
      outer_title = sprintf("Power / Type I Trade-off  |  threshold = %.4f", obj()$param_star)))
  })
  output$p2_note <- renderText({
    validate(need(ran(2), "—"))
    obj()$optimal_result$selection_note %||% "(no note)"
  })

  # ---------------- Phase 3 ----------------
  output$p3_plot <- renderPlot({
    validate(need(ran(3), "Run Phase 3 for the power vs Type I error ROC."))
    fw <- fw_rv()
    validate(need(!is.null(obj()$roc_data), "Phase 3 produced no data (check Log)."))
    with_dark(fw$power_vs_error_roc_plot(
      roc_results = obj()$roc_data, nominal_alpha = par()$test_alpha,
      optimal_result = obj()$optimal_result))
  })

  # ---------------- Phase 4 ----------------
  output$p4_plot <- renderPlot({
    validate(need(ran(4), "Run Phase 4 for power & Type I error vs sample size."))
    fw <- fw_rv()
    with_dark(fw$plot_power_type1_results(
      combined_power = obj()$combined_power_n, combined_type1 = obj()$combined_type1_n,
      ds_test_methods = par()$ds_test_methods, distributions = par()$distributions,
      sample_sizes = par()$sample_sizes, test_alpha = par()$test_alpha,
      optimal_result = obj()$optimal_result, norm_config = par()$norm_config))
  })
  output$p4_tables <- renderUI({
    validate(need(ran(4), "—"))
    tabs <- obj()$auc_tables
    do.call(tagList, lapply(names(tabs), function(d) {
      tb <- tabs[[d]]; tb[] <- lapply(tb, function(c) if (is.numeric(c)) round(c, 4) else c)
      tagList(h6(paste("Distribution:", d)),
              renderTable(tb, striped = TRUE, bordered = TRUE))
    }))
  })

  # ---------------- Phase 5 ----------------
  output$p5_plot <- renderPlot({
    validate(need(ran(5), "Run Phase 5 for power vs effect size."))
    fw <- fw_rv()
    es <- par()$effect_size_input
    eff_len <- if (length(es) > 1) length(es) else 1
    with_dark(fw$plot_power_by_effect_size(
      combined_power = obj()$combined_power_effect, fixed_n = par()$single_n,
      optimal_result = obj()$optimal_result, norm_config = par()$norm_config,
      eff_len = eff_len))
  })

  # ---------------- Log & downloads ----------------
  output$console_log   <- renderText(run_state$log %||% "Run a simulation to see output.")
  output$load_report   <- renderText({
    fw <- fw_rv(); if (is.null(fw)) return("Framework not loaded.")
    r <- attr(fw, "load_report")
    paste0("Kept expressions : ", r$kept, "\nSkipped (side-effects & examples): ",
           r$skipped, "\nLoad-time messages: ", length(r$errors),
           if (length(r$errors)) paste0("\n  - ", paste(r$errors, collapse = "\n  - ")) else "")
  })
  output$arg_dump <- renderText({
    a <- run_state$args; if (is.null(a)) return("Run a simulation first.")
    keep <- setdiff(names(a), c("gen_data","get_parameters","fn_to_get_norm_obj",
                                "fn_for_ds_test_1","fn_for_ds_test_2","norm_config"))
    lines <- vapply(keep, function(k) {
      v <- a[[k]]; paste0("  ", k, " = ",
        if (is.numeric(v)) paste(round(v, 4), collapse = ", ") else paste(v, collapse = ", "))
    }, "")
    paste0("norm_config$method = ", a$norm_config$method, "\n",
           "phases = ", paste(a$phases, collapse = ", "), "\n",
           paste(lines, collapse = "\n"))
  })
  output$downloads <- renderUI({
    r <- run_state$results
    if (is.null(r) || !("6" %in% as.character(r$params$phases_run)) && !(6 %in% r$params$phases_run))
      return(helpText("Include Phase 6 to write PDFs and RData."))
    wd <- attr(r, "wd"); files <- list.files(file.path(wd, "results"),
                                              recursive = TRUE, full.names = TRUE)
    if (!length(files)) return(helpText("No output files found."))
    do.call(tagList, lapply(seq_along(files), function(i) {
      id <- paste0("dl_", i)
      output[[id]] <- downloadHandler(
        filename = function() basename(files[i]),
        content  = function(file) file.copy(files[i], file, overwrite = TRUE))
      tags$div(downloadButton(id, basename(files[i]), class = "btn-sm mb-1"))
    }))
  })
}

shinyApp(ui, server)
