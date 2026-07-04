# pretestsim

An R package that wraps the **pretest-simulation framework**
(`user_framework_rpk_v3_0.R`) in an interactive Shiny dashboard. The framework
is bundled **unchanged** in `inst/framework/`; the package loads it at runtime
and drives its `run_simulation()` engine through all six phases.

## Install

# Install from Github
```r
# Install from GitHub
install.packages("remotes")

remotes::install_github(
  "bkkongyir1106/Academia",
  subdir = "R_pkg/pretestsim",
  dependencies = TRUE
)

library(pretestsim)
run_app()
```


App dependencies (`shiny`, `bslib`, `DT`) install automatically. The framework's
own packages are *Suggests* and are installed on first use when possible
(`MASS`, `nortest`, `moments`, `tseries`, `DescTools`, `coin`, `Rfit`,
`LaplacesDemon`, `VGAM`, `evd`); a feature that needs a missing package is simply
unavailable (e.g. the Pareto distribution needs `VGAM`).

## Run

```r
library(pretestsim)
run_app()                                  # uses the bundled framework
run_app(framework_path = "~/my_copy.R")    # or point at your own copy
```

You can also use the framework programmatically:

```r
fw  <- load_framework()
res <- do.call(fw$run_simulation, c(fw$onesample_fns, list(
  Nsim = 200, N_tradeoff = 300, test_type = "demo",
  distributions = c("exponential", "normal"),
  norm_config = list(method = "classical", config = list(norm_test = "SW")),
  threshold_grid = seq(0, 1, 0.025), tol_pos = 0.005, loss_tol = 0.01,
  test_alpha = 0.05, center_by = "median", effect_size = 0.5,
  sample_sizes = c(10, 20, 30, 40, 50), single_n = 20,
  effect_sizes_plot = c(0, 0.2, 0.5, 0.8), sig_levels = seq(0, 1, 0.005),
  norm_test = "SW", selected_tests = "SW", phases = c(1, 2, 4))))
```

## How the framework is used (unmodified)

`load_framework()` parses the bundled file and evaluates every function
definition, every `norm_config_*` preset and every `*_fns` / `*_params` bundle
exactly as written. It skips only:

- machine-specific side-effects: `setwd()`, the external ML `source()`, the
  `pacman::p_load()` block, and `load("…trained_models.RData")`;
- the example `run_simulation()` calls at the end of the file (so the app does
  not auto-launch long simulations on startup);
- `save.image()` inside Phase 6 is redirected to a no-op; `save()` of the
  compact results object still runs.

## Dashboard workflow

Everything lives in the always-visible sidebar:

- **Framework file** — confirm the auto-detected (bundled) framework, or set a
  path / upload another copy.
- **Normality approach** — *Classical* (pick the pretest + the Phase-1 ROC
  battery) or *Custom / ML*. Custom: paste an R function returning a single
  named numeric scalar. ML: give the **path** to your ML R file (sourced into
  the framework's `ml_env`) and the **path** to the `.RData` holding
  `trained_models` — loading by path avoids the upload limit; the scorers
  `ml_fn_roc` / `ml_fn_decision` are wired automatically.
- **Test family & data** — one-sample / two-sample / ANOVA / regression, the
  non-parametric alternative, the H1 distribution, `center_by`, `effect_size`.
- **Simulation parameters** — `Nsim`, `N_tradeoff`, `test_alpha`, `single_n`,
  `sample_sizes`, `effect_sizes_plot`, `tol_pos`, `loss_tol`, the
  `threshold_grid` (from/to/by), `sig_levels` step, `per_n_thresholds`.
- **Phases & run** — tick any subset of phases 1-6. If Phases 3/4/5 are
  requested while Phase 2 is skipped, a `pretest_threshold` field appears and is
  required (scalar or per-n pairs) — exactly as the framework enforces.

### Phase result tabs (rendered with the framework's own plot functions)

1. **Normality ROC** — AUROC curves comparing tests + an AUC-per-test table.
2. **Trade-off** — the trade-off plot, the optimal threshold (`param_star`),
   power gain/loss, Type I inflation.
3. **Power-Error ROC** — power vs Type I error ROC-like curves.
4. **Power & Type I vs n** — across sample sizes + AUC-Power / AUC-Type I tables.
5. **Power vs effect** — power vs effect size at the fixed `single_n`.
6. **Log & Downloads** — console output, loader report, resolved arguments, and
   download buttons for the Phase-6 PDFs and RData.

## License

MIT.
