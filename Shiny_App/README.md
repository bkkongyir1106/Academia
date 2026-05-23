# Pretest Simulation Studio — phase-driven dashboard

An interactive **R Shiny** front-end that runs your **original, unmodified**
framework (`user_framework_rpk_v3_0.R`) through the six phases of its
`run_simulation()` engine.

## Files (keep together in one folder)

| File                          | Purpose                                                              |
|-------------------------------|---------------------------------------------------------------------|
| `app.R`                       | The dashboard (UI + server).                                        |
| `framework_loader.R`          | Sources your framework as-is into an environment.                   |
| `user_framework_rpk_v3_0.R`   | **Your original file — unchanged.** Place it next to `app.R`.       |

> The app auto-detects `user_framework_rpk_v3_0.R` (also `v3.0` / `rpkg`
> spellings) in the app folder, its parent, a sibling `one_sample/` folder, and
> `/mnt/project/`. If it isn't found, use the **Framework file** panel at the top
> of the sidebar to type the full path and click **Load / Reload**, or upload the
> `.R` file directly. The panel shows the load status (definitions kept/skipped).
> Nothing in your file is edited.

## How the loader treats your code

It parses the file and evaluates **every** function definition, every
`norm_config_*` preset and every `*_fns` / `*_params` bundle exactly as written.
It skips only the parts that are machine-specific or are example runs:

- `setwd(...)`, the external ML `source(...)`, the `pacman::p_load(...)` block,
  and `load("…trained_models.RData")` — these are environment-specific.
- the 51 example `run_simulation()` invocations at the end of the file (so the
  app does not auto-launch hours of simulation on startup).
- `save.image()` inside Phase 6 is redirected to a no-op (so a long app session
  is never dumped); `save()` of the compact `results` object still runs.

## Requirements

App: `shiny`, `bslib`, `DT` (auto-installed if missing).
Framework (used via `::`): `MASS`, `nortest`, `moments`, `tseries`, `DescTools`,
`coin`, `Rfit`, `LaplacesDemon`, `VGAM`, `evd` — the loader attempts to install
any that are missing. A test/distribution simply becomes unavailable if its
package is absent (e.g. the Pareto distribution needs `VGAM`).

```r
install.packages(c("shiny","bslib","DT",
                   "nortest","moments","tseries","DescTools",
                   "coin","Rfit","LaplacesDemon","VGAM","evd"))
```

## Run

```r
shiny::runApp("path/to/folder")
```

## Workflow (mirrors the framework)

The always-visible sidebar holds every control:

0. **Framework file** — confirm the detected path (or set/upload it) and load.
1. **Normality approach** — *Classical* (pick the pretest + the Phase-1 ROC
   battery) or *Custom / ML*. In Custom you paste an R function returning a
   single named numeric scalar (the Fisher SW+AD test is pre-filled); in **ML
   mode** you upload a `.RData` containing `trained_models` and the app wires it
   into the framework's `ml_fn_roc` (compare models, Phase 1) or `ml_fn_decision`
   (single model, Phases 2-6).
2. **Test family & data** — one-sample / two-sample / ANOVA / regression, the
   non-parametric alternative, the non-normal H1 distribution, `center_by`,
   `effect_size`.
3. **Simulation parameters** — `Nsim`, `N_tradeoff`, `test_alpha`, `single_n`,
   `sample_sizes`, `effect_sizes_plot`, `tol_pos`, `loss_tol`, the
   `threshold_grid` (from/to/by), `sig_levels` step, and `per_n_thresholds`.
4. **Phases & run** — tick any subset of phases 1-6.

**Phase dependency.** If you request Phases 3, 4 or 5 but skip Phase 2, the
sidebar reveals a required `pretest_threshold` field (scalar like `0.05`, or
per-n pairs like `10=0.05, 20=0.04`) — the same rule the framework enforces.

## What each result tab shows (rendered with the framework's own plot functions)

- **Phase 1** — AUROC curves comparing normality tests + an AUC-per-test table.
- **Phase 2** — the trade-off plot, the optimal threshold (`param_star`), power
  gain/loss, Type I inflation, and the selection note.
- **Phase 3** — power vs Type I error ROC-like curves at the optimal threshold.
- **Phase 4** — power and Type I error across sample sizes + the **AUC-Power /
  AUC-Type I** tables.
- **Phase 5** — power vs effect size at the fixed `single_n`.
- **Log & Downloads** — console output, the loader report, the resolved run
  arguments, and download buttons for the Phase-6 PDFs and `*_results.RData`.
