# pretestsim

`pretestsim` is an R package for evaluating the downstream utility of adaptive
normality pretesting. It compares a fixed parametric procedure, a fixed
alternative procedure, and an adaptive procedure that chooses between them.

You can use it in two ways:

- through the interactive Shiny dashboard;
- directly from R scripts, including reproducible analyses and cluster jobs.

The Paper A edition bundles the framework with the package. Once installed,
loading it does not download framework files or install other packages.
Classical normality procedures are supported, and **Fisher SW + AD** is the
only built-in custom normality option.

## Requirements

Use R 4.1.0 or later. RStudio is optional but convenient for interactive use.

Install the dashboard dependencies and the packages used by the examples below:

```r
install.packages(c(
  "shiny", "bslib", "DT",
  "nortest", "moments", "tseries"
))
```

Some downstream procedures and distributions need additional packages. For the
additional Paper A applications, install:

```r
install.packages(c("MASS", "LaplacesDemon", "evd", "Rfit"))
```

The full list of optional dependencies is in `DESCRIPTION`. If an optional
package is missing, install the package named in the error message and rerun.

## Installation from the Publication Archive

Extract `Paper_A_publication_code.zip`. The installable source package is:

```text
Paper_A_publication/validation/pretestsim_1.0.0.tar.gz
```

In R or RStudio, run the following and select that `.tar.gz` file:

```r
install.packages(
  file.choose(),
  repos = NULL,
  type = "source"
)
```

Alternatively, supply the path directly. Replace the example path with the
location of the file on your computer:

```r
install.packages(
  "path/to/pretestsim_1.0.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

Load the package and check the version:

```r
library(pretestsim)

packageVersion("pretestsim")
# Expected for this release: '1.0.0'
```

## Installation from GitHub

The existing project uses the repository `bkkongyir1106/Academia`. **After the
publication edition has been uploaded to `user_framework_2.0/pretestsim`**, users
can install it with:

```r
install.packages("remotes")

remotes::install_github(
  "bkkongyir1106/Academia",
  subdir = "user_framework_2.0/pretestsim",
  force = TRUE,
  upgrade = "never"
)

library(pretestsim)
packageVersion("pretestsim")
```

The `subdir` must match where the package is published. For example, if the
publication folder is uploaded unchanged at the repository root, use
`subdir = "Paper_A_publication/pretestsim"` instead. These commands describe the
publication layout; preparing the local code bundle does not update GitHub.
Use the accompanying source archive to install this version independently of
the repository's current contents.

## Reinstalling or Updating

If an earlier version is already loaded, restart R first. In RStudio, choose
**Session → Restart R**. Then reinstall the new source archive:

```r
install.packages(file.choose(), repos = NULL, type = "source")

library(pretestsim)
packageVersion("pretestsim")
find.package("pretestsim")
```

Installing into the same R library replaces the earlier package version.
The bundled edition does not use the old downloaded-asset cache, so clearing
that cache is unnecessary. Do not pass the old `refresh`, `raw_base`, or
`refresh_assets` arguments to this version.

## Package Files

The package source contains the files needed to load the framework and launch
the dashboard:

```text
pretestsim/
  DESCRIPTION
  NAMESPACE
  LICENSE
  README.md
  NEWS.md
  R/
    app.R
    loader.R
    run_app.R
    pretestsim-package.R
  inst/
    framework/
      framework.R
  man/
  tests/
```

Install the whole package rather than sourcing the files in `R/` individually.
No external framework file or trained-model file is needed for normal use.

## Running the Shiny Dashboard

Open the dashboard with:

```r
library(pretestsim)

pretestsim::run_app()
```

Inside the dashboard:

1. Wait for **Loaded bundled framework. Ready to run.**
2. Choose **Classical** or **Custom** under **Normality approach**.
3. For Classical, select the adaptive-routing pretest and the Phase 1 ROC
   battery. For Custom, select **Fisher SW + AD**.
4. Choose the downstream comparison, alternative test, distribution, centering
   convention, and effect size.
5. Set the simulation parameters and select the phases to run.
6. Click **Run simulation**.
7. Inspect the phase tabs, and use **Log & Downloads** to save the results.

The app stays idle when opened and when settings are changed. The simulation
and running indicator start only after **Run simulation** is clicked. Loading
the framework again with **Load framework** does not run an analysis.

### First Test Run

Start with small settings to confirm the installation works:

| Setting | Value |
| --- | --- |
| Downstream comparison | One-sample |
| Alternative downstream test | Sign test |
| Normality approach | Classical |
| Adaptive-routing pretest | SW |
| Evaluation repetitions (`Nsim`) | 100 |
| Calibration repetitions (`N_tradeoff`) | 100 |
| Sample sizes | `10,20` |
| Threshold-grid increment | `0.1` |
| Selected phases | 1, 2 and 4 |

After the run finishes, inspect the ROC plot, selected threshold, and the power
and Type I error plots and tables. Repeat with **Custom → Fisher SW + AD** to
check that option. Small runs are software checks; their numerical results are
not precise enough for publication.

### Simulation Phases

| Phase | Result |
| --- | --- |
| 1 | Normality-test ROC curves |
| 2 | Threshold trade-off analysis and selected routing threshold |
| 3 | Downstream power versus Type I error curves |
| 4 | Power, Type I error, and AUC summaries across sample sizes |
| 5 | Power across effect sizes at a fixed sample size |

Phases 3–5 need a routing threshold. Include Phase 2 to select it automatically,
or supply a threshold when Phase 2 is omitted. The dashboard displays the
threshold input when it is needed.

The dashboard writes results into a temporary session directory. Download the
files you need before closing the R session. Use the scripted workflow below
for results saved directly to a folder you choose.

## Running Without the Shiny Dashboard

Load the bundled functions:

```r
library(pretestsim)

fw <- pretestsim::load_framework()
```

`fw` is an environment containing the framework functions. You can inspect a
function's arguments before using it:

```r
args(fw$run_simulation)
args(fw$generate_data)
```

Use `set.seed()` before a simulation to make a run reproducible under the same
R environment, settings, and package versions.

## Example: One-Sample t-Test versus Sign Test

This example uses Shapiro–Wilk normality pretesting to choose between the
one-sample t-test and the Sign test. It runs the normality ROC, threshold
selection, and sample-size evaluation phases.

```r
library(pretestsim)
fw <- pretestsim::load_framework()

set.seed(12345)

res <- fw$run_simulation(
  Nsim       = 1000,
  N_tradeoff = 1000,

  test_type     = "onesample_ttest_vs_sign_SW",
  distributions = c("exponential", "normal"),

  norm_config = list(
    method = "classical",
    config = list(norm_test = "SW")
  ),

  gen_data           = fw$onesample_data,
  get_parameters     = fw$onesample_parameters,
  fn_to_get_norm_obj = fw$raw_data,
  fn_for_ds_test_1   = fw$one_sample_t_test,
  fn_for_ds_test_2   = fw$sign_test,

  effect_size = 0.5,
  center_by   = "median",
  test_alpha  = 0.05,
  tol_pos     = 0.005,
  loss_tol    = 0.01,

  threshold_grid = seq(0, 1, by = 0.025),
  sample_sizes   = c(10, 20, 30, 40, 50),
  single_n       = 10,
  phases         = c(1, 2, 4),

  output_dir   = "onesample_ttest_vs_sign_results",
  save_results = TRUE
)
```

Inspect the returned object:

```r
names(res)
res$files
res$objects$optimal_result
res$objects$auc_tables
```

The files are saved in a subfolder named after `test_type`:

```text
onesample_ttest_vs_sign_results/
  onesample_ttest_vs_sign_SW/
```

Relative output paths are resolved against `getwd()`. Use an absolute path in
`output_dir` if you want a particular location. Reusing the same output folder
and `test_type` replaces files from the previous run.

Here, median centering defines the one-sample null. Under an asymmetric
population, the t-test's mean target can differ from that median target;
interpret its Type I error accordingly, as discussed in Paper A.

## Example: Two-Sample Welch t-Test versus Mann–Whitney U

Use the two-sample generators and downstream procedures:

```r
fw <- pretestsim::load_framework()
set.seed(22345)

res_two <- fw$run_simulation(
  Nsim = 1000,
  N_tradeoff = 1000,
  test_type = "twosample_welch_vs_mwu_SW",
  distributions = c("exponential", "normal"),
  norm_config = list(method = "classical", config = list(norm_test = "SW")),

  gen_data           = fw$two_sample_data,
  get_parameters     = fw$twosample_parameters,
  fn_to_get_norm_obj = fw$raw_data,
  fn_for_ds_test_1   = fw$twosample_t_test,
  fn_for_ds_test_2   = fw$Mann_whitney_U_test,

  effect_size = 0.5,
  center_by = "median",
  test_alpha = 0.05,
  tol_pos = 0.005,
  loss_tol = 0.01,
  threshold_grid = seq(0, 1, by = 0.025),
  sample_sizes = c(10, 20, 30, 40, 50),
  single_n = 10,
  phases = c(1, 2, 4),
  output_dir = "twosample_results",
  save_results = TRUE
)
```

Each sample size is **per group**. The parametric branch is selected only when
both samples pass the normality pretest. Under the common-shape, common-scale
location-shift model used here, the two procedures address the same median-shift
null; there is no parameter mismatch in this design.

## Example: One-Sample t-Test versus Randomized Sign Test

Change the alternative downstream function to `fw$randomized_sign_pvalue`:

```r
fw <- pretestsim::load_framework()
set.seed(32345)

res_randomized <- fw$run_simulation(
  Nsim = 1000,
  N_tradeoff = 1000,
  test_type = "onesample_ttest_vs_randomized_sign_SW",
  distributions = c("exponential", "normal"),
  norm_config = list(method = "classical", config = list(norm_test = "SW")),

  gen_data           = fw$onesample_data,
  get_parameters     = fw$onesample_parameters,
  fn_to_get_norm_obj = fw$raw_data,
  fn_for_ds_test_1   = fw$one_sample_t_test,
  fn_for_ds_test_2   = fw$randomized_sign_pvalue,

  effect_size = 0.5,
  center_by = "median",
  test_alpha = 0.05,
  tol_pos = 0.005,
  loss_tol = 0.01,
  threshold_grid = seq(0, 1, by = 0.025),
  sample_sizes = c(10, 20, 30, 40, 50),
  phases = c(2, 4),
  output_dir = "randomized_sign_results",
  save_results = TRUE
)
```

This is an additional usage example; it does not replace the ordinary Sign test
in the reported Paper A comparison.

## Example: Fisher SW + AD without Shiny

Use the bundled Fisher combination through the custom interface:

```r
fw <- pretestsim::load_framework()
set.seed(42345)

fisher_config <- list(
  method = "custom",
  config = list(
    fn = fw$fisher_combined,
    label = "Fisher SW + AD"
  )
)

res_fisher <- fw$run_simulation(
  Nsim = 1000,
  N_tradeoff = 1000,
  test_type = "onesample_ttest_vs_sign_Fisher",
  distributions = c("exponential", "normal"),
  norm_config = fisher_config,

  gen_data           = fw$onesample_data,
  get_parameters     = fw$onesample_parameters,
  fn_to_get_norm_obj = fw$raw_data,
  fn_for_ds_test_1   = fw$one_sample_t_test,
  fn_for_ds_test_2   = fw$sign_test,

  effect_size = 0.5,
  center_by = "median",
  test_alpha = 0.05,
  tol_pos = 0.005,
  loss_tol = 0.01,
  threshold_grid = seq(0, 1, by = 0.025),
  sample_sizes = c(10, 20, 30, 40, 50),
  single_n = 10,
  phases = c(1, 2, 4),
  output_dir = "fisher_results",
  save_results = TRUE
)
```

The SW and AD p-values come from the same sample and are dependent. Treat their
Fisher combination as a routing score whose performance is evaluated through
the simulations; its chi-square reference is not guaranteed to yield an exactly
uniform null p-value.

## Choosing Simulation Sizes

`N_tradeoff` controls threshold calibration. `Nsim` controls the other requested
evaluation phases. These are different from the sample size of each simulated
dataset.

| Purpose | `Nsim` | `N_tradeoff` |
| --- | ---: | ---: |
| Installation check | 100 | 100 |
| Exploratory example | 1,000 | 1,000 |
| Paper A representative applications | 100,000 | 1,000,000 |

For the publication settings, also use the manuscript threshold grid with
increments of `0.005`. Larger runs take longer; permutation procedures also
perform internal resampling within each repetition. For the exact Paper A
settings, use the dedicated publication scripts instead of adapting a quick
example by changing only the repetition counts.

## Reproducing Paper A Analyses

The complete publication archive contains scripts outside the installed package:

```text
Paper_A_publication/
  R/
  config/
  scripts/
    run_main.R
    run_supplement.R
    run_broad.R
    run_assessment.R
    recalculate_saved_assessment.R
    export_tables.R
  docs/
    result_index.csv
    analysis_notes.md
  pretestsim/
```

From that folder, run these commands in a terminal:

```sh
# List main-paper jobs without running simulations.
Rscript scripts/run_main.R plan

# Run a short software check.
Rscript scripts/run_main.R smoke

# Run the main-paper applications at publication effort.
Rscript scripts/run_main.R publication
```

The publication folder's README covers the supplementary, broader-distribution,
and strict assessment analyses. `docs/result_index.csv` maps manuscript tables
and figures to the corresponding scripts. These scripts specify the manuscript
settings explicitly; the general package defaults serve exploratory use.

## Common Troubleshooting

### The old dashboard still appears

Restart R, reinstall the publication source archive, and check:

```r
packageVersion("pretestsim")
find.package("pretestsim")
.libPaths()
```

If several R libraries contain `pretestsim`, verify that R is loading the
intended installation. Version `1.0.0` uses the bundled framework and does not
show the old remote-asset or model-loading controls.

### An unused-argument error mentions refresh or a GitHub URL

Those arguments belonged to the previous remote-loading version. Use:

```r
fw <- pretestsim::load_framework()
pretestsim::run_app()
```

### A phase tab is empty

Check that the phase was selected, then click **Run simulation**. Inspect
**Log & Downloads** for errors. If Phase 2 is omitted while running Phases 3–5,
supply a routing threshold.

### Results vary between runs

Use the same `set.seed()`, settings, R version and dependency versions. Small
Monte Carlo runs can produce unstable thresholds and performance estimates.
Increasing repetitions improves precision but does not increase the sample
size of each dataset.

### Results cannot be found after closing the dashboard

Dashboard files are temporary. Download them before closing the R session, or
use `fw$run_simulation()` with an explicit `output_dir` for persistent outputs.
Inspect `res$files` for the exact generated paths.

### A run is slow

First confirm the workflow with 100 repetitions and fewer sample sizes. Run only
the phases needed. Use the publication scripts for long batch runs, and keep the
final settings consistent with the manuscript.
