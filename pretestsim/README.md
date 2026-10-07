# pretestsim

**Simulation-based evaluation of normality pretesting in statistical inference**

`pretestsim` provides a framework for studying how normality pretesting affects
subsequent statistical inference. It compares fixed parametric and alternative
procedures with an adaptive procedure that selects between them using a
normality test. Simulation results describe Type I error, power, and the
conditions under which pretesting improves or reduces performance.

The package supports one-sample comparisons, two-sample comparisons, one-way
analysis of variance, and linear regression. Users can specify the generating
distribution, sample sizes, effect sizes, downstream procedures, and normality
assessment. Functions are provided for threshold selection, ROC analysis,
area-under-the-curve summaries, and graphical presentation of results.

Analyses can be conducted through an interactive Shiny application or directly
in R. The examples below illustrate common applications; the framework can be
configured for other simulation designs through its data-generation and testing
interfaces.

## Installation

`pretestsim` requires R (>= 4.1.0). RStudio is optional.

### Dependencies

Install the application dependencies and the packages used in the normality-test
examples:

```r
install.packages(c(
  "shiny", "bslib", "DT",
  "nortest", "moments", "tseries"
))
```

Additional distributions and downstream procedures may require:

```r
install.packages(c("MASS", "LaplacesDemon", "evd", "Rfit"))
```

The complete dependency list is given in `DESCRIPTION`. Optional packages are
required only for the procedures that use them.

### From a source archive

Download or locate `pretestsim_1.0.0.tar.gz`, then install it from R:

```r
install.packages(
  "path/to/pretestsim_1.0.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

Replace the path with the archive's location. In an interactive session, the
archive can also be selected using a file dialog:

```r
install.packages(file.choose(), repos = NULL, type = "source")
```

Check the installation:

```r
library(pretestsim)
packageVersion("pretestsim")
# [1] '1.0.0'
```

The framework is included in the package. An internet connection is not needed
to load the installed framework or run simulations.

### From GitHub

The repository installation command is:

```r
install.packages("remotes")

remotes::install_github(
  "bkkongyir1106/Academia",
  subdir = "user_framework_2.0/pretestsim",
  force = TRUE,
  upgrade = "never"
)
```

This command installs the version available at that repository path. Version
1.0.0 must be uploaded there before it can be installed through GitHub; use the
source archive to install the accompanying release. If the package is published
at a different location within the repository, adjust `subdir` accordingly.

### Updating an existing installation

Restart R before replacing a package that is already loaded. In RStudio, select
**Session → Restart R**, then reinstall the source archive or use the GitHub
command above. Confirm the version and library location:

```r
library(pretestsim)
packageVersion("pretestsim")
find.package("pretestsim")
```

Version 1.0.0 loads its bundled framework. The `refresh`, `raw_base`, and
`refresh_assets` arguments from earlier versions are no longer used.

## Interactive application

Launch the Shiny application with:

```r
library(pretestsim)
pretestsim::run_app()
```

To configure an analysis:

1. Select **Classical** or **Custom** under **Normality approach**. Classical
   mode provides standard normality tests; Custom mode currently offers
   **Fisher SW + AD**.
2. Select the downstream comparison, alternative procedure, generating
   distribution, centering convention, and effect size.
3. Specify the simulation parameters and analytical phases.
4. Click **Run simulation**.
5. View the results in the phase tabs and download files from **Log & Downloads**.

The framework loads automatically when the application opens. Simulations begin
only when **Run simulation** is selected.

### A short introductory run

The following settings provide a quick installation check:

| Setting | Value |
| --- | --- |
| Downstream comparison | One-sample |
| Alternative procedure | Sign test |
| Normality approach | Classical |
| Adaptive-routing pretest | SW |
| Evaluation repetitions (`Nsim`) | 100 |
| Calibration repetitions (`N_tradeoff`) | 100 |
| Sample sizes | `10,20` |
| Threshold-grid increment | `0.1` |
| Phases | 1, 2 and 4 |

This run produces normality ROC curves, a selected routing threshold, and power
and Type I error summaries. The small number of repetitions is intended for
checking the workflow. Repeat the run with **Custom → Fisher SW + AD** to examine
the custom normality option.

### Analytical phases

| Phase | Analysis |
| --- | --- |
| 1 | Normality-test ROC curves |
| 2 | Threshold selection based on Type I error and power |
| 3 | Downstream power versus Type I error curves |
| 4 | Power, Type I error, and AUC summaries across sample sizes |
| 5 | Power across effect sizes at a fixed sample size |

Phases 3–5 require a routing threshold. Include Phase 2 to select a threshold,
or supply one directly when Phase 2 is omitted. The application displays the
threshold input when needed.

Application outputs are stored in a temporary session directory. Download the
required files before closing the R session. Scripted analyses can instead write
to a persistent directory specified by `output_dir`.

## Programmatic use

Load the framework functions into an environment:

```r
library(pretestsim)
fw <- pretestsim::load_framework()
```

Inspect the available arguments with:

```r
args(fw$run_simulation)
args(fw$generate_data)
```

The main simulation function accepts functions for data generation, parameter
specification, extraction of the normality-assessment object, and the two
downstream procedures. This interface allows the same workflow to be used across
inferential settings. Use `set.seed()` to reproduce a simulation under the same
settings, R version, and dependency versions.

## Examples

### One-sample t-test and Sign test

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
the comparison therefore evaluates the consequences of using a mean-based
procedure for a median-based inferential target.

### Two-sample Welch t-test and Mann–Whitney U test

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
null. This interpretation depends on the common location-shift structure.

### Randomized Sign test

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

### Fisher combination of Shapiro–Wilk and Anderson–Darling tests

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

The Shapiro–Wilk and Anderson–Darling p-values are dependent because they are
computed from the same sample. The Fisher combination is evaluated as a routing
score; its chi-square reference does not imply exact null calibration.

## Simulation settings and reproducibility

`N_tradeoff` is the number of repetitions used to select the routing threshold.
`Nsim` controls the other requested evaluation phases. `sample_sizes` specifies
the size of each simulated dataset, rather than the number of repetitions.

| Purpose | Example `Nsim` | Example `N_tradeoff` |
| --- | ---: | ---: |
| Installation check | 100 | 100 |
| Exploratory analysis | 1,000 | 1,000 |
| Larger simulation study | 100,000 | 1,000,000 |

These values illustrate different computational budgets. Choose repetition
counts according to the precision needed for the operating characteristics and
threshold comparison. Increasing repetitions reduces Monte Carlo variability;
it does not increase the sample size within each dataset.

The threshold grid determines the resolution of threshold selection. For
example, `seq(0, 1, by = 0.005)` evaluates thresholds in increments of 0.005.
Permutation procedures require additional resampling within each repetition and
can take substantially longer than analytical tests.

For reproducible studies, retain the seed, simulation settings, result objects,
and software versions:

```r
sessionInfo()
saveRDS(res, file = "simulation_results.rds")
```

The full set of generated paths is available in `res$files`. A saved result can
be read with:

```r
res <- readRDS("simulation_results.rds")
```

## Citation

If you use `pretestsim` in your work, please cite the package:

> Kongyir, B., & Rudra, P. (2026). *pretestsim: Simulation Framework for Adaptive Normality
> Pretesting*. R package version 1.0.0.
> https://github.com/bkkongyir1106/Academia

To obtain the citation and BibTeX entry for the installed version, run:

```r
citation("pretestsim")
toBibtex(citation("pretestsim"))
```

## Authors

- **Benedict Kongyir** — author and maintainer.
- **Pratyaydipta Rudra** — author.

Maintainer contact: [kongyirbenkk@gmail.com](mailto:kongyirbenkk@gmail.com)

## Support

Report bugs and feature requests through the
[GitHub issue tracker](https://github.com/bkkongyir1106/Academia/issues).
Include a minimal reproducible example, the error message, and the output of
`sessionInfo()`.

### An earlier version is still loaded

Restart R and reinstall the intended source archive. Check the active package
and library paths:

```r
packageVersion("pretestsim")
find.package("pretestsim")
.libPaths()
```

When several libraries contain the package, ensure that R is loading the intended
installation. Version 1.0.0 includes the framework source within the package.

### The framework source cannot be found

Check the installed version and reinstall the complete source archive after
restarting R. Install the package as a whole, including `inst/framework`, rather
than copying or sourcing its individual application files.

### An argument from an earlier version is rejected

Use the current loading and application functions without remote-asset options:

```r
fw <- pretestsim::load_framework()
pretestsim::run_app()
```

### A result tab is empty

Confirm that the corresponding phase was selected and inspect **Log & Downloads**
for errors. When Phase 2 is omitted, Phases 3–5 require a supplied threshold.

### Results differ between runs

Use the same seed and settings in the same software environment. With few
repetitions, Monte Carlo variability can affect both threshold selection and
estimated power or Type I error. Increase the repetitions when greater precision
is required.

### A simulation takes too long

Begin with a smaller number of repetitions and fewer sample sizes, and select
only the required phases. For longer analyses, use an R script with a persistent
output directory.

## License

`pretestsim` is distributed under the MIT license. See [LICENSE](LICENSE) for
copyright information.
