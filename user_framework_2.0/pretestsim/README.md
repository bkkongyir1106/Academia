# pretestsim

`pretestsim` is an R package for running the adaptive normality-pretesting
simulation framework. It can be used in two ways:

- through the interactive Shiny dashboard;
- directly from R scripts without opening Shiny.

The package is intentionally lightweight. At runtime, it downloads and caches
the framework R code, machine-learning support code, and trained model file from
the GitHub repository.

## Installation from GitHub

Install the package from the `user_framework_2.0/pretestsim` folder of the
GitHub repository:

```r
install.packages("remotes")

remotes::install_github(
  "bkkongyir1106/Academia",
  subdir = "user_framework_2.0/pretestsim",
  force = TRUE,
  upgrade = "never"
)
```

Load the package and check that the correct version and GitHub asset path are
being used:

```r
library(pretestsim)

packageVersion("pretestsim")
pretestsim::asset_url("framework")
```

The framework asset URL should contain:

```text
/main/user_framework_2.0/user_framework_rpkg.R
```

It should not contain:

```text
/main/R_pkg/user_framework_2.0/
```

## Clean Reinstallation

If an older version of the package was installed previously, remove it and clear
the local cache before reinstalling:

```r
remove.packages("pretestsim")
unlink(tools::R_user_dir("pretestsim", "cache"), recursive = TRUE, force = TRUE)

remotes::install_github(
  "bkkongyir1106/Academia",
  subdir = "user_framework_2.0/pretestsim",
  force = TRUE,
  upgrade = "never"
)

library(pretestsim)
packageVersion("pretestsim")
```

## Required GitHub Files

The GitHub repository should contain these files in the top-level
`user_framework_2.0` folder:

```text
user_framework_2.0/
  pretestsim/
  user_framework_rpkg.R
  ML_framework_fun.R
  trained_models_multiple_sample_sizes.RData
  Run_framework.R
  run_ml_framework.R
```

The trained model file is large, so it may need to be pushed with GitHub Desktop
or command-line git rather than uploaded through the GitHub web interface.

## Running the Shiny Dashboard

To open the dashboard:

```r
library(pretestsim)

pretestsim::run_app(refresh_assets = TRUE)
```

Inside the dashboard:

1. Check **Refresh cached assets** when you want to force a fresh download from
   GitHub.
2. Click **Load framework**.
3. Select the normality approach, downstream test family, simulation settings,
   and phases.
4. Click **Run simulation**.

The first run may take longer because the framework files and trained ML model
are downloaded and cached locally.

## Running Without the Shiny Dashboard

You can also run the framework directly from an R script. This is useful for
reproducible runs, batch jobs, manuscript results, or running on a remote
machine.

First, load the framework from GitHub:

```r
library(pretestsim)

fw <- pretestsim::load_framework(
  refresh = TRUE,
  raw_base = "https://raw.githubusercontent.com/bkkongyir1106/Academia/main/user_framework_2.0"
)
```

The object `fw` is an environment containing the framework functions. For
example:

```r
ls(fw)
```

## Example: One-Sample t-Test vs Sign Test

This example runs the one-sample demonstration comparing the one-sample t-test
with the sign test, using the Shapiro-Wilk normality pretest for adaptive
routing.

```r
library(pretestsim)

fw <- pretestsim::load_framework(
  refresh = TRUE,
  raw_base = "https://raw.githubusercontent.com/bkkongyir1106/Academia/main/user_framework_2.0"
)

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

  threshold_grid    = seq(0, 1, by = 0.025),
  sample_sizes      = c(10, 20, 30, 40, 50),
  effect_sizes_plot = c(0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0),

  phases = 1:5,

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

The generated files are written to:

```text
onesample_ttest_vs_sign_results/
```

For a quick test run, reduce the simulation sizes:

```r
Nsim = 100
N_tradeoff = 100
```

For final analysis runs, increase them:

```r
Nsim = 10000
N_tradeoff = 10000
```

## Example: One-Sample t-Test vs Randomized Sign Test

The same framework can use another downstream alternative by changing
`fn_for_ds_test_2`:

```r
res_randomized <- fw$run_simulation(
  Nsim       = 1000,
  N_tradeoff = 1000,
  test_type  = "onesample_ttest_vs_randomized_sign_SW",

  distributions = c("exponential", "normal"),

  norm_config = list(
    method = "classical",
    config = list(norm_test = "SW")
  ),

  gen_data           = fw$onesample_data,
  get_parameters     = fw$onesample_parameters,
  fn_to_get_norm_obj = fw$raw_data,
  fn_for_ds_test_1   = fw$one_sample_t_test,
  fn_for_ds_test_2   = fw$randomized_sign_pvalue,

  effect_size = 0.5,
  center_by   = "median",
  test_alpha  = 0.05,
  sample_sizes = c(10, 20, 30, 40, 50),
  phases = 1:5,
  output_dir = "onesample_ttest_vs_randomized_sign_results",
  save_results = TRUE
)
```

## Example: Using Machine Learning Without Shiny

To use the trained ML normality model directly from R:

```r
library(pretestsim)

fw <- pretestsim::load_framework(refresh = TRUE)

ml_resources <- pretestsim::load_ml_resources_from_github(
  framework_env = fw,
  refresh = TRUE,
  required_sample_sizes = c(10, 20, 30, 40, 50)
)

ml_configs <- fw$make_ml_normality_configs(
  ml_resources = ml_resources,
  chosen_model = "SVM",
  models_for_roc = c("RF", "GBM", "ANN", "SVM"),
  use_majority_vote = TRUE
)

res_ml <- fw$run_simulation(
  Nsim       = 1000,
  N_tradeoff = 1000,
  test_type  = "onesample_ttest_vs_sign_ml_svm",

  distributions = c("exponential", "normal"),
  norm_config   = ml_configs$decision,

  gen_data           = fw$onesample_data,
  get_parameters     = fw$onesample_parameters,
  fn_to_get_norm_obj = fw$raw_data,
  fn_for_ds_test_1   = fw$one_sample_t_test,
  fn_for_ds_test_2   = fw$sign_test,

  effect_size = 0.5,
  center_by   = "median",
  test_alpha  = 0.05,
  sample_sizes = c(10, 20, 30, 40, 50),
  phases = 1:5,
  output_dir = "onesample_ttest_vs_sign_ml_svm_results",
  save_results = TRUE
)
```

## Local Cache

Downloaded GitHub assets are cached here:

```r
tools::R_user_dir("pretestsim", "cache")
```

Clear the cache when GitHub files have changed:

```r
unlink(tools::R_user_dir("pretestsim", "cache"), recursive = TRUE, force = TRUE)
```

Or request a refresh when loading the framework:

```r
fw <- pretestsim::load_framework(refresh = TRUE)
```

## Common Troubleshooting

If the old dashboard or old GitHub asset path appears, restart R, remove the old
package, clear the cache, and reinstall from GitHub.

If `load_framework()` gives an error about `.default_raw_base`, call the package
function explicitly:

```r
fw <- pretestsim::load_framework(refresh = TRUE)
```

If ML does not load, confirm that this URL exists:

```text
https://raw.githubusercontent.com/bkkongyir1106/Academia/main/user_framework_2.0/trained_models_multiple_sample_sizes.RData
```
