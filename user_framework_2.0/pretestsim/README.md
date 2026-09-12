# pretestsim

`pretestsim` is a GitHub-backed Shiny package for the adaptive normality
pretesting framework. The package itself is lightweight: it downloads and
caches the framework R code, ML support code, and trained model RData from the
package GitHub repository at runtime.

```r
library(pretestsim)
run_app()
```

By default, assets are resolved from:

```text
https://raw.githubusercontent.com/bkkongyir1106/Academia/main/R_pkg/user_framework_2.0
```

Override this when testing a branch, fork, or different subdirectory:

```r
options(pretestsim.github_raw_base =
  "https://raw.githubusercontent.com/bkkongyir1106/Academia/main/R_pkg/user_framework_2.0")
```

Main assets:

- `user_framework_rpkg.R`
- `ML_framework_fun.R`
- `trained_models_multiple_sample_sizes.RData`
- `Run_framework.R`
- `run_ml_framework.R`

The downloaded files are cached under `tools::R_user_dir("pretestsim", "cache")`
and can be refreshed from the app or with `cache_asset(..., refresh = TRUE)`.
