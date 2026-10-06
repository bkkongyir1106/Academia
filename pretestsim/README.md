# pretestsim 1.0.0 — Paper A edition

A bundled R simulation framework with an optional Shiny dashboard. No framework
code or model files are downloaded at runtime. Loading does not install packages.

From the publication folder, install the source package:

```sh
R CMD INSTALL pretestsim
```

Install the dashboard dependencies first if necessary:

```r
install.packages(c("shiny", "bslib", "DT"))
```

Launch the app:

```r
pretestsim::run_app()
```

The app opens in a ready state. Select the comparison, normality test, phases,
and simulation settings, then click **Run simulation**. Editing settings and
loading the framework do not run simulations or display simulation progress.
Custom mode offers **Fisher SW + AD**. Classical tests remain available.

Use the framework directly:

```r
fw <- pretestsim::load_framework()
set.seed(12345)
res <- fw$run_simulation(
  Nsim = 100, N_tradeoff = 100,
  sample_sizes = c(10, 20), threshold_grid = seq(0, 1, .1),
  phases = c(2, 4), output_dir = "example_results"
)
```

Optional procedure dependencies are listed in DESCRIPTION. Missing ones produce
an installation instruction when the relevant procedure is called. The main
paper runners also need nortest, moments and tseries; the supplement uses MASS,
LaplacesDemon, evd and Rfit. Other supported framework options may need additional
Suggested packages.

The framework file in `inst/framework/framework.R` is a release copy of
`../R/framework.R`. Use `../scripts/build_release.R` to synchronize and build.
The public API is `load_framework()` and `run_app()`; the former returns the
full function environment for programmatic analyses. Old remote-asset arguments
and helpers were removed in this publication edition.

Fisher combines the SW and AD scores from the same sample. They are dependent,
so its chi-square reference score is not guaranteed to be a calibrated uniform
null p-value. Evaluate its routing performance through the framework. It is not
used for the manuscript's reported classical-test results.
