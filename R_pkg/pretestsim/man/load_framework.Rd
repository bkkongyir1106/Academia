\name{load_framework}
\alias{load_framework}
\alias{default_framework_path}
\title{Load the original framework into an environment, unmodified}
\usage{
load_framework(path = default_framework_path())
}
\arguments{
\item{path}{Path to the framework \code{.R} file. Defaults to the bundled copy.}
}
\value{An environment containing the framework's objects (notably
\code{run_simulation}), with a \code{"load_report"} attribute.}
\description{
Parses the framework source and evaluates every function definition and preset
exactly as written, skipping only machine-specific side-effects (\code{setwd},
the ML \code{source()}, the \code{pacman} block, the trained-models
\code{load()}) and the example \code{run_simulation()} calls.
}
\examples{
\dontrun{
fw <- load_framework()
res <- do.call(fw$run_simulation, c(fw$onesample_fns, list(
  Nsim = 100, N_tradeoff = 200, test_type = "demo",
  distributions = c("exponential", "normal"),
  norm_config = list(method = "classical", config = list(norm_test = "SW")),
  threshold_grid = seq(0, 1, 0.05), tol_pos = 0.01, loss_tol = 0.01,
  test_alpha = 0.05, center_by = "median", effect_size = 0.5,
  sample_sizes = c(10, 20, 30), single_n = 20,
  effect_sizes_plot = c(0, 0.3, 0.6), sig_levels = seq(0, 1, 0.05),
  norm_test = "SW", selected_tests = "SW", phases = c(1, 2, 4))))
}
}
