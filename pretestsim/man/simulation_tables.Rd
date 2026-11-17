\name{simulation_tables}
\alias{simulation_tables}
\title{Stored simulation p-values and reuse}
\description{
The bundled framework calculates normality scores and the two fixed-test
p-values once per simulated dataset, then derives adaptive p-values by routing
those values at the chosen threshold.
}
\details{
Load the framework with \code{fw <- load_framework()}. Its
\code{run_simulation()} function saves one \file{pvalue_tables.rds} file under
\file{<output_dir>/<test_type>/pvalue_tables/}. The path is returned as
\code{res$files$pvalue_tables}. Set \code{pvalue_directory} to use another location.

The saved list is organized by phase, then batch. Each batch contains \code{H0}
and \code{H1} data frames with four columns:
\describe{
\item{Normality_Score}{A list column of normality scores, retaining separate group scores.}
\item{Test1_P_Value}{P-value from the parametric procedure.}
\item{Test2_P_Value}{P-value from the alternative procedure.}
\item{Adaptive_P_Value}{The fixed-test p-value selected by the routing rule.}
}
Phase 1 calculates normality scores alone; its downstream columns contain
\code{NA}. Phase 2 fills adaptive values after choosing its threshold.
Phase-specific simulation batches are retained; reuse does not substitute draws
from one analytical phase for another.

To rebuild a run, restore its seed and random-number settings, retain its
simulation counts, generation functions and phase order, and set
\code{reuse_pvalues = TRUE}. Missing or incompatible batches cause an error.
New simulations cannot silently replace existing p-value files. Use a new
\code{test_type} or output directory for a new simulation.

The default \code{save_results = FALSE} omits an additional RData workspace;
the p-value file and figures are still saved. The returned result contains
summaries and plot data. The dashboard offers a seed setting and reuse of the
previous run in the current session.
}
\examples{
\dontrun{
fw <- load_framework()
args <- list(
  Nsim = 100, N_tradeoff = 100,
  test_type = "reuse_example", phases = c(2, 4),
  sample_sizes = c(10, 20), single_n = 10,
  output_dir = "simulation_results"
)
set.seed(12345)
res <- do.call(fw$run_simulation, args)
saved <- readRDS(res$files$pvalue_tables)
fw$pvalue_table_index(saved)
head(saved$phase4_sample_sizes$batch_0001$H0)

set.seed(12345)
args$reuse_pvalues <- TRUE
res <- do.call(fw$run_simulation, args)
}
}
\seealso{\code{\link{load_framework}}, \code{\link{run_app}}}
