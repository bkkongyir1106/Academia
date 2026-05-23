\name{run_app}
\alias{run_app}
\title{Launch the Pretest Simulation Studio dashboard}
\usage{
run_app(framework_path = NULL, max_upload_mb = 2048,
        launch.browser = interactive(), ...)
}
\arguments{
\item{framework_path}{Optional path to a framework \code{.R} file to use instead
of the bundled copy. Sets the \code{pretestsim.framework_path} option.}
\item{max_upload_mb}{Maximum browser upload size in MB (default 2048). Large
trained-model RData is better loaded by path inside the ML panel.}
\item{launch.browser}{Passed to \code{\link[shiny]{runApp}}.}
\item{...}{Further arguments passed to \code{\link[shiny]{runApp}}.}
}
\value{Invisibly \code{NULL}; called for its side effect of running the app.}
\description{
Starts the Shiny app that drives the bundled, unmodified pretest-simulation
framework through its six phases.
}
\examples{
\dontrun{
library(pretestsim)
run_app()
run_app(framework_path = "~/my_framework.R")
}
}
