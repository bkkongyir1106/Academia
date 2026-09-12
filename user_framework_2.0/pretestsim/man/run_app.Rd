\name{run_app}
\alias{run_app}
\title{Launch the Pretest Simulation Studio dashboard}
\description{
Starts the Shiny dashboard that downloads framework assets from GitHub, runs
selected simulation phases, and displays plots, tables, logs, and downloads.
}
\usage{
run_app(raw_base = .default_raw_base(), refresh_assets = FALSE, max_upload_mb = 2048,
  launch.browser = interactive(), ...)
}
\arguments{
\item{raw_base}{Raw GitHub folder URL containing package assets.}
\item{refresh_assets}{Whether to refresh cached assets on launch.}
\item{max_upload_mb}{Maximum Shiny upload size in MB.}
\item{launch.browser}{Passed to \code{shiny::runApp()}.}
\item{...}{Additional arguments passed to \code{shiny::runApp()}.}
}
\value{
Invisibly \code{NULL}; called for its side effect of running the app.
}
