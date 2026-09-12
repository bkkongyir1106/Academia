\name{load_framework}
\alias{load_framework}
\title{Load the adaptive pretesting framework}
\description{
Downloads the framework R file from GitHub when needed, sources it into an
isolated environment, and returns that environment.
}
\usage{
load_framework(path = NULL, refresh = FALSE, raw_base = .default_raw_base())
}
\arguments{
\item{path}{Optional local framework R file.}
\item{refresh}{Whether to refresh the cached GitHub copy.}
\item{raw_base}{Raw GitHub folder URL.}
}
\value{
An environment containing framework functions such as \code{run_simulation()}.
}
