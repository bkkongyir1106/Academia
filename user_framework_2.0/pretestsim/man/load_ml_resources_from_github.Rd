\name{load_ml_resources_from_github}
\alias{load_ml_resources_from_github}
\title{Load ML resources from GitHub-backed assets}
\description{
Downloads the ML framework R file and trained-model RData from GitHub, then
delegates validation and loading to the user framework.
}
\usage{
load_ml_resources_from_github(framework_env, refresh = FALSE,
  raw_base = .default_raw_base(), required_sample_sizes = NULL)
}
\arguments{
\item{framework_env}{Environment returned by \code{load_framework()}.}
\item{refresh}{Whether to refresh cached GitHub copies.}
\item{raw_base}{Raw GitHub folder URL.}
\item{required_sample_sizes}{Optional sample sizes to map to trained bundles.}
}
\value{
The ML resources list returned by the user framework.
}
