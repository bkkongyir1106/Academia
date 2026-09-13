\name{assets}
\alias{available_assets}
\alias{asset_url}
\alias{cache_asset}
\title{Resolve GitHub-backed package assets}
\description{
Functions for listing, locating, downloading, and caching framework assets from
the package GitHub repository.
}
\usage{
available_assets()
asset_url(asset, raw_base = .default_raw_base())
cache_asset(asset, refresh = FALSE, raw_base = .default_raw_base(),
  cache_dir = .cache_dir())
}
\arguments{
\item{asset}{Asset key or file name.}
\item{raw_base}{Raw GitHub folder URL.}
\item{refresh}{Whether to download a fresh copy even when cached.}
\item{cache_dir}{Local cache directory.}
}
\value{
\code{available_assets()} returns a named character vector, \code{asset_url()}
returns a URL, and \code{cache_asset()} returns a local file path.
}
