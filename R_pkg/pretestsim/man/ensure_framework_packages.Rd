\name{ensure_framework_packages}
\alias{ensure_framework_packages}
\title{Ensure the framework's runtime packages are available}
\usage{
ensure_framework_packages()
}
\value{(invisibly) the character vector of packages that are available.}
\description{
Attempts to install (quietly) the packages the bundled framework reaches for
via \code{::}. Anything missing simply disables the feature that needs it.
Never errors.
}
