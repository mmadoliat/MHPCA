% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nfdMethods.R
\name{$.nfd}
\alias{$.nfd}
\title{Accessor for nfd Class}
\usage{
\method{$}{nfd}(x, name)
}
\arguments{
\item{x}{An object of class `nfd`.}

\item{name}{The name of the element to access.}
}
\value{
The requested element or `NULL` if it doesn't exist.
}
\description{
Accessor for nfd Class
}
\examples{
mat <- matrix(1:9, nrow = 3)
nfd_obj <- nfd(mat)
nfd_obj$nobs # Returns 3
}
