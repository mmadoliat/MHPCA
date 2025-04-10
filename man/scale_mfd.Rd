% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mfdMethods.R
\name{scale_mfd}
\alias{scale_mfd}
\title{Scale an `mfd` Object}
\usage{
scale_mfd(mfd_obj, mfd_eval_length = 100, weight = NULL)
}
\arguments{
\item{mfd_obj}{An object of class `mfd`.}

\item{mfd_eval_length}{The number of evaluation points to use for scaling.}

\item{weight}{An optional numeric value to use as the scaling factor. If NULL, the scaling factor is calculated automatically.}
}
\value{
A scaled mfd object.
}
\description{
This function scales an `mfd` object by calculating a scaling factor based on the variance of its evaluations 
or using a provided weight. It returns a new scaled `mfd` object.
}
\examples{
# Example usage:
# Assuming `mfd_obj` is a valid mfd object:
# scaled_mfd <- scale_mfd(mfd_obj, mfd_eval_length = 100)
# scaled_mfd <- scale_mfd(mfd_obj, mfd_eval_length = 100, weight = 0.5)
}
