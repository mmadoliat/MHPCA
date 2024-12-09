% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vdClass.R
\name{vd}
\alias{vd}
\alias{Vd}
\title{Define a Set of Multidimensional Vector Data Objects}
\usage{
Vd(data)
}
\arguments{
\item{data}{A numeric matrix of data points.}
}
\description{
The `vd` class represents a set of multidimensional vector data.
Vector data objects are constructed by specifying a matrix of data points.

Constructor for `vd` objects (same as Vd(...))
}
\examples{
# Create some vector data
data_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
vd_obj <- Vd(data = data_matrix)
print(vd_obj)
vd_obj$data
vd_obj$nobs

}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{data}}{A numeric matrix of data points.}

\item{\code{nobs}}{Number of observations.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-vd-new}{\code{vd$new()}}
\item \href{#method-vd-print}{\code{vd$print()}}
\item \href{#method-vd-clone}{\code{vd$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-vd-new"></a>}}
\if{latex}{\out{\hypertarget{method-vd-new}{}}}
\subsection{Method \code{new()}}{
Constructor for `vd` objects (same as Vd(...))
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{vd$new(data)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{A numeric matrix of data points.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-vd-print"></a>}}
\if{latex}{\out{\hypertarget{method-vd-print}{}}}
\subsection{Method \code{print()}}{
Print method for `vd` objects.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{vd$print(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Additional arguments to be passed to `print`.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-vd-clone"></a>}}
\if{latex}{\out{\hypertarget{method-vd-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{vd$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
