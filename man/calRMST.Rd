% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calRMST.R
\name{calRMST}
\alias{calRMST}
\title{Calculate RMST}
\usage{
calRMST(m = 1, shape = 1, k = 1, tau = NULL)
}
\arguments{
\item{m}{Median survival time.}

\item{shape}{The shape parameter in generalized gamma distributions
(flexsurv::pgengamma.orig) for event time.}

\item{k}{The k parameter in generalized gamma distributions
(flexsurv::pgengamma.orig) for event time.}

\item{tau}{A value to specify the truncation time point for the RMST calculation.}
}
\description{
A function to calculate RMST
}
\examples{
# calRMST(m=1, shape=1, k=1, tau=3)
}
