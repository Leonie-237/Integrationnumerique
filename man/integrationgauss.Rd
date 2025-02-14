% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integration_numerique.R
\name{integrationgauss}
\alias{integrationgauss}
\title{Title Numerical Integration using the Gauss-Legendre Method}
\usage{
integrationgauss(f, a, b, n)
}
\arguments{
\item{f}{function to integrate}

\item{a}{lower bound of integration}

\item{b}{upper bound of integration}

\item{n}{number of subdivisions}
}
\value{
The approximate value of the integral
}
\description{
Title Numerical Integration using the Gauss-Legendre Method
}
\examples{
f <- function(x) exp(-x^2)
resultat <- integrationgauss(f, -1, 1, 5)
}
