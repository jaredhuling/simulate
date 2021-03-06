% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_covariates.R
\name{gen_covariates}
\alias{gen_covariates}
\title{Generating random covariates}
\usage{
gen_covariates(n, p.contin, p.binary, rho = 0.5, bin.prob = 0.5)
}
\arguments{
\item{n}{integer. sample size}

\item{p.contin}{integer. number of continuous covariates}

\item{p.binary}{integer. number of binary covariates}

\item{rho}{numeric value between 0 and 1, correlation parameter for AR(1) covariance structure of
continuous covariates}

\item{bin.prob}{numeric value between 0 and 1. fraction of 1 values for binary covariates}
}
\description{
Generates a mix of binary and continuous (multivariate normal with AR(1) covariance) covariates
}
\examples{

x <- gen_covariates(10, p.contin = 3, p.binary = 2, rho = 0.9, bin.prob = 0.25)
x
}
