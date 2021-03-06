% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_coefficients.R
\name{gen_coefs}
\alias{gen_coefs}
\title{Generating random coefficients}
\usage{
gen_coefs(p, nonzero.frac = 1, family = c("uniform", "normal"),
  max.effect.size = 1, min.effect.size = max.effect.size * 0.5,
  frac.positive = 0.5, mean = 0, sd = 1)
}
\arguments{
\item{p}{integer. dimension of coefficients}

\item{nonzero.frac}{numeric value between 0 and 1, fraction of coefficients which are nonzero}

\item{family}{either \code{"uniform"} or \code{"normal"}. Specifying \code{"uniform"} will result in coefficients
generated uniformly on [\code{-max.effect.size},\code{-min.effect.size}]u[\code{min.effect.size},\code{max.effect.size}]
and specifying \code{"normal"} will result in coefficients generated randomly from a normal distribution with
mean \code{mean} and standard deviation \code{sd}}

\item{max.effect.size}{positive numeric value, maximum magnitude of coefficients if \code{family = "uniform"}}

\item{min.effect.size}{positive numeric value, minimum magnitude of coefficients if \code{family = "uniform"}}

\item{frac.positive}{positive numeric value between 0 and 1. If \code{family = "uniform"}, this is the
fraction of coefficients which have positive values}

\item{mean}{mean of random coefficients if \code{family = "normal"}}

\item{sd}{standard deviation of random coefficients if \code{family = "normal"}}
}
\description{
Generates random coeffients (some of which can be exactly zero)
}
\examples{

beta <- gen_coefs(20, nonzero.frac = 0.5,
                  max.effect.size = 0.5)
beta
}
