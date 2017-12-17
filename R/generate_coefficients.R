
#' Generating random coefficients
#'
#' @description Generates random coeffients (some of which can be exactly zero)
#'
#' @param p integer. dimension of coefficients
#' @param nonzero.frac numeric value between 0 and 1, fraction of coefficients which are nonzero
#' @param family either \code{"uniform"} or \code{"normal"}. Specifying \code{"uniform"} will result in coefficients
#' generated uniformly on [\code{-max.effect.size},\code{-min.effect.size}]u[\code{min.effect.size},\code{max.effect.size}]
#' and specifying \code{"normal"} will result in coefficients generated randomly from a normal distribution with
#' mean \code{mean} and standard deviation \code{sd}
#' @param max.effect.size positive numeric value, maximum magnitude of coefficients if \code{family = "uniform"}
#' @param min.effect.size positive numeric value, minimum magnitude of coefficients if \code{family = "uniform"}
#' @param frac.positive positive numeric value between 0 and 1. If \code{family = "uniform"}, this is the
#' fraction of coefficients which have positive values
#' @param mean mean of random coefficients if \code{family = "normal"}
#' @param sd standard deviation of random coefficients if \code{family = "normal"}
#' @export
#' @examples
#'
#' beta <- gen_coefs(20, nonzero.frac = 0.5,
#'                   max.effect.size = 0.5)
#' beta
gen_coefs <- function(p,
                      nonzero.frac = 1,
                      family = c("uniform", "normal"),
                      max.effect.size = 1,
                      min.effect.size = max.effect.size * 0.5,
                      frac.positive = 0.5,
                      mean = 0,
                      sd = 1)
{

    family <- match.arg(family)

    max.effect.size <- as.double(max.effect.size[1])
    min.effect.size <- as.double(min.effect.size[1])
    p               <- as.integer(p[1])
    nonzero.frac    <- as.double(nonzero.frac[1])
    sd              <- as.double(sd[1])

    stopifnot(nonzero.frac <= 1 & nonzero.frac > 0)
    stopifnot(max.effect.size > 0)
    stopifnot(min.effect.size >= 0)
    stopifnot(sd > 0)

    p.dep <- floor(p * nonzero.frac)


    stopifnot(p.dep <= p)
    stopifnot(min.effect.size < max.effect.size)

    if (family == "uniform")
    {
        mult.fact <- rep(1, p.dep)
        mult.fact[1:floor(frac.positive * p.dep)] <- -1
        mult.fact <- mult.fact[sample.int(p.dep, p.dep)]

        beta.nz <- runif(p.dep, min = min.effect.size, max = max.effect.size) * mult.fact
    } else
    {
        beta.nz <- rnorm(p.dep, mean = mean, sd = sd)
    }

    beta <- numeric(p)
    beta[sample.int(p, p.dep)] <- beta.nz

    beta
}
