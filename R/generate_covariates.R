
#' Generating random covariates
#'
#' @description Generates a mix of binary and continuous (multivariate normal with AR(1) covariance) covariates
#'
#' @param n integer. sample size
#' @param p.contin integer. number of continuous covariates
#' @param p.binary integer. number of binary covariates
#' @param rho numeric value between 0 and 1, correlation parameter for AR(1) covariance structure of
#' continuous covariates
#' @param bin.prob numeric value between 0 and 1. fraction of 1 values for binary covariates
#' @importFrom stats rbinom rnorm runif var
#' @importFrom MASS mvrnorm
#' @export
#' @examples
#'
#' x <- gen_covariates(10, p.contin = 3, p.binary = 2, rho = 0.9, bin.prob = 0.25)
#' x
gen_covariates <- function(n,
                           p.contin,
                           p.binary,
                           rho = 0.5,
                           bin.prob = 0.5)
{

    p.contin <- as.integer(p.contin[1])
    p.binary <- as.integer(p.binary[1])
    bin.prob <- as.double(bin.prob[1])
    rho      <- as.double(rho[1])
    n        <- as.integer(n[1])

    p.tot <- p.contin + p.binary

    stopifnot(p.tot > 0 & n > 0)
    stopifnot(p.contin >= 0)
    stopifnot(p.binary >= 0)

    stopifnot(rho <= 1 & rho >= 0)
    stopifnot(bin.prob <= 1 & bin.prob >= 0)


    # generate continuous covariates

    if (p.contin > 0)
    {
        Sig <- rho ^ abs(outer(1:p.contin, 1:p.contin, FUN = "-"))
        x.contin <- mvrnorm(n, mu = rep(0, p.contin), Sigma = Sig)
    }

    if (p.binary > 0)
    {
        x.binary <- matrix(rbinom(n * p.binary, 1, bin.prob), ncol = p.binary)
    }

    if (p.binary > 0 & p.contin > 0)
    { # only add them together if they are both asked for
        minp     <- min(p.contin, p.binary)
        whichmin <- which.min(c(p.contin, p.binary))

        if (whichmin == 1)
        {
            x <- cbind(x.contin, x.binary)
        } else
        {
            x <- cbind(x.binary, x.contin)
        }


        colgrid <- 1:ncol(x)

        intseq <- rep(1:minp, each = 2)
        intseq[2 * (1:minp)] <- intseq[2 * (1:minp)] + minp

        colgrid[1:(2*minp)] <- intseq

        # here we're mixing up the columns such that the first min(p.contin, p.binary)
        # columns alternate between binary and continuous. This is so we don't make
        # either all continuous or all covariate effects as the ones that only drive heterogeneity
        # of treatment effect
        x <- x[,colgrid]


    } else
    {
        if (p.binary > 0)
        {
            x <- x.binary
        } else
        {
            c <- x.contin
        }
    }

    x
}





