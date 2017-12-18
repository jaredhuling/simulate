
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





gen.data <- function(n,                               # training sample size
                     n.test,                          # testing sample size
                     p.contin,                        # total continuous covariates
                     p.binary,                        # total binary covariates
                     p.propens.dep,                   # number effects guide propensity function
                     p.trt.effects,                   # num effects guide trt int effects
                     p.main.effects,                  # num effects guide main effects
                     max.effect.size.main    = 1,     # max magnitude of main effects
                     max.effect.size.trt     = 1,     # max magnitude of trt int effects
                     max.effect.size.propens = 1,     # max magnitude of propensity func effects
                     quadratic.main.effects  = FALSE, # is there a quadratic main effect?
                     trt.int.effect.ints     = NULL,  # list of all pairwise covariate interactions
                     # among the treatment-covariate interactions
                     trt.int.quadratic.effects = NULL, # vector of quadratic covariate effects to add to covariate
                     #-trt interactions
                     propens.rct             = FALSE, # is trt assignment randomized?
                     snr                     = 1,     # total signal-to-noise ratio
                     trt.frac                = 0.5,   # if trt randomized, what's the fraction of trt == 1
                     rho                     = 0.5,   # rho param for AR(1) structure of contin covariate dependence
                     bin.prob                = 0.5)   # for binary covariates, frac samples == 1
{
    # make sure there aren't more active effects than total covariates
    p.tot <- p.contin + p.binary
    stopifnot(p.propens.dep  <= p.tot)
    stopifnot(p.trt.effects  <= p.tot)
    stopifnot(p.main.effects <= p.tot)

    # gen covariates
    x        <- gen_covariates(n,      p.contin, p.binary, rho, bin.prob)
    x.test   <- gen_covariates(n.test, p.contin, p.binary, rho, bin.prob)

    # gen beta for propensity function

    trt.int.coef <- log(trt.frac / (1 - trt.frac))

    beta.propens <- c(trt.int.coef,
                      gen_coefs(p.tot, p.propens.dep/p.tot, max.effect.size.propens))

    # gen treatment assignment
    trt      <- 2 * gen_trt(~ ., data.frame(x),
                            beta.propens, indep = propens.rct, trt.frac, nrow(x)) - 1
    trt.test <- 2 * gen_trt(~ ., data.frame(x.test),
                            beta.propens, indep = propens.rct, trt.frac, nrow(x.test)) - 1

    # gen beta for main effects
    beta.main <- gen_coefs(p.tot, p.main.effects/p.tot, max.effect.size.main)

    # gen beta for trt int effects
    beta.trt  <- gen_coefs(p.tot, p.trt.effects/p.tot, max.effect.size.trt)


    if (!is.null(trt.int.effect.ints))
    {
        stopifnot(is.list(trt.int.effect.ints))

        ids <- unlist(trt.int.effect.ints)

        stopifnot(max(ids) <= p.tot)
    }

    if (!is.null(trt.int.quadratic.effects))
    {
        stopifnot(max(trt.int.quadratic.effects) <= p.tot)
    }


    # generate responses
    if (quadratic.main.effects)
    {
        y.raw      <- drop(x      %*% beta.main) +
            0.5 * drop(x      %*% beta.main) ^ 2 +
            trt * drop(x %*% beta.trt)
        y.test.raw <- drop(x.test %*% beta.main) +
            0.5 * drop(x.test %*% beta.main) ^ 2 +
            trt.test * drop(x.test %*% beta.trt)
    } else
    {
        y.raw      <- drop(x      %*% beta.main) + trt * drop(x %*% beta.trt)
        y.test.raw <- drop(x.test %*% beta.main) + trt.test * drop(x.test %*% beta.trt)
    }

    delta.train <- drop(x %*% beta.trt)
    delta.test  <- drop(x.test %*% beta.trt)

    if (!is.null(trt.int.effect.ints))
    {
        p.intint <- length(trt.int.effect.ints)
        beta.trt.int <- gen_coefs(p.intint, 1, max.effect.size.trt)
        for (i in 1:p.intint)
        {
            delta.train <- delta.train +
                drop(beta.trt.int[i] * apply(x[,trt.int.effect.ints[[i]]], 1, prod))

            delta.test <- delta.test +
                drop(beta.trt.int[i] * apply(x.test[,trt.int.effect.ints[[i]]], 1, prod))
        }
    }

    if (!is.null(trt.int.quadratic.effects))
    {
        p.intquad     <- length(trt.int.quadratic.effects)
        beta.trt.quad <- gen_coefs(p.intquad, 1, max.effect.size.trt * 0.5)

        delta.train   <- delta.train + drop((x[,trt.int.quadratic.effects] ^ 2)      %*% beta.trt.quad)
        delta.test    <- delta.test  + drop((x.test[,trt.int.quadratic.effects] ^ 2) %*% beta.trt.quad)
    }

    trt.optimal      <- 2 * (drop(x %*% beta.trt) > 0) - 1
    trt.test.optimal <- 2 * (drop(x.test %*% beta.trt) > 0) - 1

    signal <- var(c(y.raw, y.test.raw))
    sd.use <- signal / snr

    # add error term to responses
    y      <- y.raw      + rnorm(n, sd = sd.use)
    y.test <- y.test.raw + rnorm(n.test, sd = sd.use)

    list(y = y, x = x, trt = trt,
         y.test = y.test, x.test = x.test, trt.test = trt.test,
         trt.optimal = trt.optimal,
         trt.test.optimal = trt.test.optimal,
         delta.train = delta.train,
         delta.test = delta.test,
         beta.trt = beta.trt,
         beta.main = beta.main,
         sd.error = sd.use)
}
