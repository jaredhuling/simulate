

gen_covariates <- function(n,
                           p.contin,
                           p.binary,
                           rho = 0.5,
                           bin.prob = 0.5)
{

    p.contin <- as.integer(p.contin[1])
    p.binary <- as.integer(p.binary[1])
    n <- as.integer(n[1])

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

gen.effects <- function(p, p.dep, max.effect.size = 1)
{
    stopifnot(p.dep <= p)
    beta.nz <- runif(p.dep, min = max.effect.size * 0.5, max = max.effect.size) *
        (2 * rbinom(p.dep, 1, 0.5) - 1)

    beta <- numeric(p)
    beta[sample.int(p, p.dep)] <- beta.nz

    beta
}

gen.trt.assignment <- function(x, indep = FALSE, beta, trt.frac = 0.5)
{
    if (indep)
    {
        trt <- rbinom(NROW(x), 1, trt.frac)
    } else
    {
        p <- NCOL(x)
        stopifnot(length(beta) == p)

        prob <- 1 / (1 + exp(-drop(x %*% beta)))

        trt <- rbinom(NROW(x), 1, prob)
    }

    trt
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
    x        <- gen.cov(n,      p.contin, p.binary, rho, bin.prob)
    x.test   <- gen.cov(n.test, p.contin, p.binary, rho, bin.prob)

    # gen beta for propensity function
    beta.propens <- gen.effects(p.tot, p.propens.dep, max.effect.size.propens)

    # gen treatment assignment
    trt      <- 2 * gen.trt.assignment(x,      indep = propens.rct, beta.propens, trt.frac) - 1
    trt.test <- 2 * gen.trt.assignment(x.test, indep = propens.rct, beta.propens, trt.frac) - 1

    # gen beta for main effects
    beta.main <- gen.effects(p.tot, p.main.effects, max.effect.size.main)

    # gen beta for trt int effects
    beta.trt  <- gen.effects(p.tot, p.trt.effects, max.effect.size.trt)


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
        beta.trt.int <- gen.effects(p.intint, p.intint, max.effect.size.trt)
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
        p.intquad <- length(trt.int.quadratic.effects)
        beta.trt.quad <- gen.effects(p.intquad, p.intquad, max.effect.size.trt * 0.5)

        delta.train <- delta.train + drop((x[,trt.int.quadratic.effects] ^ 2) %*%  beta.trt.quad)
        delta.test  <- delta.test + drop((x.test[,trt.int.quadratic.effects] ^ 2) %*%  beta.trt.quad)
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
