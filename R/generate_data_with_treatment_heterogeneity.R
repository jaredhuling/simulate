

#' Generating data with heterogeneity of treatment effect
#'
#' @description Generates random data with heterogeneity of treatment effect
#'
#' @param n integer sample size
#' @param n.test integer sample size for test dataset
#' @param p.contin number of continuous covariates
#' @param p.binary number of binary covariates
#' @param family type of outcome to generate. \code{"gaussian"} available
#' @param p.propens.dep integer number of variables that have nonzero effects in the the propensity score model
#' @param p.trt.effects integer number of variables that drive treatment heterogeneity
#' @param p.main.effects integer number of variables that have nonzero main effects for the outcome (ie number
#' of prognostic variables)
#' @param max.effect.size.main positive value for maximum magnitude of coefficients for main effects
#' @param max.effect.size.trt positive value for maximum magnitude of coefficients for treatment-covariate
#' interaction effects
#' @param max.effect.size.propens positive value for maximum magnitude of coefficients for propensity
#' score model effects
#' @param quadratic.main.effects logical value. Should quadratic terms be added to the covariate main
#' effects? These terms will not impact treatment heterogeneity
#' @param trt.int.effect.ints either \code{NULL} for no added interactions or a list of vectors of indices of variables. Each vector
#' as an element in the list will represent covariate-covariate interactions to be added to the treatment-covariate
#' interactions. A pairwise interaction of X_1 and X_2 here is in effect a 3-way interaction between Trt, X_1, and X_2
#' @param trt.int.quadratic.effects either \code{NULL} for no added quadratic terms of a vector of indices
#' for variables which have quadratic terms added to the covariate-treatment interactions, e.g.
#' \code{trt.int.quadratic.effects = c(1, 4)} will add \code{Trt:X1^2 + Trt:X4^2} to the model
#' @param propens.rct is the propensity model from a randomized controlled trial? If \code{TRUE}, treatments will
#' be assigned randomly taking value 1 with probability \code{trt.frac}
#' @param snr positive value for the overall signal-to-noise ratio of the model
#' @param trt.frac value between 0 and 1. If \code{propens.rct = TRUE}, \code{trt.frac} is the fraction of observations with trt = 1
#' @param rho value between 0 and 1 for correlation term in AR(1) structure of continuous covariates
#' @param bin.prob value between 0 and 1 for the fraction of 1's in binary covariates
#' @importFrom stats lm
#' @export
#' @examples
#'
#' set.seed(1)
#' quad.eff <- c(8, 9, 10)
#' int.eff  <- list(c(1,2), c(3,4), c(2,3), c(4, 6), c(5,7))
#'
#' data <- gen_hte_data(n = 500, n.test = 10000,
#'                      p.contin = 10,
#'                      p.binary = 5,
#'                      p.propens.dep = 5,
#'                      p.trt.effects = 4,
#'                      p.main.effects = 10,
#'                      max.effect.size.main = 1,
#'                      max.effect.size.trt  = 1,
#'                      max.effect.size.propens = 0.5,
#'                      quadratic.main.effects = TRUE,
#'                      propens.rct = FALSE,
#'                      trt.int.effect.ints = int.eff,
#'                      trt.int.quadratic.effects = quad.eff,
#'                      snr = 1,
#'                      trt.frac = 0.5,
#'                      rho = 0.5,
#'                      bin.prob = 0.25)
#'
#' df <- data.frame(y = data$y,
#'                  trt = data$trt,
#'                  data$x)
#' form <- paste("y ~ trt*(", paste(colnames(df)[-c(1:2)], collapse = "+"), ")")
#' summary(lmf <- lm(as.formula(form), data = df))
#'
#' df0 <- df1 <- data.frame(trt = 1, data$x.test)
#' df0$trt <- -1
#'
#' pred1 <- predict(lmf, df1)
#' pred0 <- predict(lmf, df0)
#'
#' est.opt.trt <- apply(cbind(pred0, pred1), 1, function(rr) c(-1,1)[which.max(rr)])
#'
#' table(est.opt.trt, data$trt.test.optimal)
#'
#' mean(est.opt.trt == data$trt.test.optimal)
#'
gen_hte_data <- function(n,                               # training sample size
                         n.test,                          # testing sample size
                         p.contin,                        # total continuous covariates
                         p.binary,                        # total binary covariates
                         family = c("gaussian"),
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

    family <- match.arg(family)

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
                      gen_coefs(p = p.tot, nonzero.frac = p.propens.dep/p.tot,
                                max.effect.size = max.effect.size.propens))

    # gen treatment assignment
    trt      <- 2 * gen_trt(~ ., data.frame(x),
                            beta.propens, indep = propens.rct, trt.frac, nrow(x)) - 1
    trt.test <- 2 * gen_trt(~ ., data.frame(x.test),
                            beta.propens, indep = propens.rct, trt.frac, nrow(x.test)) - 1

    # gen beta for main effects
    beta.main <- gen_coefs(p = p.tot, nonzero.frac = p.main.effects/p.tot,
                           max.effect.size = max.effect.size.main)

    # gen beta for trt int effects
    beta.trt  <- gen_coefs(p = p.tot, nonzero.frac = p.trt.effects/p.tot,
                           max.effect.size = max.effect.size.trt)


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
        beta.trt.int <- gen_coefs(p = p.intint, nonzero.frac = 1,
                                  max.effect.size = max.effect.size.trt)
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
        beta.trt.quad <- gen_coefs(p = p.intquad, nonzero.frac = 1,
                                   max.effect.size = max.effect.size.trt * 0.5)

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
