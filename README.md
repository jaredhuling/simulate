
[![Build Status](https://travis-ci.org/jaredhuling/simulate.svg?branch=master)](https://travis-ci.org/jaredhuling/simulate)

Installing the 'simulate' package
=================================

Install the development version using the **devtools** package:

``` r
devtools::install_github("jaredhuling/simulate")
```

or by cloning and building using `R CMD INSTALL`

Quick Usage Overview
====================

Load the package:

``` r
library(simulate)
```

Simulate some data!

``` r
set.seed(1)
n.obs  <- 100
x <- gen_covariates(n.obs, p.contin = 5, p.binary = 5, bin.prob = 0.25)

head(x)
```

    ##            [,1] [,2]        [,3] [,4]      [,5] [,6]       [,7] [,8]
    ## [1,] -0.3462054    0  0.19671927    0  1.090017    0  0.1220082    0
    ## [2,] -0.2855640    0 -0.94962126    0  1.381942    0 -0.2009052    0
    ## [3,] -1.5147703    0  1.22075585    0  1.219996    1  0.6142724    0
    ## [4,] -0.5552653    1 -1.13214016    0 -1.565413    0 -0.9955410    1
    ## [5,] -0.0633063    0  0.01139449    0 -1.538331    0 -0.6034691    1
    ## [6,] -0.1249445    0  2.33351270    0  1.740155    1 -0.5123285    1
    ##            [,9] [,10]
    ## [1,]  0.9467321     0
    ## [2,] -0.9539595     0
    ## [3,]  0.7716022     1
    ## [4,] -0.9764151     0
    ## [5,]  1.7194675     0
    ## [6,] -1.4200947     1

Simulate random coefficients:

``` r
beta <- gen_coefs(10, nonzero.frac = 0.5, max.effect.size = 0.5,
                  min.effect.size = 0.05)

beta
```

    ##  [1] -0.2227978  0.0000000  0.0000000  0.0000000  0.0000000  0.3864677
    ##  [7]  0.3251893  0.3631059 -0.0942763  0.0000000
