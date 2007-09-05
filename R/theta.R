## theta.R (2002-08-28)

##   Population Parameter THETA

## theta.h: using homozigosity
## theta.k: using expected number of alleles
## theta.s: using segregating sites in DNA sequences

## Copyright 2002 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

theta.h <- function(x, standard.error = FALSE)
{
    HE <- H(x, variance = TRUE)
    sdH <- HE[2]
    HE <- HE[1]
    f <- function(th) HE - th * (1 + (2 * (1 + th)) / ((2 + th) * (3 + th)))
    th <- uniroot(f, interval = c(0, 1))$root
    if (standard.error) {
        SE <- (2 + th)^2 * (2 + th)^3 * sdH /
          HE^2 * (1 + th) * ((2 + th) * (3 + th) * (4 + th) + 10 * (2 + th) + 4)
        return(c(th, SE))
    }
    else return(th)
}

theta.k <- function(x, n = NULL, k = NULL)
{
    if (is.null(n)) {
        if (!is.factor(x)) {
            if (is.numeric(x)) {
                n <- sum(x)
                k <- length(x)
            }
            else x <- factor(x)
        }
        if (is.factor(x)) { # ne pas remplacer par `else'...
            n <- length(x)
            k <- nlevels(x)
        }
    }
    f <- function(th) th * sum(1 / (th + (0:(n - 1)))) - k
    th <- uniroot(f, interval = c(1e-8, 100))$root
    return(th)
}

theta.s <- function(s, n, variance = FALSE)
{
    a1 <- sum(1 / (1:(n - 1)))
    th <- s / a1
    if (variance) {
        a2 <- sum(1 / (1:(n - 1))^2)
        var.th <- (a1^2 * s + a2 * s^2) / (a1^2 * (a1^2 + a2))
        return(c(th, var.th))
    }
    else return(th)
}
