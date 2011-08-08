## mantel.test.R (2011-06-22)

##   Mantel Test for Similarity of Two Matrices

## Copyright 2002-2006 Ben Bolker and Julien Claude

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

perm.rowscols <- function(m1, n)
{
    s <- sample(1:n)
    m1[s, s]
}

### calculate the Mantel z-statistic for two square matrices m1 and m2
mant.zstat <- function(m1, m2) sum(lower.triang(m1 * m2))

lower.triang <- function(m)
{
    d <- dim(m)
    if (d[1] != d[2]) print("Warning: non-square matrix")
    m[col(m) <= row(m)]
}

mantel.test <- function (m1, m2, nperm = 1000, graph = FALSE,
                         alternative = "two.sided", ...)
{
    alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
    n <- nrow(m1)
    realz <- mant.zstat(m1, m2)
    nullstats <- replicate(nperm, mant.zstat(m1, perm.rowscols(m2, n)))
    pval <- switch(alternative,
                   "two.sided" = 2 * sum(abs(nullstats) > abs(realz)),
                   "less" = sum(nullstats < realz),
                   "greater" = sum(nullstats > realz))
    pval <- pval / nperm
    if (graph) {
        plot(density(nullstats), type = "l", ...)
        abline(v = realz)
    }
    list(z.stat = realz, p = pval, alternative = alternative)
}
