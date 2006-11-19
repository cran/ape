### gammaStat.R (2006-10-04)
###
###    Gamma-Statistic of Pybus and Harvey
###
### Copyright 2002-2006 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

gammaStat <- function(phy)
{
    if (class(phy) != "phylo") stop('object "phy" is not of class "phylo"')
    N <- length(phy$tip.label)
    bt <- sort(branching.times(phy))
    g <- rev(c(bt[1], diff(bt))) # internode intervals are from past to present
    ST <- sum((2:N) * g)
    stat <- sum(cumsum((2:(N - 1)) * g[-(N - 1)]))/(N - 2)
    m <- ST/2
    s <- ST * sqrt(1/(12 * (N - 2)))
    (stat - m)/s
}
