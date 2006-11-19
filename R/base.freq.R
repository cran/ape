### base.freq.R (2003-08-13)
###
###    Base frequencies from DNA Sequences
###
### Copyright 2002-2003 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

base.freq <- function(x)
{
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.matrix(x)) dim(x) <- NULL
    if (is.list(x)) x <- unlist(x)
    x <- x[grep("[acgt]", x)] # get only the known bases
    n <- length(x)
    table(factor(x, levels = c("a", "c", "g", "t"))) / n
}

GC.content <- function(x)
{
    BF <- base.freq(x)
    sum(BF[2:3])
}
