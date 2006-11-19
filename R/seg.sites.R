### seg.sites.R (2004-03-19)
###
###    Find Segregating Sites in DNA Sequences
###
### Copyright 2004 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

seg.sites <- function(X)
{
    if (is.list(X)) {
        if (length(unique(unlist(lapply(X, length)))) > 1)
          stop("sequences in list must have the same lengths")
        X <- matrix(unlist(X), nrow = length(X), byrow = TRUE)
    }
    if (is.data.frame(X)) X <- as.matrix(X)
    which(apply(X, 2, function(x) length(unique(x)) > 1))
}
