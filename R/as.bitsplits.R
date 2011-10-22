## as.bitsplits.R (2011-10-19)

##   Conversion Among Split Classes

## Copyright 2011 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

as.bitsplits <- function(x) UseMethod("as.bitsplits")

as.bitsplits.prop.part <- function(x)
{
    foo <- function(vect, RAWVECT) {
        res <- RAWVECT
        for (y in vect) {
            i <- ceiling(y/8)
            res[i] <- res[i] | as.raw(2^(8 - ((y - 1) %% 8) - 1))
        }
        res
    }

    N <- length(x) # number of splits
    n <- length(x[[1]]) # number of tips

    nr <- ceiling(n/8)
    mat <- raw(N * nr)
    dim(mat) <- c(nr, N)

    RAWVECT <- raw(nr)

    for (i in 1:N) mat[, i] <- foo(x[[i]], RAWVECT)

    ## add the n trivial splits of size 1... :
    mat.bis <- raw(n * nr)
    dim(mat.bis) <- c(nr, n)
    for (i in 1:n) mat.bis[, i] <- foo(i, RAWVECT)

    ## ... drop the trivial split of size n... :
    mat <- cbind(mat.bis, mat[, -1, drop = FALSE])

    ## ... update the split frequencies... :
    freq <- attr(x, "number")
    freq <- c(rep(freq[1L], n), freq[-1L])

    ## ... and numbers:
    N <- N + n - 1L

    structure(list(matsplit = mat, labels = attr(x, "labels"),
                   freq =  freq), class = "bitsplits")
}

print.bitsplits <- function(x, ...)
{
    cat('Object of class "bitsplits"\n')
    cat('   ', length(x$labels), 'tips\n')
    cat('   ', length(x$freq), 'partitions\n\n')
}
