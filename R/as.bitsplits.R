## as.bitsplits.R (2021-12-27)

##   Conversion Among Split Classes

## Copyright 2011-2021 Emmanuel Paradis, 2019 Klaus Schliep

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
    n <- length(x$freq)
    cat("Object of class \"bitsplits\"\n")
    cat("   ", length(x$labels), "tips\n")
    cat("   ", n, "partition")
    if (n > 1) cat("s")
    cat("\n")
}

sort.bitsplits <- function(x, decreasing = FALSE, ...)
{
    o <- order(x$freq, decreasing = decreasing)
    x$matsplit <- x$matsplit[, o]
    x$freq <- x$freq[o]
    x
}

as.prop.part <- function(x, ...) UseMethod("as.prop.part")

as.prop.part.bitsplits <- function(x, include.trivial = FALSE, ...)
{
    decodeBitsplits <- function(x) {
        f <- function(y) rev(rawToBits(y)) == as.raw(1)
        which(unlist(lapply(x, f)))
    }
    N <- ncol(x$matsplit) # nb of splits
    n <- length(x$labels) # nb of tips
    Nres <- if (include.trivial) N + 1L else N
    res <- vector("list", Nres)
    if (include.trivial) res[[1]] <- 1:n
    j <- if (include.trivial) 2L else 1L
    for (i in 1:N) {
        res[[j]] <- decodeBitsplits(x$matsplit[, i])
        j <- j + 1L
    }
    attr(res, "number") <- if (include.trivial) c(N, x$freq) else x$freq
    attr(res, "labels") <- x$labels
    class(res) <- "prop.part"
    res
}

bitsplits <- function(x)
{
    if (inherits(x, "phylo")) {
        x <- list(x)
        class(x) <- "multiPhylo"
    } else {
        if (!inherits(x, "multiPhylo"))
            stop('x is not of class "phylo" or "multiPhylo"')
    }
    if (any(is.rooted(x)))
        stop("bitsplits() accepts only unrooted trees")
    x <- .compressTipLabel(x)
    labs <- attr(x, "TipLabel")
    n <- length(labs)
    nr <- ceiling(n/8)
    ans <- .Call(bitsplits_multiPhylo, x, n, nr)
    nc <- ans[[3]]
    if (nc) {
        o <- ans[[1]][1:(nr * nc)]
        freq <- ans[[2]][1:nc]
    } else {
        o <- raw()
        freq <- integer()
    }
    dim(o) <- c(nr, nc)
    structure(list(matsplit = o, labels = labs, freq = freq),
              class = "bitsplits")
}

countBipartitions <- function(phy, X)
{
    split <- bitsplits(phy)
    SPLIT <- bitsplits(X)
    .Call("CountBipartitionsFromSplits", split, SPLIT)
}
