## sh.test.R (2008-07-11)

##   Shimodaira-Hasegawa Test

## Copyright 2006-2008 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

sh.test <- function(..., x, model = DNAmodel(), B = 100)
{
    ## Prepare the list of trees:
    phy <- list(...)
    if (length(phy) == 1 && class(phy[[1]]) != "phylo")
      phy <- unlist(phy, recursive = FALSE)
    ntree <- length(phy)

    ## Arrange the sequences as a matrix:
    if (is.list(x)) {
        nm <- names(x)
        n <- length(x)
        x <- unlist(x)
        nL <- length(x)
        x <- matrix(x, n, nL/n, byrow = TRUE)
        rownames(x) <- nm
    }

    ## Step 1:
    foo <- function(PHY)
      attr(mlphylo(x, PHY, model, search.tree = FALSE, quiet = TRUE), "loglik")
    Talpha <- sapply(phy, foo)
    Talpha <- max(Talpha) - Talpha

    ## Do the bootstrap resampling (Step 2):
    M <- matrix(NA, ntree, B)
    for (i in 1:B) {
        boot.samp <- x[, sample(ncol(x), replace = TRUE)]
        for (j in 1:ntree)
          M[j, i] <- attr(mlphylo(boot.samp, phy[[j]], model,
                                  search.tree = FALSE, quiet = TRUE),
                          "loglik")
    }
    M <- M - rowMeans(M) # Step 3
    ## Step 4: <FIXME> This can greatly simplified </FIXME>
    for (i in 1:B)
      for (j in 1:ntree)
        M[j, i] <- max(M[j, i] - M[, i])
    ## Step 5:
    count <- numeric(ntree)
    for (j in 1:ntree)
      count[j] <- sum(M[j, ] > Talpha[j])
    count <- count/B
    names(count) <- names(phy)
    count
}
