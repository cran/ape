### dist.topo.R  (2005-05-31)
###
###     Topological Distances, Tree Bipartition, and
###              Bootstrapping Phylogenies
###
### Copyright 2005 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
###
### This file is part of the `ape' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

dist.topo <- function(x, y)
{
    bp1 <- .Call("bipartition", as.integer(x$edge[, 1]),
                 as.integer(x$edge[, 2]), PACKAGE = "ape")
    ## we remove the 2 first elements of 'bp1' since:
    ##  - the 1st one contains all tips (maybe need to be changed!)
    ##  - the 2nd one is redundant with the other node linked to the root
    bp1 <- lapply(bp1[-(1:2)], function(xx) sort(x$tip.label[xx]))
    bp2 <- .Call("bipartition", as.integer(y$edge[, 1]),
                 as.integer(y$edge[, 2]), PACKAGE = "ape")
    bp2 <- lapply(bp2[-(1:2)], function(xx) sort(y$tip.label[xx]))
    q1 <- length(bp1)
    q2 <- length(bp2)
    p <- 0
    for (i in 1:length(bp1)) {
        for (j in 1:length(bp2)) {
            if (identical(all.equal(bp1[[i]], bp2[[j]]), TRUE)) {
                p <- p + 1
                break
            }
        }
    }
    dT <- if (q1 == q2) 2 * (q1 - p) else 2 * (min(q1, q2) - p) + abs(q1 - q2)
    dT
}

prop.part <- function(...)
{
    obj <- list(...)
    if (length(obj) == 1 & class(obj[[1]]) != "phylo")
      obj <- unlist(obj, recursive = FALSE)
    ntree <- length(obj)
    bp <- .Call("bipartition", as.integer(obj[[1]]$edge[, 1]),
                as.integer(obj[[1]]$edge[, 2]), PACKAGE = "ape")
    clades <- lapply(bp, function(xx) sort(obj[[1]]$tip.label[xx]))
    no <- rep(1, length(clades))

    if (ntree > 1) {
        for (k in 2:ntree) {
            bp <- .Call("bipartition", as.integer(obj[[k]]$edge[, 1]),
                        as.integer(obj[[k]]$edge[, 2]), PACKAGE = "ape")
            bp <- lapply(bp, function(xx) sort(obj[[k]]$tip.label[xx]))
            for (i in 1:length(bp)) {
                done <- FALSE
                for (j in 1:length(clades)) {
                    if (identical(all.equal(bp[[i]], clades[[j]]), TRUE)) {
                        no[j] <- no[j] + 1
                        done <- TRUE
                        break
                    }
                }
                if (!done) {
                    clades <- c(clades, bp[i])
                    no <- c(no, 1)
                }
            }
        }
    }
    attr(clades, "number") <- no
    clades
}

prop.clades <- function(phy, ..., part = NULL)
{
    if (is.null(part)) {
        obj <- list(...)
        if (length(obj) == 1 & class(obj[[1]]) != "phylo")
          obj <- unlist(obj, recursive = FALSE)
        part <- prop.part(obj)
    }
    bp <- .Call("bipartition", as.integer(phy$edge[, 1]),
                as.integer(phy$edge[, 2]), PACKAGE = "ape")
    bp <- lapply(bp, function(xx) sort(phy$tip.label[xx]))
    n <- numeric(length(bp))
    for (i in 1:length(bp)) {
        for (j in 1:length(part)) {
            if (identical(all.equal(bp[[i]], part[[j]]), TRUE)) {
                n[i] <- attr(part, "number")[j]
                done <-  TRUE
                break
            }
        }
    }
    n
}

boot.phylo <- function(phy, x, FUN, B = 100, block = 1)
{
    if (is.list(x)) {
        nm <- names(x)
        n <- length(x)
        x <- unlist(x)
        nL <- length(x)
        x <- matrix(x, n, nL/n, byrow = TRUE)
        rownames(x) <- nm
    }
    boot.tree <- list()
    length(boot.tree) <- B
    for (i in 1:B) {
        if (block > 1) {
            y <- seq(block, ncol(x), block)
            boot.i <- sample(y, replace = TRUE)
            boot.samp <- numeric(ncol(x))
            boot.samp[y] <- boot.i
            for (j in 1:(block - 1))
              boot.samp[y - j] <- boot.i - j
        } else boot.samp <- sample(ncol(x), replace = TRUE)
        boot.tree[[i]] <- FUN(x[, boot.samp])
    }
    prop.clades(phy, boot.tree)
}
