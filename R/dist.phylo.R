### dist.phylo.R  (2005-09-14)
###
###     Pairwise Distances from a Phylogenetic Tree
###
### Copyright 2002-2005 Emmanuel Paradis
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

dist.phylo <- function(phy, full = FALSE)
{
    if (class(phy) != "phylo") stop('object "phy" is not of class "phylo"')
    if (is.null(phy$edge.length)) stop("your tree has no branch lengths defined")
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    N <- nb.tip + nb.node

    x <- matrix(NA, N, N)
    mode(x) <- "numeric"
    x[cbind(1:N, 1:N)] <- 0
    dimnames(x)[1:2] <- list(as.character(c(1:nb.tip, -(1:nb.node))))

    for (i in 1:dim(phy$edge)[1])
      x[phy$edge[i, 2], phy$edge[i, 1]] <-
        x[phy$edge[i, 1], phy$edge[i, 2]] <- phy$edge.length[i]

    ok <- c(rep(TRUE, nb.tip), rep(FALSE, nb.node))
    names(ok) <- dimnames(x)[[1]]

    basal <- ""
    while(!identical(basal, "-1")) {
        term <- names(ok[ok])
        ind <- phy$edge[, 2] %in% term
        basal <- names(which(table(phy$edge[ind, 1]) > 1))
        for (nod in basal) {
            i.anc <- which(phy$edge[, 2] == nod)
            l <- phy$edge.length[i.anc]
            anc <- phy$edge[i.anc, 1]

            desc <- phy$edge[which(phy$edge[, 1] == nod), 2]
            ## Here we need to check that all the branches found in the next
            ## few lines just above are `available' for `clustering'; this may
            ## not be the case if other sister-branches have daughter-branches
            ## which are not yet defined, for instance if there is a multichotomy.
            if (all(desc %in% term)) {
                ## compute the distances ...
                for (i in 1:(length(desc) - 1)) {
                    if (as.numeric(desc[i]) > 0) d1 <- desc[i] else {
                        ## lin <- x[desc[i], 1:nb.tip]
                        lin <- x[desc[i], ]
                        d1 <- names(lin[!is.na(lin)])
                    }
                    for (j in (i + 1):length(desc)) {
                        if (as.numeric(desc[j]) > 0) d2 <- desc[j] else {
                            ## lin <- x[desc[j], 1:nb.tip]
                            lin <- x[desc[j], ]
                            d2 <- names(lin[!is.na(lin)])
                        }
                        for (y in d1)
                          for (z in d2)
                            x[y, z] <- x[z, y] <- x[nod, y] + x[nod, z]
                        ## compute the distances between the tips in `d2'
                        ## and the ancestor of the current node
                        x[d2, anc] <- x[anc, d2] <- x[nod, d2] + l
                    }
                    ## compute the distances between the tips in `d1'
                    ## and the ancestor of the current node
                    x[d1, anc] <- x[anc, d1] <- x[nod, d1] + l
                }
                ok[desc] <- FALSE
                ok[nod] <- TRUE
            }
        }
    }
    if (!full) {
        x <- x[1:nb.tip, 1:nb.tip]
        dimnames(x)[1:2] <- list(phy$tip.label)
    }
    x
}
