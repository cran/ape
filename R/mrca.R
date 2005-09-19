### mrca.R  (2005-09-17)
###
###    Find Most Recent Common Ancestors Between Pairs
###
### Copyright 2005 Emmanuel Paradis
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

mrca <- function(phy, full = FALSE, as.numeric = FALSE)
{
    if (class(phy) != "phylo") stop('object "phy" is not of class "phylo"')
###    if (!is.rooted(phy)) stop("the tree must be rooted.")
    ## Get all clades:
    BP <- .Call("bipartition", as.integer(phy$edge[, 1]),
                as.integer(phy$edge[, 2]), PACKAGE = "ape")
    mode(phy$edge) <- "numeric"
    nb.tip <- max(phy$edge)
    nb.node <- -min(phy$edge)
    N <- nb.tip + nb.node
    ROOT <- nb.tip + 1
    ## Here we change all node numbers {-1, -2, ...} to
    ## {nb.tip + 1, nb.tip + 2, ...}:
    phy$edge[phy$edge < 0] <- nb.tip - phy$edge[phy$edge < 0]
    ## In the following matrix, numeric indexing will be used:
    M <- numeric(N * N)
    dim(M) <- c(N, N)

    ## We start at the root:
    next.node <- ROOT
    while (length(next.node)) {
        tmp <- numeric(0)
        for (anc in next.node) {
            ## Find the branches which `anc' is the ancestor...:
            id <- which(phy$edge[, 1] == anc)
            ## ... and get their descendants:
            desc <- phy$edge[id, 2]
            ## `anc' is itself the MRCA of its direct descendants:
            M[anc, desc] <- M[desc, anc] <- anc
            ## Find all 2-by-2 combinations of `desc': `anc'
            ## is their MRCA:
            for (i in 1:length(desc))
              M[cbind(desc[i], desc[-i])] <- anc
            ## If one element of `desc' is a node, then the tips it
            ## leads to and the other elements of `desc' have also
            ## `anc' as MRCA!
            for (i in 1:length(desc)) {
                if (desc[i] < ROOT) next
                ## (get the tips:)
                tips <- BP[[desc[i] - nb.tip]]
                ## Same thing for the nodes...
                node.desc <- numeric(0)
                for (k in 1:nb.node) {
                    if (k == desc[i] - nb.tip) next
                    ## If the clade of the current node is a
                    ## subset of desc[i], then it is one of its
                    ## descendants:
                    if (all(BP[[k]] %in% tips))
                      node.desc <- c(node.desc, k)
                }
                ## all nodes and tips which are descendants of
                ## `desc[i]':
                ALLDESC <- c(tips, node.desc + nb.tip)
                M[ALLDESC, desc[-i]] <- M[desc[-i], ALLDESC] <- anc
                for (j in 1:length(desc)) {
                    if (j == i || desc[j] < ROOT) next
                    tips2 <- BP[[desc[j] - nb.tip]]
                    node.desc <- numeric(0)
                    for (k in 1:nb.node) {
                        if (k == desc[j] - nb.tip) next
                        if (all(BP[[k]] %in% tips2))
                          node.desc <- c(node.desc, k)
                    }
                    ALLDESC2 <- c(tips2, node.desc + nb.tip)
                    M[ALLDESC, ALLDESC2] <- M[ALLDESC2, ALLDESC] <- anc
                }
                ## `anc' is also the MRCA of itself and its descendants:
                M[ALLDESC, anc] <- M[anc, ALLDESC] <- anc
            }
            ## When it is done, `desc' i stored to become
            ## the new `next.node', if they are nodes:
            tmp <- c(tmp, desc[desc > nb.tip])
        }
        next.node <- tmp
    }
    M[cbind(1:N, 1:N)] <- 1:N
    M[M > nb.tip] <- -(M[M > nb.tip] - nb.tip)
    if (!as.numeric) mode(M) <- "character"
    if (full)
      dimnames(M)[1:2] <- list(as.character(c(1:nb.tip, -(1:nb.node))))
    else {
        M <- M[1:nb.tip, 1:nb.tip]
        dimnames(M)[1:2] <- list(phy$tip.label)
    }
    M
}
