### balance.R  (2002-08-28)
###
###     Balance of a Dichotomous Phylogenetic Tree
###
### Copyright 2002 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

balance <- function(phy)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    if (nb.node != nb.tip - 1) stop("\"phy\" is not fully dichotomous")
    balance <- matrix(NA, nb.tip + nb.node, 2)
    balance[1:nb.tip, ] <- 0.5

    rownames(balance) <- as.character(c(1:nb.tip, -(1:nb.node)))
    ## `unused' says if the tip or node has NOT been used
    ## to compute balance at a lower level...
    unused <- rep(TRUE, nb.tip + nb.node)
    names(unused) <- rownames(balance)

    while(sum(unused) > 1) {
        term <- names(balance[, 1][!is.na(balance[, 1]) & unused])
        ind <- as.logical(match(phy$edge[, 2], term))
        ind[is.na(ind)] <- FALSE
        term.br <- matrix(phy$edge[ind], length(term), 2)
        ## extract the nodes with 2 branches above
        basal <- names(which(table(term.br[, 1]) == 2))
        for (nod in basal) {
            pair.ind <- which(phy$edge[, 1] == nod)
            pair <- phy$edge[pair.ind, 2]
            balance[nod, ] <- c(sum(balance[pair[1], ]), sum(balance[pair[2], ]))
            unused[pair] <- FALSE
        }
    }
    balance <- balance[-(1:nb.tip), ]
    if (!is.null(phy$node.label)) rownames(balance) <- phy$node.label
    balance
}
