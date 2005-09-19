### branching.times.R  (2004-04-23)
###
###     Branching Times of a Phylogenetic Tree
###
### Copyright 2002-2004 Emmanuel Paradis
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

branching.times <- function(phy)
{
    if (class(phy) != "phylo")
      stop("object \"phy\" is not of class \"phylo\"")
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    xx <- as.numeric(rep(NA, 1 + nb.node))
    names(xx) <- as.character(c(-(1:nb.node), 1))
    xx["-1"] <- 0
    for (i in 2:length(xx)) {
        nod <- names(xx[i])
        ind <- which(phy$edge[, 2] == nod)
        base <- phy$edge[ind, 1]
        xx[i] <- xx[base] + phy$edge.length[ind]
    }
    depth <- xx[length(xx)]
    branching.times <- depth - xx[-length(xx)]
    if (!is.null(phy$node.label))
      names(branching.times) <- phy$node.label
    branching.times
}
