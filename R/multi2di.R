### multi2di.R  (2005-04-15)
###
###     Collapse and Resolve Multichotomies
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

multi2di <- function(phy)
{
    X <- table(phy$edge[, 1])
    target <- names(which(X > 2))
    if (is.null(target)) return(phy)
    next.node <- min(as.numeric(phy$edge)) - 1
    for (node in target) {
        ind <- which(phy$edge[, 1] == node)
        new.nodes <- as.character(next.node:(next.node - length(ind) + 3))
        N <- length(new.nodes)
        phy$edge[ind[-1], 1] <- c(new.nodes, new.nodes[N])
        mat <- matrix(c(node, new.nodes[-N], new.nodes), ncol = 2)
        phy$edge <- rbind(phy$edge, mat)
        phy$edge.length <- c(phy$edge.length, rep(0, N))
        next.node <- next.node - N
    }
    read.tree(text = write.tree(phy, multi.line = FALSE))
}

di2multi <- function(phy, tol = 1e-8)
{
    ## We select only the internal branches which are
    ## significantly small:
    ind <- which(phy$edge.length < tol & as.numeric(phy$edge[, 2]) <  0)
    if (!length(ind)) return(phy)
    for (i in ind) {
        des <- phy$edge[i, 2]
        phy$edge[which(phy$edge[, 1] == des), 1] <- phy$edge[i, 1]
    }
    phy$edge <- phy$edge[-ind, ]
    phy$edge.length <- phy$edge.length[-ind]
    read.tree(text = write.tree(phy, multi.line = FALSE))
}
