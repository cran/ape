### nj.R (2005-12-06)
###
###        Neighbor-Joining Tree Estimation
###
### Copyright 2004-2005 Emmanuel Paradis
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

nj <- function(X)
{
    if (is.matrix(X)) X <- as.dist(X)
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    edge1 <- edge2 <- numeric(2 * N - 3)
    edge.length <- numeric(2 * N - 3)
    ans <- .C("nj", as.double(X), as.integer(N), as.integer(edge1),
              as.integer(edge2), as.double(edge.length), PACKAGE = "ape")
    edge <- cbind(ans[[3]], ans[[4]])
    mode(edge) <- "character"
    obj <- list(edge = edge, edge.length = ans[[5]], tip.label = labels)
    class(obj) <- "phylo"
    read.tree(text = write.tree(obj, multi.line = FALSE))
}
