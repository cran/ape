### zoom.R  (2004-12-17)
###
###     Zoom on a Portion of a Phylogeny
###
### Copyright 2004 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

zoom <- function(phy, focus, subtree = FALSE, col = rainbow, ...)
{
    if (!is.list(focus)) focus <- list(focus)
    n <- length(focus)
    for (i in 1:n)
      if (is.character(focus[[i]]))
        focus[[i]] <- which(phy$tip.label == focus[[i]])
    if (is.function(col))
      if (deparse(substitute(col)) == "grey")
        col <- grey(1:n/n) else col <- col(n)
    ext <- list()
    length(ext) <- n
    for (i in 1:n)
      ext[[i]] <- drop.tip(phy, phy$tip.label[-focus[[i]]],
                           subtree = subtree)
    nc <- round(sqrt(n)) + 1
    nr <- ceiling(sqrt(n))
    M <- matrix(0, nr, nc)
    x <- c(rep(1, nr), 2:(n + 1))
    M[1:length(x)] <- x
    layout(M, c(1, rep(3 / (nc - 1), nc - 1)))
    phy$tip.label <- rep("", length(phy$tip.label))
    colo <- rep("black", dim(phy$edge)[1])
    for (i in 1:n)
      colo[which.edge(phy, focus[[i]])] <- col[i]
    plot.phylo(phy, edge.color = colo, ...)
    for (i in 1:n)
      plot.phylo(ext[[i]], edge.color = col[i], ...)
}
