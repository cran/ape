### rtree.R (2005-06-14)
###
###          Generates Random Trees
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

rtree <- function(n, rooted = TRUE, tip.label = NULL, br = runif, ...)
{
    nbr <- 2 * n - 2
    if (!rooted) nbr <- nbr - 1
    edge <- matrix(NA, nbr, 2)
    if (n == 2) {
        if (rooted) {
            edge[, 1] <- -1
            edge[, 2] <- 1:2
        } else edge[1, ] <- 1:2
    } else if (n == 3) {
        if (rooted) {
            edge[, 1] <- c(-1, -2, -2, -1)
            edge[, 2] <- c(-2, 1:3)
        } else {
            edge[, 1] <- -1
            edge[, 2] <- 1:3
        }
    } else if (n == 4 && !rooted) {
        edge[, 1] <- c(-1, -2, -2, -1, -1)
        edge[, 2] <- c(-2, 1:4)
    } else {
        if (rooted) {
            edge[1:2, 1] <- 1
            rnd <- ceiling(runif(n - 2, 0, 2:(n - 1)))
            for (i in 2:(n - 1)) {
                edge[which(is.na(edge[1:((i - 1) * 2), 2]))[rnd[i - 1]], 2] <- i
                edge[(i - 1) * 2 + 1, 1] <- i
                edge[(i - 1) * 2 + 2, 1] <- i
            }
        } else {
            edge[1:3, 1] <- 1                          # != the rooted case
            rnd <- ceiling(runif(n - 3, 0, 2:(n - 3))) # != the rooted case
            for (i in 2:(n - 2)) {
                edge[which(is.na(edge[1:((i - 1) * 2), 2]))[rnd[i - 1]], 2] <- i
                edge[(i - 1) * 2 + 2, 1] <- i # these two lines differ
                edge[(i - 1) * 2 + 3, 1] <- i # from the rooted case
            }
        }
        edge <- -edge
        edge[, 2][is.na(edge[, 2])] <- 1:n
    }
    mode(edge) <- "character"
    phy <- list(edge = edge)
    phy$tip.label <-
      if (is.null(tip.label)) paste("t", 1:n, sep = "")
      else sample(tip.label)
    class(phy) <- "phylo"
    if (is.function(br)) phy$edge.length <- br(nbr, ...)
    read.tree(text = write.tree(phy, multi.line = FALSE))
}

rcoal <- function(n, tip.label = NULL, br = rexp, ...)
{
    nbr <- 2 * n - 2
    edge <- matrix(NA, nbr, 2)
    x <- br(n - 1, ...) # coalescence times
    if (n == 2) {
        edge[, 1] <- -1
        edge[, 2] <- 1:2
        edge.length <- rep(x, 2)
    } else if (n == 3) {
        edge[, 1] <- c(-1, -2, -2, -1)
        edge[, 2] <- c(-2, 1:3)
        edge.length <- c(x[2], x[1], x[1], sum(x))
    } else {
        edge.length <- numeric(nbr)
        h <- numeric(2 * n - 1)
        names(h) <- as.character(c(1:n, -(1:(n - 1))))
        node.height <- cumsum(x)
        pool <- 1:n
        nextnode <- 1 - n
        for (i in 1:(n - 1)) {
            y <- sample(pool, size = 2)
            edge[(i - 1) * 2 + 1:2, 2] <- y
            edge[(i - 1) * 2 + 1:2, 1] <- nextnode
            edge.length[(i - 1) * 2 + 1:2] <- node.height[i] - h[as.character(y)]
            h[as.character(nextnode)] <- node.height[i]
            pool <- c(pool[! pool %in% y], nextnode)
            nextnode <- nextnode + 1
        }
    }
    mode(edge) <- "character"
    phy <- list(edge = edge, edge.length = edge.length)
    phy$tip.label <-
      if (is.null(tip.label)) paste("t", 1:n, sep = "")
      else tip.label
    class(phy) <- "phylo"
    read.tree(text = write.tree(phy, multi.line = FALSE))
}
