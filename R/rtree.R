### rtree.R  (2004-05-23)
###
###                  Generates Random Trees
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

rtree <- function(n, br = TRUE)
{
    edge <- rep(NA, 4 * n - 4)
    dim(edge) <- c(2 * n - 2, 2)
    mode(edge) <- "numeric"
    edge[1:2, 1] <- 1

    rnd <- ceiling(runif(n - 2, 0, 2:(n - 1)))

    for (i in 2:(n - 1)) {
        edge[which(is.na(edge[1:((i - 1) * 2), 2]))[rnd[i - 1]], 2] <- i
        edge[(i - 1) * 2 + 1, 1] <- i
        edge[(i - 1) * 2 + 2, 1] <- i
    }
    edge <- -edge
    edge[, 2][is.na(edge[, 2])] <- 1:n
    mode(edge) <- "character"
    phy <- list(edge = edge, tip.label = paste("t", 1:n, sep = ""))
    class(phy) <- "phylo"
    if (br) phy$edge.length <- runif(2 * n - 2)
    read.tree(text = write.tree(phy))
}
