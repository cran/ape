### dist.gene.pairwise.R  (2002-08-28)
###
###     Pairwise Distances from Genetic Data
###
### Copyright 2002 Emmanuel Paradis
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

dist.gene.pairwise <- function(x, variance = FALSE)
{
    if (is.data.frame(x)) x <- as.matrix(x)
    L <- ncol(x)
    n <- nrow(x)
    D <- matrix(NA, n, n)
    diag(D) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            D[i, j] <- D[j, i] <- L - sum(x[i, ] == x[j, ])
        }
    }
    if (!is.null(rownames(x))) rownames(D) <- colnames(D) <- rownames(x)
    if (variance) {
        var.D <- D * (L - D) / L
        return(list(distances = D, variance = var.D))
    }
    else return(D)
}

dist.gene.percentage <- function(x, variance = FALSE)
{
    L <- ncol(x)
    D <- dist.gene.pairwise(x) / L
    if (variance) {
        var.D <- D * (1 - D) / L
        return(list(pairwise.distances = D, variance = var.D))
    }
    else return(D)
}

dist.gene <- function(x, method = "pairwise", variance = FALSE)
{
    if (method == "pairwise")
      return(dist.gene.pairwise(x, variance = variance))
    if (method == "percentage")
      return(dist.gene.percentage(x, variance = variance))
}
