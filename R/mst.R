### mst.R  (2003-02-20)
###
###     Minimum Spanning Tree
###
### Copyright 2002 Yvonnick Noel <noel@univ-lille3.fr>,
###     Julien Claude <claude@isem.univ-montp2.fr> and
###     Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

mst <- function(X)
{
    if (class(X) == "dist") X <- as.matrix(X)
    n <- dim(X)[1]
    N <- matrix(0, n, n)
    tree <- NULL
    large.value <- max(X) + 1
    diag(X) <- large.value
    index.i <- 1

    for (i in 1:(n - 1)) {
        tree <- c(tree, index.i)
        m <- apply(as.matrix(X[, tree]), 2, min)  #calcul les minimum par colonne
        a <- sort.index(X[, tree])[1, ]
        b <- sort.index(m)[1]
        index.j <- tree[b]
        index.i <- a[b]
        
        N[index.i, index.j] <- 1
        N[index.j, index.i] <- 1
        
        for (j in tree) {
            X[index.i, j] <- large.value
            X[j, index.i] <- large.value
        }
    }
    dimnames(N) <- dimnames(X)
    class(N) <- "mst"
    return(N)
}

### Function returning an index matrix for an increasing sort
sort.index <- function(X)
{
    if(length(X) == 1) return(1)                  # sorting a scalar?
    if(!is.matrix(X)) X <- as.matrix(X)           # force vector into matrix
    n <- nrow(X)
    apply(X, 2, function(v) order(rank(v)))       # find the permutation
}

plot.mst <- function(x, graph = "circle", x1 = NULL, x2 = NULL, ...)
{
    m <- x
    rm(x)
    if (is.null(x1) | is.null(x2)) {
        if (graph == "circle") {
            ang <- seq(0, 2 * pi, length = n + 1)
            x <- cos(ang)
            y <- sin(ang)
            plot(x, y, type = "n", xlab = "", ylab = "",
                 xaxt = "n", yaxt = "n", bty = "n", ...)
        }
        if (graph == "nsca") {
            XY <- nsca(m)
            x <- XY[, 1]
            y <- XY[, 2]
            plot(XY, type = "n", xlab = "\"nsca\" - axis 1",
                 ylab = "\"nsca\" - axis 2", ...)
        }
    }
    else {
        x <- x1
        y <- x2
        plot(x, y, type = "n", xlab = deparse(substitute(x1)),
             ylab = deparse(substitute(x2)), ...)
    }
    for (i in 1:n) {
        w1 <- which(m[i, ] == 1)
        segments(x[i], y[i], x[w1], y[w1])
    }
    points(x, y, pch = 21, col = "black", bg = "white", cex = 3)
    text(x, y, 1:n, cex = 0.8)
}

nsca <- function(A)
{
    Dr <- apply(A, 1, sum)
    Dc <- apply(A, 2, sum)

    eig.res <- eigen(diag(1 / sqrt(Dr)) %*% A %*% diag(1 / sqrt(Dc)))
    r <- diag(1 / Dr) %*% (eig.res$vectors)[, 2:4]
    ## The next line has been changed by EP (20-02-2003)
    ## dimnames(r)[[1]] <- dimnames(A)[[1]]
    rownames(r) <- rownames(A)
    r
}
