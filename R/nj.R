### nj.R  (2004-12-07)
###
###        Neighbor-Joining Tree Estimation
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

nj <- function(X)
{
    foo <- function(i, j) {
        ## it is assumed that i != j (always)
        if (length(i) > length(j)) j <- rep(j, length = length(i))
        if (length(i) < length(j)) i <- rep(i, length = length(j))
        comp <- i > j
        if (any(comp)) {
            tmp <- i
            i[comp] <- j[comp]
            j[comp] <- tmp[comp]
        }
        n * (i - 1) - i * (i - 1) / 2 + j - i
    }
    X <- as.dist(X)
    N <- n <- attr(X, "Size") # n will be updated
    labels <- attr(X, "Labels")
    OTU <- 1:n
    edge <- matrix(NA, 2 * n - 3, 2) # the NJ tree is unrooted!
    edge.length <- numeric(2 * n - 3)
    cur.nod <- -n + 2 # the NJ tree is unrooted! (bis)
    DI <- numeric(n - 2) # store the 'internal' distances of a OTU to correct the internal branch lengths
    j <- 1
    for (i in 1:(N - 3)) {
        rownum <- numeric(0)
        for (i in 2:n) rownum <- c(rownum, i:n)
        colnum <- rep(1:(n - 1), (n - 1):1)
        SUMD <- sum(X)
        pX <- matrix(NA, n, n - 1) # D_ij with j != i
        pX[1, ] <- X[1:(n - 1)] # for j = 2 TO n
        for (k in 2:(n - 1)) {
            sel <- c(n * (1:(k-1) - 1) - 1:(k-1) * (1:(k-1) - 1) / 2 + k - 1:(k-1),
                     n * (k - 1) - k * (k - 1) / 2 + (k+1):n - k)
            pX[k, ] <- X[sel]
        }
        pX[n, ] <- X[n * (1:(n-1) - 1) - 1:(n-1) * (1:(n-1) - 1) / 2 + n - 1:(n-1)]
        S <- numeric(length(X))
        for (ki in 1:(n - 1)) {
            for (kj in (ki + 1):n) {
                ind <- n * (ki - 1) - ki * (ki - 1) / 2 + kj - ki
                A <- sum(pX[ki, ]) - X[ind]
                B <- sum(pX[kj, ]) - X[ind]
                S[ind] <- (A + B) / (2 * n - 4) + 0.5 * X[ind] +
                  (SUMD - A - B - X[ind]) / (n - 2)
            }
        }
        smallest <- which.min(S)
        is <- rownum[smallest]
        js <- colnum[smallest]
        edge[j:(j + 1), 1] <- cur.nod
        edge[j, 2] <- OTU[is]
        edge[j + 1, 2] <- OTU[js]
        if (is < js) {
            xi <- pX[is, ][-(js - 1)]
            xj <- pX[js, ][-is]
        } else {
            xi <- pX[is, ][-js]
            xj <- pX[js, ][-(is - 1)]
        }
        A <- sum(xi) / (n - 2)
        B <- sum(xj) / (n - 2)
        edge.length[j] <- (X[smallest] + A - B) / 2
        edge.length[j + 1] <- (X[smallest] + B - A) / 2
        ## now update the distance matrix
        new <- (xi + xj) / 2
        DI[-cur.nod] <- X[smallest]
        del <- c(foo(is, (1:n)[-is]),
                 foo(js, (1:n)[-c(is, js)]))
        X <- X[-del]
        OTU <- OTU[-c(is, js)]
        X <- c(new, X)
        OTU <- c(cur.nod, OTU)
        cur.nod <- cur.nod + 1
        n <- n - 1
        j <- j + 2
    }
    roots <- (length(edge.length) - 2):length(edge.length)
    edge[roots, 1] <- -1
    edge[roots, 2] <- OTU
    edge.length[roots] <- c((X[1] + X[2] - X[3]) / 2,
                            (X[1] + X[3] - X[2]) / 2,
                            (X[3] + X[2] - X[1]) / 2)
    ## estimate internal branch lengths
    int <- which(edge[, 2] < 0)
    edge.length[int] <- edge.length[int] - DI[-edge[int, 2]] / 2
    ## finalize the "phylo" object
    mode(edge) <- "character"
    obj <- list(edge = edge, edge.length = edge.length, tip.label = labels)
    class(obj) <- "phylo"
    obj
    read.tree(text = write.tree(obj, multi.line = FALSE))
}
