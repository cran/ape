## SDM.R (2011-10-11)

## Construction of Consensus Distance Matrix With SDM

## Copyright 2011 Andrei-Alin Popescu

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

SDM <- function(...)
{
    st <- list(...) # first half contains matrices, second half s_p
    k <- length(st)/2
    labels <- attr(as.dist(st[[1]]), "Labels")
    tot <- length(rownames(st[[1]]))
    for (i in 2:k) {
        labels <- union(labels, attr(as.dist(st[[i]]), "Labels"))
        tot <- tot + length(rownames(st[[i]]))
    }
    sp <- mat.or.vec(1,k)
    for (i in c(k+1:k))
        sp[i - k] <- st[[i]]

    astart <- mat.or.vec(1, tot) # start of aip, astart[p] is start of aip
    astart[1] <- k
    for (i in 2:k)
        astart[i] <- astart[i - 1] + length(rownames(st[[i - 1]]))

    a <- mat.or.vec(1, k + tot + k + length(labels))
    ## first k are alphas, subseqeunt ones aip
    ## each matrix p starting at astart[p], next are
    ## Lagrange multipliers, miu, niu, lambda in that order
    n <- length(labels)
    miustart <- k + tot
    niustart <- miustart + n
    lambstart <- niustart + k - 1

    X <- mat.or.vec(n, n)
    V <- mat.or.vec(n, n)
    w <- mat.or.vec(n, n)

    col <- mat.or.vec(k + tot + k + length(labels), 1) # free terms of system

    dimnames(X) <- list(labels, labels)
    dimnames(V) <- list(labels, labels)
    dimnames(W) <- list(labels, labels)

    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            for (p in 1:k) {
                d <- st[[p]]
                if (is.element(rownames(X)[i], rownames(d)) & is.element(colnames(X)[j], colnames(d))) {
                    w[i, j] <- w[j, i] <- w[i, j] + sp[p]
                }
            }
        }
    }

    Q <- mat.or.vec(length(labels) + k + k + tot, length(labels) + k + k + tot)
    ## first decompose first sum in paper
    for (p in 1:k) {
        d_p <- st[[p]]
        for (l in 1:k) { # first compute coeficients of alphas
            sum <- 0
            dijp <- -1
            if (l == p) { # calculate alpha_p
                for (i in 1:n) {
                    for (j in 1:n) { #check if {i,j}\subset L_l
                        d <- st[[l]]
                        if (i != j & is.element(labels[i],rownames(d)) & is.element(labels[j],colnames(d))) {
                            dij <- d[rownames(d) == labels[i], colnames(d) == labels[j]]
                            sum <- sum + ((dij*dij) - sp[l]*dij*dij/w[i,j])
                            ipos <- which(rownames(d) == labels[i])
                            jpos <- which(rownames(d) == labels[j])
                            Q[p, astart[l] + ipos] <- Q[p, astart[l] + ipos] + (dij - (sp[l]*dij/w[i,j]))
                            Q[p, astart[l] + jpos] <- Q[p, astart[l] + jpos] + (dij - (sp[l]*dij/w[i,j]))
                        }
                    }
                }
            } else {
                for (i in 1:n) {
                    for (j in 1:n) { #check if {i,j}\subset L_l
                        d <- st[[l]]
                        if (i != j & is.element(labels[i], rownames(d)) & is.element(labels[j], colnames(d)) & is.element(labels[i], rownames(d_p)) & is.element(labels[j], colnames(d_p))) {
                            dij <- d[rownames(d) == labels[i], colnames(d) == labels[j]]
                            dijp <- d_p[rownames(d_p) == labels[i], colnames(d_p) == labels[j]]
                            sum <- sum - sp[l]*dij*dijp/w[i, j]
                            ipos <- which(rownames(d) == labels[i])
                            jpos <- which(rownames(d) == labels[j])
                            Q[p,astart[l] + ipos] <- Q[p, astart[l] + ipos] - sp[l]*dijp/w[i, j]
                            Q[p,astart[l] + jpos] <- Q[p, astart[l] + jpos] - sp[l]*dijp/w[i, j]
                        }
                    }
                }
            }
            Q[p, l] <- sum
        }
        Q[p, lambstart + 1] <- 1
    }
    r <- k
    for (p in 1:k) {
        dp <- st[[p]]
        for (i in 1:n) {
            if (is.element(labels[i], rownames(dp))) {
                r <- r + 1
                for (l in 1:k) {
                    d <- st[[l]]
                    if (l == p) {
                        for (j in 1:n) {
                            if (i != j & is.element(labels[j], rownames(dp))) {
                                dij <- d[rownames(d) == labels[i], colnames(d) == labels[j]]
                                Q[r, l] <- Q[r, l] + (dij - sp[l]*dij/w[i, j])
                                ipos <- which(rownames(d) == labels[i])
                                jpos <- which(rownames(d) == labels[j])
                                Q[r, astart[l] + ipos] <- Q[r, astart[l] + ipos] + (1 - sp[l]/w[i, j])
                                Q[r, astart[l] + jpos] <- Q[r, astart[l] + jpos] + (1 - sp[l]/w[i, j])
                            }
                        }
                    } else {
                        for (j in 1:n) {
                            if (i != j & is.element(labels[j], rownames(dp)) & is.element(labels[i], rownames(d)) & is.element(labels[j], colnames(d))) {
                                dij <- d[rownames(d) == labels[i], colnames(d) == labels[j]]
                                Q[r,l] <- Q[r,l] - sp[l]*dij/w[i, j]
                                ipos <- which(rownames(d) == labels[i])
                                jpos <- which(rownames(d) == labels[j])
                                Q[r, astart[l] + ipos] <- Q[r, astart[l] + ipos] - sp[l]/w[i, j]
                                Q[r, astart[l] + jpos] <- Q[r, astart[l] + jpos] - sp[l]/w[i, j]
                            }
                        }
                    }
                }
                if (p < k) Q[r, ] <- Q[r, ] * sp[p]
                Q[r, miustart + i] <- 1
                if (p < k) Q[r, niustart + p] <- 1
            }
        }
    }
    r <- r + 1
    col[r] <- k
    for (i in 1:k) Q[r,i] <- 1

    for (i in 1:n) {
        r <- r + 1
        for (p in 1:k) {
            d <- st[[p]]
            if (is.element(labels[i], rownames(d))) {
                ipos <- which(rownames(d) == labels[i])
                Q[r, astart[p] + ipos] <- 1
            }
        }
    }
    for (p in 1:(k - 1)) {
        r <- r + 1
        for (i in 1:n) {
            d <- st[[p]]
            if (is.element(labels[i], rownames(d))) {
                ipos <- which(rownames(d) == labels[i])
                Q[r, astart[p] + ipos] <- 1
            }
        }
    }
    a <- solve(Q, col, 1e-19)
    for(i in 1:n) {
        for(j in 1:n) {
            sum <- 0
            sumv <- 0
            for(p in 1:k) {
                d <- st[[p]]
                if (is.element(labels[i], rownames(d)) & is.element(labels[j], rownames(d))) {
                    ipos <- which(rownames(d) == labels[i])
                    jpos <- which(rownames(d) == labels[j])
                    sum <- sum + sp[p]*(a[p]*d[ipos, jpos] + a[astart[p] + ipos] + a[astart[p] + jpos])
                    sumv <- sumv + sp[p]*(a[p]*d[ipos, jpos])*(a[p]*d[ipos, jpos])
                }
            }
            X[i, j] <- sum/w[i, j]
            V[i, j] <- sumv/(w[i, j] * w[i, j])
            if (i == j)
                X[i, j] <- V[i, j] <- 0
        }
    }
    list(X, V)
}
