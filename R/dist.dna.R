### dist.dna.R  (2004-08-31)
###
###     Pairwise Distances from DNA Sequences
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

dist.dna <- function(x, y = NULL, variance = FALSE, gamma = NULL,
                     method = "Kimura", basefreq = NULL, GCcontent = NULL)
{
    if (!is.null(y)) x <- cbind(x, y)
    if (is.list(x)) {
        nm <- names(x)
        n <- length(x)
        if (length(unique(unlist(lapply(x, length)))) != 1)
            stop("DNA sequences in list not of the same length")
        x <- unlist(x)
        nL <- length(x)
        dim(x) <- c(nL / n, n)
        colnames(x) <- nm
    }
    n <- ncol(x)
    foo <- paste("dist.dna.", method, sep = "")
    D <- matrix(NA, n, n)
    diag(D) <- 0
    if (variance) var.D <- D
    if (method == "JukesCantor" | method == "Kimura") {
        expr <- parse(text = paste(foo, "(x[, i], x[, j], ", "variance = ", variance,
                                   ", gamma = ", gamma, ")", sep = ""))
    }
    if (method == "TajimaNei") {
        g <- if (is.null(basefreq)) base.freq(x) else basefreq
        expr <- parse(text = paste(foo, "(x[, i], x[, j], ", "variance = ", variance,
                                   ", basefreq = g)", sep = ""))
    }
    if (method == "Tamura") {
        GC <- if (is.null(GCcontent)) GC.content(x) else GCcontent
        expr <- parse(text = paste(foo, "(x[, i], x[, j], ", "variance = ", variance,
                                   ", GCcontent = GC)", sep = ""))
    }
    if (method == "TamuraNei") {
        g <- if (is.null(basefreq)) base.freq(x) else basefreq
        expr <- parse(text = paste(foo, "(x[, i], x[, j], ", "variance = ", variance,
                                   ", gamma = ", gamma, ", basefreq = g)", sep = ""))
    }
    for (i in 1:(n - 1)) for (j in (i + 1):n) {
        ## if both sequences are identical, do not call the function...
        d <- if (all(x[, i] == x[, j])) c(0, 0) else eval(expr)
        D[i, j] <- D[j, i] <- d[1]
        if (variance) var.D[i, j] <- var.D[j, i] <- d[2]
    }
    if (variance) {
        if (!is.null(colnames(x)))
          rownames(D) <- colnames(D) <- rownames(var.D) <- colnames(var.D) <- colnames(x)
        return(list(D = D, var.D = var.D))
    }
    else {
        if (!is.null(colnames(x))) rownames(D) <- colnames(D) <- colnames(x)
        return(D)
    }
}

dist.dna.JukesCantor <- function(x, y, variance = FALSE, gamma = NULL)
{
    L <- length(x)
    Nd <- sum(x != y)
    p <- Nd / L
    D <-if (is.null(gamma)) -0.75 * log(1 - 4 * p / 3) else 0.75 * gamma * ((1 - 4 * p / 3)^(1 / gamma) - 1)
    if (variance) {
        var.D <- if (is.null(gamma)) p * (1 - p) / ((1 - 4 * p / 3)^2 * L) else p * (1 - p) / ((1 - 4 * p / 3)^(-2 / (gamma + 1)) * L)
        return(c(D, var.D))
    }
    else return(D)
}

dist.dna.TajimaNei <- function(x, y, variance = FALSE, basefreq = NULL)
{
    Nd <- sum(x != y)
    if (Nd == 0) {
        if (variance) return(c(0, 0)) else return(0)
    } else {
        g <- if (is.null(basefreq)) base.freq(c(x, y)) else basefreq
        L <- length(x)
        p <- Nd / L
        X <- table(x, y) / L
        sel <- col(X) > row(X)
        c <- sum((X[sel] + t(X)[sel])^2 / (2 * c(g[1] * g[2:4], g[2] * g[3:4], g[3] * g[4])))
        b <- (1 - sum(g^2) + p^2 / c) / 2
        D <- -b * log(1 - p / b)
        if (variance) {
            var.D <- b^2 * p * (1 - p) / ((b - p)^2 * L)
            return(c(D, var.D))
        }
        else return(D)
    }
}

dist.dna.Kimura <- function(x, y, variance = FALSE, gamma = NULL)
{
    d <- x != y
    Nd <- sum(d)
    if (Nd == 0) {
        if (variance) return(c(0, 0)) else return(0)
    } else {
        L <- length(x)
        pw.diff <- cbind(x[d], y[d])
        PuPy <- ifelse(pw.diff == "a" | pw.diff == "g", "R", "Y")
        Nv <- sum(PuPy[, 1] != PuPy[, 2])
        Ns <- Nd - Nv
        P <- Ns / L
        Q <- Nv / L
        a1 <- 1 - 2 * P - Q
        a2 <- 1 - 2 * Q
        if (is.null(gamma)) D <- -0.5 * log(a1 * sqrt(a2))
        else {
            b <- -1 / gamma
            D <- gamma * (a1^b + 0.5 * a2^b - 1.5) / 2
        }
        if (variance) {
            if (is.null(gamma)) {
                c1 <- 1 / a1
                c2 <- 1 / a2
                c3 <- (c1 + c2) / 2
            }
            else {
                b <- -(1 / gamma + 1)
                c1 <- a1^b
                c2 <- a2^b
                c3 <- (c1 + c2) / 2            
            }
            var.D <- (c1^2 * P + c3^2 * Q - (c1 * P + c3 * Q)^2) / L
            return(c(D, var.D))
        }
        else return(D)
    }
}

dist.dna.Tamura <- function(x, y, variance = FALSE, GCcontent = NULL)
{
    d <- x != y
    Nd <- sum(d)
    if (Nd == 0) {
        if (variance) return(c(0, 0)) else return(0)
    } else {
        GC <- if (is.null(GCcontent)) GC.content(c(x, y)) else GCcontent
        L <- length(x)
        pw.diff <- cbind(x[d], y[d])
        PuPy <- ifelse(pw.diff == "a" | pw.diff == "g", "R", "Y")
        Nv <- sum(PuPy[, 1] != PuPy[, 2])
        Ns <- Nd - Nv
        P <- Ns / L
        Q <- Nv / L
        wg <- 2 * GC * (1 - GC)
        a1 <- 1 - P / wg - Q
        a2 <- 1 - 2 * Q
        D <- -wg * log(a1) - 0.5 * (1 - wg) * log(a2)
        if (variance) {
            c1 <- 1 / a1
            c2 <- 1 / a2
            c3 <- wg * (c1 - c2) + c2
            var.D <- (c1^2 * P + c3^2 * Q - (c1 * P + c3 * Q)^2) / L
            return(c(D, var.D))
        }
        else return(D)
    }
}

dist.dna.TamuraNei <- function(x, y, variance = FALSE, basefreq = NULL, gamma = NULL)
{
    d <- x != y
    Nd <- sum(d)
    if (Nd == 0) {
        if (variance) return(c(0, 0)) else return(0)
    } else {
        g <- if (is.null(basefreq)) base.freq(c(x, y)) else basefreq
        L <- length(x)
        gR <- g[1] + g[3]
        gY <- g[2] + g[4]
        k1 <- 2 * g[1] * g[3] / gR
        k2 <- 2 * g[2] * g[4] / gY
        k3 <- 2 * (gR * gY - g[1] * g[3] * gY / gR - g[2] * g[4] * gR / gY)
        pw.diff <- cbind(x[d], y[d])
        PuPy <- ifelse(pw.diff == "a" | pw.diff == "g", "R", "Y")
        Nv <- sum(PuPy[, 1] != PuPy[, 2])
        Ns <- Nd - Nv
        P <- Ns / L
        Q <- Nv / L
        P1 <- (sum(pw.diff[, 1] == "a" & pw.diff[, 2] == "g") + sum(pw.diff[, 1] == "g" & pw.diff[, 2] == "a")) / L
        P2 <- P - P1
        w1 <- 1 - P1 / k1 - Q / (2 * gR)
        w2 <- 1 - P2 / k2 - Q / (2 * gY)
        w3 <- 1 - Q / (2 * gR * gY)
        if (is.null(gamma)) {
            k4 <- 2 * ((g[1]^2 + g[3]^2) / (2 * gR^2) + (g[3]^2 + g[4]^2) / (2 * gY^2))
            c1 <- 1 / w1
            c2 <- 1 / w2
            c3 <- 1 / w3
            c4 <- k1 * c1 / (2 * gR) + k2 * c2 / (2 * gY) + k4 * c3
            D <- -k1 * log(w1) - k2 * log(w2) - k3 * log(w3)
        }
        else {
            k4 <- 2 * (g[1] * g[3] + g[2] * g[4] + gR * gY)
            b <- -1 / gamma
            c1 <- w1^b
            c2 <- w2^b
            c3 <- w3^b
            c4 <- k1 * c1 / (2 * gR) + k2 * c2 / (2 * gY) + k3 * c3 / (2 * gR * gY)
            D <- gamma * (k1 * w1^b + k2 * w2^b + k3 * w3^b - k4)
        }
        if (variance) {
            var.D <- (c1^2 * P1 + c2^2 * P2 + c4^2 * Q - (c1 * P1 + c2 * P2 + c4 * Q)^2) / L
            return(c(D, var.D))
        }
        else return(D)
    }
}

