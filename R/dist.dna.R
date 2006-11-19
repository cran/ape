### dist.dna.R (2005-11-13)
###
###    Pairwise Distances from DNA Sequences
###
### Copyright 2002-2006 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

dist.dna <- function(x, model = "K80", variance = FALSE,
                     gamma = FALSE, pairwise.deletion = FALSE,
                     base.freq = NULL, as.matrix = FALSE)
{
    MODELS <- c("raw", "JC69", "K80", "F81", "K81",
                "F84", "T92", "TN93", "GG95", "BH87",
                "logdet", "paralin")
    imod <- which(MODELS == model)
    if (gamma && imod %in% c(1, 5:7, 9:12)) {
        warning(paste("gamma-correction not available for model", model))
        gamma <- FALSE
    }
    ## if the data are in matrix-form, transpose them to have
    ## individuals as cols and sites as rows
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.matrix(x)) x <- t(x)
    ## if the data are given as a list, we transform it as a matrix
    ## with individuals as cols and sites as rows
    if (is.list(x)) {
        nm <- names(x)
        n <- length(x)
        if (length(unique(unlist(lapply(x, length)))) != 1)
            stop("DNA sequences in list not of the same length.")
        x <- unlist(x)
        nL <- length(x)
        dim(x) <- c(nL / n, n)
        colnames(x) <- nm
    }
    BF <- if (is.null(base.freq)) base.freq(x) else base.freq
    if (!pairwise.deletion) {
        sel <- !apply(x, 1, function(x) any(x %in% c("n", "-")))
        x <- x[sel, ]
    }
    s <- nrow(x) # the number of sites
    n <- ncol(x) # the number of individuals
    var <- if (variance) numeric(n*(n - 1)/2) else 0
    if (imod <= 9) {
        if (!gamma) gamma <- alpha <- 0
        else {
            alpha <- gamma
            gamma <- 1
        }
        d <- .C("dist_dna", as.character(x), as.integer(n), as.integer(s),
                as.integer(imod), as.double(numeric(n*(n - 1)/2)),
                as.double(BF), as.integer(pairwise.deletion),
                as.integer(variance), as.double(var),
                as.integer(gamma), as.double(alpha),
                NAOK = TRUE, PACKAGE = "ape")
        if (variance) var <- d[[9]]
        d <- d[[5]]
    } else {
        d <- switch(imod - 9, dist.dna.BH87(x, n, s, variance),
                    dist.dna.logdet(x, n, s, variance),
                    dist.dna.paralin(x, n, s, variance))
        if (variance) {
            var <- d[[2]]
            d <- d[[1]]
        }
    }
    if (imod != 10) {
        attr(d, "Size") <- n
        attr(d, "Labels") <- colnames(x)
        attr(d, "Diag") <- attr(d, "Upper") <- FALSE
        attr(d, "call") <- match.call()
        attr(d, "method") <- model
        class(d) <- "dist"
        if (as.matrix) d <- as.matrix(d)
    }
    if (variance) attr(d, "variance") <- var
    d
}

dist.dna.BH87 <- function(x, n, s, variance)
{
    d <- matrix(0, n, n)
    if (variance) var <- d
    baz <- c("a", "c", "g", "t")
    for (i in 1:(n - 1)) {
        f1 <- base.freq(x[, i])
        for (j in (i + 1):n) {
            if (i == j) next
            f2 <- base.freq(x[, j])
            y <- table(x[, i], x[, j])[baz, baz]
            y <<- y
            P12 <- y/rowSums(y)#apply(y, 2, "*", f1)
            P12 <<- P12
            y <- t(y)
            P21 <- y/rowSums(y)#apply(t(y), 2, "*", f2)
            d[i, j] <- -log(det(P12))/4
            d[j, i] <- -log(det(P21))/4
            if (variance) {
                M12 <- solve(P12)
                M21 <- solve(P21)
                var[i, j] <- rep(f1, 4)*(t(M12)^2*P12 - 1)/(16*s)
                var[j, i] <- rep(f2, 4)*(t(M21)^2*P21 - 1)/(16*s)
            }
        }
    }
    if (variance) d <- list(d, var)
    d
}

dist.dna.logdet <- function(x, n, s, variance)
{
    d <- numeric(n*(n - 1)/2)
    if (variance) var <- d
    baz <- c("a", "c", "g", "t")
    k <- 1
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            P <- table(x[, i], x[, j])[baz, baz]/s
            d[k] <- (-log(det(P))/4 - log(4))/4
            if (variance)
              var[k] <- sum(solve(P)^2 * t(P) - 1)/(16*s)
            k <- k + 1
        }
    }
    if (variance) d <- list(d, var)
    d
}

dist.dna.paralin <- function(x, n, s, variance)
{
    d <- numeric(n*(n - 1)/2)
    if (variance) var <- d
    baz <- c("a", "c", "g", "t")
    k <- 1
    for (i in 1:(n - 1)) {
        f1 <- base.freq(x[, i])
        for (j in (i + 1):n) {
            f2 <- base.freq(x[, j])
            P <- table(x[, i], x[, j])[baz, baz]/s
            d[k] <- -0.25*log(det(P)/sqrt(prod(f1)*prod(f2)))
            if (variance)
              var[k] <- sum(solve(P)^2 * t(P) - 1/
                            sqrt(rep(f1, each = 4)*rep(f2, each = 4)))/(16*s)
            k <- k + 1
        }
    }
    if (variance) d <- list(d, var)
    d
}
