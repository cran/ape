## DNA.R (2007-05-01)
##
##    Comparisons and Manipulations of DNA Sequences
##
## Copyright 2002-2007 Emmanuel Paradis
##
## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

as.alignment <- function(x)
{
    if (is.list(x)) n <- length(x)
    if (is.matrix(x)) n <- dim(x)[1]
    seq <- character(n)
    if (is.list(x)) {
        nam <- names(x)
        for (i in 1:n)
          seq[i] <- paste(x[[i]], collapse = "")
    }
    if (is.matrix(x)) {
        nam <- dimnames(x)[[1]]
        for (i in 1:n)
          seq[i] <- paste(x[i, ], collapse = "")
    }
    obj <- list(nb = n, seq = seq, nam = nam, com = NA)
    class(obj) <- "alignment"
    obj
}

"[.DNAbin" <- function(x, i, j, drop = TRUE)
{
    class(x) <- NULL
    if (is.matrix(x)) {
        if (nargs() == 2 && !missing(i)) ans <- x[i]
        else {
            nd <- dim(x)
            if (missing(i)) i <- 1:nd[1]
            if (missing(j)) j <- 1:nd[2]
            ans <- x[i, j, drop = drop]
        }
    } else {
        if (missing(i)) i <- 1:length(x)
        ans <- x[i]
    }
    structure(ans, class = "DNAbin")
}

rbind.DNAbin <- function(...)
### works only with matrices for the moment
{
    obj <- list(...)
    nobj <- length(obj)
    if (nobj == 1) stop("only one matrix to bind.")
    NC <- ncol(obj[[1]])
    for (i in 2:nobj)
      if(ncol(obj[[i]]) != NC)
        stop("matrices do not have the same number of columns.")
    for (i in 1:nobj) class(obj[[i]]) <- NULL
    ans <- obj[[1]]
    for (i in 2:nobj) ans <- rbind(ans, obj[[i]])
    structure(ans, class = "DNAbin")
}

cbind.DNAbin <- function(..., check.names = TRUE)
### works only with matrices for the moment
{
    obj <- list(...)
    nobj <- length(obj)
    if (nobj == 1) stop("only one matrix to bind.")
    NR <- nrow(obj[[1]])
    for (i in 2:nobj)
      if(nrow(obj[[i]]) != NR)
        stop("matrices do not have the same number of rows.")
    for (i in 1:nobj) class(obj[[i]]) <- NULL
    nms <- rownames(obj[[1]])
    if (check.names) {
        for (i in 2:nobj)
          if (all(rownames(obj[[i]]) %in% nms))
            obj[[i]] <- obj[[i]][nms, ]
        else stop("row names do not match among matrices.")
    }
    ans <- matrix(unlist(obj), NR)
    rownames(ans) <- nms
    structure(ans, class = "DNAbin")
}

print.DNAbin <- function(x, ...)
{
    if (is.list(x)) cat(length(x), "DNA sequences in binary format.\n")
    else if (is.matrix(x)) cat(dim(x)[1], "DNA sequences in binary format.\n")
    else cat("1 DNA sequence in binary format.\n")
}

summary.DNAbin <- function(object, printlen = 6, digits = 3, ...)
{
    if (is.list(object)) {
        n <- length(object)
        nms <- names(object)
        if (n == 1) {
            cat("1 DNA sequence in binary format stored in a list.\n\n")
            cat("Sequence length:", length(object[[1]]), "\n\n")
            cat("Label:", nms, "\n\n")
        } else {
            cat(n, "DNA sequences in binary format stored in a list.\n\n")
            cat("Summary of sequence lengths:\n")
            print(summary(unlist(lapply(object, length))))
            TAIL <- "\n\n"
            if (printlen < n) {
                nms <- nms[1:printlen]
                TAIL <- "...\n\n"
            }
            cat("\nLabels:", paste(nms, collapse = " "), TAIL)
        }
    } else if (is.matrix(object)) {
        nd <- dim(object)
        nms <- rownames(object)
        cat(nd[1], "DNA sequences in binary format stored in a matrix.\n\n")
        cat("All sequences of same length:", nd[2], "\n")
        TAIL <- "\n\n"
        if (printlen < nd[1]) {
            nms <- nms[1:printlen]
            TAIL <- "...\n\n"
        }
        cat("\nLabels:", paste(nms, collapse = " "), TAIL)
    } else {
        cat("1 DNA sequence in binary format stored in a vector.\n\n")
        cat("Sequence length:", length(object), "\n\n")
    }
    cat("Base composition:\n")
    print(round(base.freq(object), digits))
}

as.DNAbin <- function(x, ...) UseMethod("as.DNAbin")

._cs_<- letters[c(1, 7, 3, 20, 18, 13, 23, 19, 11, 25, 22, 8, 4, 2, 14)]

._bs_<- c(136, 72, 40, 24, 192, 160, 144, 96, 80, 48, 224, 176, 208, 112, 240)

as.DNAbin.character <- function(x, ...)
{
    n <- length(x)
    ans <- raw(n)
    for (i in 1:15)
      ans[which(x == ._cs_[i])] <- as.raw(._bs_[i])
    ans[which(x == "-")] <- as.raw(4)
    ans[which(x == "?")] <- as.raw(2)
    if (is.matrix(x)) {
        dim(ans) <- dim(x)
        dimnames(ans) <- dimnames(x)
    }
    class(ans) <- "DNAbin"
    ans
}

as.DNAbin.list <- function(x, ...)
{
    obj <- lapply(x, as.DNAbin)
    class(obj) <- "DNAbin"
    obj
}

as.character.DNAbin <- function(x, ...)
{
    f <- function(xx) {
        ans <- character(length(xx))
        for (i in 1:15)
          ans[which(xx == ._bs_[i])] <- ._cs_[i]
        ans[which(xx == 4)] <- "-"
        ans[which(xx == 2)] <- "?"
        if (is.matrix(xx)) {
            dim(ans) <- dim(xx)
            dimnames(ans) <- dimnames(xx)
        }
        ans
    }
    if (is.list(x)) lapply(x, f) else f(x)
}

base.freq <- function(x)
{
    if (is.list(x)) x <- unlist(x)
    if (is.matrix(x)) dim(x) <- NULL
    n <- length(x)
    BF <- .C("BaseProportion", as.raw(x), as.integer(n),
             double(4), PACKAGE = "ape")[[3]]
    names(BF) <- letters[c(1, 3, 7, 20)]
    BF
}

GC.content <- function(x)
{
    BF <- base.freq(x)
    sum(BF[2:3])
}

seg.sites <- function(x)
{
    n <- dim(x)
    s <- n[2]
    n <- n[1]
    ans <- .C("SegSites", x, as.integer(n), as.integer(s),
              integer(s), PACKAGE = "ape")
    which(as.logical(ans[[4]]))
}

nuc.div <- function(x, variance = FALSE, pairwise.deletion = FALSE)
{
    if (pairwise.deletion && variance)
      warning("cannot compute the variance of nucleotidic diversity\nwith pairwise deletion: try 'pairwise.deletion = FALSE' instead.")

    n <- dim(x)
    s <- n[2]
    n <- n[1]

    ## <FIXME> this should be safely deleted
    if (!pairwise.deletion) {
        keep <- .C("GlobalDeletionDNA", x, as.integer(n),
                   as.integer(s), as.integer(rep(1, s)),
                   PACKAGE = "ape")[[4]]
        x <- x[,  as.logical(keep)]
        s <- dim(x)[2]
    }
    ## </FIXME>

    ans <- .C("NuclearDiversity", x, as.integer(n), as.integer(s),
              as.integer(pairwise.deletion), double(1), PACKAGE = "ape")[[5]]

    if (variance) {
        var <- (n + 1)*ans/(3*(n + 1)*s) + 2*(n^2 + n + 3)*ans/(9*n*(n - 1))
        ans <- c(ans, var)
    }
    ans
}

dist.dna <- function(x, model = "K80", variance = FALSE, gamma = FALSE,
                     pairwise.deletion = FALSE, base.freq = NULL,
                     as.matrix = FALSE)
{
    MODELS <- c("RAW", "JC69", "K80", "F81", "K81", "F84", "T92", "TN93",
                "GG95", "LOGDET", "BH87", "PARALIN")
    imod <- which(MODELS == toupper(model))

    if (imod == 11 && variance) {
        warning("computing variance temporarily not available for model BH87.")
        variance <- FALSE
    }

    if (gamma && imod %in% c(1, 5:7, 9:12)) {
        warning(paste("gamma-correction not available for model", model))
        gamma <- FALSE
    }

    if (is.matrix(x)) {
        nms <- dimnames(x)[[1]]
        n <- dim(x)
        s <- n[2]
        n <- n[1]
    }
    if (is.list(x)) {
        if (length(unique(unlist(lapply(x, length)))) != 1)
          stop("DNA sequences in list not of the same length.")
        nms <- names(x)
        n <- length(x)
        s <- length(x[[1]])
        x <- matrix(unlist(x), n, s, byrow = TRUE)
    }
    BF <- if (is.null(base.freq)) base.freq(x) else base.freq
    if (!pairwise.deletion) {
        keep <- .C("GlobalDeletionDNA", x, as.integer(n),
                   as.integer(s), as.integer(rep(1, s)),
                   PACKAGE = "ape")[[4]]
        x <- x[,  as.logical(keep)]
        s <- dim(x)[2]
    }
    Ndist <- if (imod == 11) n*n else n*(n - 1)/2
    var <- if (variance) double(Ndist) else 0
    if (!gamma) gamma <- alpha <- 0
    else alpha <- gamma <- 1
    d <- .C("dist_dna", x, as.integer(n), as.integer(s),
            as.integer(imod), double(Ndist), BF,
            as.integer(pairwise.deletion), as.integer(variance),
            var, as.integer(gamma), alpha, PACKAGE = "ape")
    if (variance) var <- d[[9]]
    d <- d[[5]]
    if (imod == 11) {
        dim(d) <- c(n, n)
        dimnames(d) <- list(nms, nms)
    } else {
        attr(d, "Size") <- n
        attr(d, "Labels") <- nms
        attr(d, "Diag") <- attr(d, "Upper") <- FALSE
        attr(d, "call") <- match.call()
        attr(d, "method") <- model
        class(d) <- "dist"
        if (as.matrix) d <- as.matrix(d)
    }
    if (variance) attr(d, "variance") <- var
    d
}
