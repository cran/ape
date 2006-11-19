### write.dna.R (2003-11-19)
###
###     Write DNA Sequences in a File
###
### Copyright 2003-2006 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

write.dna <- function(x, file, format = "interleaved", append = FALSE,
                      nbcol = 6, colsep = " ", colw = 10, indent = NULL,
                      blocksep = 1)
{
    format <- match.arg(format, c("interleaved", "sequential", "fasta"))
    if (is.matrix(x)) {
        N <- dim(x)[1]
        xx <- vector("list", N)
        for (i in 1:N) xx[[i]] <- x[i, ]
        names(xx) <- rownames(x)
        x <- xx
        rm(xx)
    } else N <- length(x)
    if (is.null(names(x))) names(x) <- as.character(1:N)
    if (is.null(indent)) {
        indent <- if (format %in% c("interleaved", "sequential")) 10 else  0
    }
    if (indent == "") indent <- 0
    if (is.numeric(indent)) indent <- paste(rep(" ", indent), collapse = "")
    if (format == "interleaved") {
        if (blocksep) {
            blockseparation <- TRUE
            blocksep <- paste(rep("\n", blocksep), collapse = "")
        } else {
            blockseparation <- FALSE
        }
        if (nbcol < 0) format <- "sequential"
    }
    zz <- if (append) file(file, "a") else file(file, "w")
    if (format %in% c("interleaved", "sequential")) {
        S <- unique(unlist(lapply(x, length)))
        ## check that all sequences have the same length
        if (length(S) != 1)
          stop("sequences must have the same length for interleaved or sequential format.")
        ## truncate names if necessary
        if (any(nchar(names(x)) > 10)) {
            warning("at least one name was longer than 10 characters;\nthey will be truncated which may lead to some redundancy.\n")
            names(x) <- substr(names(x), 1, 10)
        }
        for (i in 1:N) {
            nam <- names(x)[i]
            nch <- nchar(nam)
            if (nch < 10)
              names(x)[i] <- paste(nam, paste(rep(" ", 10 - nch), collapse = ""), sep = "")
        }
        cat(N, S, "\n", file = zz)
        if (nbcol < 0) {
            nb.block <- 1
            nbcol <- totalcol <- ceiling(S / colw)
        } else {
            nb.block <- ceiling(S / (colw * nbcol))
            totalcol <- ceiling(S / colw)
        }
        ## Prepare the sequences in a matrix which elements are
        ## strings with `colw' characters.
        SEQ <- matrix(NA, N, totalcol)
        mode(SEQ) <- "character"
        for (i in 1:N) {
            X <- paste(x[[i]], collapse= "")
            for (j in 1:totalcol) SEQ[i, j] <- substr(X, 1 + (j - 1)*colw, colw + (j - 1)*colw)
        }
    }
    if (format == "interleaved") {
        ## Write the first block with the taxon names
        if (nb.block == 1) {
            for (i in 1:N) {
                cat(names(x)[i], file = zz)
                cat(SEQ[i, ], sep = colsep, file = zz)
                cat("\n", file = zz)
            }
        } else {
            for (i in 1:N) {
                cat(names(x)[i], file = zz)
                cat(SEQ[i, 1:nbcol], sep = colsep, file = zz)
                cat("\n", file = zz)
            }
        }
        ## Write the other blocks
        if (nb.block > 1) {
            for (k in 2:nb.block) {
                if (blockseparation) cat(blocksep, file = zz)
                if (k == nb.block) {
                    for (i in 1:N) {
                        cat(indent, file = zz)
                        cat(SEQ[i, (1 + (k - 1)*nbcol):ncol(SEQ)], sep = colsep, file = zz)
                        cat("\n", file = zz)
                    }
                } else {
                    for (i in 1:N) {
                        cat(indent, file = zz)
                        cat(SEQ[i, (1 + (k - 1)*nbcol):(nbcol + (k - 1)*nbcol)], sep = colsep, file = zz)
                        cat("\n", file = zz)
                    }
                }
            }
        }
    }
    if (format == "sequential") {
        if (nb.block == 1) {
            for (i in 1:N) {
               cat(names(x)[i], file = zz)
               cat(SEQ[i, 1:nbcol], sep = colsep, file = zz)
               cat("\n", file = zz)
           }
        } else {
            for (i in 1:N) {
                cat(names(x)[i], file = zz)
                cat(SEQ[i, 1:nbcol], sep = colsep, file = zz)
                cat("\n", file = zz)
                for (k in 2:nb.block) {
                    if (k == nb.block) {
                        cat(indent, file = zz)
                        cat(SEQ[i, (1 + (k - 1)*nbcol):ncol(SEQ)], sep = colsep, file = zz)
                        cat("\n", file = zz)
                    } else {
                        cat(indent, file = zz)
                        cat(SEQ[i, (1 + (k - 1)*nbcol):(nbcol + (k - 1)*nbcol)], sep = colsep, file = zz)
                        cat("\n", file = zz)
                    }
                }
            }
        }
    }
    if (format == "fasta") {
        for (i in 1:N) {
            cat(">", names(x)[i], file = zz)
            cat("\n", file = zz)
            X <- paste(x[[i]], collapse= "")
            S <- length(x[[i]])
            if (nbcol < 0) {
                nb.block <- 1
                nbcol <- totalcol <- ceiling(S / colw)
            } else {
                totalcol <- ceiling(S / colw)
                nb.block <- ceiling(totalcol / nbcol)
            }
            SEQ <- character(totalcol)
            for (j in 1:totalcol) SEQ[j] <- substr(X, 1 + (j - 1)*colw, colw + (j - 1)*colw)
            for (k in 1:nb.block) {
                if (k == nb.block) {
                    cat(indent, file = zz)
                    cat(SEQ[(1 + (k - 1)*nbcol):length(SEQ)], sep = colsep, file = zz)
                    cat("\n", file = zz)
                } else {
                    cat(indent, file = zz)
                    cat(SEQ[(1 + (k - 1)*nbcol):(nbcol + (k - 1)*nbcol)], sep = colsep, file = zz)
                    cat("\n", file = zz)
                }
            }
        }
    }
    close(zz)
}
