### read.dna.R (2006-12-29)
###
###     Read DNA Sequences in a File
###
### Copyright 2003-2006 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

read.dna <- function(file, format = "interleaved", skip = 0,
                     nlines = 0, comment.char = "#", seq.names = NULL)
{
    getTaxaNames <- function(x) {
        x <- sub("^ +", "", x) # remove the leading spaces
        x <- sub(" +$", "", x) # remove the trailing spaces
        x <- sub("^['\"]", "", x) # remove the leading quotes
        x <- sub("['\"]$", "", x) # remove the trailing quotes
        x
    }
    format <- match.arg(format, c("interleaved", "sequential", "fasta"))
    phylip <- if (format %in% c("interleaved", "sequential")) TRUE else FALSE
    X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE,
              skip = skip, nlines = nlines, comment.char = comment.char)
    if (phylip) {
        fl <- X[1]
        oop <- options(warn = -1)
        ## need to remove the possible leading spaces in the first line
        fl.num <- as.numeric(unlist(strsplit(gsub("^ +", "", fl), " +")))
        options(oop)
        if (all(is.na(fl.num)))
          stop("the first line of the file must contain the dimensions of the data")
        if (length(fl.num) != 2)
          stop("the first line of the file must contain TWO numbers")
        else {
            n <- fl.num[1]
            s <- fl.num[2]
        }
        X <- X[-1]
        obj <- vector("character", n*s)
        dim(obj) <- c(n, s)
    }
    if (format == "interleaved") {
        fl <- X[1]
        fl <- unlist(strsplit(fl, NULL))
        bases <- grep("[-AaCcGgTtUuMmRrWwSsYyKkVvHhDdBbNn]", fl)
        z <- diff(bases)
        for (i in 1:length(z)) if (all(z[i:(i + 8)] == 1)) break
        start.seq <- bases[i]
        if (is.null(seq.names))
          seq.names <- getTaxaNames(substr(X[1:n], 1, start.seq - 1))
        X[1:n] <- substr(X[1:n], start.seq, nchar(X[1:n]))
        X <- gsub(" ", "", X)
        nl <- length(X)
        for (i in 1:n)
          obj[i, ] <- unlist(strsplit(X[seq(i, nl, n)], NULL))
    }
    if (format == "sequential") {
        fl <- X[1]
        taxa <- character(n)
        j <- 1
        for (i in 1:n) {
            bases <- grep("[-AaCcGgTtUuMmRrWwSsYyKkVvHhDdBbNn]",
                          unlist(strsplit(X[j], NULL)))
            z <- diff(bases)
            for (k in 1:length(z)) if (all(z[k:(k + 8)] == 1)) break
            start.seq <- bases[k]
            taxa[i] <- substr(X[j], 1, start.seq - 1)
            sequ <- substr(X[j], start.seq, nchar(X[j]))
            sequ <- gsub(" ", "", sequ)
            j <- j + 1
            while (nchar(sequ) < s) {
                sequ <- paste(sequ, gsub(" " , "", X[j]), sep = "")
                j <- j + 1
            }
            obj[i, ] <- unlist(strsplit(sequ, NULL))
        }
        if (is.null(seq.names)) seq.names <- getTaxaNames(taxa)
    }
    if (format == "fasta") {
        start <- grep("^ {0,}>", X)
        taxa <- X[start]
        n <- length(taxa)
        obj <- vector("list", n)
        if (is.null(seq.names)) {
            taxa <- sub("^ {0,}>", "", taxa) # remove the hook and the spaces before
            seq.names <- getTaxaNames(taxa)
        }
        start <- c(start, length(X) + 1) # this avoids the following to crash when `i = n'
        for (i in 1:n)
          obj[[i]] <- unlist(strsplit(gsub(" ", "",
                                           X[(start[i] + 1):(start[i + 1] - 1)]),
                                      NULL))
    }
    if (phylip) {
        rownames(obj) <- seq.names
        obj <- tolower(obj)
    } else {
        names(obj) <- seq.names
        obj <- lapply(obj, tolower)
    }
    obj
}
