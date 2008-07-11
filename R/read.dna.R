## read.dna.R (2008-07-03)

##   Read DNA Sequences in a File

## Copyright 2003-2008 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

read.dna <- function(file, format = "interleaved", skip = 0,
                     nlines = 0, comment.char = "#", seq.names = NULL,
                     as.character = FALSE)
{
    getTaxaNames <- function(x) {
        x <- sub("^['\" ]+", "", x) # remove the leading quotes and spaces
        x <- sub("['\" ]+$", "", x) #   "     "  trailing  "     "    "
        x
    }
    getNucleotide <- function(x) {
        x <- gsub(" ", "", x)
        x <- strsplit(x, NULL)
        unlist(x)
    }
    formats <- c("interleaved", "sequential", "fasta", "clustal")
    format <- match.arg(format, formats)
    phylip <- if (format %in% formats[1:2]) TRUE else FALSE
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE,
              skip = skip, nlines = nlines, comment.char = comment.char)
    pat.base <- "[-AaCcGgTtUuMmRrWwSsYyKkVvHhDdBbNn?]{10}"
    if (phylip) {
        ## need to remove the possible leading spaces in the first line
        fl <- gsub("^ +", "", X[1])
        fl <- as.numeric(unlist(strsplit(fl, " +")))
        if (length(fl) != 2 || any(is.na(fl)))
            stop("the first line of the file must contain the dimensions of the data")
        n <- fl[1]
        s <- fl[2]
        obj <- matrix("", n, s)
        X <- X[-1]
    }
    if (format == "interleaved") {
        start.seq <- regexpr(pat.base, X[1])[1]
        if (is.null(seq.names))
            seq.names <- getTaxaNames(substr(X[1:n], 1, start.seq - 1))
        X[1:n] <- substr(X[1:n], start.seq, nchar(X[1:n]))
        nl <- length(X)
        for (i in 1:n)
            obj[i, ] <- getNucleotide(X[seq(i, nl, n)])
    }
    if (format == "sequential") {
        taxa <- character(n)
        j <- 1
        for (i in 1:n) {
            start.seq <- regexpr(pat.base, X[j])[1]
            taxa[i] <- substr(X[j], 1, start.seq - 1)
            sequ <- getNucleotide(substr(X[j], start.seq, nchar(X[j])))
            j <- j + 1
            while (length(sequ) < s) {
                sequ <- c(sequ, getNucleotide(X[j]))
                j <- j + 1
            }
            obj[i, ] <- sequ
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
            obj[[i]] <- getNucleotide(X[(start[i] + 1):(start[i + 1] - 1)])
    }
    if (format == "clustal") {
        X <- X[-1]
        ## find where the 1st sequence starts
        start.seq <- regexpr(pat.base, X[1])[1]
        ## find the lines with *********....
        nspaces <- paste("^ {", start.seq - 1, "}", sep = "", collapse = "")
        stars <- grep(nspaces, X)
        ## we now know how many sequences in the file:
        n <- stars[1] - 1
        ## get the sequence names in the same way than "interleaved":
        if (is.null(seq.names))
            seq.names <- getTaxaNames(substr(X[1:n], 1, start.seq - 1))
        ## need to remove the sequence names before getting the sequences:
        X <- substr(X, start.seq, nchar(X))
        nl <- length(X)
        ## find the length of the 1st sequence:
        tmp <- getNucleotide(X[seq(1, nl, n + 1)])
        s <- length(tmp)
        obj <- matrix("", n, s)
        obj[1, ] <- tmp
        for (i in 2:n)
            obj[i, ] <- getNucleotide(X[seq(i, nl, n + 1)])
    }
    if (format != "fasta") {
        rownames(obj) <- seq.names
        obj <- tolower(obj)
    } else {
        names(obj) <- seq.names
        obj <- lapply(obj, tolower)
    }
    if (!as.character) obj <- as.DNAbin(obj)
    obj
}
