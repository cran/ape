### read.GenBank.R (2007-06-27)
###
###    Read DNA Sequences from GenBank via Internet
###
### Copyright 2002-2007 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

read.GenBank <- function(access.nb, seq.names = access.nb,
                         species.names = TRUE, as.character = FALSE)
{
    N <- length(access.nb)
    ## If there are more than 400 sequences, we need to break down the
    ## requests, otherwise there is a segmentation fault.
    nrequest <- N %/% 400 + as.logical(N %% 400)
    X <- character(0)
    for (i in 1:nrequest) {
        a <- (i - 1) * 400 + 1
        b <- 400 * i
        if (i == nrequest) b <- N
        URL <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                     paste(access.nb[a:b], collapse = ","),
                     "&rettype=gb", sep = "")
        X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
    }
    FI <- grep("^ {0,}ORIGIN", X) + 1
    LA <- which(X == "//") - 1
    obj <- list()
    length(obj) <- N
    for (i in 1:N) {
        ## remove all spaces and digits
        tmp <- gsub("[[:digit:] ]", "", X[FI[i]:LA[i]])
        obj[[i]] <- unlist(strsplit(tmp, NULL))
    }
    names(obj) <- seq.names
    if (!as.character) obj <- as.DNAbin(obj)
    if (species.names) {
        tmp <- character(N)
        sp <- grep("ORGANISM", X)
        for (i in 1:N)
          tmp[i] <- unlist(strsplit(X[sp[i]], " +ORGANISM +"))[2]
        attr(obj, "species") <- gsub(" ", "_", tmp)
    }
    obj
}
