## clustal.R (2016-11-29)

##   Multiple Sequence Alignment with External Applications

## Copyright 2011-2016 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.errorAlignment <- function(exec, prog)
{
    dirs <- strsplit(Sys.getenv("PATH"), .Platform$path.sep)[[1]]
    paste0("\n   cannot find executable ", sQuote(exec), " on your computer.\n",
           "  It is recommended that you place the executable of ", prog, "\n",
           "  in a directory on the PATH of your computer which is:\n",
           paste(sort(dirs), collapse = "\n"))
}

clustalomega <- function(x, exec = NULL, MoreArgs = "", quiet = TRUE)
{
    os <- Sys.info()[1]
    if (is.null(exec)) {
        exec <- switch(os, "Linux" = "clustalo", "Darwin" = "clustalo",
                       "Windows" = "clustalo.exe")
    }

    if (missing(x)) {
        out <- system(paste(exec, "-h"))
        if (out == 127) stop(.errorAlignment(exec, "Clustal-Omega"))
        return(invisible(NULL))
    }

    x <- as.list(x)
    labels.bak <- names(x)
    names(x) <- paste0("Id", 1:length(x))

    d <- tempdir()
    inf <- paste(d, "input_clustalo.fas", sep = "/")
    outf <- paste(d, "output_clustalo.fas", sep = "/")
    write.dna(x, inf, "fasta")
    opts <- paste("-i", inf, "-o", outf, "--force")
    if (!quiet) opts <- paste(opts, "-v")
    opts <- paste(opts, MoreArgs)
    out <- system(paste(exec, opts), ignore.stdout = quiet)
    if (out == 127) stop(.errorAlignment(exec, "Clustal-Omega"))
    res <- read.dna(outf, "fasta")
    rownames(res) <- labels.bak
    res
}

clustal <- function(x, pw.gapopen = 10, pw.gapext = 0.1,
                    gapopen = 10, gapext = 0.2, exec = NULL,
                    MoreArgs = "", quiet = TRUE, original.ordering = TRUE)
{
    os <- Sys.info()[1]
    if (is.null(exec)) {
        exec <- switch(os, "Linux" = "clustalw", "Darwin" = "clustalw2",
                       "Windows" = "clustalw2.exe")
    }

    if (missing(x)) {
        out <- system(paste(exec, "-help"))
        if (out == 127) stop(.errorAlignment(exec, "Clustal"))
        return(invisible(NULL))
    }

    x <- as.list(x)
    labels.bak <- names(x)
    names(x) <- paste0("Id", 1:length(x))

    d <- tempdir()
    inf <- paste(d, "input_clustal.fas", sep = "/")
    outf <- paste(d, "input_clustal.aln", sep = "/")
    write.dna(x, inf, "fasta")
    prefix <- c("-INFILE", "-PWGAPOPEN", "-PWGAPEXT", "-GAPOPEN", "-GAPEXT")
    suffix <- c(inf, pw.gapopen, pw.gapext, gapopen, gapext)
    opts <- paste(prefix, suffix, sep = "=", collapse = " ")
    opts <- paste(opts, MoreArgs)
    out <- system(paste(exec, opts), ignore.stdout = quiet)
    if (out == 127) stop(.errorAlignment(exec, "Clustal"))
    res <- read.dna(outf, "clustal")
    if (original.ordering) res <- res[labels(x), ]
    rownames(res) <- labels.bak
    res
}

muscle <- function(x, exec = "muscle", MoreArgs = "", quiet = TRUE, original.ordering = TRUE)
{
    if (missing(x)) {
        out <- system(exec)
        if (out == 127) stop(.errorAlignment(exec, "MUSCLE"))
        return(invisible(NULL))
    }

    x <- as.list(x)
    labels.bak <- names(x)
    names(x) <- paste0("Id", 1:length(x))

    d <- tempdir()
    inf <- paste(d, "input_muscle.fas", sep = "/")
    outf <- paste(d, "output_muscle.fas", sep = "/")
    write.dna(x, inf, "fasta")
    opts <- paste("-in", inf, "-out", outf)
    if (quiet) opts <- paste(opts, "-quiet")
    opts <- paste(opts, MoreArgs)
    out <- system(paste(exec, opts))
    if (out == 127) stop(.errorAlignment(exec, "MUSCLE"))
    res <- read.dna(outf, "fasta")
    if (original.ordering) res <- res[labels(x), ]
    rownames(res) <- labels.bak
    res
}

tcoffee <- function(x, exec = "t_coffee", MoreArgs = "", quiet = TRUE, original.ordering = TRUE)
{
    if (missing(x)) {
        out <- system(exec)
        if (out == 127) stop(.errorAlignment(exec, "T-Coffee"))
        return(invisible(NULL))
    }

    x <- as.list(x)
    labels.bak <- names(x)
    names(x) <- paste0("Id", 1:length(x))

    d <- tempdir()
    od <- setwd(d)
    on.exit(setwd(od))
    inf <- "input_tcoffee.fas"
    write.dna(x, inf, "fasta")
    opts <- paste(inf, MoreArgs)
    if (quiet) opts <- paste(opts, "-quiet=nothing")
    out <- system(paste(exec, opts))
    if (out == 127) stop(.errorAlignment(exec, "T-Coffee"))
    res <- read.dna("input_tcoffee.aln", "clustal")
    if (original.ordering) res <- res[labels(x), ]
    rownames(res) <- labels.bak
    res
}
