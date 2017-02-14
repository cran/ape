## apetools.R (2017-02-03)

##   APE Tools

## Copyright 2005-2017 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.file.extensions <-
    list(clustal = "aln", fasta = c("fas", "fasta"),
         fastq = c("fq", "fastq"), newick = c("nwk", "newick", "tre", "tree"),
         nexus = c("nex", "nexus"), phylip = "phy")

Xplorefiles <- function(from = "HOME", recursive = TRUE, ignore.case = TRUE)
{
    if (from == "HOME") from <- Sys.getenv("HOME")
    FILES <- list.files(path = from, recursive = recursive, full.names = TRUE)
    ext <- if (exists(".file.extensions", envir = .PlotPhyloEnv))
               get(".file.extensions", envir = .PlotPhyloEnv)
           else .file.extensions
    res <- vector("list", length(ext))
    names(res) <- names(ext)
    for (i in seq_along(res)) {
        e <- paste0("\\.", ext[[i]], "$")
        if (length(e) > 1) e <- paste(e, collapse = "|")
        x <- grep(e, FILES, ignore.case = ignore.case, value = TRUE)
        res[[i]] <- data.frame(File = x, Size = file.size(x),
                               stringsAsFactors = FALSE)
    }
    res
}

editFileExtensions <- function()
{
    foo <- function(x) {
        n <- length(x)
        if (n < m) x[(n + 1):m] <- NA
        x
    }
    res <- if (exists(".file.extensions", envir = .PlotPhyloEnv))
               get(".file.extensions", envir = .PlotPhyloEnv)
           else .file.extensions
    m <- max(lengths(res, FALSE))
    res <- lapply(res, foo)
    res <- as.data.frame(res, stringsAsFactors = FALSE)
    res <- edit(res)
    res <- lapply(res, function(x) x[!is.na(x)])
    assign(".file.extensions", res, envir = .PlotPhyloEnv)
}

bydir <- function(x)
{
    nofile <- which(sapply(x, nrow) == 0)
    if (length(nofile)) x <- x[-nofile]
    if (!length(x)) {
        cat("No file\n")
        return(invisible(NULL))
    }
    for (i in seq_along(x)) x[[i]]$Type <- names(x)[i]
    x <- do.call(rbind, x)
    x <- x[order(x$File), ]
    SPLIT <- strsplit(x$File, "/")
    LL <- lengths(SPLIT)
    foo <- function(i, PATH) {
        K <- grep(paste0("^", PATH, "/"), x$File)
        sel <- intersect(K, which(LL == i + 1L))
        if (length(sel)) {
            y <- x[sel, ]
            y$File <- gsub(".*/", "", y$File)
            cat("\n", PATH, "/\n", sep = "")
            print(y, row.names = FALSE)
        }
        if (length(sel) < length(K)) {
            d <- setdiff(K, sel)
            subdir <- unlist(lapply(SPLIT[d], "[", i + 1L))
            for (z in unique(subdir))
                foo(i + 1L, paste(PATH, z, sep = "/"))
        }
    }
    top <- unlist(lapply(SPLIT, "[", 1L))
    for (z in unique(top)) foo(1L, z)
}
