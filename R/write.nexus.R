## write.nexus.R (2011-03-26)

##   Write Tree File in Nexus Format

## Copyright 2003-2011 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

write.nexus <- function(..., file = "", translate = TRUE, original.data = TRUE)
{
    obj <- list(...)
    ## We insure that all trees are in a list, even if there is a single one:
    if (length(obj) == 1) {
        if (class(obj[[1]]) == "phylo") ntree <- 1
        else {
            obj <- obj[[1]] # NOT use unlist()
            ntree <- length(obj)
        }
    } else ntree <- length(obj)
    cat("#NEXUS\n", file = file)
    cat(paste("[R-package APE, ", date(), "]\n\n", sep = ""),
        file = file, append = TRUE)
    if (original.data) {
        if (!is.null(attr(obj[[1]], "origin"))) {
            if (!file.exists(attr(obj[[1]], "origin"))) {
                warning(paste("the file", attr(obj[[1]], "origin"),
                              "cannot be found,
the original data won't be written with the tree."))
                original.data <- FALSE
            }
            else {
                ORI <- scan(file = attr(obj[[1]], "origin"), what = character(),
                            sep = "\n", skip = 1)
                start <- grep("BEGIN TAXA;", ORI)
                ORI <- ORI[-(1:(start - 1))]
                ORI <- gsub("ENDBLOCK;", "END;", ORI)
                endblock <- grep("END;", ORI)
                start <- grep("BEGIN TREES;", ORI)
                end <- endblock[endblock > start][1]
                cat(ORI[1:(start - 1)], file = file, append = TRUE, sep = "\n")
                ORI <- ORI[-(1:end)]
            }
        }
        else original.data <- FALSE
    }
    N <- length(obj[[1]]$tip.label)
    if (!original.data) {
        cat("BEGIN TAXA;\n", file = file, append = TRUE)
        cat(paste("\tDIMENSIONS NTAX = ", N, ";\n", sep = ""),
            file = file, append = TRUE)
        cat("\tTAXLABELS\n", file = file, append = TRUE)
        cat(paste("\t\t", obj[[1]]$tip.label, sep = ""),
            sep = "\n", file = file, append = TRUE)
        cat("\t;\n", file = file, append = TRUE)
        cat("END;\n", file = file, append = TRUE)
    }
    cat("BEGIN TREES;\n", file = file, append = TRUE)
    if (translate) {
        cat("\tTRANSLATE\n", file = file, append = TRUE)
        obj <- .compressTipLabel(obj)
        X <- paste("\t\t", 1:N, "\t", attr(obj, "TipLabel"), ",", sep = "")
        ## We remove the last comma:
        X[length(X)] <- gsub(",", "", X[length(X)])
        cat(X, file = file, append = TRUE, sep = "\n")
        cat("\t;\n", file = file, append = TRUE)
        class(obj) <- NULL
        for (i in 1:ntree)
            obj[[i]]$tip.label <- as.character(1:N)
    } else {
        if (is.null(attr(obj, "TipLabel"))) {
            for (i in 1:ntree)
                obj[[i]]$tip.label <- checkLabel(obj[[i]]$tip.label)
        } else {
            attr(obj, "TipLabel") <- checkLabel(attr(obj, "TipLabel"))
            obj <- .uncompressTipLabel(obj)
        }
    }

    title <- names(obj)
    if (is.null(title))
        title <- rep("UNTITLED", ntree)
    else {
        if (any(s <- title == "")) title[s] <- "UNTITLED"
    }

    for (i in 1:ntree) {
        if (class(obj[[i]]) != "phylo") next
        if (is.rooted(obj[[i]]))
          cat("\tTREE *,", title[i], "= [&R] ", file = file, append = TRUE)
        else cat("\tTREE *", title[i], "= [&U] ", file = file, append = TRUE)
        cat(write.tree(obj[[i]], file = ""),
            "\n", sep = "", file = file, append = TRUE)
    }
    cat("END;\n", file = file, append = TRUE)
    if(original.data) cat(ORI, file = file, append = TRUE, sep = "\n")
}
