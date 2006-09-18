### write.nexus.R (2006-09-09)
###
###          Write Tree File in Nexus Format
###
### Copyright 2003-2006 Emmanuel Paradis
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

write.nexus <- function(..., file = "", translate = TRUE, original.data = TRUE)
{
    obj <- list(...)
    ## We insure that all trees are in a list, even if there is a single one:
    if (length(obj) == 1) {
        if (class(obj[[1]]) == "phylo") ntree <- 1
        else {
            obj <- unlist(obj, recursive = FALSE)
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
        ## We take arbitrarily the labels of the first tree, and
        ## translate them as "1", "2", "3", ...
        cat("\tTRANSLATE\n", file = file, append = TRUE)
        tmp <- checkLabel(obj[[1]]$tip.label)
        X <- paste("\t\t", 1:N, "\t", tmp, ",", sep = "")
        ## We remove the last comma:
        X[length(X)] <- gsub(",", "", X[length(X)])
        cat(X, file = file, append = TRUE, sep = "\n")
        cat("\t;\n", file = file, append = TRUE)
        token <- as.character(1:N)
        names(token) <- obj[[1]]$tip.label
        obj[[1]]$tip.label <- token
        if (ntree > 1)
          for (i in 2:ntree)
            obj[[i]]$tip.label <- token[obj[[i]]$tip.label]
    } else {
        for (i in 1:ntree)
          obj[[i]]$tip.label <- checkLabel(obj[[i]]$tip.label)
    }
    for (i in 1:ntree) {
        if (class(obj[[i]]) != "phylo") next
        if (is.rooted(obj[[i]]))
          cat("\tTREE * UNTITLED = [&R] ", file = file, append = TRUE)
        else cat("\tTREE * UNTITLED = [&U] ", file = file, append = TRUE)
        cat(write.tree(obj[[i]], file = "", multi.line = FALSE),
            "\n", sep = "",file = file, append = TRUE)
    }
    cat("END;\n", file = file, append = TRUE)
    if(original.data) cat(ORI, file = file, append = TRUE, sep = "\n")
}
