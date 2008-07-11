## makeLabel.R (2008-07-03)

##   Label Management

## Copyright 2008 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

makeLabel <- function(x, ...) UseMethod("makeLabel")

makeLabel.character <- function(x, len = 99, space = "_",
          make.unique = TRUE, illegal = "():;,[]", quote = FALSE, ...)
{
    x <- gsub("[[:space:]]", space, x)
    if (illegal != "") {
        illegal <- unlist(strsplit(illegal, NULL))
        for (i in illegal) x <- gsub(i, "", x, fixed = TRUE)
    }
    if (quote) len <- len - 2
    nc <- nchar(x) > len
    if (any(nc)) x[nc] <- substr(x[nc], 1, len)
    tab <- table(x)
    if (all(tab == 1)) make.unique <- FALSE
    if (make.unique) {
        dup <- tab[which(tab > 1)]
        nms <- names(dup)
        for (i in 1:length(dup)) {
            j <- which(x == nms[i])
            end <- nchar(x[j][1])
            ## w: number of characters to be added as suffix
            w <- floor(log10(dup[i])) + 1
            suffix <- formatC(1:dup[i], width = w, flag = "0")
            if (end + w > len) {
                start <- end - w + 1
                substr(x[j], start, end) <- suffix
            } else x[j] <- paste(x[j], suffix, sep = "")
        }
    }
    if (quote) x <- paste('"', x, '"', sep = "")
    x
}

makeLabel.phylo <- function(x, tips = TRUE, nodes = TRUE, ...)
{
    if (tips)
        x$tip.label <- makeLabel.character(x$tip.label, ...)
    if (!is.null(x$node.label) && nodes)
        x$node.label <- makeLabel.character(x$node.label, ...)
    x
}

makeLabel.multiPhylo <- function(x, tips = TRUE, nodes = TRUE, ...)
{
    y <- attr(x, "TipLabel")
    if (is.null(y)) {
        for (i in 1:length(x))
            x[[i]] <- makeLabel.phylo(x[[i]], tips = tips, nodes = nodes, ...)
    } else {
        attr(x, "TipLabel") <- makeLabel.character(y, ...)
    }
    x
}

makeLabel.DNAbin <- function(x, ...)
{
    if (is.vector(x) || is.list(x))
        names(x) <- makeLabel.character(names(x), ...)
    else rownames(x) <- makeLabel.character(rownames(x), ...)
    x
}
