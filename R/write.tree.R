## write.tree.R (2009-03-23)

##   Write Tree File in Parenthetic Format

## Copyright 2002-2009 Emmanuel Paradis and Daniel Lawson

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

checkLabel <- function(x, ...)
{
    ## delete all leading and trailing spaces and tabs, and
    ## the leading left and trailing right parentheses:
    ## (the syntax will work with any mix of these characters,
    ##  e.g., "    ( ( ((  " will correctly be deleted)
    x <- gsub("^[[:space:]\\(]+", "", x)
    x <- gsub("[[:space:]\\)]+$", "", x)
    ## replace all spaces and tabs by underscores:
    x <- gsub("[[:space:]]", "_", x)
    ## remove all commas, colons, and semicolons
    x <- gsub("[,:;]", "", x)
    ## replace left and right parentheses with dashes:
    x <- gsub("[\\(\\)]", "-", x)
    ## delete extra underscores and extra dashes:
    x <- gsub("_{2,}", "_", x)
    x <- gsub("-{2,}", "-", x)
    x
}

write.tree <-
    function (phy, file = "", append = FALSE, digits = 10, tree.names = FALSE)
{
    if (is.logical(tree.names)) {
        output.tree.names <- tree.names
        tree.names <- NULL
    } else if (is.character(tree.names)) {
        output.tree.names <- TRUE
        names(tree) <- tree.names
    }
    if (output.tree.names)
        names(tree) <- checkLabel(names(tree))
    if (class(phy) == "multiPhylo") {
        write.tree(phy[[1]], file = file, append = append,
                   digits = digits, tree.names = names[1])
        if (length(phy) > 1)
            for (i in 2:length(phy)) write.tree(phy[[i]], file = file,
                append = TRUE, digits = digits, tree.names = names(phy)[i])
        return(invisible(NULL))
    }
    if (class(phy) != "phylo")
        stop("object \"phy\" is not of class \"phylo\"")
    brl <- !is.null(phy$edge.length)
    nodelab <- !is.null(phy$node.label)
    phy$tip.label <- checkLabel(phy$tip.label)
    if (nodelab) phy$node.label <- checkLabel(phy$node.label)
    f.d <- paste("%.", digits, "g", sep = "")
    cp <- function(s) STRING <<- paste(STRING, s, sep = "")
    add.internal <- function(i) {
        cp("(")
        br <- which(phy$edge[, 1] == i)
        for (j in br) {
            desc <- phy$edge[j, 2]
            if (desc > n)
                add.internal(desc)
            else add.terminal(j)
            if (j != br[length(br)])
                cp(",")
        }
        cp(")")
        if (nodelab)
            cp(phy$node.label[i - n])
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[which(phy$edge[,
                2] == i)]))
        }
    }
    add.terminal <- function(i) {
        cp(phy$tip.label[phy$edge[i, 2]])
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[i]))
        }
    }
    n <- length(phy$tip.label)
    STRING <-
        if (output.tree.names) paste(tree.names, "(", sep = "") else "("
    br <- which(phy$edge[, 1] == n + 1)
    for (j in br) {
        desc <- phy$edge[j, 2]
        if (desc > n)
            add.internal(desc)
        else add.terminal(j)
        if (j != br[length(br)])
            cp(",")
    }
    if (is.null(phy$root.edge)) {
        cp(")")
        if (nodelab)
            cp(phy$node.label[1])
        cp(";")
    }
    else {
        cp(")")
        if (nodelab)
            cp(phy$node.label[1])
        cp(":")
        cp(sprintf(f.d, phy$root.edge))
        cp(";")
    }
    if (file == "")
        return(STRING)
    else cat(STRING, file = file, append = append, sep = "\n")
}

