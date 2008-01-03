## write.tree.R (2007-12-22)

##   Write Tree File in Parenthetic Format

## Copyright 2002-2007 Emmanuel Paradis

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

write.tree <- function(phy, file = "", append = FALSE,
                       digits = 10)
{
    if (class(phy) == "multiPhylo") {
        write.tree(phy[[1]], file = file,
                   append = append, digits = digits)
        if (length(phy) > 1)
            for (i in 2:length(phy))
                write.tree(phy[[i]], file = file,
                           append = TRUE, digits = digits)
        return(invisible(NULL))
    }

    if (class(phy) != "phylo")
      stop('object "phy" is not of class "phylo"')

    brl <- !is.null(phy$edge.length)

### Ne serait-il pas plus efficace de créer des node labels vides
### "" et d'éviter l'évaluation if (nodelab) ????
### Autre possibilité : créer plusieurs variants de ces fonctions
### (au moins deux avec/sans edge.length)

### Encore autre chose: les appels à which ne peuvent-ils pas
### être évités ??? surtout si l'arbre est en cladewise order...

    nodelab <- !is.null(phy$node.label)
    phy$tip.label <- checkLabel(phy$tip.label)
    if (nodelab)
      phy$node.label <- checkLabel(phy$node.label)

    f.d <- paste("%.", digits, "g", sep = "")

    cp <- function(s) STRING <<- paste(STRING, s, sep = "")
    add.internal <- function(i) {
        cp("(")
        br <- which(phy$edge[, 1] == i)
        for (j in br) {
            desc <- phy$edge[j, 2]
            if (desc > n) add.internal(desc) else add.terminal(j)
            if (j != br[length(br)]) cp(",")
        }
        cp(")")
        if (nodelab) cp(phy$node.label[i - n])
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[which(phy$edge[, 2] == i)]))
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
    STRING <- "("
    br <- which(phy$edge[, 1] == n + 1)
    for (j in br) {
        desc <- phy$edge[j, 2]
        if (desc > n) add.internal(desc) else add.terminal(j)
        if (j != br[length(br)]) cp(",")
    }
    if (is.null(phy$root.edge)) {
        cp(")")
        if (nodelab) cp(phy$node.label[1])
        cp(";")
    } else {
        cp(")")
        if (nodelab) cp(phy$node.label[1])
        cp(":")
        cp(sprintf(f.d, phy$root.edge))
        cp(";")
    }
    if (file == "") return(STRING)
    else cat(STRING, file = file, append = append, sep = "\n")
}
