## summary.phylo.R (2008-04-22)

##   Print Summary of a Phylogeny

## Copyright 2003-2008 Emmanuel Paradis, and 2006 Ben Bolker

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

Ntip <- function(phy)
{
    if (class(phy) != "phylo")
      stop('object "phy" is not of class "phylo"')
    length(phy$tip.label)
}

Nnode <- function(phy, internal.only = TRUE)
{
    if (class(phy) != "phylo")
      stop('object "phy" is not of class "phylo"')
    if (internal.only) return(phy$Nnode)
    phy$Nnode + length(phy$tip.label)
}

Nedge <- function(phy)
{
    if (class(phy) != "phylo")
      stop('object "phy" is not of class "phylo"')
    dim(phy$edge)[1]
}

summary.phylo <- function(object, ...)
{
    cat("\nPhylogenetic tree:", deparse(substitute(object)), "\n\n")
    nb.tip <- length(object$tip.label)
    nb.node <- object$Nnode
    cat("  Number of tips:", nb.tip, "\n")
    cat("  Number of nodes:", nb.node, "\n")
    if (is.null(object$edge.length))
      cat("  No branch lengths.\n")
    else {
        cat("  Branch lengths:\n")
        cat("    mean:", mean(object$edge.length), "\n")
        cat("    variance:", var(object$edge.length), "\n")
        cat("    distribution summary:\n")
        print(summary(object$edge.length)[-4])
    }
    if (is.null(object$root.edge))
      cat("  No root edge.\n")
    else
      cat("  Root edge:", object$root.edge, "\n")
    if (nb.tip <= 10) {
        cat("  Tip labels:", object$tip.label[1], "\n")
        cat(paste("             ", object$tip.label[-1]), sep = "\n")
    }
    else {
        cat("  First ten tip labels:", object$tip.label[1], "\n")
        cat(paste("                       ", object$tip.label[2:10]), sep = "\n")
    }
    if (is.null(object$node.label))
      cat("  No node labels.\n")
    else {
        if (nb.node <= 10) {
            cat("  Node labels:", object$node.label[1], "\n")
            cat(paste("              ", object$node.label[-1]), sep = "\n")
        }
        else {
            cat("  First ten node labels:", object$node.label[1], "\n")
            cat(paste("                        ", object$node.label[2:10]), sep = "\n")

        }
    }
    if (!is.null(attr(object, "loglik"))) {
        cat("Phylogeny estimated by maximum likelihood.\n")
        cat("  log-likelihood:", attr(object, "loglik"), "\n\n")
        npart <- length(attr(object, "para"))
        for (i in 1:npart) {
            cat("partition ", i, ":\n", sep = "")
            print(attr(object, "para")[[i]])
            if (i == 1) next
            else cat("  contrast parameter (xi):",
                     attr(object, "xi")[i - 1], "\n")
        }
    }
}

### by BB:
print.phylo <- function(x, printlen = 6,...)
{
    nb.tip <- length(x$tip.label)
    nb.node <- x$Nnode
    cat(paste("\nPhylogenetic tree with", nb.tip, "tips and", nb.node,
              "internal nodes.\n\n"))
    cat("Tip labels:\n")
    if (nb.tip > printlen) {
        cat(paste("\t", paste(x$tip.label[1:printlen],
                              collapse=", "), ", ...\n", sep = ""))
    } else print(x$tip.label)
    if (!is.null(x$node.label)) {
        cat("\tNode labels:\n")
        if (nb.node > printlen) {
            cat(paste("\t", paste(x$node.label[1:printlen],
                                 collapse=", "), ",...\n", sep = ""))
        } else print(x$node.label)
    }
    rlab <- if (is.rooted(x)) "Rooted" else "Unrooted"
    cat("\n", rlab, "; ", sep="")

    blen <- if (is.null(x$edge.length)) "no branch lengths." else
    "includes branch lengths."
    cat(blen, "\n", sep = "")
}

print.multiPhylo <- function(x, details = FALSE, ...)
{
    N <- length(x)
    cat(N, "phylogenetic trees\n")
    if (details)
      for (i in 1:N)
        cat("tree", i, ":", length(x[[i]]$tip.label), "tips\n")
    cat("\n")
}

"[[.multiPhylo" <- function(x, i)
{
    class(x) <- NULL
    phy <- x[[i]]
    if (!is.null(attr(x, "TipLabel")))
        phy$tip.label <- attr(x, "TipLabel")
    phy
}

`$.multiPhylo` <- function(x, name) x[[name]]

"[.multiPhylo" <- function(x, i)
{
    class(x) <- NULL
    structure(x[i], TipLabel = attr(x, "TipLabel"),
              class = "multiPhylo")
}

str.multiPhylo <- function(object, ...)
{
    class(object) <- NULL
    cat('Class "multiPhylo"\n')
    str(object, ...)
}
