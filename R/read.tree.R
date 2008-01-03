## read.tree.R (2007-12-22)

##   Read Tree Files in Parenthetic Format

## Copyright 2002-2007 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

tree.build <- function(tp)
{
    add.internal <- function() {
        edge[j, 1] <<- current.node
        edge[j, 2] <<- current.node <<- node <<- node + 1
        j <<- j + 1
    }
    add.terminal <- function() {
        edge[j, 1] <<- current.node
        edge[j, 2] <<- tip
        X <- unlist(strsplit(tpc[k], ":"))
        tip.label[tip] <<- X[1]
        edge.length[j] <<- as.numeric(X[2])
        k <<- k + 1
        tip <<- tip + 1
        j <<- j + 1
    }
    go.down <- function() {
        l <- which(edge[, 2] == current.node)
        X <- unlist(strsplit(tpc[k], ":"))
        node.label[current.node - nb.tip] <<- X[1]
        edge.length[l] <<- as.numeric(X[2])
        k <<- k + 1
        current.node <<- edge[l, 1]
    }
    if (!length(grep(",", tp))) {
        obj <- list(edge = matrix(c(2, 1), 1, 2))
        tp <- unlist(strsplit(tp, "[\\(\\):;]"))
        obj$edge.length <- as.numeric(tp[3])
        obj$Nnode <- 1
        obj$tip.label <- tp[2]
        if (length(tp) == 4) obj$node.label <- tp[4]
        class(obj) <- "phylo"
        return(obj)
    }
    tsp <- unlist(strsplit(tp, NULL))
    tpc <- unlist(strsplit(tp, "[\\(\\),;]"))
    tpc <- tpc[tpc != ""]
    skeleton <- tsp[tsp == "(" | tsp == ")" | tsp == "," | tsp == ";"]
    nsk <- length(skeleton)
    nb.node <- sum(skeleton == ")")
    nb.tip <- sum(skeleton == ",") + 1
    ## We will assume there is an edge at the root;
    ## if so, it will be removed and put into a vector
    nb.edge <- nb.node + nb.tip
    node.label <- character(nb.node)
    tip.label <- character(nb.tip)

    edge.length <- numeric(nb.edge)
    edge <- matrix(NA, nb.edge, 2)
    current.node <- node <- nb.tip + 1 # node number
    edge[nb.edge, 1] <- 0    # see comment above
    edge[nb.edge, 2] <- node #

    ## j: index of the line number of edge
    ## k: index of the line number of tpc
    ## tip: tip number
    j <- k <- tip <- 1

    for (i in 2:nsk) {
        if (skeleton[i] == "(") add.internal() # add an internal branch (on top)
        if (skeleton[i] == ",") {
            if (skeleton[i - 1] != ")") add.terminal() # add a terminal branch
        }
        if (skeleton[i] == ")") {
            if (skeleton[i - 1] == ",") { # add a terminal branch and go down one level
                add.terminal()
                go.down()
            }
            if (skeleton[i - 1] == ")") go.down() # go down one level
        }
    }
    if (is.na(node.label[1])) node.label[1] <- ""
    edge <- edge[-nb.edge, ]
    root.edge <- edge.length[nb.edge]
    edge.length <- edge.length[-nb.edge]
    obj <- list(edge = edge, edge.length = edge.length, Nnode = nb.node,
                tip.label = tip.label, node.label = node.label,
                root.edge = root.edge)
    if (all(obj$node.label == "")) obj$node.label <- NULL
    if (is.na(obj$root.edge)) obj$root.edge <- NULL
    if (all(is.na(obj$edge.length))) obj$edge.length <- NULL # added 2005-08-18
    class(obj) <- "phylo"
    obj
}

read.tree <- function(file = "", text = NULL, tree.names = NULL,
                      skip = 0, comment.char = "#", ...)
{
    if (!is.null(text)) {
        if (!is.character(text))
          stop("argument `text' must be of mode character")
        tree <- text
    } else {
        tree <- scan(file = file, what = character(), sep = "\n", quiet = TRUE,
                     skip = skip, comment.char = comment.char, ...)
    }
    ## Suggestion from Eric Durand and Nicolas Bortolussi (added 2005-08-17):
    if (identical(tree, character(0))) {
        warning("empty character string.")
        return(NULL)
    }
    tree <- gsub("[ \t]", "", tree)
    tsp <- unlist(strsplit(tree, NULL))
    ind <- which(tsp == ";")
    nb.tree <- length(ind)
    x <- c(1, ind[-nb.tree] + 1)
    y <- ind - 1
    ## Suggestion from Olivier François (added 2006-07-15):
    if (is.na(y[1])) return(NULL)
    else {
        STRING <- vector("list", nb.tree)
        for (i in 1:nb.tree)
          STRING[[i]] <- paste(tsp[x[i]:y[i]], sep = "", collapse = "")
    }
    obj <- vector("list", nb.tree)
    for (i in 1:nb.tree) {
        obj[[i]] <- if (length(grep(":", STRING[[i]]))) tree.build(STRING[[i]]) else clado.build(STRING[[i]])
        ## Check here that the root edge is not incorrectly represented
        ## in the object of class "phylo" by simply checking that there
        ## is a bifurcation at the root
        ROOT <- length(obj[[i]]$tip.label) + 1
        if(sum(obj[[i]]$edge[, 1] == ROOT) == 1 && dim(obj[[i]]$edge)[1] > 1) {
            stop(paste("There is apparently two root edges in your file: cannot read tree file.\n  Reading Newick file aborted at tree no.", i, sep = ""))
        }
    }
    if (nb.tree == 1) obj <- obj[[1]] else {
        if (is.null(tree.names))
          tree.names <- paste("tree", 1:nb.tree, sep = "")
        names(obj) <- tree.names
        class(obj) <- "multiPhylo"
    }
    obj
}
