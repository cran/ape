### read.nexus.R (2006-11-24)
###
###     Read Tree File in Nexus Format
###
### Copyright 2003-2006 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

clado.build <- function(tp) {
    add.internal <- function() {
        edge[j, 1] <<- current.node
        node <<- node + 1
        edge[j, 2] <<- current.node <<- node
        j <<- j + 1
    }
    add.terminal <- function() {
        edge[j, 1] <<- current.node
        edge[j, 2] <<- tip
        tip.label[tip] <<- tpc[k]
        k <<- k + 1
        tip <<- tip + 1
        j <<- j + 1
    }
    go.down <- function() {
        l <- which(edge[, 2] == current.node)
        node.label[current.node - nb.tip] <<- tpc[k]
        k <<- k + 1
        current.node <<- edge[l, 1]
    }
    if (!length(grep(",", tp))) {
        obj <- list(edge = matrix(c(2, 1), 1, 2), Nnode = 1)
        tp <- unlist(strsplit(tp, "[\\(\\);]"))
        obj$tip.label <- tp[2]
        if (length(tp) == 3) obj$node.label <- tp[3]
        class(obj) <- "phylo"
        return(obj)
    }
    tsp <- unlist(strsplit(tp, NULL))
    tp <- gsub(")", ")NA", tp)
    tp <- gsub(" ", "", tp)
    tpc <- unlist(strsplit(tp, "[\\(\\),;]"))
    tpc <- tpc[tpc != ""]
    skeleton <- tsp[tsp == "(" | tsp == ")" | tsp == "," | tsp == ";"]
    nsk <- length(skeleton)
    nb.node <- length(skeleton[skeleton == ")"])
    nb.tip <- length(skeleton[skeleton == ","]) + 1
    ## We will assume there is an edge at the root;
    ## if so, it will be removed and put in a vector
    nb.edge <- nb.node + nb.tip
    node.label <- character(nb.node)
    tip.label <- character(nb.tip)

    edge <- matrix(NA, nb.edge, 2)
    current.node <- node <- nb.tip + 1 # node number
    edge[nb.edge, 1] <- 0    # see comment above
    edge[nb.edge, 2] <- node #

    ## j: index of the line number of edge
    ## k: index of the line number of tpc
    ## tip: tip number
    j <- k <- tip <- 1

    for (i in 2:nsk) {
        if (skeleton[i] == "(") add.internal()      # add an internal branch (on top)
        if (skeleton[i] == ",") {
            if (skeleton[i - 1] != ")") add.terminal()   # add a terminal branch
        }
        if (skeleton[i] == ")") {
            if (skeleton[i - 1] == ",") {   # add a terminal branch and go down one level
                add.terminal()
                go.down()
            }
            if (skeleton[i - 1] == ")") go.down()   # go down one level
        }
    }
#    if(node.label[1] == "NA") node.label[1] <- ""
    edge <- edge[-nb.edge, ]
    obj <- list(edge = edge, tip.label = tip.label,
                Nnode = nb.node, node.label = node.label)
    obj$node.label <- if (all(obj$node.label == "NA")) NULL else gsub("^NA", "", obj$node.label)
    class(obj) <- "phylo"
    return(obj)
}

read.nexus <- function(file, tree.names = NULL)
{
    X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE)
    ## first remove all the comments
    LEFT <- grep("\\[", X)
    RIGHT <- grep("\\]", X)
    if (length(LEFT)) {
        for (i in length(LEFT):1) {
            if (LEFT[i] == RIGHT[i]) {
                X[LEFT[i]] <- gsub("\\[.*\\]", "", X[LEFT[i]])
            } else {
                X[LEFT[i]] <- gsub("\\[.*", "", X[LEFT[i]])
                X[RIGHT[i]] <- gsub(".*\\]", "", X[RIGHT[i]])
                if (LEFT[i] < RIGHT[i] - 1) X <- X[-((LEFT[i] + 1):(RIGHT[i] - 1))]
            }
        }
    }
    X <- gsub("ENDBLOCK;", "END;", X, ignore.case = TRUE)
    endblock <- grep("END;", X, ignore.case = TRUE)
    semico <- grep(";", X)
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    translation <- FALSE
    if (length(i2) == 1) if (i2 > i1) translation <- TRUE
    if (translation) {
        end <- semico[semico > i2][1]
        x <- paste(X[i2:end], sep = "", collapse = "")
        x <- gsub("TRANSLATE", "", x, ignore.case = TRUE)
        x <- unlist(strsplit(x, "[,; \t]"))
        x <- x[x != ""]
        TRANS <- matrix(x, ncol = 2, byrow = TRUE)
        TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
    }
    start <- if (translation)  semico[semico > i2][1] + 1 else semico[semico > i1][1]
    end <- endblock[endblock > i1][1] - 1
    tree <- paste(X[start:end], sep = "", collapse = "")
    tree <- gsub(" ", "", tree)
    tree <- unlist(strsplit(tree, "[=;]"))
    tree <- tree[grep("[\\(\\)]", tree)]
    nb.tree <- length(tree)
    STRING <- as.list(tree)
    trees <- list()
    for (i in 1:nb.tree) {
        obj <- if (length(grep(":", STRING[[i]]))) tree.build(STRING[[i]]) else clado.build(STRING[[i]])
        if (translation) {
            for (j in 1:length(obj$tip.label)) {
                ind <- which(obj$tip.label[j] == TRANS[, 1])
                obj$tip.label[j] <- TRANS[ind, 2]
            }
            if (!is.null(obj$node.label)) {
                for (j in 1:length(obj$node.label)) {
                    ind <- which(obj$node.label[j] == TRANS[, 1])
                    obj$node.label[j] <- TRANS[ind, 2]
                }
            }
        }
        trees[[i]] <- obj
        ## Check here that the root edge is not incorrectly represented
        ## in the object of class "phylo" by simply checking that there
        ## is a bifurcation at the root (node "-1")
        if (sum(trees[[i]]$edge[, 1] == "-1") == 1 && dim(trees[[i]]$edge)[1] > 1) {
            warning("The root edge is apparently not correctly represented\nin your tree: this may be due to an extra pair of\nparentheses in your file; the returned object has been\ncorrected but your file may not be in a valid Newick\nformat")
            ind <- which(trees[[i]]$edge[, 1] == "-1")
            trees[[i]]$root.edge <- trees[[i]]$edge.length[ind]
            trees[[i]]$edge.length <- trees[[i]]$edge.length[-ind]
            trees[[i]]$edge <- trees[[i]]$edge[-ind, ]
            for (j in 1:length(trees[[i]]$edge))
              if (as.numeric(trees[[i]]$edge[j]) < 0)
                trees[[i]]$edge[j] <- as.character(as.numeric(trees[[i]]$edge[j]) + 1)
            ## Check a second time and if there is still a problem...!!!
            if(sum(trees[[i]]$edge[, 1] == "-1") == 1)
              stop("There is apparently two root edges in your file: cannot read tree file")
        }
    }
    if (nb.tree == 1) trees <- trees[[1]] else {
        names(trees) <- if (is.null(tree.names)) paste("tree", 1:nb.tree, sep = "") else tree.names
        class(trees) <- c("phylo", "multi.tree")
    }
    if (length(grep("[\\/]", file)) == 1) attr(trees, "origin") <- file
    else attr(trees, "origin") <- paste(getwd(), file, sep = "/")
    trees
}
