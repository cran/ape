## read.nexus.R (2011-03-26)

##   Read Tree File in Nexus Format

## Copyright 2003-2011 Emmanuel Paradis and 2010 Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.treeBuildWithTokens <- function(x)
{
    phy <- .Call("treeBuildWithTokens", x, PACKAGE = "ape")
    dim(phy[[1]]) <- c(length(phy[[1]])/2, 2)
    nms <- c("edge", "edge.length", "Nnode", "node.label", "root.edge")
    if (length(phy) == 4) nms <- nms[-5]
    names(phy) <- nms
    if (all(phy$node.label == "")) phy$node.label <- NULL
    class(phy) <- "phylo"
    phy
}

clado.build <- function(tp)
{
    add.internal <- function() {
        edge[j, 1] <<- current.node
        node <<- node + 1
        edge[j, 2] <<- current.node <<- node
        index[node] <<- j # set index
        j <<- j + 1
    }
    add.terminal <- function() {
        edge[j, 1] <<- current.node
        edge[j, 2] <<- tip
        index[tip] <<- j # set index
        tip.label[tip] <<- tpc[k]
        k <<- k + 1
        tip <<- tip + 1
        j <<- j + 1
    }
    go.down <- function() {
        l <- index[current.node]
        node.label[current.node - nb.tip] <<- tpc[k]
        k <<- k + 1
        current.node <<- edge[l, 1]
    }
    if (!length(grep(",", tp))) {
        obj <- list(edge = matrix(c(2, 1), 1, 2), Nnode = 1)
        tp <- unlist(strsplit(tp, "[\\(\\);]"))
        obj$tip.label <- tp[2]
        if (tp[3] != "") obj$node.label <- tp[3]
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

    index <- numeric(nb.edge + 1)
    index[node] <- nb.edge
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
    edge <- edge[-nb.edge, ]
    obj <- list(edge = edge, tip.label = tip.label,
                Nnode = nb.node, node.label = node.label)
    obj$node.label <-
        if (all(obj$node.label == "NA")) NULL
        else gsub("^NA", "", obj$node.label)
    class(obj) <- "phylo"
    obj
}

read.nexus <- function(file, tree.names = NULL)
{
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    ## remove all comments
    ## (this might not work if there are square brackets within the comments)
    LEFT <- grep("\\[", X)
    RIGHT <- grep("\\]", X)
    if (length(LEFT)) { # in case there are no comments at all
        w <- LEFT == RIGHT
        if (any(w)) { # in case all comments use at least 2 lines
            s <- LEFT[w]
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
            ## The above regexp was quite tough to find: it makes
            ## possible to delete series of comments on the same line:
            ##       ...[...]xxx[...]...
            ## without deleting the "xxx". This regexp is in three parts:
            ##       \\[      [^]]*       \\]
            ## where [^]]* means "any character, except "]", repeated zero
            ## or more times" (note that the ']' is not escaped here).
            ## The previous version was:
            ##       X[s] <- gsub("\\[.*\\]", "", X[s])
            ## which deleted the "xxx". (EP  2008-06-24)
        }
        w <- !w
        if (any(w)) {
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1))
                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
        }
    }
    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    semico <- grep(";", X)
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    translation <- if (length(i2) == 1 && i2 > i1) TRUE else FALSE
    if (translation) {
        end <- semico[semico > i2][1]
        x <- X[(i2 + 1):end] # assumes there's a 'new line' after "TRANSLATE"
        ## x <- gsub("TRANSLATE", "", x, ignore.case = TRUE)
        x <- unlist(strsplit(x, "[,; \t]"))
        x <- x[nzchar(x)]
        TRANS <- matrix(x, ncol = 2, byrow = TRUE)
        TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
        n <- dim(TRANS)[1]
    }
    start <-
        if (translation) semico[semico > i2][1] + 1
        else semico[semico > i1][1]
    end <- endblock[endblock > i1][1] - 1
    tree <- X[start:end]
    rm(X)
###    tree <- gsub("^.*= *", "", tree)
    ## check whether there are empty lines from the above manips:
    tree <- tree[tree != ""]
    semico <- grep(";", tree)
    Ntree <- length(semico) # provisional -- some ";" may actually mark end of commands
    ## are some trees on several lines?
    ## -- this actually 'packs' all characters ending with a ";" in a single string --
    if (Ntree == 1 && length(tree) > 1) STRING <- paste(tree, collapse = "") else {
        if (any(diff(semico) != 1)) {
            STRING <- character(Ntree)
            s <- c(1, semico[-Ntree] + 1)
            j <- mapply(":", s, semico)
            if (is.list(j)) {
                for (i in 1:Ntree)
                    STRING[i] <- paste(tree[j[[i]]], collapse = "")
            } else {
                for (i in 1:Ntree)
                    STRING[i] <- paste(tree[j[, i]], collapse = "")
            }
        } else STRING <- tree
    }
    rm(tree)
    ## exclude the possible command lines ending with ";":
    STRING <- STRING[grep("^[[:blank:]]*tree.*= *", STRING, ignore.case = TRUE)]
    Ntree <- length(STRING) # update Ntree
    ## get the tree names:
    nms.trees <- sub(" *= *.*", "", STRING) # only the first occurence of "="
    nms.trees <- sub("^ *tree *", "", nms.trees, ignore.case = TRUE)
    STRING <- sub("^.*= *", "", STRING) # delete title and 'TREE' command with 'sub'
    STRING <- gsub(" ", "", STRING) # delete all white spaces
    colon <- grep(":", STRING)
    if (!length(colon)) {
        trees <- lapply(STRING, clado.build)
    } else if (length(colon) == Ntree) {
        trees <-
            if (translation) lapply(STRING, .treeBuildWithTokens)
            else lapply(STRING, tree.build)
    } else {
        trees <- vector("list", Ntree)
        trees[colon] <- lapply(STRING[colon], tree.build)
        nocolon <- (1:Ntree)[!1:Ntree %in% colon]
        trees[nocolon] <- lapply(STRING[nocolon], clado.build)
        if (translation) {
            for (i in 1:Ntree) {
                tr <- trees[[i]]
                for (j in 1:n) {
                    ind <- which(tr$tip.label[j] == TRANS[, 1])
                    tr$tip.label[j] <- TRANS[ind, 2]
                }
                if (!is.null(tr$node.label)) {
                    for (j in 1:length(tr$node.label)) {
                        ind <- which(tr$node.label[j] == TRANS[, 1])
                        tr$node.label[j] <- TRANS[ind, 2]
                    }
                }
                trees[[i]] <- tr
            }
            translation <- FALSE
        }
    }
    for (i in 1:Ntree) {
        tr <- trees[[i]]
        ## Check here that the root edge is not incorrectly represented
        ## in the object of class "phylo" by simply checking that there
        ## is a bifurcation at the root
        if (!translation) n <- length(tr$tip.label)
        ROOT <- n + 1
        if (sum(tr$edge[, 1] == ROOT) == 1 && dim(tr$edge)[1] > 1) {
            stop(paste("There is apparently two root edges in your file: cannot read tree file.\n  Reading NEXUS file aborted at tree no.", i, sep = ""))
        }
    }
    if (Ntree == 1) {
        trees <- trees[[1]]
        if (translation) {
            trees$tip.label <-
                if (length(colon)) TRANS[, 2] else
                TRANS[, 2][as.numeric(trees$tip.label)]
        }
    } else {
        if (!is.null(tree.names)) names(trees) <- tree.names
        if (translation) {
            if (length(colon) == Ntree) # .treeBuildWithTokens() was used
                attr(trees, "TipLabel") <- TRANS[, 2]
            else { # reassign the tip labels then compress
                for (i in 1:Ntree)
                    trees[[i]]$tip.label <-
                        TRANS[, 2][as.numeric(trees[[i]]$tip.label)]
                trees <- .compressTipLabel(trees)
            }
        }
        class(trees) <- "multiPhylo"
        if (!all(nms.trees == "")) names(trees) <- nms.trees
    }
    if (length(grep("[\\/]", file)) == 1)
        if (!file.exists(file)) # suggestion by Francois Michonneau
            file <- paste(getwd(), file, sep = "/")
    attr(trees, "origin") <- file
    trees
}
