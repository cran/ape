## read.tree.R (2017-09-21)

##   Read Tree Files in Parenthetic Format

## Copyright 2002-2017 Emmanuel Paradis, Daniel Lawson and Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

##tree.build <- function(tp)
##{
##    add.internal <- function() {
##        edge[j, 1] <<- current.node
##        edge[j, 2] <<- current.node <<- node <<- node + 1L
##        index[node] <<- j # set index
##        j <<- j + 1L
##    }
##    add.terminal <- function() {
##        edge[j, 1] <<- current.node
##        edge[j, 2] <<- tip
##        index[tip] <<- j # set index
##        X <- unlist(strsplit(tpc[k], ":"))
##        tip.label[tip] <<- X[1]
##        edge.length[j] <<- as.numeric(X[2])
##        k <<- k + 1L
##        tip <<- tip + 1L
##        j <<- j + 1L
##    }
##    go.down <- function() {
##        l <- index[current.node]
##        X <- unlist(strsplit(tpc[k], ":"))
##        node.label[current.node - nb.tip] <<- X[1]
##        edge.length[l] <<- as.numeric(X[2])
##        k <<- k + 1L
##        current.node <<- edge[l, 1]
##    }
##    if (!length(grep(",", tp))) {
##        obj <- list(edge = matrix(c(2L, 1L), 1, 2))
##        tp <- unlist(strsplit(tp, "[\\(\\):;]"))
##        obj$edge.length <- as.numeric(tp[3])
##        obj$Nnode <- 1L
##        obj$tip.label <- tp[2]
##        if (tp[4] != "") obj$node.label <- tp[4]
##        class(obj) <- "phylo"
##        return(obj)
##    }
##
##    tpc <- unlist(strsplit(tp, "[\\(\\),;]"))
##    tpc <- tpc[nzchar(tpc)]
##    ## the following 2 lines are (slightly) faster than using gsub()
##    tsp <- unlist(strsplit(tp, NULL))
##    skeleton <- tsp[tsp %in% c("(", ")", ",", ";")]
##    nsk <- length(skeleton)
##    nb.node <- sum(skeleton == ")")
##    nb.tip <- sum(skeleton == ",") + 1
##    ## We will assume there is an edge at the root;
##    ## if so, it will be removed and put into a vector
##    nb.edge <- nb.node + nb.tip
##    node.label <- character(nb.node)
##    tip.label <- character(nb.tip)
##
##    edge.length <- numeric(nb.edge)
##    edge <- matrix(0L, nb.edge, 2)
##    current.node <- node <- as.integer(nb.tip + 1) # node number
##    edge[nb.edge, 2] <- node
##    index <- numeric(nb.edge + 1) # hash index to avoid which
##    index[node] <- nb.edge
##
##    ## j: index of the line number of edge
##    ## k: index of the line number of tpc
##    ## tip: tip number
##    j <- k <- tip <- 1L
##
##    for (i in 2:nsk) {
##        if (skeleton[i] == "(") add.internal() # add an internal branch (on top)
##        if (skeleton[i] == ",") {
##            if (skeleton[i - 1] != ")") add.terminal() # add a terminal branch
##        }
##        if (skeleton[i] == ")") {
##            if (skeleton[i - 1] == ",") { # add a terminal branch and go down one level
##                add.terminal()
##                go.down()
##            }
##            ## added by Klaus to allow singleton nodes (2017-05-26):
##            if (skeleton[i - 1] == "(") {
##                add.terminal()
##                go.down()
##            }
##            ## end
##            if (skeleton[i - 1] == ")") go.down() # go down one level
##        }
##    }
##
##    edge <- edge[-nb.edge, ]
##    obj <- list(edge = edge, Nnode = nb.node, tip.label = tip.label)
##    root.edge <- edge.length[nb.edge]
##    edge.length <- edge.length[-nb.edge]
##    if (!all(is.na(edge.length))) # added 2005-08-18
##        obj$edge.length <- edge.length
##    if (is.na(node.label[1])) node.label[1] <- ""
##    if (any(nzchar(node.label))) obj$node.label <- node.label
##    if (!is.na(root.edge)) obj$root.edge <- root.edge
##    class(obj) <- "phylo"
##    attr(obj, "order") <- "cladewise"
##    obj
##}

read.tree <- function(file = "", text = NULL, tree.names = NULL, skip = 0,
    comment.char = "", keep.multi = FALSE, ...)
{
    if (!is.null(text)) {
        if (!is.character(text))
            stop("argument `text' must be of mode character")
        tree <- text
    } else {
        tree <- scan(file = file, what = "", sep = "\n", quiet = TRUE,
                     skip = skip, comment.char = comment.char, ...)
    }

    ## Suggestion from Eric Durand and Nicolas Bortolussi (added 2005-08-17):
    if (identical(tree, character(0))) {
        warning("empty character string.")
        return(NULL)
    }

    tree <- gsub("[ \t]", "", tree)
    tree <- gsub("''", "", tree)

    single_quotes <- function(x, start = 1L) {
        z <- unlist(gregexpr("'", x))
        if (length(z) %% 2)
            stop("wrong number of single quotes around labels")
        l <- length(z) / 2
        tmp <- strsplit(x, "'")[[1]]
        ind_orig <- 2L * (1L:l)
        tmp_label <- paste0("@_", start:(start + l - 1), "_@")
        orig_label <- tmp[ind_orig]  #paste0("'", tmp[ind_orig], "'")
        names(orig_label) <- tmp_label
        for (i in 1:l) tmp[2 * i] <- tmp_label[i]
        tmp <- paste0(tmp, collapse = "")
        list(tmp, orig_label)
    }

    ## replace labels with single quotes
    z <- grepl("'", tree)
    if (any(z)) {
        Ntree <- length(tree)
        tmp_label <- vector("list", Ntree)
        for (i in 1:Ntree) {
            if (z[i]) {
                TMP <- single_quotes(tree[i])
                tree[i] <- TMP[[1]]
                tmp_label[[i]] <- TMP[[2]]
            }
        }
    }

    y <- unlist(gregexpr(";", tree))

### replace comments may handle them different later on

    ## if one tree per line much faster
    if (identical(y, nchar(tree))) { # check if always one tree per line
        Ntree <- length(y)
        STRING <- character(Ntree)
        for (i in 1:Ntree) {
            STRING[i] <- gsub("\\[[^]]*\\]", "", tree[i]) # delete comments (fix 2015-01-12)
        }
    } else {
        ## tree <- paste0(tree)
        ## tree <- unlist(strsplit(tree, ";"))
        ## tree <- paste0(tree, ";")
        tree <- unlist(strsplit(tree, NULL))
        y <- which(tree == ";")
        Ntree <- length(y)
        x <- c(1, y[-Ntree] + 1)
        ## Suggestion from Olivier Francois (added 2006-07-15):
        if (is.na(y[1])) return(NULL)
        STRING <- character(Ntree)
        for (i in 1:Ntree) {
            tmp <- paste0(tree[x[i]:y[i]], collapse = "")
            STRING[i] <- gsub("\\[[^]]*\\]", "", tmp) # delete comments (fix 2015-01-12)
        }
    }

    ## remove possible leading and trailing underscores
    STRING <- gsub("^_+", "", STRING)
    STRING <- gsub("_+$", "", STRING)

    getTreeName <- function(x) {
        res <- rep("", length(x))
        i <- regexpr("\\(", x)
        s <- i > 1
        if (any(s)) res[s] <- substr(x[s], 1, i[s] - 1)
        res
    }

    tmpnames <- getTreeName(STRING)
    if (is.null(tree.names) && any(nzchar(tmpnames))) tree.names <- tmpnames

    colon <- grep(":", STRING)
    if (!length(colon)) {
        obj <- lapply(STRING, .cladoBuild)
    } else if (length(colon) == Ntree) {
        obj <- lapply(STRING, .treeBuild)
    } else {
        obj <- vector("list", Ntree)
        obj[colon] <- lapply(STRING[colon], .treeBuild)
        nocolon <- (1:Ntree)[!1:Ntree %in% colon]
        obj[nocolon] <- lapply(STRING[nocolon], .cladoBuild)
    }

    for (i in 1:Ntree) {
        if (z[i]) {
            tmp_lab <- tmp_label[[i]]
            tip.label <- obj[[i]]$tip.label
            node.label <- obj[[i]]$node.label
            ind <- match(tip.label, names(tmp_lab))
            ind2 <- which(!is.na(ind))
            tip.label[ind2] <- tmp_lab[ind[ind2]]
            tmp_lab <- tmp_lab[-ind[ind2]]

            ind <- match(node.label, names(tmp_lab))
            ind2 <- which(!is.na(ind))
            if (length(ind2)) {
                node.label[ind2] <- tmp_lab[ind[ind2]]
                tmp_lab <- tmp_lab[-ind[ind2]]
            }
            if (length(tmp_lab)) {
                for (j in 1:length(tmp_lab)) {
                    node.label <- gsub(names(tmp_lab)[j], tmp_lab[j], node.label)
                    tip.label <- gsub(names(tmp_lab)[j], tmp_lab[j], tip.label)
                }
            }
            obj[[i]]$tip.label <- tip.label
            obj[[i]]$node.label <- node.label
        }
    }
    if (Ntree == 1 && !keep.multi) obj <- obj[[1]] else {
        if (!is.null(tree.names)) names(obj) <- tree.names
        class(obj) <- "multiPhylo"
    }
    obj
}
