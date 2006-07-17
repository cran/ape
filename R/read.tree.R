### read.tree.R (2005-08-18)
###
###     Read Tree File in Parenthetic Format
###
### Copyright 2002-2005 Emmanuel Paradis
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

tree.build <- function(tp) {
    add.internal <- function() {
        edge[j, 1] <<- current.node
        node <<- node - 1
        edge[j, 2] <<- current.node <<- node
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
        node.label[-current.node] <<- X[1]
        edge.length[l] <<- as.numeric(X[2])
        k <<- k + 1
        current.node <<- edge[l, 1]
    }
    tsp <- unlist(strsplit(tp, NULL))
    tpc <- unlist(strsplit(tp, "[\\(\\),;]"))
    tpc <- tpc[tpc != ""]
    skeleton <- tsp[tsp == "(" | tsp == ")" | tsp == "," | tsp == ";"]
    nsk <- length(skeleton)
    nb.node <- length(skeleton[skeleton == ")"])
    nb.tip <- length(skeleton[skeleton == ","]) + 1
    ## We will assume there is an edge at the root;
    ## if so, it will be removed and put into a vector
    nb.edge <- nb.node + nb.tip
    node.label <- character(nb.node)
    tip.label <- character(nb.tip)

    edge.length <- numeric(nb.edge)
    edge <- matrix(NA, nb.edge, 2)
    edge[nb.edge, 1] <- 0  # see comment above
    edge[nb.edge, 2] <- -1 #
    node <- -1                              # node number
    current.node <- node
    j <- 1                                  # index of the line number of edge
    k <- 1                                  # index of the line number of tpc
    tip <- 1                                # tip number

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
    if (is.na(node.label[1])) node.label[1] <- ""
    edge <- edge[-nb.edge, ]
    mode(edge) <- "character"
    root.edge <- edge.length[nb.edge]
    edge.length <- edge.length[-nb.edge]
    obj <- list(edge = edge,
                edge.length = edge.length,
                tip.label = tip.label,
                node.label = node.label,
                root.edge = root.edge)
    if (all(obj$node.label == "")) obj$node.label <- NULL
    if (is.na(obj$root.edge)) obj$root.edge <- NULL
    if (all(is.na(obj$edge.length))) obj$edge.length <- NULL # added 2005-08-18
    class(obj) <- "phylo"
    obj
}

read.tree <- function(file = "", format = "Newick", rooted = TRUE, text = NULL,
                      tree.names = NULL, skip = 0, comment.char = "#", ...)
{
    if (!is.null(text)) {
        if (!is.character(text))
          stop("argument `text' must be of mode character")
        tree <- text
    }
    else {
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
        STRING <- list()
        for (i in 1:nb.tree)
          STRING[[i]] <- paste(tsp[x[i]:y[i]], sep = "", collapse = "")
    }
    list.obj <- list()
    for (i in 1:nb.tree) {
        list.obj[[i]] <- if (length(grep(":", STRING[[i]]))) tree.build(STRING[[i]]) else clado.build(STRING[[i]])
        ## Check here that the root edge is not incorrectly represented
        ## in the object of class "phylo" by simply checking that there
        ## is a bifurcation at the root (node "-1")
        if(sum(list.obj[[i]]$edge[, 1] == "-1") == 1) {
            warning("The root edge is apparently not correctly represented\nin your tree: this may be due to an extra pair of\nparentheses in your file; the returned object has been\ncorrected but your file may not be in a valid Newick\nformat")
            ind <- which(list.obj[[i]]$edge[, 1] == "-1")
            list.obj[[i]]$root.edge <- list.obj[[i]]$edge.length[ind]
            list.obj[[i]]$edge.length <- list.obj[[i]]$edge.length[-ind]
            list.obj[[i]]$edge <- list.obj[[i]]$edge[-ind, ]
            for (j in 1:length(list.obj[[i]]$edge))
              if (as.numeric(list.obj[[i]]$edge[j]) < 0)
                list.obj[[i]]$edge[j] <- as.character(as.numeric(list.obj[[i]]$edge[j]) + 1)
            ## Check a second time and if there is still a problem...!!!
            if(sum(list.obj[[i]]$edge[, 1] == "-1") == 1)
              stop("There is apparently two root edges in your file:\ncannot read tree file")
        }
    }
    if (nb.tree == 1) list.obj <- list.obj[[1]] else {
        if (is.null(tree.names)) names(list.obj) <- paste("tree", 1:nb.tree, sep = "")
        else names(list.obj) <- tree.names
        class(list.obj) <- c("multi.tree", "phylo")
    }
    list.obj
}
