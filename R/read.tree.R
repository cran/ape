### read.tree.R  (2003-06-06)
###
###     Read Tree File in Parenthetic Format
###
### Copyright 2003 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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
    ## if so, it will be removed and put in a vector
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
    ## The following lines are for the transition between R 1.7.0
    ## and R 1.8.0 (EP 2003-05-29)
    if (as.numeric(R.Version()$minor) >= 8) {
        if (is.na(node.label[1])) node.label[1] <- ""
    }
    else {
        if (node.label[1] == "NA") node.label[1] <- ""
    }
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
    tree <- gsub("[ \t]", "", tree)
    tsp <- unlist(strsplit(tree, NULL))
    ind <- which(tsp == ";")
    nb.tree <- length(ind)
    x <- c(1, ind[-nb.tree] + 1)
    y <- ind - 1
    STRING <- list()
    for (i in 1:nb.tree) STRING[[i]] <- paste(tsp[x[i]:y[i]], sep = "", collapse = "")
    list.obj <- list()
    for (i in 1:nb.tree) {
        if (length(grep(":", STRING[[i]]) > 0)) list.obj[[i]] <- tree.build(STRING[[i]])
        else list.obj[[i]] <- clado.build(STRING[[i]])
    }
    if (nb.tree == 1) list.obj <- list.obj[[1]]
    else {
        if (is.null(tree.names)) names(list.obj) <- paste("tree", 1:nb.tree, sep = "")
        else names(list.obj) <- tree.names
        class(list.obj) <- c("phylo", "multi.tree")
    }
    list.obj
}
