### read.nexus.R  (2003-08-15)
###
###     Read Tree File in Nexus Format
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

clado.build <- function(tp) {
    add.internal <- function() {
        edge[j, 1] <<- current.node
        node <<- node - 1
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
        node.label[-current.node] <<- tpc[k]
        k <<- k + 1
        current.node <<- edge[l, 1]
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
#    if(node.label[1] == "NA") node.label[1] <- ""
    edge <- edge[-nb.edge, ]
    mode(edge) <- "character"
    obj <- list(edge = edge,
                tip.label = tip.label,
                node.label = node.label)
    if (all(obj$node.label == "NA")) obj$node.label <- NULL
    else obj$node.label <- gsub("^NA", "", obj$node.label)
    class(obj) <- "phylo"
    return(obj)
}

read.nexus <- function(file, tree.names = NULL)
{
    X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE)
    X <- gsub("[Ee][Nn][Dd][Bb][Ll][Oo][Cc][Kk];", "END;", X)
    endblock <- grep("[Ee][Nn][Dd];", X)
    semico <- grep(";", X)
    i1 <- grep("[Bb][Ee][Gg][Ii][Nn] [Tt][Rr][Ee][Ee][Ss];", X)
    i2 <- grep("[Tt][Rr][Aa][Nn][Ss][Ll][Aa][Tt][Ee]", X)
    translation <- FALSE
    if (length(i2) == 1) if (i2 > i1) translation <- TRUE
    if (translation) {
        end <- semico[semico > i2][1]
        x <- paste(X[i2:end], sep = "", collapse = "")
        x <- gsub("[Tt][Rr][Aa][Nn][Ss][Ll][Aa][Tt][Ee]", "", x)
        x <- unlist(strsplit(x, "[,; \t]"))
        x <- x[x != ""]
        TRANS <- matrix(x, ncol = 2, byrow = TRUE)
        TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
    }
    if (translation) start <- semico[semico > i2][1] + 1 else start <- semico[semico > i1][1]
    end <- endblock[endblock > i1][1] - 1
    tree <- paste(X[start:end], sep = "", collapse = "")
    if (length(grep("\\[&U\\]", tree)) > 0) {
        warning("at least one tree was unrooted\n(the current version of ape does not distinguish\nrooted and unrooted trees)")
        tree <- gsub("\\[&[Uu]\\]", "", tree)
    }
    tree <- gsub("\\[&[Rr]\\]", "", tree)
    tree <- gsub(" ", "", tree)
    tree <- unlist(strsplit(tree, "[=;]"))
    tree <- tree[grep("[\\(\\)]", tree)]
    nb.tree <- length(tree)
    STRING <- as.list(tree)
    trees <- list()
    for (i in 1:nb.tree) {
        if (length(grep(":", STRING[[i]]))) obj <- tree.build(STRING[[i]]) else obj <- clado.build(STRING[[i]])
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
    }
    if (nb.tree == 1) trees <- trees[[1]] else {
        if (is.null(tree.names)) names(trees) <- paste("tree", 1:nb.tree, sep = "")
        else names(trees) <- tree.names
        class(trees) <- c("phylo", "multi.tree")
    }
    if (length(grep("[\\/]", file)) == 1) attr(trees, "origin") <- file else attr(trees, "origin") <- paste(getwd(), file, sep = "/")
    return(trees)
}
