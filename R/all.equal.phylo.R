### all.equal.phylo.R  (2002-08-28)
###
###     Global Comparison of two Phylogenies
###
### Copyright 2002 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

all.equal.phylo <- function(target, current, ...)
{
    phy1 <- target
    phy2 <- current
    rm(target, current)
    ## compare topologies using numbers of branches from each tip to the root
    nb.br.tip2root <- function(phy, nb.tip)
    {
        seq.nod <- list()
        for (i in as.character(1:nb.tip)) {
            vec <- i
            j <- i
            while (j != "-1") {
                ind <- which(phy$edge[, 2] == j)
                j <- phy$edge[ind, 1]
                vec <- c(vec, j)
            }
            seq.nod[[i]] <- vec
        }
        return(unlist(lapply(seq.nod, length)) - 1)
    }
    if (class(phy1) != "phylo") stop("object \"target\" is not of class \"phylo\"")
    if (class(phy2) != "phylo") stop("object \"current\" is not of class \"phylo\"")
    ## first, do a global comparison of both phylogenies
    if (identical(all.equal.list(phy1, phy2, ...), TRUE)) return(TRUE)
    else {
        tmp <- as.numeric(phy1$edge)
        nb.tip1 <- max(tmp)
        nb.node1 <- -min(tmp)
        tmp <- as.numeric(phy2$edge)
        nb.tip2 <- max(tmp)
        nb.node2 <- -min(tmp)
        ## compare numbers of tips and of nodes
        equal.nb.tip <- equal.nb.node <-  FALSE
        if (nb.tip1 == nb.tip2) equal.nb.tip <- TRUE
        if (nb.node1 == nb.node2) equal.nb.node <- TRUE
        if (!equal.nb.tip) {
            cat("Number of tips differ:", nb.tip1, "and", nb.tip2, "\n")
            cat("---Comparison stopped here.---\n\n")
        }
        else {
            cat("Number of tips are equal:", nb.tip1, "\n")
            if (!equal.nb.node) {
                cat("Number of nodes differ:", nb.node1, "and", nb.node2, "\n")
                cat("---Comparison stopped here.---\n\n")
            }
            else {
                cat("Number of nodes are equal:", nb.node1, "\n")
                ## compare tip labels
                equal.tip.label <- FALSE
                if (all.equal(sort(phy1$tip.label), sort(phy2$tip.label)))
                  equal.tip.label <- TRUE
                if (!equal.tip.label) cat("Tip labels differ.\n")
                else cat("Tip labels are the same.\n")
                ## compare topologies
                nb.br.tip2root1 <- sort(nb.br.tip2root(phy1, nb.tip1))
                nb.br.tip2root2 <- sort(nb.br.tip2root(phy2, nb.tip2))
                equal.nb.br.tip2root <- all.equal(nb.br.tip2root1, nb.br.tip2root2, ...)
                if (identical(equal.nb.br.tip2root, TRUE)) {
                    cat("Topologies are the same.\n")
                    cat("The labeled trees are the same.\n")
                    cat("Both phylogenies seem identical.\n")
                }
                else {
                    if (equal.nb.br.tip2root[2] == "TRUE") cat("Topologies are the same\n")
                    else cat("Topologies differ.\n")
                    lab1 <- t1$tip.label[as.numeric(names(nb.br.tip2root1))]
                    lab2 <- t2$tip.label[as.numeric(names(nb.br.tip2root2))]
                    if (identical(all.equal(lab1, lab2), TRUE)) {
                        cat("The labeled trees are the same.\n")
                        cat("Both phylogenies seem identical.\n")
                    }
                    else cat("The labeled trees differ.\n")
                }
            }
        }
    }
}
