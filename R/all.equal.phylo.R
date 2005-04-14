### all.equal.phylo.R  (2005-03-25)
###
###     Global Comparison of two Phylogenies
###
### Copyright 2005 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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
    ## first, do a global comparison of both phylogenies
    if (identical(all.equal.list(target, current, ...), TRUE)) return(TRUE)
    else {
        tmp <- as.numeric(target$edge)
        nb.tip1 <- max(tmp)
        nb.node1 <- -min(tmp)
        tmp <- as.numeric(current$edge)
        nb.tip2 <- max(tmp)
        nb.node2 <- -min(tmp)
        ## compare numbers of tips and of nodes
        equal.nb.tip <- equal.nb.node <-  FALSE
        if (nb.tip1 == nb.tip2) equal.nb.tip <- TRUE
        if (nb.node1 == nb.node2) equal.nb.node <- TRUE
        if (!equal.nb.tip) {
            msg <- c(paste("Number of tips differ:",
                           nb.tip1, "and", nb.tip2),
                     "---Comparison stopped here.---")
        }
        else {
            msg <- paste("Number of tips are equal:", nb.tip1)
            if (!equal.nb.node)
              msg <- c(msg, paste("Number of nodes differ:",
                                  nb.node1, "and", nb.node2),
                       "---Comparison stopped here.---")
            else {
                msg <- c(msg, paste("Number of nodes are equal:",
                                    nb.node1))
                ## compare tip labels
                equal.tip.label <- FALSE
                if (all.equal(sort(target$tip.label), sort(current$tip.label)))
                  equal.tip.label <- TRUE
                if (!equal.tip.label)
                  msg <- c(msg, "Tip labels differ.")
                else msg <- c(msg, "Tip labels are the same.")
                ## compare topologies
                tip2root1 <- sort(nb.br.tip2root(target, nb.tip1))
                tip2root2 <- sort(nb.br.tip2root(current, nb.tip2))
                equal.nb.br.tip2root <- all.equal(tip2root1, tip2root2, ...) # not sure about the dot-dot-dot here (EP 7-12-2002)
                if (identical(equal.nb.br.tip2root, TRUE)) {
                    msg <- c(msg, "Topologies are the same.")
                    names(tip2root1) <-
                      target$tip.label[as.numeric(names(tip2root1))]
                    names(tip2root2) <-
                      current$tip.label[as.numeric(names(tip2root2))]
                    tip2root1 <- tip2root1[sort(names(tip2root1))]
                    tip2root2 <- tip2root2[sort(names(tip2root2))]
                    if (identical(all.equal(tip2root1, tip2root2), TRUE)) return(TRUE)
                    else
                      msg <- c(msg, "The labeled trees differ.")
                }
                else
                  msg <- c(msg, "Topologies differ.",
                           "---Comparison stopped here.---")
            }
        }
        msg
    }
}
