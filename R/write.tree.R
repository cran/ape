### write.tree.R  (2002-09-26)
###
###     Write Tree File in Parenthetic Format
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

write.tree <- function(phy, file = "", append = FALSE, format = "Newick")
{
    cp <- function(s) STRING <<- paste(STRING, s, sep = "")
    add.internal <- function(i)
    {
        cp("(")
        br <- which(phy$edge[, 1] == i)
        for (j in br) {
            desc <- phy$edge[j, 2]
            if (desc < 0) add.internal(desc) else add.terminal(j)
            if (j != br[length(br)]) cp(",")
        }
        cp(")")
        if (!is.null(phy$node.label)) cp(phy$node.label[-i])
        cp(":")
        cp(as.character(phy$edge.length[which(phy$edge[, 2] == i)]))
    }
    add.terminal <- function(i)
    {
        cp(phy$tip.label[phy$edge[i, 2]])
        cp(":")
        cp(as.character(phy$edge.length[i]))
    }
    if (class(phy) != "phylo")
      stop(paste("object \"", deparse(substitute(phy)),
                 "\" is not of class \"phylo\""), sep = "")
    nb.node <- -min(as.numeric(phy$edge))
    mode(phy$edge) <- "numeric"
    STRING <- "("
    br <- which(phy$edge[, 1] == -1)
    for (j in br) {
        desc <- phy$edge[j, 2]
        if (desc < 0) add.internal(desc) else add.terminal(j)
        if (j != br[length(br)]) cp(",")
    }
    if (is.null(phy$root.edge)) {
        cp(")")
        if (!is.null(phy$node.label)) cp(phy$node.label[1])
        cp(";")
    }
    else {
        cp(")")
        if (!is.null(phy$node.label)) cp(phy$node.label[1])
        cp(":")
        cp(as.character(phy$root.edge))
        cp(";")
    }
    cat(STRING, "\n", file = file, append = append, sep = "")
}
