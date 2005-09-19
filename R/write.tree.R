### write.tree.R (2005-08-17)
###
###     Write Tree File in Parenthetic Format
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

write.tree <- function(phy, file = "", append = FALSE, format = "Newick",
                       multi.line = TRUE, digits = 10)
{
    brl <- !is.null(phy$edge.length)
    nodelab <- !is.null(phy$node.label)

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
        if (nodelab) cp(phy$node.label[-i])
        if (brl) {
            cp(":")
            cp(formatC(phy$edge.length[which(phy$edge[, 2] == i)],
                       format = "f", digits = digits))
        }
    }
    add.terminal <- function(i)
    {
        cp(phy$tip.label[phy$edge[i, 2]])
        if (brl) {
            cp(":")
            cp(formatC(phy$edge.length[i], format = "f", digits = digits))
        }
    }
    if (class(phy) != "phylo")
      stop(paste("object \"", deparse(substitute(phy)),
                 "\" is not of class \"phylo\""), sep = "")
    mode(phy$edge) <- "numeric"
    nb.node <- -min(phy$edge)
    STRING <- "("
    br <- which(phy$edge[, 1] == -1)
    for (j in br) {
        desc <- phy$edge[j, 2]
        if (desc < 0) add.internal(desc) else add.terminal(j)
        if (j != br[length(br)]) cp(",")
    }
    if (is.null(phy$root.edge)) {
        cp(")")
        if (nodelab) cp(phy$node.label[1])
        cp(";")
    } else {
        cp(")")
        if (nodelab) cp(phy$node.label[1])
        cp(":")
        cp(formatC(phy$root.edge, format = "f", digits = digits))
        cp(";")
    }
    if (nchar(STRING) < 80) multi.line <- FALSE
    if (multi.line) {
        tmp <- unlist(strsplit(STRING, NULL))
        wh <- grep("[,:]", tmp)
        fin <- seq(7, length(wh), by = 7)
        fin <- if (fin[length(fin)] == length(wh))
          wh[fin] else c(wh[fin], length(tmp))
        debut <- c(1, fin[-length(fin)] + 1)
        STRING <- character(length(fin))
        for (i in 1:length(STRING))
          STRING[i] <- paste(tmp[debut[i]:fin[i]], collapse = "")
    }
    if (file == "") return(STRING)
    else cat(STRING, file = file, append = append, sep = "\n")
}
