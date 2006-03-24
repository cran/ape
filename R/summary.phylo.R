### summary.phylo.R (2006-03-24)
###
###       Print Summary of a Phylogeny
###
### Copyright 2003-2006 Emmanuel Paradis
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

summary.phylo <- function(object, ...)
{
    cat("\nPhylogenetic tree:", deparse(substitute(object)), "\n\n")
    tmp <- as.numeric(object$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    cat("  Number of tips:", nb.tip, "\n")
    cat("  Number of nodes:", nb.node, "\n")
    if (is.null(object$edge.length))
      cat("  No branch lengths.\n")
    else {
        cat("  Branch lengths:\n")
        cat("    mean:", mean(object$edge.length), "\n")
        cat("    variance:", var(object$edge.length), "\n")
        cat("    distribution summary:\n")
        print(summary(object$edge.length)[-4])
    }
    if (is.null(object$root.edge))
      cat("  No root edge.\n")
    else
      cat("  Root edge:", object$root.edge, "\n")
    if (nb.tip <= 10) {
        cat("  Tip labels:", object$tip.label[1], "\n")
        cat(paste("             ", object$tip.label[-1]), sep = "\n")
    }
    else {
        cat("  First ten tip labels:", object$tip.label[1], "\n")
        cat(paste("                       ", object$tip.label[2:10]), sep = "\n")
    }
    if (is.null(object$node.label))
      cat("  No node labels.\n")
    else {
        if (nb.node <= 10) {
            cat("  Node labels:", object$node.label[1], "\n")
            cat(paste("              ", object$node.label[-1]), sep = "\n")
        }
        else {
            cat("  First ten node labels:", object$node.label[1], "\n")
            cat(paste("                        ", object$node.label[2:10]), sep = "\n")

        }
    }
    if (!is.null(attr(object, "loglik"))) {
        cat("Phylogeny estimated by maximum likelihood.\n")
        cat("  log-likelihood:", attr(object, "loglik"), "\n\n")
        npart <- length(attr(object, "para"))
        for (i in 1:npart) {
            cat("partition ", i, ":\n", sep = "")
            print(attr(object, "para")[[i]])
            if (i == 1) next
            else cat("  contrast parameter (xi):",
                     attr(object, "xi")[i - 1], "\n")
        }
    }
}
