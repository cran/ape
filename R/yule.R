### yule.R  (2003-12-24)
###
###     Fits Yule Model to a Phylogenetic Tree
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

yule <- function(phy)
{
    if (!is.binary.tree(phy))
      stop("tree must be dichotomous to fit the Yule model.")
    bt <- rev(sort(branching.times(phy))) # branching times from past to present
    ni <- cumsum(rev(table(bt))) + 1
    X <- sum(phy$edge.length)
    nb.node <- -min(as.numeric(phy$edge))
    if (is.null(phy$root.edge)) {
        nb.node <- nb.node - 1
    } else {
        X <- X + phy$root.edge
        ni <- c(1, ni)
    }
    lambda <- nb.node / X
    se <- lambda / sqrt(nb.node)
    loglik <- -lambda * X + sum(log(ni[-length(ni)])) + nb.node * log(lambda)
    obj <- list(lambda = lambda, se = se, loglik = loglik)
    class(obj) <- "yule"
    obj
}
