### yule.cov.R  (2004-10-15)
###
###           Fits the `Yule Model With Covariates'
###
### Copyright 2004 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

yule.cov <- function(phy, formula, data = NULL)
{
    if (is.null(data)) data <- parent.frame()
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    if (!is.null(phy$node.label)) phy$node.label <- NULL
    bt <- sort(branching.times(phy)) # branching times (from present to past)
    bt <- rev(bt) # branching times from past to present
    ni <- cumsum(rev(table(bt))) + 1
    X <- model.matrix(formula, data)
    rownames(X) <- as.character(c(1:nb.tip, -(1:nb.node)))
    Xi <- X[phy$edge[, 1], ]
    Xj <- X[phy$edge[, 2], ]
    dev <- function(b) {
        2 * sum(((1 / (1 + exp(-(Xi %*% b)))) +
                 (1 / (1 + exp(-(Xj %*% b)))))
                * phy$edge.length / 2) -
         2 * (sum(log(ni[-length(ni)])) +
              sum(log((1 / (1 + exp(-(X[as.numeric(rownames(X)) < -1, ] %*% b)))))))
    }
    out <- nlm(function(p) dev(p),
               p = c(rep(0, ncol(X) - 1), -1),
               hessian = TRUE)
    Dev <- out$minimum
    para <- matrix(NA, ncol(X), 2)
    para[, 1] <- out$estimate
    para[, 2] <- sqrt(diag(solve(out$hessian)))
    rownames(para) <- colnames(X)
    colnames(para) <- c("Estimate", "StdErr")
    cat("\n---- Yule Model With Covariates ----\n\n")
    cat("    Phylogenetic tree:", deparse(substitute(phy)), "\n")
    cat("       Number of tips:", nb.tip, "\n")
    cat("      Number of nodes:", nb.node, "\n")
    cat("             Deviance:", Dev, "\n")
    cat("       Log-likelihood:", -Dev/2, "\n\n")
    cat("  Parameter estimates:\n")
    print(para)
    cat("\n")
}
