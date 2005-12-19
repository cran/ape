### compar.ou.R  (2005-09-15)
###
###    Ornstein--Uhlenbeck Model for Continuous Characters
###
### Copyright 2005 Emmanuel Paradis
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

compar.ou <- function(x, phy, node = NULL, alpha = NULL)
{
    if (class(phy) != "phylo")
      stop('object "phy" is not of class "phylo".')
    if (!is.numeric(x)) stop("'x' must be numeric.")
    if (!is.null(names(x))) {
        if(!any(is.na(match(names(x), phy$tip.label))))
          x <- x[phy$tip.label]
        else
          warning('the names of argument "x" and the names of the tip labels
did not match: the former were ignored in the analysis.')
    }
    if (is.numeric(node)) node <- as.character(node)
    if (is.null(node)) node <- character(0)
    if ("-1" %in% node) node <- node[-which(node == "-1")]
    nb.tip <- max(as.numeric(phy$edge))
    bt <- branching.times(phy)
    Tmax <- bt["-1"]
    Wend <- matrix(0, nb.tip, length(node) + 1)
    colnames(Wend) <- c(names(sort(bt[node])), "-1")
    Wstart <- Wend
    Wstart[, ncol(Wstart)] <- Tmax
    root2tip <- .Call("seq_root2tip", phy$edge[, 1],
                      phy$edge[, 2], PACKAGE = "ape")
    for (i in 1:nb.tip) {
        last.change <- names(Tmax)
        for (j in root2tip[[i]]) {#[-1]) {# don't need to look at the root
            if (j %in% node) {
                Wend[i, last.change] <-
                  Wstart[i, as.character(j)] <- bt[as.character(j)]
                last.change <- as.character(j)
            }
        }
    }
    W <- cophenetic.phylo(phy)
    dev <- function(p) {
        M <- rowSums(exp(-p[1] * Wstart) - exp(-p[1] * Wend) * p[-(1:2)])
        V <- exp(-p[1]*W) * (1 - exp(-2*p[1]*(Tmax - W/2)))
        nb.tip*log(2*pi*p[2]) + log(det(V)) +
          (t(x - M) %*% chol2inv(V) %*% (x - M)) / p[2]
    }
    if (is.null(alpha))
      out <- nlm(function(p) dev(p),
                 p = c(0.1, 1, rep(mean(x), ncol(Wstart))),
                 hessian = TRUE)
    else
      out <- nlm(function(p) dev(c(alpha, p)),
                 p = c(1, rep(mean(x), ncol(Wstart))),
                 hessian = TRUE)
    para <- cbind(out$estimate, sqrt(diag(solve(out$hessian))))
    nms <- c("sigma2", paste("theta", 1:ncol(Wstart), sep = ""))
    if (is.null(alpha)) nms <- c("alpha", nms)
    dimnames(para) <- list(nms, c("estimate", "stderr"))
    obj <- list(deviance = out$minimum, para = para, call = match.call())
    class(obj) <- "compar.ou"
    obj
}
