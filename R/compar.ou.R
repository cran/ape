### compar.ou.R (2006-10-05)
###
###   Ornstein--Uhlenbeck Model for Continuous Characters
###
### Copyright 2005-2006 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

compar.ou <- function(x, phy, node = NULL, alpha = NULL)
{
    if (class(phy) != "phylo")
      stop('object "phy" is not of class "phylo".')
    if (!is.numeric(x)) stop("'x' must be numeric.")
    if (!is.null(names(x))) {
        if (all(names(x) %in% phy$tip.label)) x <- x[phy$tip.label]
        else warning('the names of argument "x" and the names of the tip labels did not match: the former were ignored in the analysis.')
    }
    nb.tip <- length(phy$tip.label)
    root <- nb.tip + 1
    if (is.null(node)) node <- numeric(0)
    if (root %in% node) node <- node[node != root]
    bt <- branching.times(phy)
    Tmax <- bt[1]
    Wend <- matrix(0, nb.tip, length(node) + 1)
    colnames(Wend) <- c(names(sort(bt[node])), as.character(root))
    Wstart <- Wend
    Wstart[, ncol(Wstart)] <- Tmax
    root2tip <- .Call("seq_root2tip", phy$edge, nb.tip,
                      phy$Nnode, PACKAGE = "ape")
    for (i in 1:nb.tip) {
        last.change <- names(Tmax)
        for (j in root2tip[[i]]) {#[-1]) {# don't need to look at the root
            if (j %in% node) {
                jb <- as.character(j)
                Wend[i, last.change] <- Wstart[i, jb] <- bt[jb]
                last.change <- jb
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
