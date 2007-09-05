## yule.R (2006-10-04)

##     Fits Yule Model to a Phylogenetic Tree

## yule: standard Yule model (constant birth rate)
## yule.cov: Yule model with covariates

## Copyright 2003-2006 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

yule <- function(phy)
{
    if (!is.binary.tree(phy))
      stop("tree must be dichotomous to fit the Yule model.")
    bt <- rev(sort(branching.times(phy))) # branching times from past to present
    ni <- cumsum(rev(table(bt))) + 1
    X <- sum(phy$edge.length)
    nb.node <- phy$Nnode
    if (is.null(phy$root.edge)) {
        nb.node <- nb.node - 1
    } else {
        X <- X + phy$root.edge
        ni <- c(1, ni)
    }
    lambda <- nb.node/X
    se <- lambda/sqrt(nb.node)
    loglik <- -lambda*X + sum(log(ni[-length(ni)])) + nb.node*log(lambda)
    obj <- list(lambda = lambda, se = se, loglik = loglik)
    class(obj) <- "yule"
    obj
}

yule.cov <- function(phy, formula, data = NULL)
{
    if (is.null(data)) data <- parent.frame()
    n <- length(phy$tip.label)
    nb.node <- phy$Nnode
    if (!is.null(phy$node.label)) phy$node.label <- NULL
    bt <- sort(branching.times(phy)) # branching times (from present to past)
    bt <- rev(bt) # branching times from past to present
    ni <- cumsum(rev(table(bt))) + 1
    X <- model.matrix(formula, data)
    Xi <- X[phy$edge[, 1], ]
    Xj <- X[phy$edge[, 2], ]
    dev <- function(b) {
        2 * sum(((1/(1 + exp(-(Xi %*% b)))) +
                 (1/(1 + exp(-(Xj %*% b)))))
                * phy$edge.length/2) -
         2 * (sum(log(ni[-length(ni)])) +
              sum(log((1/(1 + exp(-(X[-(1:(n + 1)), ] %*% b)))))))
    }
    out <- nlm(function(p) dev(p), p = c(rep(0, ncol(X) - 1), -1),
               hessian = TRUE)
    Dev <- out$minimum
    para <- matrix(NA, ncol(X), 2)
    para[, 1] <- out$estimate
    if (any(out$gradient == 0))
      warning("The likelihood gradient seems flat in at least one dimension (null gradient):\ncannot compute the standard-errors of the transition rates.\n")
    else para[, 2] <- sqrt(diag(solve(out$hessian)))
    rownames(para) <- colnames(X)
    colnames(para) <- c("Estimate", "StdErr")
    cat("\n---- Yule Model with Covariates ----\n\n")
    cat("    Phylogenetic tree:", deparse(substitute(phy)), "\n")
    cat("       Number of tips:", n, "\n")
    cat("      Number of nodes:", nb.node, "\n")
    cat("             Deviance:", Dev, "\n")
    cat("       Log-likelihood:", -Dev/2, "\n\n")
    cat("  Parameter estimates:\n")
    print(para)
    cat("\n")
}
