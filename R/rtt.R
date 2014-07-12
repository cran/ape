## rtt.R (2014-06-16)

##   Root a tree by root-to-tip regression

## Copyright (c) 2014, Rosemary McCloskey, BC Centre for Excellence in HIV/AIDS

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

if (getRversion() >= "2.15.1") utils::globalVariables("mclapply")

rtt <- function(t, tip.dates, ncpu = 1, objective = "correlation", opt.tol = .Machine$double.eps^0.25)
{
    if (ncpu > 1) library(parallel)

    ## These are objective functions which can be used to evaluate the "goodness" of
    ## a regression fit.
    if (objective == "correlation")
        objective <- function(x, y) cor.test(y, x)$estimate
    else if (objective == "rsquared")
        objective <- function(x, y) summary(lm(y ~ x))$r.squared
    else if (objective == "rms")
        objective <- function(x, y) -summary(lm(y ~ x))$sigma^2
    else
        stop('objective must be one of "correlation", "rsquared", or "rms"')

    t <- unroot(t)
    dist <- dist.nodes(t)[, 1:(t$Nnode + 2)]

    ## Do root-to-tip regressions for every possible choice of root.
    fits <- if (ncpu > 1)
        unlist(mclapply(1:nrow(dist), function(row) objective(tip.dates, dist[row, ]),
                                      mc.cores = ncpu))
        else unlist(lapply(1:nrow(dist), function(row) objective(tip.dates, dist[row, ])))

    ## Find the best one (highest value of objective function).
    fit.edge <- apply(t$edge, 2, function(e) fits[e])
    obj.edge <- apply(fit.edge, 1, mean)
    ## Compatibility with Path-O-Gen: find the last maximum, not the first.
    best.edge <- length(obj.edge) - which.max(rev(obj.edge)) + 1
    best.edge.parent <- t$edge[best.edge, 1]
    best.edge.child <- t$edge[best.edge, 2]
    best.edge.length <- t$edge.length[best.edge]

    ## Find the best location on that edge.
    f <- function(x) {
        dist <- x * dist[best.edge.parent, ] + (1 - x) * dist[best.edge.child, ]
        objective(tip.dates, dist)
    }
    best.pos <- optimize(f, c(0, 1), maximum = TRUE, tol = opt.tol)$maximum

    ## Reroot the tree at the optimal location
    new.root <- list(edge = matrix(c(2L, 1L), 1, 2), tip.label = "new.root", edge.length = 1, Nnode = 1L, root.edge = 1)
    class(new.root) <- "phylo"
    t <- bind.tree(t, new.root, where = best.edge.child, position = best.pos * best.edge.length)
    t <- collapse.singles(t)
    t <- root(t, "new.root")
    drop.tip(t, "new.root")
}
