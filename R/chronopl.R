### chronopl.R (2007-01-17)
###
###    Molecular Dating With Penalized Likelihood
###
### Copyright 2005-2007 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

chronopl <- function(phy, lambda, node.age = 1, node = "root",
                     CV = FALSE)
{
    n <- length(phy$tip.label)
    n.node <- phy$Nnode
    if (n != n.node + 1)
      stop("the tree must be rooted AND dichotomous")
    N <- dim(phy$edge)[1]
    ROOT <- n + 1
    if (node == "root") node <- ROOT
    ini.rate <- phy$edge.length
    ## `known.ages' contains the index of all nodes (internal and
    ## terminal) of known age:
    known.ages <- c(1:n, node)
    ## `unknown.ages' contains the index of the nodes of unknown age:
    unknown.ages <- ((n + 1):(n + n.node))[-(node - n)]
    ## `basal' contains the indices of the basal edges (ie, linked to the root):
    basal <- which(phy$edge[, 1] == ROOT)

    ## `ind' contains in its 1st column the index of all nonbasal
    ## edges, and in its second column the index of the edges
    ## where these edges come from (ie, this matrix contains pairs
    ## of contiguous edges), eg:

    ##         ___b___    ind:
    ##        |           |   |   |
    ## ___a___|           | b | a |
    ##        |           | c | a |
    ##        |___c___    |   |   |

    ind <- matrix(NA, N - 2, 2)
    j <- 1
    for (i in 1:N) {
        if (phy$edge[i, 1] == ROOT) next
        ind[j, 1] <- i
        ind[j, 2] <- which(phy$edge[, 2] == phy$edge[i, 1])
        j <- j + 1
    }

    age <- rep(0, 2*n - 1)
    age[node] <- node.age

    tmp <- reorder(phy, "pruningwise")
    ini.time <- .C("node_depth", as.integer(n), as.integer(n.node),
                   as.integer(tmp$edge[, 1]), as.integer(tmp$edge[, 2]),
                   as.integer(N), double(n + n.node), DUP = FALSE,
                   PACKAGE = "ape")[[6]][-(1:n)] - 1
    ini.time <- ini.time/max(ini.time)
    ini.time <- ini.time*node.age/ini.time[known.ages[-(1:n)] - n]
    ## check that there are no negative branch lengths:
    ini.time[known.ages[-(1:n)] - n] <- node.age
    it <- c(age[1:n], ini.time)
    ibl <- it[phy$edge[, 1]] - it[phy$edge[, 2]]
    if (any(ibl < 0)) {
        for (i in which(ibl < 0))
          if (phy$edge[i, 1] %in% node)
            ini.time[phy$edge[i, 2]] <- ini.time[phy$edge[i, 1]] - 1e-3
          else
            ini.time[phy$edge[i, 1]] <- ini.time[phy$edge[i, 2]] + 1e-3
    }

    ploglik <- function(rate, node.time) {
        age[unknown.ages] <- node.time
        real.edge.length <- age[phy$edge[, 1]] - age[phy$edge[, 2]]
        B <- rate*real.edge.length
        loglik <- sum(-B + phy$edge.length*log(B) -
                      lfactorial(phy$edge.length))
        loglik - lambda * (sum((rate[ind[, 1]] - rate[ind[, 2]])^2)
                           + var(rate[basal]))
    }

    out <- nlm(function(p) -ploglik(p[1:N], p[-(1:N)]),
               p = c(ini.rate, ini.time[unknown.ages - n]),
               iterlim = 500)

    attr(phy, "ploglik") <- -out$minimum
    attr(phy, "rates") <- out$estimate[1:N]
    age[unknown.ages] <- out$estimate[-(1:N)]
    if (CV) ophy <- phy
    phy$edge.length <- age[phy$edge[, 1]] - age[phy$edge[, 2]]
    if (CV)
      attr(phy, "D2") <-
        chronopl.cv(ophy, lambda, node.age, node, n)
    phy
}

chronopl.cv <- function(ophy, lambda, node.age, nodes, n)
### ophy: the original phylogeny
### n: number of tips
### Note that we assume here that the order of the nodes
### in node.label are not modified by the drop.tip operation
{
    cat("Doing cross-validation\n")
    BT <- branching.times(ophy)
    D2 <- numeric(n)

    for (i in 1:n) {
        cat("  dropping tip", i, "\n")
        tr <- drop.tip(ophy, i)
        j <- which(ophy$edge[, 2] == i)
        if (ophy$edge[j, 1] %in% nodes) {
            k <- which(nodes == ophy$edge[j, 1])
            nodes <- nodes[-k]
            node.age <- node.age[-k]
        }
        if (length(nodes)) {
            chr <- chronopl(tr, lambda, node.age, nodes)
            ## <FIXME> à vérifier:
            ## tmp <- BT[as.character(ophy$edge[j, 1])]
            tmp <- BT[-(ophy$edge[j, 1] - n)]
            ## </FIXME>
            D2[i] <- sum((tmp - branching.times(chr))^2 / tmp)
        } else D2[i] <- 0
    }
    D2
}
