### chronopl.R  (2005-10-20)
###
###     Molecular Dating With Penalized Likelihood
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

chronopl <- function(phy, lambda, node.age = NULL, nodes = NULL,
                     CV = FALSE)
{
    if (!is.binary.tree(phy)) stop("the tree must be dichotomous")
    if (is.null(node.age)) {
        node.age <- 1
        names(node.age) <- nodes <- "-1"
    } else names(node.age) <- as.character(nodes)
    if (CV) nodes.save <- nodes
    ini.rate <- phy$edge.length
    N <- dim(phy$edge)[1] # number of branches
    n <- (N + 2) / 2      # number of tips
    nodes <- as.character(-(1:(n - 1)))
    nb.basal.br <- table(phy$edge[, 1])["-1"]
    basal <- numeric(0)
    ind <- matrix(NA, N - nb.basal.br, 2)
    j <- 1
    for (i in 1:N) {
        if (phy$edge[i, 1] == "-1") {
            basal <- c(basal, i)
            next
        }
        ind[j, 1] <- i
        ind[j, 2] <- which(phy$edge[, 2] == phy$edge[i, 1])
        j <- j + 1
    }

    tip.age <- rep(0, n)
    names(tip.age) <- as.character(1:n)
    known.ages <- c(node.age, tip.age)
    nms.unknown.ages <- nodes[!nodes %in% names(node.age)]

    ini.time <- node.depth(phy$edge)[1:(n - 1)] - 1
    ini.time <- ini.time / max(ini.time) # 'ini.time' has names!
    ini.time <- ini.time * node.age / ini.time[names(node.age)]
    ## check that there is no negative branch lengths:
    ini.time[names(node.age)] <- node.age
    it <- c(ini.time, tip.age)
    ibl <- it[phy$edge[, 1]] - it[phy$edge[, 2]]
    if (any(ibl < 0)) {
        for (i in which(ibl < 0))
          if (phy$edge[i, 1] %in% names(node.age))
            ini.time[phy$edge[i, 2]] <- ini.time[phy$edge[i, 1]] - 1e-3
          else
            ini.time[phy$edge[i, 1]] <- ini.time[phy$edge[i, 2]] + 1e-3
    }

    ploglik <- function(rate, node.time) {
        ## since the names are not passed by nlm, they must be
        ## re-inserted here:
        names(node.time) <- nms.unknown.ages
        time <- c(node.time, known.ages)
        real.edge.length <- time[phy$edge[, 1]] - time[phy$edge[, 2]]
        B <- rate * real.edge.length
        loglik <- sum(-B + phy$edge.length * log(B) -
                      lfactorial(phy$edge.length))
        loglik - lambda * (sum((rate[ind[, 1]] - rate[ind[, 2]])^2)
                           + var(rate[basal]))
    }

    out <- nlm(function(p) -ploglik(p[1:N], p[-(1:N)]),
               p = c(ini.rate, ini.time[nms.unknown.ages]),
               iterlim = 500)

    attr(phy, "ploglik") <- -out$minimum
    attr(phy, "rates") <- out$estimate[1:N]
    node.ages <- out$estimate[-(1:N)]
    names(node.ages) <- nms.unknown.ages
    node.ages <- c(node.ages, known.ages)
    tmp <- node.ages[phy$edge[, 1]] - node.ages[phy$edge[, 2]]
    names(tmp) <- NULL
    ophy <- phy
    phy$edge.length <- tmp

    if (CV) {
        BT <- branching.times(phy)
        attr(phy, "D2") <-
          chronopl.cv(ophy, lambda, node.age, nodes.save, n, BT)
    }

    phy
}

chronopl.cv <- function(ophy, lambda, node.age, nodes, n, BT)
### ophy: the original phylogeny
### n: number of tips
### BT: branching times of the estimated chronogram
###
### Note that we assume here that the order of the nodes
### in node.label are not modified by the drop.tip operation
{
    cat("Doing cross-validation\n")
    D2 <- numeric(n)
    nodes <- as.character(nodes)

    if (is.null(node.age)) {
        node.age <- 1
        names(node.age) <- nodes <- "-1"
    } else names(node.age) <- nodes

    for (i in 1:n) {
        cat("  dropping tip", i, "\n")
        tr <- drop.tip(ophy, i)
        j <- which(ophy$edge[, 2] == as.character(i))
        if (ophy$edge[j, 1] %in% names(node.age)) {
            k <- which(names(node.age) == ophy$edge[j, 1])
            nodes <- nodes[-k]
            node.age <- node.age[-k]
        }
        if (length(nodes)) {
            chr <- chronopl(tr, lambda, node.age, nodes)
            tmp <- BT[-as.numeric(ophy$edge[j, 1])]
            D2[i] <- sum((tmp - branching.times(chr))^2 / tmp)
        } else D2[i] <- 0
    }
    D2
}
