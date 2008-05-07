## chronopl.R (2008-03-26)

##   Molecular Dating With Penalized Likelihood

## Copyright 2005-2008 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

chronopl <- function(phy, lambda, age.min = 1, age.max = NULL,
                     node = "root", S = 1, tol = 1e-8,
                     CV = FALSE, eval.max = 500, iter.max = 500, ...)
{
    n <- length(phy$tip.label)
    ROOT <- n + 1
    if (length(node) == 1 && node == "root") node <- ROOT
    if (any(node <= n))
        stop("node numbers should be greater than the number of tips")
    zerobl <- which(phy$edge.length <= 0)
    if (length(zerobl)) {
        if (any(phy$edge[zerobl, 2] <= n))
            stop("at least one terminal branch is of length zero:
  you should remove it to have a meaningful estimation.")
        else {
            warning("at least one internal branch is of length zero:
  it was collapsed and some nodes have been deleted.")
            if (length(node) == 1 && node == ROOT)
                phy <- di2multi(phy)
            else {
                tmp <- FALSE
                if (is.null(phy$node.label)) {
                    tmp <- !tmp
                    phy$node.label <- paste("node", 1:phy$Nnode)
                }
                node.lab <- phy$node.label[node - n]
                phy <- di2multi(phy)
                node <- match(node.lab, phy$node.label) + n
                if (tmp) phy$node.label <- NULL
            }
        }
    }
    m <- phy$Nnode
    el <- phy$edge.length
    e <- phy$edge
    N <- dim(e)[1]
    TIPS <- 1:n
    EDGES <- 1:N

    ini.rate <- el
    el <- el/S

    ## `basal' contains the indices of the basal edges
    ## (ie, linked to the root):
    basal <- which(e[, 1] == ROOT)
    Nbasal <- length(basal)

    ## `ind' contains in its 1st column the index of all nonbasal
    ## edges, and in its second column the index of the edges
    ## where these edges come from (ie, this matrix contains pairs
    ## of contiguous edges), eg:

    ##         ___b___    ind:
    ##        |           |   |   |
    ## ___a___|           | b | a |
    ##        |           | c | a |
    ##        |___c___    |   |   |

    ind <- matrix(0L, N - Nbasal, 2)
    ind[, 1] <- EDGES[-basal]
    ind[, 2] <- match(e[EDGES[-basal], 1], e[, 2])

    age <- numeric(n + m)

    ##ini.time <- node.depth(phy)[-TIPS] - 1
    ini.time <- node.depth(phy) - 1

    ## first, rescale all times with respect to
    ## the age of the 1st element of `node':
    ratio <- age.min[1]/ini.time[node[1]]
    ini.time <- ini.time*ratio

    if (length(node) > 1) {
        ini.time[node] <- age.min
        real.edge.length <- ini.time[e[, 1]] - ini.time[e[, 2]]
        while (any(real.edge.length <= 0)) {
            for (i in EDGES) {
                if (real.edge.length[i] <= 0) {
                    if (e[i, 1] %in% node) {
                        ini.time[e[i, 2]] <-
                            ini.time[e[, 2]] - 2*real.edge.length[i]
                        next
                    }
                    if (e[i, 2] %in% node) {
                        ini.time[e[i, 1]] <-
                            ini.time[e[, 1]] + 2*real.edge.length[i]
                        next
                    }
                    ini.time[e[i, 2]] <-
                        ini.time[e[, 2]] - real.edge.length[i]
                    ini.time[e[i, 1]] <-
                        ini.time[e[, 1]] + real.edge.length[i]
                }
            }
            real.edge.length <- ini.time[e[, 1]] - ini.time[e[, 2]]
        }
    }

    ## `unknown.ages' will contain the index of the nodes of unknown age:
    unknown.ages <- n + 1:m

    ## define the bounds for the node ages:
    lower <- rep(tol, length(unknown.ages))
    upper <- rep(1/tol, length(unknown.ages))

    if (!is.null(age.max)) { # are some nodes known within some intervals?
        lower[node - n] <- age.min
        upper[node - n] <- age.max
        interv <- which(age.min != age.max)
        node <- node[-interv]
        if (length(node)) age[node] <- age.min[-interv]
    } else age[node] <- age.min

    if (length(node)) {
        unknown.ages <- unknown.ages[n - node]
        lower <- lower[n - node]
        upper <- upper[n - node]
    }

    ## `known.ages' contains the index of all nodes (internal and
    ## terminal) of known age:
    known.ages <- c(TIPS, node)

    ## concatenate the bounds for the rates:
    lower <- c(rep(tol, N), lower)
    upper <- c(rep(1 - tol, N), upper)

    minusploglik.gr <- function(rate, node.time) {
        grad <- numeric(N + length(unknown.ages))
        age[unknown.ages] <- node.time
        real.edge.length <- age[e[, 1]] - age[e[, 2]]
        if (any(real.edge.length < 0)) {
            grad[] <- 0
            return(grad)
        }
        ## gradient for the rates:
        ## the parametric part can be calculated without a loop:
        grad[EDGES] <- real.edge.length - el/rate
        if (Nbasal == 2) { # the simpler formulae if there's a basal dichotomy
            grad[basal[1]] <-
                grad[basal[1]] + lambda*(rate[basal[1]] - rate[basal[2]])
            grad[basal[2]] <-
                grad[basal[2]] + lambda*(rate[basal[2]] - rate[basal[1]])
        } else { # the general case
            for (i in 1:Nbasal)
                grad[basal[i]] <- grad[basal[i]] +
                    lambda*(2*rate[basal[i]]*(1 - 1/Nbasal) -
                            2*sum(rate[basal[-i]])/Nbasal)/(Nbasal - 1)
        }

        for (i in EDGES) {
            ii <- c(which(e[, 2] == e[i, 1]), which(e[, 1] == e[i, 2]))
            if (!length(ii)) next
            grad[i] <- grad[i] + lambda*(2*length(ii)*rate[i] - 2*sum(rate[ii]))
        }

        ## gradient for the 'node times'
        for (i in 1:length(unknown.ages)) {
            nd <- unknown.ages[i]
            ii <- which(e[, 1] == nd)
            grad[i + N] <-
                sum(rate[ii] - el[ii]/real.edge.length[ii])#, na.rm = TRUE)
            if (nd != ROOT) {
                ii <- which(e[, 2] == nd)
                grad[i + N] <- grad[i + N] -
                    rate[ii] + el[ii]/real.edge.length[ii]
            }
        }
        grad
    }

    minusploglik <- function(rate, node.time) {
        age[unknown.ages] <- node.time
        real.edge.length <- age[e[, 1]] - age[e[, 2]]
        if (any(real.edge.length < 0)) return(1e50)
        B <- rate*real.edge.length
        loglik <- sum(-B + el*log(B) - lfactorial(el))
        -(loglik - lambda*(sum((rate[ind[, 1]] - rate[ind[, 2]])^2)
                           + var(rate[basal])))
    }

    out <- nlminb(c(ini.rate, ini.time[unknown.ages]),
                  function(p) minusploglik(p[EDGES], p[-EDGES]),
                  function(p) minusploglik.gr(p[EDGES], p[-EDGES]),
                  control = list(eval.max = eval.max, iter.max = iter.max, ...),
                  lower = lower, upper = upper)

    attr(phy, "ploglik") <- -out$objective
    attr(phy, "rates") <- out$par[EDGES]
    attr(phy, "message") <- out$message
    age[unknown.ages] <- out$par[-EDGES]
    if (CV) ophy <- phy
    phy$edge.length <- age[e[, 1]] - age[e[, 2]]
    if (CV) attr(phy, "D2") <-
        chronopl.cv(ophy, lambda, age.min, age.max, node,
                    n, S, tol, eval.max, iter.max, ...)
    phy
}

chronopl.cv <- function(ophy, lambda, age.min, age.max, nodes,
                        n, S, tol, eval.max, iter.max, ...)
### ophy: the original phylogeny
### n: number of tips
### Note that we assume here that the order of the nodes
### in node.label are not modified by the drop.tip operation
{
    cat("Doing cross-validation\n")
    BT <- branching.times(ophy)
    D2 <- numeric(n)

    cat("  dropping tip")
    for (i in 1:n) {
        cat(" ", i, sep = "")
        tr <- drop.tip(ophy, i)
        j <- which(ophy$edge[, 2] == i)
        if (ophy$edge[j, 1] %in% nodes) {
            k <- which(nodes == ophy$edge[j, 1])
            node <- nodes[-k]
            agemin <- age.min[-k]
            agemax <- age.max[-k]
        } else node <- nodes
        if (length(node)) {
            chr <- chronopl(tr, lambda, age.min, age.max, node,
                            S, tol, FALSE, eval.max, iter.max, ...)
            D2[i] <- sum((tmp - branching.times(chr))^2 / tmp)
        } else D2[i] <- 0
    }
    cat("\n")
    D2
}
