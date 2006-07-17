### ace.R  (2006-06-27)
###
###            Ancestral Character Estimation
###
### Copyright 2005-2006 Emmanuel Paradis and Ben Bolker
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

ace <- function(x, phy, type = "continuous", method = "ML", CI = TRUE,
                model = if (type == "continuous") "BM" else "ER",
                scaled = TRUE, kappa = 1, corStruct = NULL, ip = 0.1)
{
    if (class(phy) != "phylo")
      stop('object "phy" is not of class "phylo".')
    type <- match.arg(type, c("continuous", "discrete"))
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    if (nb.node != nb.tip - 1)      stop('"phy" is not rooted AND fully dichotomous.')
    if (length(x) != nb.tip)
      stop("length of phenotypic and of phylogenetic data do not match.")
    if (!is.null(names(x))) {
        if(!any(is.na(match(names(x), phy$tip.label))))
          x <- x[phy$tip.label]
        else warning('the names of argument "x" and the names of the tip labels
did not match: the former were ignored in the analysis.')
    }
    obj <- list()
    if (kappa != 1) phy$edge.length <- phy$edge.length^kappa
    if (type == "continuous") {
        if (method == "pic") {
            if (model != "BM")
              stop('the "pic" method can be used only with model = "BM".')
            phenotype <- as.numeric(rep(NA, nb.tip + nb.node))
            names(phenotype) <- as.character(c(1:nb.tip, -(1:nb.node)))
            if (is.null(names(x))) phenotype[1:nb.tip] <- x
            else for (i in 1:nb.tip) phenotype[i] <- x[phy$tip.label[i]]
            bl <- phy$edge.length
            unused <- rep(TRUE, nb.tip + nb.node)
            names(unused) <- names(phenotype)
            contr <- numeric(nb.node)
            names(contr) <- as.character(-(1:nb.node))
            if (CI) var.con <- contr
            while (sum(unused) > 1) {
                term <- names(phenotype[!is.na(phenotype) & unused])
                ind <- as.logical(match(phy$edge[, 2], term))
                ind[is.na(ind)] <- FALSE
                term.br <- matrix(phy$edge[ind], length(term), 2)
                basal <- names(which(table(term.br[, 1]) == 2))
                for (nod in basal) {
                    pair.ind <- which(phy$edge[, 1] == nod)
                    i <- pair.ind[1]
                    j <- pair.ind[2]
                    pair <- phy$edge[pair.ind, 2]
                    a <- pair[1]
                    b <- pair[2]
                    if (scaled)
                      contr[nod] <- (phenotype[a] - phenotype[b])/sqrt(bl[i] + bl[j])
                    else contr[nod] <- phenotype[a] - phenotype[b]
                    if (CI) var.con[nod] <- bl[i] + bl[j]
                    unused[pair] <- FALSE
                    phenotype[nod] <-
                      (phenotype[a] * bl[j] + phenotype[b] * bl[i])/(bl[i] + bl[j])
                    k <- which(phy$edge[, 2] == nod)
                    bl[k] <- bl[k] + bl[i] * bl[j]/(bl[i] + bl[j])
                }
            }
            obj$ace <- phenotype[-(1:nb.tip)]
            se <- sqrt(var.con)
            if (CI) {
                CI95 <- matrix(NA, nb.node, 2)
                CI95[, 1] <- obj$ace + se * qnorm(0.025)
                CI95[, 2] <- obj$ace - se * qnorm(0.025)
                obj$CI95 <- CI95
            }
        }
        if (method == "ML") {
            if (model == "BM") {
                dev <- function(p) {
                    X <- c(x, p[-1])
                    names(X) <- as.character(c(1:nb.tip, -(1:nb.node)))
                    x1 <- X[phy$edge[, 1]]
                    x2 <- X[phy$edge[, 2]]
                    -2 * (-sum((x1 - x2)^2/phy$edge.length)/(2*p[1]) -
                          nb.node * log(p[1]))
                }
                out <- nlm(function(p) dev(p),
                           p = c(1, rep(mean(x), nb.node)), hessian = TRUE)
                obj$loglik <- -out$minimum / 2
                obj$ace <- out$estimate[-1]
                se <- sqrt(diag(solve(out$hessian)))
                obj$sigma2 <- c(out$estimate[1], se[1])
                se <- se[-1]
                if (CI) {
                    CI95 <- matrix(NA, nb.node, 2)
                    CI95[, 1] <- obj$ace + se * qt(0.025, nb.node)
                    CI95[, 2] <- obj$ace - se * qt(0.025, nb.node)
                    obj$CI95 <- CI95
                }
            }
        }
        if (method == "GLS") {
            if (is.null(corStruct))
              stop('you must give a correlation structure if method = "GLS".')
            if (class(corStruct)[1] == "corMartins")
              M <- corStruct[1] * cophenetic.phylo(phy, full = TRUE)
            if (class(corStruct)[1] == "corGrafen")
              phy <- compute.brlen(attr(object, "tree"),
                                   method = "Grafen",
                                   power = exp(corStruct[1]))
            if (class(corStruct)[1] %in% c("corBrownian", "corGrafen")) {
                dis <- cophenetic.phylo(attr(corStruct, "tree"), full = TRUE)
                MRCA <- mrca(attr(corStruct, "tree"), full = TRUE)
                M <- dis["-1", MRCA]
                dim(M) <- rep(sqrt(length(M)), 2)
            }
            varAY <- M[-(1:nb.tip), 1:nb.tip]
            varA <- M[-(1:nb.tip), -(1:nb.tip)]
            V <- corMatrix(Initialize(corStruct, data.frame(x)),
                           corr = FALSE)
            tmp <- solve(V)
            obj$ace <- varAY %*% tmp %*% x
            if (CI) {
                CI95 <- matrix(NA, nb.node, 2)
                se <- sqrt((varA - varAY %*% tmp %*% t(varAY))[cbind(1:nb.node, 1:nb.node)])
                CI95[, 1] <- obj$ace + se * qnorm(0.025)
                CI95[, 2] <- obj$ace - se * qnorm(0.025)
                obj$CI95 <- CI95
            }
        }
    } else { # type == "discrete"
        if (method != "ML")
          stop("only ML estimation is possible for discrete characters.")
        if (!is.factor(x)) x <- factor(x)
        nl <- nlevels(x)
        lvls <- levels(x)
        x <- as.integer(x)
        if (is.character(model)) {
            rate <- matrix(NA, nl, nl)
            if (model == "ER") np <- rate[] <- 1
            if (model == "ARD") {
                np <- nl * (nl - 1)
                rate[col(rate) != row(rate)] <- 1:np
            }
            if (model == "SYM") {
                np <- nl * (nl - 1)/2
                rate[col(rate) < row(rate)] <- 1:np
                rate <- t(rate)
                rate[col(rate) < row(rate)] <- 1:np
            }
        } else {
            rate <- model
            np <- max(rate)
        }
        rate[cbind(1:nl, 1:nl)] <- 0
        rate[rate == 0] <- np + 1

        liks <- matrix(0, nb.tip + nb.node, nl)
        rownames(liks) <- as.character(c(1:nb.tip, -(1:nb.node)))
        for (i in 1:nb.tip) liks[i, x[i]] <- 1

        Q <- matrix(0, nl, nl)
        dev <- function(p, output.liks = FALSE) {
            Q[] <- c(p, 0)[rate]
            diag(Q) <- -rowSums(Q)
            unused <- rep(TRUE, nb.tip + nb.node)
            done <- c(rep(TRUE, nb.tip), rep(FALSE, nb.node))
            names(unused) <- names(done) <- rownames(liks)
            while (sum(unused) > 1) {
                term <- names(which(done & unused))
                ind <- as.logical(match(phy$edge[, 2], term))
                ind[is.na(ind)] <- FALSE
                term.br <- matrix(phy$edge[ind], length(term), 2)
                basal <- names(which(table(term.br[, 1]) == 2))
                for (nod in basal) {
                    pair.ind <- which(phy$edge[, 1] == nod)
                    i <- pair.ind[1]
                    j <- pair.ind[2]
                    pair <- phy$edge[pair.ind, 2]
                    a <- pair[1]
                    b <- pair[2]
                    tmp <- eigen(Q * phy$edge.length[i], symmetric = FALSE)
                    P1 <- tmp$vectors %*% diag(exp(tmp$values)) %*% solve(tmp$vectors)
                    tmp <- eigen(Q * phy$edge.length[j], symmetric = FALSE)
                    P2 <- tmp$vectors %*% diag(exp(tmp$values)) %*% solve(tmp$vectors)
                    liks[nod, ] <- P1 %*% liks[a, ] * P2 %*% liks[b, ]
                    unused[pair] <- FALSE
                }
                done[basal] <- TRUE
            }
            if (output.liks) return(liks[-(1:nb.tip), ])
            - 2 * log(sum(liks["-1", ]))
        }
        out <- nlm(function(p) dev(p),
                   p = rep(ip, length.out = np),
                   hessian = TRUE)
        obj$loglik <- -out$minimum / 2
        obj$rates <- out$estimate
        if (any(out$gradient == 0))
          warning("The likelihood gradient seems flat in at least one dimension (null gradient):\ncannot compute the standard-errors of the transition rates.\n")
        else obj$se <- sqrt(diag(solve(out$hessian)))
        obj$index.matrix <- rate
        diag(obj$index.matrix) <- NA
        if (CI) {
            lik.anc <- dev(obj$rates, TRUE)
            lik.anc <- lik.anc / rowSums(lik.anc)
            colnames(lik.anc) <- lvls
            obj$lik.anc <- lik.anc
        }
    }
    obj$call <- match.call()
    class(obj) <- "ace"
    obj
}

logLik.ace <- function(object, ...) object$loglik

deviance.ace <- function(object, ...) -2*object$loglik

AIC.ace <- function(object, ..., k = 2)
{
    if (is.null(object$loglik)) return(NULL)
    ## Trivial test of "type"; may need to be improved
    ## if other models are included in ace(type = "c")
    np <- if (!is.null(object$sigma2)) 1 else length(object$rates)
    -2*object$loglik + np*k
}

### by BB:
anova.ace <- function(object, ...)
{
    X <- c(list(object), list(...))
    df <- sapply(lapply(X, "[[", "rates"), length)
    ll <- sapply(X, "[[", "loglik")
    ## check if models are in correct order?
    dev <- c(NA, 2*diff(ll))
    ddf <- c(NA, diff(df))
    table <- data.frame(ll, df, ddf, dev,
                        pchisq(dev, ddf, lower.tail = FALSE))
    dimnames(table) <- list(1:length(X), c("Log lik.", "Df",
                                           "Df change", "Deviance",
                                           "Pr(>|Chi|)"))
    structure(table, heading = "Likelihood Ratio Test Table",
              class = c("anova", "data.frame"))
}
