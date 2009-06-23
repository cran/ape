## ace.R (2009-06-10)

##   Ancestral Character Estimation

## Copyright 2005-2009 Emmanuel Paradis and Ben Bolker

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

ace <- function(x, phy, type = "continuous", method = "ML", CI = TRUE,
                model = if (type == "continuous") "BM" else "ER",
                scaled = TRUE, kappa = 1, corStruct = NULL, ip = 0.1)
{
    if (!inherits(phy, "phylo"))
      stop('object "phy" is not of class "phylo".')
    if (is.null(phy$edge.length))
        stop("tree has no branch lengths")
    type <- match.arg(type, c("continuous", "discrete"))
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    if (nb.node != nb.tip - 1)
      stop('"phy" is not rooted AND fully dichotomous.')
    if (length(x) != nb.tip)
      stop("length of phenotypic and of phylogenetic data do not match.")
    if (!is.null(names(x))) {
        if(all(names(x) %in% phy$tip.label))
          x <- x[phy$tip.label]
        else warning('the names of argument "x" and the tip labels of the tree
did not match: the former were ignored in the analysis.')
    }
    obj <- list()
    if (kappa != 1) phy$edge.length <- phy$edge.length^kappa
    if (type == "continuous") {
        if (method == "pic") {
            if (model != "BM")
              stop('the "pic" method can be used only with model = "BM".')
            ## See pic.R for some annotations.
            phy <- reorder(phy, "pruningwise")
            phenotype <- numeric(nb.tip + nb.node)
            phenotype[1:nb.tip] <- if (is.null(names(x))) x else x[phy$tip.label]
            contr <- var.con <- numeric(nb.node)
            ans <- .C("pic", as.integer(nb.tip), as.integer(nb.node),
                      as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]),
                      as.double(phy$edge.length), as.double(phenotype),
                      as.double(contr), as.double(var.con),
                      as.integer(CI), as.integer(scaled),
                      PACKAGE = "ape")
            obj$ace <- ans[[6]][-(1:nb.tip)]
            names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
            if (CI) {
                se <- sqrt(ans[[8]])
                CI95 <- matrix(NA, nb.node, 2)
                CI95[, 1] <- obj$ace + se * qnorm(0.025)
                CI95[, 2] <- obj$ace - se * qnorm(0.025)
                obj$CI95 <- CI95
            }
        }
        if (method == "ML") {
            if (model == "BM") {
                tip <- phy$edge[, 2] <= nb.tip
                dev.BM <- function(p) {
                    x1 <- p[-1][phy$edge[, 1] - nb.tip]
                    x2 <- numeric(length(x1))
                    x2[tip] <- x[phy$edge[tip, 2]]
                    x2[!tip] <- p[-1][phy$edge[!tip, 2] - nb.tip]
                    -2 * (-sum((x1 - x2)^2/phy$edge.length)/(2*p[1]) -
                          nb.node * log(p[1]))
                }
                out <- nlm(function(p) dev.BM(p),
                           p = c(1, rep(mean(x), nb.node)), hessian = TRUE)
                obj$loglik <- -out$minimum / 2
                obj$ace <- out$estimate[-1]
                names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
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
              M <- corStruct[1] * dist.nodes(phy)
            if (class(corStruct)[1] == "corGrafen")
              phy <- compute.brlen(attr(corStruct, "tree"),
                                   method = "Grafen",
                                   power = exp(corStruct[1]))
            if (class(corStruct)[1] %in% c("corBrownian", "corGrafen")) {
                dis <- dist.nodes(attr(corStruct, "tree"))
                MRCA <- mrca(attr(corStruct, "tree"), full = TRUE)
                M <- dis[as.character(nb.tip + 1), MRCA]
                dim(M) <- rep(sqrt(length(M)), 2)
            }
            varAY <- M[-(1:nb.tip), 1:nb.tip]
            varA <- M[-(1:nb.tip), -(1:nb.tip)]
            V <- corMatrix(Initialize(corStruct, data.frame(x)),
                           corr = FALSE)
            invV <- solve(V)
            o <- gls(x ~ 1, data.frame(x), correlation = corStruct)
            GM <- o$coefficients
            obj$ace <- drop(varAY %*% invV %*% (x - GM) + GM)
            names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
            if (CI) {
                CI95 <- matrix(NA, nb.node, 2)
                se <- sqrt((varA - varAY %*% invV %*% t(varAY))[cbind(1:nb.node, 1:nb.node)])
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
                np <- nl*(nl - 1)
                rate[col(rate) != row(rate)] <- 1:np
            }
            if (model == "SYM") {
                np <- nl * (nl - 1)/2
                rate[col(rate) < row(rate)] <- 1:np
                rate <- t(rate)
                rate[col(rate) < row(rate)] <- 1:np
            }
        } else {
            if (ncol(model) != nrow(model))
              stop("the matrix given as `model' is not square")
            if (ncol(model) != nl)
              stop("the matrix `model' must have as many rows
as the number of categories in `x'")
            rate <- model
            np <- max(rate)
        }
        index.matrix <- rate
        index.matrix[cbind(1:nl, 1:nl)] <- NA
        rate[cbind(1:nl, 1:nl)] <- 0
        rate[rate == 0] <- np + 1 # to avoid 0's since we will use this an numeric indexing

        liks <- matrix(0, nb.tip + nb.node, nl)
        TIPS <- 1:nb.tip
        liks[cbind(TIPS, x)] <- 1
        phy <- reorder(phy, "pruningwise")

        Q <- matrix(0, nl, nl)
        ## from Rich FitzJohn:
        comp <- numeric(nb.tip + nb.node) # Storage...
        dev <- function(p, output.liks = FALSE) {
            Q[] <- c(p, 0)[rate]
            diag(Q) <- -rowSums(Q)
            for (i  in seq(from = 1, by = 2, length.out = nb.node)) {
                j <- i + 1
                anc <- phy$edge[i, 1]
                des1 <- phy$edge[i, 2]
                des2 <- phy$edge[j, 2]
                v.l <- matexpo(Q * phy$edge.length[i]) %*% liks[des1, ]
                v.r <- matexpo(Q * phy$edge.length[j]) %*% liks[des2, ]
                v <- v.l * v.r
                comp[anc] <- sum(v)
                liks[anc, ] <- v/comp[anc]
            }
            if (output.liks) return(liks[-TIPS, ])
            -2 * sum(log(comp[-TIPS]))
        }
        out <- nlminb(rep(ip, length.out = np), function(p) dev(p),
                      lower = rep(0, np), upper = rep(Inf, np))
        obj$loglik <- -out$objective/2
        obj$rates <- out$par
        oldwarn <- options("warn")
        options(warn = -1)
        h <- nlm(function(p) dev(p), p = obj$rates, iterlim = 1,
                 stepmax = 0, hessian = TRUE)$hessian
        options(oldwarn)
        if (any(h == 0))
          warning("The likelihood gradient seems flat in at least one dimension (gradient null):\ncannot compute the standard-errors of the transition rates.\n")
        else obj$se <- sqrt(diag(solve(h)))
        obj$index.matrix <- index.matrix
        if (CI) {
            obj$lik.anc <- dev(obj$rates, TRUE)
            colnames(obj$lik.anc) <- lvls
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
