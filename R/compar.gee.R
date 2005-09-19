### compar.gee.R (2004-10-12)
###
###     Comparative Analysis with GEEs
###
### Copyright 2002-2004 Emmanuel Paradis
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

compar.gee <- function(formula, data = NULL, family = "gaussian", phy,
                       scale.fix = FALSE, scale.value = 1)
{
    if (is.null(data)) data <- parent.frame() else {
        if(!any(is.na(match(rownames(data), phy$tip.label))))
          data <- data[phy$tip.label, ]
        else warning("the rownames of the data.frame and the names of the tip labels
do not match: the former were ignored in the analysis.")
    }
    if (is.null(phy$edge.length))
      stop("the tree has no branch lengths.")
    R <- vcv.phylo(phy, cor = TRUE)
    id <- rep(1, dim(R)[1])
    geemod <- do.call("gee", list(formula, id, data = data, family = family, R = R,
                                  corstr = "fixed", scale.fix = scale.fix,
                                  scale.value = scale.value))
    W <- geemod$naive.variance
    if (family == "binomial")
      W <- summary(glm(formula, family = quasibinomial, data = data))$cov.scaled
    N <- geemod$nobs
    nb.node <- -min(as.numeric(phy$edge))
    ## xx: vecteur donnant la distance d'un noeud ou tip à partir de la racine
    xx <- as.numeric(rep(NA, N + nb.node))
    names(xx) <- as.character(c(-(1:nb.node), 1:N))
    xx["-1"] <- 0
    for (i in 2:length(xx)) {
        nod <- names(xx[i])
        ind <- which(phy$edge[, 2] == nod)
        base <- phy$edge[ind, 1]
        xx[i] <- xx[base] + phy$edge.length[ind]
    }
    dfP <- sum(phy$edge.length) * N / sum(xx[as.character(1:N)])
    obj <- list(call = geemod$call,
                nobs = N,
                coefficients = geemod$coefficients,
                residuals = geemod$residuals,
                family = geemod$family$family,
                link = geemod$family$link,
                scale = geemod$scale,
                W = W,
                dfP = dfP)
    class(obj) <- "compar.gee"
    obj
}

print.compar.gee <- function(x, ...)
{
    nas <- is.na(x$coef)
    coef <- x$coef[!nas]
    cnames <- names(coef)
    coef <- matrix(rep(coef, 4), ncol = 4)
    dimnames(coef) <- list(cnames,
                           c("Estimate", "S.E.", "t", "Pr(T > |t|)"))
    df <- x$dfP - dim(coef)[1]
    coef[, 2] <- sqrt(diag(x$W))
    coef[, 3] <- coef[, 1]/coef[, 2]
    coef[, 4] <- 2 * (1 -  pt(abs(coef[, 3]), df))
    residu <- quantile(as.vector(x$residuals))
    names(residu) <- c("Min", "1Q", "Median", "3Q", "Max")
    cat("\nCall:\n")
    cat("  formula: ")
    print(x$call$formula)
    cat("\nNumber of observations: ", x$nobs, "\n")
    cat("\nModel:\n")
    cat(" Link:                     ", x$link, "\n")
    cat(" Variance to Mean Relation:", x$family, "\n")
    cat("\nSummary of Residuals:\n")
    print(residu)
    if (any(nas))
        cat("\n\nCoefficients: (", sum(nas), " not defined because of singularities)\n",
            sep = "")
    else cat("\n\nCoefficients:\n")
    print(coef)
    cat("\nEstimated Scale Parameter: ", x$scale)
    cat("\n\"Phylogenetic\" df (dfP): ", x$dfP, "\n")
}
