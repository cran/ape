### bd.ext.R  (2004-08-20)
###
###     Extended Version of the Birth-Death Models to
###       Estimate Speciation and Extinction Rates
###
### Copyright 2004 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

bd.ext <- function(phy, S)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    if (!is.null(names(S))) S <- S[phy$tip.label]
    N <- length(S)
    x <- branching.times(phy)
    x <- c(x[1], x)
    trm.br <- phy$edge.length[phy$edge[, 2] > 0]
    dev <- function(a, r)
    {
        -2 * (sum(log((N - 1):1))
              + (N - 2) * log(r)
              + (3 * N) * log(1 - a)
              + 2 * r * sum(x[2:N])
              - 2 * sum(log(exp(r * x[2:N]) - a))
              + r * sum(trm.br)
              + sum((S - 1) * log(exp(r * trm.br) - 1))
              - sum((S + 1) * log(exp(r * trm.br) - a)))
    }
    out <- nlm(function(p) dev(p[1], p[2]), c(0, 0.2), hessian = TRUE)
    if (out$estimate[1] < 0) {
        out <- nlm(function(p) dev(0, p), 0.2, hessian = TRUE)
        para <- c(0, out$estimate)
        se <- c(0, sqrt(diag(solve(out$hessian))))
    }
    else {
        para <- out$estimate
        se <- sqrt(diag(solve(out$hessian)))
    }
    Dev <- out$minimum
###     Extended Version of the Birth-Death Models to
###       Estimate Speciation and Extinction Rates

    cat("\nExtended Version of the Birth-Death Models to\n")
    cat("    Estimate Speciation and Extinction Rates\n\n")
    cat("    Data: phylogenetic:", deparse(substitute(phy)), "\n")
    cat("             taxonomic:", deparse(substitute(S)), "\n")
    cat("        Number of tips:", N, "\n")
    cat("              Deviance:", Dev, "\n")
    cat("        Log-likelihood:", -Dev/2, "\n")
    cat("   Parameter estimates:\n")
    cat("      d / b =", para[1], "  StdErr =", se[1], "\n")
    cat("      b - d =", para[2], "  StdErr =", se[2], "\n")
    cat("   (b: speciation rate, d: extinction rate)\n")
}
