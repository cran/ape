### birthdeath.R  (2003-05-29)
###
###       Estimation of Speciation and Extinction Rates
###                 With Birth-Death Models
### Copyright 2003 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

birthdeath <- function(phy)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    N <- max(as.numeric(phy$edge))
    x <- c(NA, branching.times(phy))
    dev <- function(a, r) {
        -2 * (sum(log((N - 1):1))
              + (N - 2) * log(r)
              + r * sum(x[3:N])
              + N * log(1 - a)
              - 2 * sum(log(exp(r * x[2:N]) - a)))
    }
    out <- nlm(function(p) dev(p[1], p[2]), c(0.1, 0.2), hessian = TRUE)
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
    ## compute the 95 % profile likelihood CIs
    ## (not very clean... but seems to work -- EP 2003-03-29)
    CI <- matrix(NA, 2, 2)
    foo <- function(p) dev(p, para[2]) - 3.84 - Dev
    inc <- 1e-2
    lo <- para[1] - inc
    up <- para[1] + inc
    while (foo(lo) < 0) lo <- lo - inc
    while (foo(up) < 0) up <- up + inc
    CI[1, 1] <- uniroot(foo, lower = lo, upper = para[1])$root
    if (CI[1, 1] < 0) CI[1, 1] <- 0
    CI[1, 2] <- uniroot(foo, lower = para[1], upper = up)$root
    foo <- function(p) dev(para[1], p) - 3.84 - Dev
    lo <- para[2] - inc
    up <- para[2] + inc
    while (foo(lo) < 0) lo <- lo - inc
    while (foo(up) < 0) up <- up + inc
    CI[2, 1] <- uniroot(foo, lower = lo, upper = para[2])$root
    CI[2, 2] <- uniroot(foo, lower = para[2], upper = up)$root
    cat("\nEstimation of Speciation and Extinction Rates\n")
    cat("            With Birth-Death Models\n\n")
    cat("     Phylogenetic tree:", deparse(substitute(phy)), "\n")
    cat("        Number of tips:", N, "\n")
    cat("              Deviance:", Dev, "\n")
    cat("        Log-likelihood:", -Dev/2, "\n")
    cat("   Parameter estimates:\n")
    cat("      d / b =", para[1], "  StdErr =", se[1], "\n")
    cat("      b - d =", para[2], "  StdErr =", se[2], "\n")
    cat("   (b: speciation rate, d: extinction rate)\n")
    cat("   Profile likelihood 95 % confidence intervals:\n")
    cat("      d / b: [", CI[1, 1], ", ", CI[1, 2], "]", "\n", sep = "")
    cat("      b - d: [", CI[2, 1], ", ", CI[2, 2], "]", "\n\n", sep = "")
}
