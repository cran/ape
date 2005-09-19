### mantel.test.R (2005-06-13)
###
###     Mantel Test for Similarity of Two Matrices
###
### Copyright 2002-2005 Ben Bolker <bolker@zoo.ufl.edu>,
###      and Julien Claude <claude@isem.univ-montp2.fr>
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

perm.rowscols <- function(m1, n)
{
    s <- sample(1:n)
    m1[s, s]
}

### calculate the Mantel z-statistic for two square matrices m1 and m2
mant.zstat <- function(m1, m2) sum(lower.triang(m1 * m2))

lower.triang <- function(m)
{
    d <- dim(m)
    if (d[1] != d[2]) print("Warning: non-square matrix")
    m[col(m) <= row(m)]
}

mantel.test <- function (m1, m2, nperm = 1000, graph = FALSE, ...)
{
    n <- nrow(m1)
    realz <- mant.zstat(m1, m2)
    nullstats <- replicate(nperm, mant.zstat(m1, perm.rowscols(m2, n)))
    pval <- (sum(nullstats > realz) - 1)/nperm
    if (graph) {
        plot(density(nullstats), type = "l", ...)
        abline(v = realz)
    }
    list(z.stat = realz, p = pval)
}
