### gammaStat.R (2002-08-28)
###
###     Gamma-Statistic of Pybus and Harvey
###
### Copyright 2002 Emmanuel Paradis
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

gammaStat <- function(phy)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    N <- max(as.numeric(phy$edge))
    bt <- sort(branching.times(phy))
    g <- rev(c(bt[1], diff(bt))) # internode intervals are from past to present
    ST <- sum((2:N) * g)
    stat <- sum(cumsum((2:(N - 1)) * g[-(N - 1)])) / (N - 2)
    m <- ST / 2
    s <- ST * sqrt(1 / (12 * (N - 2)))
    (stat - m) / s
}
