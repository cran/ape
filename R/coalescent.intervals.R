### coalescent.intervals.R  (2002-09-12)
###
###     Constructs objects with information on coalescent intervals
###
### Copyright 2002 Korbinian Strimmer <strimmer@stat.uni-muenchen.de>
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


coalescent.intervals <- function(x) UseMethod("coalescent.intervals")

# set up coalescent interval object (from NH tree)
coalescent.intervals.phylo <- function(x)
{
    if (class(x) != "phylo") stop("object \"x\" is not of class \"phylo\"")

    # ensure we have a BINARY tree
    if (!is.binary.tree(x)) stop("object \"x\" is not a binary tree")
    # ordered branching times
    t <- sort(branching.times(x))
    lt <- length(t)

    # interval widths
    w <- numeric(lt)
    w[1] <- t[1]
    for (i in 2:lt) w[i] <- t[i] - t[i - 1]

    l <- (lt+1):2       # number of lineages

    obj <- list(
     lineages=l,
     interval.length=w,
     interval.count=lt,
     total.depth =sum(w))
    class(obj) <- "coalescentIntervals"
    return(obj)
}


# set up coalescent interval object from vector of interval length
coalescent.intervals.default <- function(x)
{
  if (!is.vector(x)) stop("argument \"x\" is not a vector of interval lengths")

  # x = list of the widths of each interval
  lt <- length(x)
  l <- (lt+1):2           # number of lineages at the beginning of each interval

  obj <- list(
     lineages=l,
     interval.length=x,
     interval.count=lt,
     total.depth =sum(x))
    class(obj) <- "coalescentIntervals"
    return(obj)
}
