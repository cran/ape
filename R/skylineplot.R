### skylineplot.R  (2002-09-12)
###
###     Various methods to plot skyline objects (= skyline plots)
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


# plot skyline 
plot.skyline <- function(x, reverse.time=FALSE, ...)
{
  if (class(x) != "skyline")
    stop("object \"x\" is not of class \"skyline\"")
  t <- x$time
  m <- x$population.size

  if (reverse.time)
    plot(c(0,t),c(m,m[[length(m)]]),type="s",
     xlab="time (present to past)",ylab="effective population size",log="y", ...)
  else
    plot(-c(0,t),c(m,m[[length(m)]]),type="s",
     xlab="time (past to present)",ylab="effective population size",log="y", ...)
}


# plot another skyline plot on top
lines.skyline <- function(x, reverse.time=FALSE, ...)
{
  if (class(x) != "skyline")
    stop("object \"x\" is not of class \"skyline\"")
  t <- x$time
  m <- x$population.size

  if (reverse.time)
    lines(c(0,t),c(m,m[[length(m)]]),type="s", ...)
  else
    lines(-c(0,t),c(m,m[[length(m)]]),type="s", ...)
}


# convenience short cut (almost compatible with APE 0.1)
skylineplot <- function(z, ...) plot(skyline(z, ...))


#input: phylogenetic tree
skylineplot.deluxe <- function(tree, ...)
{
  if (class(tree) != "phylo")
    stop("object \"tree\" is not of class \"phylo\"")

  ci <- coalescent.intervals(tree)
  classic <- skyline(ci)
  generalized <- skyline(ci, -1)
  plot(classic,col=grey(.8), ...)
  lines(generalized, ...)
  return(generalized)
}



