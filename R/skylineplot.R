### skylineplot.R  (2004-07-4)
###
###     Various methods to plot skyline objects (= skyline plots)
###
### Copyright 2002-2004 Korbinian Strimmer <strimmer@stat.uni-muenchen.de>
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
plot.skyline <- function(x, show.years=FALSE, subst.rate, present.year, ...)
{
  if (class(x) != "skyline")
    stop("object \"x\" is not of class \"skyline\"")
  t <- x$time
  m <- x$population.size
  lm <- length(m)

  if (show.years)
  {
    plot((-c(0,t))/subst.rate+present.year,c(m,m[lm]),type="s",
     xlab="time (years)",ylab="effective population size",log="y", ...)
  
  }
  else
  {
    plot(c(0,t),c(m,m[lm]),type="s", xlim=c(t[lm],0),
     xlab="time (past to present in units of substitutions)",ylab="effective population size",log="y", ...)
  }

}

# plot another skyline plot on top
lines.skyline <- function(x, show.years=FALSE, subst.rate, present.year, ...)
{
  if (class(x) != "skyline")
    stop("object \"x\" is not of class \"skyline\"")
  t <- x$time
  m <- x$population.size
  lm <- length(m)


  if (show.years)
  {
    lines((-c(0,t))/subst.rate+present.year,c(m,m[lm]),type="s", ...)
  }
  else
  {
    lines(c(0,t),c(m,m[lm]),type="s", ...)
  }
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



