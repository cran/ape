### plot.popsize.R  (2004-07-4)
###
###     Plot population size in dependence of time
###
### Copyright 2004 Rainer Opgen-Rhein <opgen@stat.uni-muenchen.de> and
###                Korbinian Strimmer <strimmer@stat.uni-muenchen.de> 
###
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




plot.popsize <- function(x, show.median=TRUE,
    show.years=FALSE, subst.rate, present.year, ...)
{
  if (class(x) != "popsize")
    stop("object \"x\" is not of class \"popsize\"")

  ylim <- c(min(popsize[,2:5]),max(popsize[,2:5]))
  if (show.years)
  {
    x1 <- -x[,1]/subst.rate+present.year
    xlab <- "time (years)"
    xlim <- c(min(x1),max(x1))
  }
  else 
  {
    x1 <- x[,1]
    xlab <- "time (past to present in units of substitutions)"
    xlim <- c(max(x1),min(x1))
  }

  if (show.median)
    plot(x1,x[,3],type="s", xlim=xlim, ylim=ylim, xlab=xlab,ylab="effective population size",log="y", lwd=2.5, ...) #median
  else
    plot(x1,x[,2],type="s", xlim=xlim, ylim=ylim, xlab=xlab,ylab="effective population size",log="y", lwd=2.5, ...) #median 

  lines(x1,x[,4], ...)
  lines(x1,x[,5], ...)
}



lines.popsize <- function(x, show.median=TRUE,
    show.years=FALSE, subst.rate, present.year, ...)
{
  if (class(x) != "popsize")
    stop("object \"x\" is not of class \"popsize\"")
  
  
  if (show.years)
  {
    x1 <- -x[,1]/subst.rate+present.year
  }
  else 
  {
    x1 <- x[,1]
  }


  if (show.median)
    lines(x1,x[,3], lwd=2.5, ...) #median
  else
    lines(x1,x[,2], lwd=2.5, ...) #median 

    
  lines(x1,x[,4], ...)
  lines(x1,x[,5], ...)
}

