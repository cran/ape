### collapsed.intervals.R  (2002-09-12)
###
###     Collapsed coalescent intervals (e.g. for the skyline plot)
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

# construct collapsed intervals from coalescent intervals
collapsed.intervals <- function(ci, epsilon=0.0)
{
  if (class(ci) != "coalescentIntervals")
    stop("object \"ci\" is not of class \"coalescentIntervals\"")

  sz <- ci$interval.length
  lsz <- length(sz)
  idx <- c <- 1:lsz

  p <- 1
  w <- 0
  
  # starting from tips collapes intervals
  # until total size is >= epsilon
  for (i in 1:lsz)
  {
    idx[[i]] <- p
    w <- w + sz[[i]]
    if (w >= epsilon)
    {
      p <- p+1
      w <- 0
    }
  }

  # if last interval is smaller than epsilon merge
  # with second last interval
  lastInterval <- idx==p
  if ( sum(sz[lastInterval]) < epsilon )
  {
    p <- p-1
    idx[lastInterval] <- p
  }

  obj <- list(
     lineages=ci$lineages,
     interval.length=ci$interval.length,
     collapsed.interval=idx, # collapsed intervals (via reference)
     interval.count=ci$interval.count,
     collapsed.interval.count = idx[[ci$interval.count]],
     total.depth =ci$total.depth,
     epsilon = epsilon
    )
  class(obj) <- "collapsedIntervals"
  
  return(obj)
}
