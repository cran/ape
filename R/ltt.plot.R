### ltt.plot.R  (2002-08-28)
###
###     Lineages Through Time Plot
###
### Copyright 2002 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

ltt.plot <- function(phy, ...)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    time <- sort(branching.times(phy))
    N <- 1:(length(time) + 1)
    plot(-c(rev(time), 0), N, xlab = "Time", ylab = "N",
         xaxs = "r", yaxs = "r", type = "S", ...)
}
