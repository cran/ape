### skyline.plot.R  (2002-08-28)
###
###     Skyline Plot of Estimated Effective Population Size
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

skyline.plot <- function(phy, ...)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    bt <- sort(branching.times(phy))
    g <- numeric(length(bt))
    g[1] <- bt[1]
    for (i in 2:length(bt)) g[i] <- bt[i] - bt[i - 1]
    g <- rev(g) # internode intervals are now from past to present
    M <- numeric(length(bt))
    for (i in 1:length(bt)) M[i] <- i * (i + 1) * g[i] / 2
    plot(-c(rev(bt), 0), c(M[1], M), xlab = "Time", ylab = expression(italic(M[i])),
         log = "y", xaxs = "r", yaxs = "r", type = "S", ...)
}
