### plot.ancestral.R  (2005-12-04)
###
###     Plotting Ancestral Characters on a Tree
###
### Copyright 2005 Julien Dutheil <julien.dutheil@univ-montp2.fr>
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

plot.ancestral <- function(x,
    which=names(x$node.character),
    n.col=10, col.fun=function(n) rainbow(n, start=0.4, end=0),
    plot.node.values=FALSE,
    ask = prod(par("mfcol")) < length(which) && dev.interactive(),
    ...)
{
  if (!("ancestral" %in% class(x)))
      stop("object \"phy\" is not of class \"ancestral\"")
  states <- rbind(x$node.character, x$tip.character)
  cols <- col.fun(n.col)
  if(ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  for(state in which) {
    a <- states[x$edge[,2],state]
    b <- round((n.col-1)*(a-min(a))/(max(a)-min(a)))+1
    if(plot.node.values) {
      x$node.label <- x$node.character[,state]
      plot.phylo(x, edge.color=cols[b], show.node.label=TRUE, sub=state, ...)
    } else {
      plot.phylo(x, edge.color=cols[b], sub=state, ...)
    }
  }
}
