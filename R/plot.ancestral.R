## plot.ancestral.R (2005-12-04)

##   Plotting Ancestral Characters on a Tree

## Copyright 2005 Julien Dutheil

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

plot.ancestral <- function(x, which=names(x$node.character),
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
