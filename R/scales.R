## scales.R (2008-02-28)

##   Add a Scale Bar or Axis to a Phylogeny Plot

## add.scale.bar: add a scale bar to a phylogeny plot
## axisPhylo: add a scale axis on the side of a phylogeny plot

## Copyright 2002-2008 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

add.scale.bar <- function(x = 0, y = 1, length = NULL, ...)
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (is.null(length)) {
        nb.digit <- ceiling(log10(mean(lastPP$xx))) - 2
        length <- eval(parse(text = paste("1e", nb.digit, sep = "")))
    }
    segments(x, y, x + length, y)
    text(x + length * 1.1, y, as.character(length), adj = c(0, 0.5), ...)
}

axisPhylo <- function(side = 1, ...)
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (lastPP$type %in% c("phylogram", "cladogram")) {
        if (lastPP$direction %in% c("rightwards", "leftwards")) {
            x <- pretty(lastPP$xx)
            if (lastPP$direction == "rightwards") maxi <- max(lastPP$xx)
            else {
                maxi <- min(lastPP$xx)
                x <- -x
            }
        } else {
            x <- pretty(lastPP$yy)
            if (lastPP$direction == "upwards") maxi <- max(lastPP$yy)
            else {
                maxi <- min(lastPP$yy)
                x <- -x
            }
        }
    }
    axis(side = side, at = c(maxi - x), labels = abs(x), ...)
}
