## scales.R (2009-07-23)

##   Add a Scale Bar or Axis to a Phylogeny Plot

## add.scale.bar: add a scale bar to a phylogeny plot
## axisPhylo: add a scale axis on the side of a phylogeny plot

## Copyright 2002-2009 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

add.scale.bar <- function(x, y, length = NULL, ...)
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (is.null(length)) {
        nb.digit <-
          if (lastPP$direction %in% c("rightwards", "leftwards")) diff(range(lastPP$xx))
          else diff(range(lastPP$yy))
        nb.digit <- ceiling(log10(nb.digit)) - 2
        length <- eval(parse(text = paste("1e", nb.digit, sep = "")))
    }
    if (missing(x) || missing(y)) {
        switch(lastPP$direction,
               "rightwards" = {
                   x <- 0
                   y <- 1
                   segments(x, y, x + length, y)
                   text(x + length * 1.1, y, as.character(length), adj = c(0, 0.5), ...)
               },
               "leftwards" = {
                   x <- max(lastPP$xx)
                   y <- 1
                   segments(x - length, y, x, y)
                   text(x - length * 1.1, y, as.character(length), adj = c(1, 0.5), ...)
               },
               "upwards" = {
                   x <- max(lastPP$xx)
                   y <- 0
                   segments(x, y, x, y + length)
                   text(x, y + length * 1.1, as.character(length), adj = c(0, 0.5), srt = 90, ...)
               },
               "downwards" = {
                   x <- 1
                   y <- max(lastPP$yy)
                   segments(x, y - length, x, y)
                   text(x, y - length * 1.1, as.character(length), adj = c(0, 0.5), srt = 270, ...)
               }
         )
    }
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
