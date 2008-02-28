## zzz.R (2008-02-08)

##   Library Loading

## Copyright 2003-2008 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.PlotPhyloEnv <- new.env()

.First.lib <- function(lib, pkg) {
    require(nlme, quietly = TRUE)
    library.dynam("ape", pkg, lib)
}
