## zzz.R (2003-05-05)

##   Library Loading

## Copyright 2003 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.First.lib <- function(lib, pkg) {
    require(gee)
    require(nlme)
    require(lattice)
    library.dynam("ape", pkg, lib)
}
