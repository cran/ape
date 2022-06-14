## C_wrappers.R
## 
## Allow C functions to be called by external packages

#' @useDynLib ape node_depth
#' @keywords internal
#' @export
C_node_depth <- function (nTip, nNode, parent, child, nEdge) {
  .C("node_depth", as.integer(nTip), as.integer(nNode), as.integer(parent), 
     as.integer(child), as.integer(nEdge), double(nTip + nNode))[[6]]
}

#' @useDynLib ape neworder_phylo
C_neworder_phylo <- function (nTaxa, parent, child, nb.edge, whichwise) {
  .C('neworder_phylo', as.integer(nTaxa), as.integer(parent), as.integer(child), 
     as.integer(nb.edge), integer(nb.edge), as.integer(whichwise), NAOK = TRUE)[[5]]
}

#' @useDynLib ape neworder_pruningwise
C_neworder_pruningwise <- function (nTaxa, nb.node, parent, child, nb.edge) {
  .C('neworder_pruningwise', as.integer(nTaxa), as.integer(nb.node), as.integer(parent), 
     as.integer(child), as.integer(nb.edge), integer(nb.edge))[[6]]
}

