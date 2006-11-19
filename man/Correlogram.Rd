\name{Correlogram}
\alias{Correlogram}
\alias{correlogram.formula}
\alias{correlogram.phylo}
\title{Compute a correlogram}
\usage{
  correlogram.formula(formula, data = NULL, use="all.obs")
  correlogram.phylo(x, phy, nclass = NULL, breaks = NULL)
}
\arguments{
  \item{x}{a vector containing the data to plot}
  \item{formula}{a formula of the kind \code{y1+..+yn ~ g1/../gn}, where \code{y} are the data to plot and all \code{g} are the nested levels.}
  \item{data}{a dataframe containing all variables in the formula.
		If NULL, variables are taken from the environment.}
	\item{use}{an optional character string giving a method for computing 
		Moran's indices in the presence of missing values.  This must be
    one of the strings '"all.obs"', '"complete.obs"' or '"pairwise.complete.obs"'.
	  In the first case, 'NA' values produce an error.
	  If the seond option is used, all lines including 'NA' will be ignored.
	  In the later case, 'NA' will be removed separately for each correlogram.
	  If only one 'y' is given in the formula, the two last options will send the same result.
	  Finally, unless 'use' is set to '"all.obs'", lines width 'NA' values in grouping variables will be removed.}
  \item{phy}{an object of class "phylo" width edge lengths}
  \item{breaks,nclass}{The class bounds or the number of classes to use. If both are NULL, compute Moran's I on the whole distance matrix.}
}
\description{
  Compute a correlogram from taxonomic variables or from a phylogenetic tree with branch lengths.
  The first method relies on the \code{\link{weight.taxo}} function, and the second relies on the \code{\link{discrete.dist}} function.
  Both methods send an object of class \code{correlogram} to be plotted by
  the \code{\link{plot.correlogram}} method.
  For the \code{correlogram.formula} function, if several \code{y} are specified, an object of class \code{correlogramList}
  (a list of \code{correlogram} objects) is sent.
}
\details{
  See example of the \code{\link{Moran.I}} function to see how the \code{correlogram.formula} function works.
  To deal with phylogenies, the \code{correlogram.phylo} function creates classes according to distances intervals.
  Such intervals may be specified using the \code{breaks} argument or by giving a number of classes (\code{nclass} argument).
}
\value{
  An object of class \code{correlogram}, containing:
  \item{obs}{all measured Moran's I}
  \item{p.values}{the p-values of each I}
  \item{labels}{the names of each level}
  or an object of class \code{correlogramList} containing a list of correlogram objects. 
}
\section{Warning}{
  \code{correlogram.phylo} will return NAs if void classes are used.
  This may happen if \code{breaks} if not properly defined, or sometimes with the \code{nclass=argument}, depending on the tree used.
  Usually, you'll have to pull classes.
}
\author{Julien Dutheil \email{julien.dutheil@univ-montp2.fr}}
\seealso{
  \code{\link{plot.correlogram}, \link{plot.correlogramList}}
}
\examples{
  data(carnivora)
  ### Using the formula interface:
  co <- correlogram.formula(
			log10(SW) ~ Order/SuperFamily/Family/Genus,
			data=carnivora)
  co
  plot(co)
	### Several correlograms on the same plot:
  cos <- correlogram.formula(
			log10(SW) + log10(FW) ~ Order/SuperFamily/Family/Genus,
			data=carnivora)
  names(cos)
  plot(cos)

  ### Using the phylo interface:
  ### (the same analysis than in help(pic)...)

  cat("((((Homo:0.21,Pongo:0.21):0.28,",
     "Macaca:0.49):0.13,Ateles:0.62):0.38,Galago:1.00);",
     file = "ex.tre", sep = "\n")
  tree.primates <- read.tree("ex.tre")
  X <- c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968)
  Y <- c(4.74493, 3.33220, 3.36730, 2.89037, 2.30259)
  ### Since this is a small tree, 2 classes is a reasonable number:
  coX <- correlogram.phylo(X, tree.primates, nclass=2)
  coY <- correlogram.phylo(Y, tree.primates, nclass=2)
  plot(coX)
  plot(coY)
  ### Nothing significant...
  ### Computing Moran's I on the whole matrix:
  coX2 <- correlogram.phylo(X, tree.primates); coX2
  ### Significant at the 5% level
  coY2 <- correlogram.phylo(Y, tree.primates); coY2
  ### Not significant
  unlink("ex.tre") # delete the file "ex.tre"	
}
\keyword{regression}
\keyword{dplot}
