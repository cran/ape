\name{drop.tip}
\alias{drop.tip}
\title{Remove Tips in a Phylogenetic Tree}
\usage{
drop.tip(phy, tip, trim.internal = TRUE)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{tip}{a vector of mode numeric or character specifying the tips
    to delete.}
  \item{trim.internal}{a logical specifying whether to delete the
    corresponding internal branches.}
}
\description{
  This function removes the terminal branches of a phylogenetic tree,
  possibly removing the corresponding internal branches.
}
\details{
  The argument \code{tip} can be either character or numeric. In the
  first case, it gives the labels of the tips to be deleted; in the
  second case the numbers of these labels in the vector
  \code{phy$tip.label} are given.

  If \code{trim.internal = FALSE}, the new tips are given \code{"NA"} as
  labels, unless there are node labels in the tree in which case they
  are used.
}
\value{
  an object of class \code{"phylo"}.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{bind.tree}}, \code{\link{root}}, \code{\link{rotate}}
}
\examples{
data(bird.families)
tip <- c(
"Eopsaltriidae", "Acanthisittidae", "Pittidae", "Eurylaimidae",
"Philepittidae", "Tyrannidae", "Thamnophilidae", "Furnariidae",
"Formicariidae", "Conopophagidae", "Rhinocryptidae", "Climacteridae",
"Menuridae", "Ptilonorhynchidae", "Maluridae", "Meliphagidae",
"Pardalotidae", "Petroicidae", "Irenidae", "Orthonychidae",
"Pomatostomidae", "Laniidae", "Vireonidae", "Corvidae",
"Callaeatidae", "Picathartidae", "Bombycillidae", "Cinclidae",
"Muscicapidae", "Sturnidae", "Sittidae", "Certhiidae",
"Paridae", "Aegithalidae", "Hirundinidae", "Regulidae",
"Pycnonotidae", "Hypocoliidae", "Cisticolidae", "Zosteropidae",
"Sylviidae", "Alaudidae", "Nectariniidae", "Melanocharitidae",
"Paramythiidae","Passeridae", "Fringillidae")
plot(drop.tip(bird.families, tip))
plot(drop.tip(bird.families, tip, trim.internal = FALSE))
}
\keyword{manip}
