\name{carnivora}
\docType{data}
\alias{carnivora}
\title{Carnivora body sizes and life history traits}
\description{
  Dataset adapted from Gittleman (1986), including 2 morphological variables (body and brain sizes),  8 life history traits variables and 4 taxonomic variables.
}
\usage{data(carnivora)}
\format{
  A data frame with 112 observations on 17 variables.

  \tabular{rlll}{
    [,1]  \tab Order       \tab factor  \tab Carnivora order \cr
    [,2]  \tab SuperFamily \tab factor  \tab Super family (Caniformia or Feliformia) \cr
    [,3]  \tab Family      \tab factor  \tab Carnivora family \cr
    [,4]  \tab Genus       \tab factor  \tab Carnivora genus \cr
    [,5]  \tab Species     \tab factor  \tab Carnivora species \cr
    [,6]  \tab FW          \tab numeric \tab Female body weight (kg) \cr
    [,7]  \tab SW          \tab numeric \tab Average body weight of adult male and adult female (kg) \cr
    [,8]  \tab FB          \tab numeric \tab Female brain weight (g) \cr
    [,9]  \tab SB          \tab numeric \tab Average brain weight of adult male and adult female (g) \cr
    [,10] \tab LS          \tab numeric \tab Litter size \cr
    [,11] \tab GL          \tab numeric \tab Gestation length (days) \cr
    [,12] \tab BW          \tab numeric \tab Birth weigth (g) \cr
    [,13] \tab WA          \tab numeric \tab Weaning age (days) \cr
    [,14] \tab AI          \tab numeric \tab Age of independance (days) \cr
    [,15] \tab LY          \tab numeric \tab Longevity (months) \cr
    [,16] \tab AM          \tab numeric \tab Age of sexual maturity (days) \cr
    [,17] \tab IB          \tab numeric \tab Inter-birth interval (months) \cr
  }
}
\source{
  Gittleman, J. L. (1986) Carnivore life history patterns: allometric,
  phylogenetic and ecological associations. \emph{American Naturalist},
  \bold{127}: 744--771.
}
\examples{
  data(carnivora);
  # This is figure 1 of Gittleman 1986:
  library(lattice)
  trellis.device(color=FALSE)
  xyplot(BW ~ FW, groups=Family, data=carnivora, auto.key=TRUE, xlog=TRUE,
      scale=list(log=TRUE), ylim=c(1, 2000))
  trellis.device(color=FALSE)
  xyplot(BW ~ FB, groups=Family, data=carnivora, auto.key=TRUE, xlog=TRUE,
      scale=list(log=TRUE), ylim=c(1, 2000))
}
\keyword{datasets}

