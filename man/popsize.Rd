\name{popsize}
\alias{mcmc.popsize}
\alias{extract.popsize}
\alias{plot.popsize}
\alias{lines.popsize}

\title{Reversible Jump MCMC to Infer Demographic History}
\usage{
mcmc.popsize(tree,nstep, thinning=1, burn.in=0,prior.lambda=0.5, max.nodes=30,progress.bar=TRUE)
extract.popsize(mcmc.out, credible.interval=0.95, time.points=200, thinning=1, burn.in=0)
plot.popsize(x, show.median=TRUE, show.years=FALSE, subst.rate, present.year, ...)
lines.popsize(x, show.median=TRUE,show.years=FALSE, subst.rate, present.year, ...)

}
\arguments{
  \item{tree}{Either an ultrametric tree (i.e. an object of class \code{"phylo"}), 
           or coalescent intervals (i.e. an object of class \code{"coalescentIntervals"}) }
  
  \item{nstep}{Number of MCMC steps, i.e. length of the Markov chain (suggested value: 10,000-50,000)}
  
  \item{thinning}{Thinning factor (suggest value: 10-100).}
  
  \item{burn.in}{Number of steps dropped from the chain to allow for a burn-in phase (suggest value: 1000)}
  
  \item{prior.lambda}{The lambda parameter of the prior distribution for the number of internal nodes of the approximating
  spline for the demographic function (suggested value: 0.1-1.0).}
  
  \item{max.nodes}{Upper limit for the number of internal nodes of the approximating spline (default: 30).}
  
  \item{progress.bar}{Show progress bar during the MCMC run.}	  
  
  \item{mcmc.out}{output from \code{mcmc.popsize} - this is needed as input for \code{extract.popsize}}
  
 \item{credible.interval}{"width" of the confidence band (default: 0.95)}

 \item{time.points}{number of discrete time points in the table output by \code{extract.popsize}}

 \item{x}{table with population size versus time, as computed by \code{extract.popsize} }

 \item{show.median}{plot median rather than mean as point estimate for demographic function (default: TRUE)}

  \item{show.years}{option that determines whether the time is plotted in units of
        of substitutions (default) or in years (requires specification of substution rate
	and year of present).}
	
\item{subst.rate}{substitution rate (see option show.years).}	
\item{present.year}{present year (see option show.years).}	

  \item{\dots}{further arguments to be passed on  to \code{plot}.}  
}
\description{
 These functions implement a reversible jump MCMC framework to infer the demographic history,
 as well as corresponding confidence bands,
 from a genealogical tree. The computed demographic history is a continous
 and smooth function in time.
 
 
 \code{mcmc.popsize} runs the actual MCMC chain and outputs information about the
 sampling steps, \code{extract.popsize} generates from this MCMC
 output a table of population size in time, and  \code{plot.popsize} and \code{lines.popsize}
 provide utility functions to plot the corresponding demographic functions.
}

\details{
 Please refer to Opgen-Rhein et al. (2004) for methodological details, and the help page of
 \code{\link{skyline}} for information on a related approach.
}


\author{Rainer Opgen-Rhein (\url{http://www.stat.uni-muenchen.de/~opgen/}) and 
        Korbinian Strimmer (\url{http://www.stat.uni-muenchen.de/~strimmer/}).
        Parts of the rjMCMC sampling procedure are adapted from R code by Karl Browman
	 (\url{http://www.biostat.jhsph.edu/~kbroman/})}

\seealso{
\code{\link{skyline}} and \code{\link{skylineplot}}. }



\references{
  Opgen-Rhein, R., L. Fahrmeir, and K. Strimmer, K. 2003. Inference of demographic history
  from genealogical trees using reversible jump Markov chain Monte Carlo. 
  Submitted to \emph{BMC Evolutionary Biology}.
}


\examples{
library(ape)

# get tree
data("hivtree.newick") # example tree in NH format
tree.hiv <- read.tree(text = hivtree.newick) # load tree

# run mcmc chain
mcmc.out <- mcmc.popsize(tree.hiv, nstep=100, thinning=1, burn.in=0) # toy run
#mcmc.out <- mcmc.popsize(tree.hiv, nstep=10000, thinning=5, burn.in=500) # remove comments!!

# make list of population size versus time
popsize  <- extract.popsize(mcmc.out)

# plot and compare with skyline plot
sk <- skyline(tree.hiv)
plot(sk, lwd=1, lty=3, show.years=TRUE, subst.rate=0.0023, present.year = 1997)
lines(popsize, show.years=TRUE, subst.rate=0.0023, present.year = 1997)
}
\keyword{manip}
