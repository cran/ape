### read.GenBank.R  (2005-05-31)
###
###     Read DNA Sequences from GenBank via Internet
###
### Copyright 2005 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
###
### This file is part of the `ape' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

read.GenBank <- function(access.nb, seq.names = access.nb,
                         species.names = TRUE)
{
    N <- length(access.nb)
    ## If there are more than 400 sequences, we need to break down the
    ## requests, otherwise there is a segmentation fault.
    nrequest <- N %/% 400 + as.logical(N %% 400)
    X <- character(0)
    for (i in 1:nrequest) {
        a <- (i - 1) * 400 + 1
        b <- 400 * i
        if (i == nrequest) b <- N
        URL <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                     paste(access.nb[a:b], collapse = ","),
                     "&rettype=gb", sep = "")
        X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
    }
    FI <- grep("ORIGIN", X) + 1
    LA <- which(X == "//") - 1
    obj <- list()
    length(obj) <- N
    for (i in 1:N) {
        ## remove all spaces and digits
        tmp <- gsub("[[:digit:] ]", "", X[FI[i]:LA[i]])
        obj[[i]] <- unlist(strsplit(tmp, NULL))
    }
    names(obj) <- seq.names
    if (species.names) {
        tmp <- character(N)
        sp <- grep("ORGANISM", X)
        for (i in 1:N)
          tmp[i] <- unlist(strsplit(X[sp[i]], " +ORGANISM +"))[2]
        attr(obj, "species") <- gsub(" ", "_", tmp)
    }
    obj
}
