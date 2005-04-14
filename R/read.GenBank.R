### read.GenBank.R  (2005-03-29)
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
    obj <- lapply(as.list(access.nb), get.dna)
    if (species.names) species <- names(obj)
    names(obj) <- seq.names
    if (species.names) {
        species <- lapply(as.list(access.nb), get.species.names)
        attr(obj, "species") <- gsub(" ", "_", unlist(species))
    }
    obj
}

get.dna <- function(access.nb)
{
    URL <- paste("http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?val=",
                 access.nb, sep = "")
    tmp <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
    FI <- grep("ORIGIN", tmp) + 1
    LA <- which(tmp == "//") - 1
    tmp <- gsub("[[:digit:] ]", "", tmp[FI:LA]) # remove all spaces and digits
    tmp <- unlist(strsplit(tmp, NULL))
    tmp
}

get.species.names <- function(access.nb)
{
    URL <- paste("http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?val=",
                 access.nb, sep = "")
    tmp <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
    unlist(strsplit(tmp[grep("ORGANISM", tmp)], "[<>]"))[3]
}
