### read.GenBank.R  (2002-08-28)
###
###     Read DNA Sequences from GenBank via Internet
###
### Copyright 2002 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

read.GenBank <- function(access.nb, seq.names = access.nb)
{
    obj <- lapply(as.list(access.nb), get.dna)
    names(obj) <- seq.names
    return(obj)
}

get.dna <- function(access.nb)
{    
    URL <- paste("http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?val=",
                 access.nb, sep = "")
    tmp <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
    FI <- grep("ORIGIN", tmp) + 1
    LA <- which(tmp == "//") - 1
    tmp <- unlist(strsplit(tmp[FI:LA], NULL))
    tmp <- tmp[grep("[acgt]", tmp)]
    return(tmp)
}
