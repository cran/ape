### read.dna.R  (2003-01-31)
###
###     Read DNA Sequences in a File
###
### Copyright 2003 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

read.dna <- function(file, format = "interleaved", skip = 0,
                     nlines = 0, comment.char = "#")
{
    X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE,
              skip = skip, nlines = nlines, comment.char = comment.char,
              blank.lines.skip = FALSE)
    fl <- X[1]
    fl.num <- as.numeric(unlist(strsplit(fl, " +")))
    if (all(!is.na(fl.num))) {
        n <- fl.num[1]
        X <- X[-1]
        fl <- X[1]
    }
    else {
        if (format == "interleaved") n <- which(X == "")[1] - 1
        else n <- length(X)
    }
    fl <- unlist(strsplit(fl, NULL))
    bases <- grep("[AaCcGgTt]", fl)
    z <- diff(bases)
    for (i in 1:length(z)) if (all(z[i:(i + 9)] == 1)) break
    start.seq <- bases[i]
    end.seq <- length(fl)
    taxa <- substr(X[1:n], 1, start.seq - 1)
    taxa <- sub(" +$", "", taxa)
    taxa <- sub("^['\"]", "", taxa)
    taxa <- sub("['\"]$", "", taxa)
    taxa <- gsub("_", " ", taxa)
    X <- X[X != ""]
    X <- substr(X, start.seq, end.seq)
    X <- gsub(" ", "", X)
    n.line <- length(X)
    obj <- list()
    length(obj) <- n
    if (format == "interleaved") {
        for (i in 1:n) {
            sequ <- paste(X[seq(i, n.line, n)], collapse = "")
            obj[[i]] <- unlist(strsplit(sequ, NULL))
        }
    }
    else for (i in 1:n) obj[[i]] <- unlist(strsplit(X[i], NULL))
    names(obj) <- taxa
    obj
}
