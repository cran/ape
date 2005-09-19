### read.dna.R (2004-04-15)
###
###     Read DNA Sequences in a File
###
### Copyright 2003-2005 Emmanuel Paradis
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
                     nlines = 0, comment.char = "#", seq.names = NULL)
{
    format <- match.arg(format, c("interleaved", "sequential", "fasta"))
    X <- scan(file = file, what = character(), sep = "\n", quiet = TRUE,
              skip = skip, nlines = nlines, comment.char = comment.char)
    fl <- X[1]
    old.warn <- options("warn")
    options(warn = -1)
    ## need to remove the possible leading spaces in the first line
    fl.num <- as.numeric(unlist(strsplit(gsub("^ +", "", fl), " +")))
    options(warn = unlist(old.warn))
    if (format == "interleaved") {
        if (all(is.na(fl.num))) stop("the first line of the file must contain the dimensions of the data\n  if the format is interleaved.")
        if (length(fl.num) != 2) {
            stop("the first line of the file must contain TWO numbers\n  if the format is interleaved.")
        }
        else {
            n <- fl.num[1]
            s <- fl.num[2]
        }
        X <- X[-1]
        fl <- X[1]
        fl <- unlist(strsplit(fl, NULL))
        bases <- grep("[-AaCcGgTtUuMmRrWwSsYyKkVvHhDdBbNn]", fl)
        z <- diff(bases)
        for (i in 1:length(z)) if (all(z[i:(i + 8)] == 1)) break
        start.seq <- bases[i]
        taxa <- substr(X[1:n], 1, start.seq - 1)
        taxa <- sub("^ +", "", taxa) # remove the leading spaces
        taxa <- sub(" +$", "", taxa) # remove the trailing spaces
        taxa <- sub("^['\"]", "", taxa) # remove the leading quotes
        taxa <- sub("['\"]$", "", taxa) # remove the trailing quotes
        X[1:n] <- substr(X[1:n], start.seq, nchar(X[1:n]))
        X <- gsub(" ", "", X)
        nl <- length(X)
        obj <- list()
        length(obj) <- n
        for (i in 1:n) {
            sequ <- paste(X[seq(i, nl, n)], collapse = "")
            obj[[i]] <- unlist(strsplit(sequ, NULL))
        }
    }
    if (format == "sequential") {
        if (all(is.na(fl.num))) stop("the first line of the file must contain the dimensions of the data\n  if the format is sequential.")
        if (length(fl.num) != 2) {
            stop("the first line of the file must contain TWO numbers\n  if the format is sequential.")
        }
        else {
            n <- fl.num[1]
            s <- fl.num[2]
        }
        X <- X[-1]
        fl <- X[1]
        obj <- list()
        length(obj) <- n
        taxa <- character(n)
        j <- 1
        for (i in 1:n) {
            bases <- grep("[AaCcGgTtUuMmRrWwSsYyKkVvHhDdBbNn]", unlist(strsplit(X[j], NULL)))
            z <- diff(bases)
            for (k in 1:length(z)) if (all(z[k:(k + 8)] == 1)) break
            start.seq <- bases[k]
            taxa[i] <- substr(X[j], 1, start.seq - 1)
            sequ <- substr(X[j], start.seq, nchar(X[j]))
            sequ <- gsub(" ", "", sequ)
            j <- j + 1
            while (nchar(sequ) < s) {
                sequ <- paste(sequ, gsub(" " , "", X[j]), sep = "")
                j <- j + 1
            }
            obj[[i]] <- unlist(strsplit(sequ, NULL))
        }
        taxa <- sub("^ +", "", taxa) # remove the leading spaces
        taxa <- sub(" +$", "", taxa) # remove the trailing spaces
        taxa <- sub("^['\"]", "", taxa) # remove the leading quotes
        taxa <- sub("['\"]$", "", taxa) # remove the trailing quotes
    }
    if (format == "fasta") {
        start <- grep("^ {0,}>", X)
        taxa <- X[start]
        n <- length(taxa)
        obj <- list()
        length(obj) <- n
        taxa <- sub("^ {0,}> {0,}", "", taxa) # remove the hook and the spaces before and after
        taxa <- sub(" +$", "", taxa) # remove the trailing spaces
        taxa <- sub("^['\"]", "", taxa) # remove the leading quotes
        taxa <- sub("['\"]$", "", taxa) # remove the trailing quotes
        start <- c(start, length(X) + 1) # this avoids the following to crash when `i = n'
        for (i in 1:n) obj[[i]] <- unlist(strsplit(gsub(" ", "",
                                                        paste(X[(start[i] + 1):(start[i + 1] - 1)],
                                                              collapse = "")),
                                                   NULL))
    }
    names(obj) <- if (is.null(seq.names)) taxa else seq.names
    obj <- lapply(obj, tolower)
    obj
}
