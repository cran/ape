### write.dna.R  (2003-07-03)
###
###     Write DNA Sequences in a File
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

write.dna <- function(x, file = "", format = "interleaved", append = FALSE)
{
    if (is.null(names(x))) names(x) <- as.character(1:length(x))
    margin <- max(nchar(names(x))) + 1

    if (format == "interleaved") {
        N <- length(x)
        S <- length(x[[1]])
        if (append) cat(N, S, "\n", file = file, append = TRUE)
        else cat(N, S, "\n", file = file)
        for (i in 1:length(x)) {
            n <- length(x[[i]])
            ny <- ceiling(n / 10)
            if (ny == 1) y <- paste(x[[i]], collapse = " ")
            else {
                y <- numeric(ny)
                rest <- n %% 10
                if (rest == 0) {
                    for (j in 1:ny) y[j] <- paste(x[[i]][(j*10 - 9):(j*10)], collapse = "")
                }
                else {
                    for (j in 1:(ny - 1)) y[j] <- paste(x[[i]][(j*10 - 9):(j*10)], collapse = "")
                    y[ny] <- paste(x[[i]][(n - rest + 1):n], collapse = "")
                }
            }
            x[[i]] <- y
        }
        for (i in 1:length(x)) {
            n <- length(x[[i]])
            ny <- ceiling(n / 6)
            if (ny == 1) y <- paste(x[[i]], collapse = " ")
            else {
                y <- numeric(ny)
                rest <- n %% 6
                if (rest == 0) {
                    for (j in 1:ny) y[j] <- paste(x[[i]][(j*6 - 5):(j*6)], collapse = " ")
                }
                else {
                    for (j in 1:(ny - 1)) y[j] <- paste(x[[i]][(j*6 - 5):(j*6)], collapse = " ")
                    y[ny] <- paste(x[[i]][(n - rest + 1):n], collapse = " 	")
                }
            }
            x[[i]] <- y
        }
        k <- max(unlist(lapply(x, length)))
        for (i in 1:length(x)) {
            cat(names(x)[i], rep(" ", margin - nchar(names(x)[i])), x[[i]][1], "\n",
                sep = "", file = file, append = TRUE)
        }
        if (k > 1) {
            for (j in 2:k) {
                for (i in 1:length(x)) {
                    if (length(x[[i]]) < k) cat("\n", file = file, append = TRUE)
                    else cat(rep(" ", margin), x[[i]][k], "\n", sep = "", file = file, append = TRUE)
                }
            }
        }
    }
    if (format == "sequential") {
        N <- length(x)
        S <- length(x[[1]])
        if (append) cat(N, S, "\n", file = file, append = TRUE)
        else cat(N, S, "\n", file = file)
         for (i in 1:length(x)) {
            cat(names(x)[i], rep(" ", margin - nchar(names(x)[i])),
                paste(x[[i]], collapse = ""), "\n",
                sep = "", file = file, append = TRUE)
        }
    }
}
