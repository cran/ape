### write.nexus.R  (2004-12-20)
###
###           Write Tree File in Nexus Format
###
### Copyright 2004 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

write.nexus <- function(phy, file = "", translate = TRUE, original.data = TRUE)
{
    if (class(phy) != "phylo")
      stop(paste("object \"", deparse(substitute(phy)),
                 "\" is not of class \"phylo\""), sep = "")
    cat("#NEXUS\n", file = file)
    cat(paste("[R-package APE, ", date(), "]\n\n", sep = ""), file = file, append = TRUE)
    if(original.data){
        if (!is.null(attr(phy, "origin"))) {
            if (!file.exists(attr(phy, "origin"))) {
                warning(paste("the file", attr(phy, "origin"),
                              "cannot be found,\nthe original data won't be written with the tree."))
                original.data <- FALSE
            }
            else {
                ORI <- scan(file = attr(phy, "origin"), what = character(),
                            sep = "\n", skip = 1)
                start <- grep("BEGIN TAXA;", ORI)
                ORI <- ORI[-(1:(start - 1))]
                ORI <- gsub("ENDBLOCK;", "END;", ORI)
                endblock <- grep("END;", ORI)
                start <- grep("BEGIN TREES;", ORI)
                end <- endblock[endblock > start][1]
                cat(ORI[1:(start - 1)], file = file, append = TRUE, sep = "\n")
                ORI <- ORI[-(1:end)]
            }
        }
        else original.data <- FALSE
    }
    if (!original.data) {
        cat("BEGIN TAXA;\n", file = file, append = TRUE)
        cat(paste("\tDIMENSIONS NTAX = ", length(phy$tip.label), ";\n", sep = ""),
            file = file, append = TRUE)
        cat("\tTAXLABELS\n", file = file, append = TRUE)
        cat(paste("\t\t", gsub(" ", "_", phy$tip.label), sep = ""),
            sep = "\n", file = file, append = TRUE)
        cat("\t;\n", file = file, append = TRUE)
        cat("END;\n", file = file, append = TRUE)
    }

    cat("BEGIN TREES;\n", file = file, append = TRUE)
    phy$tip.label <- gsub(" ", "_", phy$tip.label)
    if (translate) {
        cat("\tTRANSLATE\n", file = file, append = TRUE)
        X <- paste("\t\t", 1:length(phy$tip.label), "\t", phy$tip.label, ",", sep = "")
        X[length(X)] <- gsub(",", "", X[length(X)])
        cat(X, file = file, append = TRUE, sep = "\n")
        cat("\t;\n", file = file, append = TRUE)
        phy$tip.label <- as.character(1:length(phy$tip.label))
    }
    cat("\tTREE * UNTITLED = [&R] ", file = file, append = TRUE)
    cat(write.tree(phy, file = "", multi.line = F), file = file, append = TRUE)
    cat("\nEND;\n", file = file, append = TRUE)
    if(original.data) cat(ORI, file = file, append = TRUE, sep = "\n")
}
