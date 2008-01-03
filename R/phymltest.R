## phymltest.R (2005-11-10)

##   Fits a Bunch of Models with PHYML

## Copyright 2004-2005 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.phymltest.model <- c("JC69", "JC69+I", "JC69+G", "JC69+I+G",
                      "K80", "K80+I", "K80+G", "K80+I+G",
                      "F81", "F81+I", "F81+G", "F81+I+G",
                      "F84", "F84+I", "F84+G", "F84+I+G",
                      "HKY85", "HKY85+I", "HKY85+G", "HKY85+I+G",
                      "TN93", "TN93+I", "TN93+G", "TN93+I+G",
                      "GTR", "GTR+I", "GTR+G", "GTR+I+G")

.phymltest.nfp <- c(1, 2, 2, 3, 2, 3, 3, 4, 4, 5, 5, 6, 5, 6, 6, 7,
                    5, 6, 6, 7, 6, 7, 7, 8, 9, 10, 10, 11)

phymltest <- function(seqfile, format = "interleaved", itree = NULL,
                      exclude = NULL, execname, path2exec = NULL)
{
    windoz <- .Platform$OS.type == "windows"
    if (missing(execname)) {
        if (windoz) execname <- "phyml_w32"
        else stop("you must give an executable file name for PHYML")
    }
    outfile <- paste(seqfile, "_phyml_stat.txt", sep = "")
    inp <- seqfile
    if (file.exists(outfile)) inp <- c(inp, "A")
    if (file.exists(paste(seqfile, "_phyml_tree.txt", sep = "")))
      inp <- c(inp, "A")
    if (format != "interleaved") inp <- c(inp, "I")
    if (!is.null(itree)) inp <- c(inp, "U", itree)
    N <- length(.phymltest.model)
    input.model <- list(c(rep("M", 5), "Y"),
                        c(rep("M", 5), "V", rep("Y", 2)),
                        c(rep("M", 5), "R", "A", rep("Y", 2)),
                        c(rep("M", 5), "R", "A", "Y", "V", rep("Y", 2)),
                        c(rep("M", 6), "T", rep("Y", 2)),
                        c(rep("M", 6), "T", "Y", "V", rep("Y", 2)),
                        c(rep("M", 6), "T", "Y", "R", "A", rep("Y", 2)),
                        c(rep("M", 6), "T", "Y", "R", "A", "Y", "V", rep("Y", 2)),
                        c(rep("M", 7), "Y"),
                        c(rep("M", 7), "V", rep("Y", 2)),
                        c(rep("M", 7), "R", "A", rep("Y", 2)),
                        c(rep("M", 7), "V", "Y", "R", "A", rep("Y", 2)),
                        c("M", "T", rep("Y", 2)),
                        c("M", "T", "Y", "V", rep("Y", 2)),
                        c("M", "T", "Y", "R", "A", rep("Y", 2)),
                        c("M", "T", "Y", "V", "Y", "R", "A", rep("Y", 2)),
                        c("T", rep("Y", 2)),
                        c("T", "Y", "V", rep("Y", 2)),
                        c("T", "Y", "R", "A", rep("Y", 2)),
                        c("T", "Y", "V", "Y", "R", "A", rep("Y", 2)),
                        c(rep("M", 2), "T", rep("Y", 2)),
                        c(rep("M", 2), "T", "Y", "V", rep("Y", 2)),
                        c(rep("M", 2), "T", "Y", "R", "A", rep("Y", 2)),
                        c(rep("M", 2), "T", "Y", "R", "A", "Y", "V", rep("Y", 2)),
                        c(rep("M", 3), "Y"),
                        c(rep("M", 3), "V", rep("Y", 2)),
                        c(rep("M", 3), "R", "A", rep("Y", 2)),
                        c(rep("M", 3), "V", "Y", "R", "A", rep("Y", 2)))
    loglik <- numeric(N)
    names(input.model) <- names(loglik) <- .phymltest.model
    if (is.null(path2exec)) exec <- execname
    else exec <- paste(path2exec, execname, sep = "/")
    imod <- if (is.null(exclude)) 1:N else (1:N)[!.phymltest.model %in% exclude]
    for (i in imod) {
        if (i == 2) {
            if (length(inp) == 1) inp <- c(inp, rep("A", 2))
            else if (inp[2] != "A") inp <- c(inp[1], rep("A", 2), inp[2:length(inp)])
        }
        if (windoz) system(exec, input = c(inp, input.model[[i]]))
        else {
            cat(c(inp, input.model[[i]]), file = "f", sep = "\n")
            system(paste(exec, "f", sep = " < "))
        }
        loglik[i] <- scan(paste(seqfile, "_phyml_lk.txt", sep = ""), quiet = TRUE)
    }
    unlink("f")
    loglik <- loglik[imod]
    class(loglik) <- "phymltest"
    loglik
}

print.phymltest <- function(x, ...)
{
    nfp <- .phymltest.nfp[.phymltest.model %in% names(x)]
    X <- cbind(nfp, x, 2 * (nfp - x))
    rownames(X) <- names(x)
    colnames(X) <- c("nb.free.para", "loglik", "AIC")
    print(X)
}

summary.phymltest <- function(object, ...)
{
    nfp <- .phymltest.nfp[.phymltest.model %in% names(object)]
    N <- length(object)
    model1 <- model2 <- character(0)
    chi2 <- df <- P.val <- numeric(0)
    for (i in 1:(N - 1)) {
        for (j in (i + 1):N) {
            if (nfp[i] >= nfp[j]) next
            m1 <- unlist(strsplit(names(object)[i], "\\+"))
            m2 <- unlist(strsplit(names(object)[j], "\\+"))
            if (m1[1] == "K80" && m2[1] == "F81") next
            ## à vérifier que ds les 2 lignes suivantes les conversions
            ## se font bien correctement!!!!
            if (length(grep("\\+I", names(object)[i])) > 0 && length(grep("\\+I", names(object)[j])) == 0) next
            if (length(grep("\\+G", names(object)[i])) > 0 && length(grep("\\+G", names(object)[j])) == 0) next
            ## Now we should be sure that m1 is nested in m2.
            chi2 <- c(chi2, 2 * (object[j] - object[i]))
            df <- c(df, nfp[j] - nfp[i])
            P.val <- c(P.val, 1 - pchisq(2 * (object[j] - object[i]), nfp[j] - nfp[i]))
            model1 <- c(model1, names(object)[i])
            model2 <- c(model2, names(object)[j])
        }
    }
    data.frame(model1, model2, chi2, df, P.val = round(P.val, 4))
}

plot.phymltest <- function(x, main = NULL, col = "blue", ...)
{
    nfp <- .phymltest.nfp[.phymltest.model %in% names(x)]
    N <- length(x)
    aic <- 2 * (nfp - x)
    if (is.null(main))
      main <- paste("Akaike information criterion for",
                    deparse(substitute(x)))
    plot(rep(1, N), aic, bty = "n", xaxt = "n", yaxt = "n",
         type = "n", xlab = "", ylab = "", main = main, ...)
    axis(side = 2, pos = 0.85, las = 2)
    abline(v = 0.85)
    y.lab <- seq(min(aic), max(aic), length = N)
    segments(0.85, sort(aic), 1.1, y.lab, col = col)
    text(1.1, y.lab,
         parse(text = sub("\\+G", "\\+Gamma", names(sort(aic)))),
         adj = 0)
}
