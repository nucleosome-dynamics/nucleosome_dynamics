#!/usr/bin/env Rscript

###############################################################################

library(IRanges)
library(plyr)

###############################################################################

sd2stiffness <- function (sd, t=310.15,  rt=8.3144598*0.239005736)
   # cal/mol * bp^2
    (rt * t) / (sd**2)

###############################################################################

dyadPos <- function (x)
    (start(x) + end(x)) / 2

fitNuc <- function (nuc, reads)
{
    nuc.reads <- selectReads(reads, nuc)
    dyads <- dyadPos(nuc.reads)
    c(k  = max(coverage(nuc.reads)),
      sd = sd(dyads),
      m  = mean(dyads))
}

fitChrNucs <- function (x, reads)
    adply(x, 1, fitNuc, reads, .progress="text")

fitIt <- function (calls, reads)
{
    if (class(reads) == "RangedData") {
        f <- function (x, reads) {
            chr <- x[1, "seqname"]
            message(chr)
            fitChrNucs(x, reads[[chr]])
        }
        ddply(calls, "seqname", f, ranges(reads))
    } else {
        fitChrNucs(calls, reads)
    }
}

###############################################################################
