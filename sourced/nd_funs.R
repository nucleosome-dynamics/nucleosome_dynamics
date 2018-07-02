#!/usr/bin/Rscript

###############################################################################

library(IRanges)
library(nucleR)
library(htSeqTools)

source(paste(SOURCE.DIR,
             "fp.R",
             sep="/"))

###############################################################################

hsSorter <- function (df) {
    df <- df[order(df$coord), ]
    df <- df[order(df$chr), ]
    unrowname(df)
}

###############################################################################

newApplyThreshold <- function (x, indel.threshs, shift.threshs) {
    threshs <- switch(x[1, "type"],
                      "INCLUSION" = indel.threshs,
                      "EVICTION"  = indel.threshs,
                      "SHIFT +"   = shift.threshs,
                      "SHIFT -"   = shift.threshs)
    x[x$nreads >= threshs[1] & x$score >= threshs[2], ]
}

###############################################################################

nucleRCall <- function (reads, mc.cores=1) {
    f.reads <- filterDuplReads(reads,
                               fdrOverAmp=0.05,
                               components=1,
                               mc.cores=mc.cores)
    prep <- processReads(f.reads,
                         type="paired",
                         fragmentLen=170)
    cover <- coverage.rpm(prep)
    fft <- filterFFT(cover, pcKeepComp=0.02, mc.cores=mc.cores)
    peaks <- peakDetection(fft,
                           width=125,
                           threshold="35%",
                           score=FALSE,
                           mc.cores=mc.cores)
    scores <- peakScoring(peaks,
                          fft,
                          threshold="35%",
                          dyad.length=50,
                          mc.cores=mc.cores)
    merged <- mergeCalls(scores, min.overlap=50)
    merged
}

chrAssignNucs <- function (df, calls, threshold=0.3) {
    q <- IRanges(start=df$start, end=df$end)
    ovs <- findOverlaps(query=q, subject=calls)
    ds <- width(ranges(ovs, q, calls))

    ovlps <- as.data.frame(ovs)
    ovlps$rel.ovlp <- ds / width(q[ovlps$queryHits])
    ovlps <- ddply(ovlps, "queryHits", function (x) x[which.max(x$rel.ovlp), ])
    ovlps <- ovlps[ovlps$rel.ovlp > threshold, ]

    missings <- which(!seq_along(q) %in% ovlps$queryHits)

    nucs <- c(ovlps$subjectHits, rep(0, length(missings)))
    ordering <- order(c(ovlps$queryHits, missings))

    nucs[ordering]
}

assignNucs <- function (hs, calls, threshold=0.3) {
    message("assigning nucleosomes")
    chr.calls <- ranges(calls)
    f <- function (x) {
        chr <- x[1, "chr"]
        calls <- chr.calls[[chr]]
        x$nuc <- chrAssignNucs(x, calls, 0.4)
        x
    }
    ddply(hs, "chr", f, .progress="text")
}

###############################################################################

combineShiftTypes <- function (x, y) {
    if (x == "SHIFT +" && y == "SHIFT +") {
        return("SHIFT +")
    } else if (x == "SHIFT -" && y == "SHIFT -") {
        return("SHIFT -")
    } else if (x == "SHIFT +" && y == "SHIFT -") {
        return("DECREASED FUZZINESS")
    } else if (x == "SHIFT -" && y == "SHIFT +") {
        return("INCREASED FUZZINESS")
    }
}

combineThem <- function (x) {
    nreads <- sum(x$nreads)
    start <- min(x$start)
    end <- max(x$end)
    type <- do.call(combineShiftTypes, as.list(x$type[order(x$start)]))
    changedArea <- sum(x$changedArea)
    involvedArea <- sum(x$involvedArea)
    score <- changedArea / involvedArea

    chr <- x[1, "chr"]
    nuc <- x[1, "nuc"]

    data.frame(nreads       = nreads,
               start        = start,
               end          = end,
               type         = type,
               changedArea  = changedArea,
               involvedArea = involvedArea,
               score        = score,
               chr          = chr,
               nuc          = nuc)
}

shiftCombiner <- function (df) {
    f <- function (x) {
        if (nrow(x) != 2) {
            return(x)
        } else {
            combineThem(x)
        }
    }
    ddply(df, "chr", ddply, "nuc", f)
}

###############################################################################
