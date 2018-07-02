#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(plyr))
source(paste(SOURCE.DIR,
             "helperfuns.R",
             sep="/"))

###############################################################################
# Helpers #####################################################################

.mid <- function(x)
    floor((start(x)+end(x)) / 2)

.gimmeDist <- function(x, open.thresh)
    ifelse(is.na(x),        "-",
    ifelse(x > open.thresh, "open",
    ifelse(x < 120,         "overlap",
                            "close")))

detectNuc <- function (xs, pos, margin, strand, nucpos, position="tss")
{
    a <- strand == "+" && nucpos == "p1"
    b <- strand == "-" && nucpos == "m1"
    c <- strand == "-" && nucpos == "p1"
    d <- strand == "+" && nucpos == "m1"

    flipper <- ifelse(position == "tss", `(`, `!`)

    after <- flipper(a || b)
    before <- flipper(c || d)

    if (after) {
        shift <- `-`
        comp <- `>`
        closest <- which.min
    } else if (before) {
        shift <- `+`
        comp <- `<`
        closest <- which.max
    }

    subxs <- xs[comp(.mid(xs), shift(pos, margin)), ]
    x <- subxs[closest(.mid(subxs)), ]
}

getDescr <- function (dist, m1.class, p1.class, open.thresh)
{
    if (!is.null(m1.class) & !is.null(p1.class)) {
        dist.class <- .gimmeDist(dist, open.thresh)
        descr <- paste(m1.class, dist.class, p1.class, sep="-")
    } else if (is.null(p1.class)) {
        descr <- "+1_missing"
    } else if (is.null(m1.class)) {
        descr <- "-1_missing"
    } else {
        descr <- NA
    }
    descr
}

checkPos <- function (nuc.pos, pos, window)
    nuc.pos > (pos - window) && nuc.pos < (pos + window)

###############################################################################
# Top level wrapper ###########################################################

patternsByChrDF <- function (calls, df, col.id="name", col.chrom="chrom",
                             col.pos="pos", col.strand="strand", ...,
                             mc.cores=1)
{   # Works by first splitting by chromosomes to save some redundant subsetting
    iterFun <- function (chr.genes) {
        chrom <- chr.genes[1, col.chrom]
        cat(chrom,"\n")
        chr.nucs <- calls[space(calls) == chrom, ]
        if (nrow(chr.nucs) == 0) {
            emptyChrom(df         = chr.genes,
                       chrom      = chrom,
                       col.id     = col.id,
                       col.pos    = col.pos,
                       col.strand = col.strand)
        } else {
            res <- nucleosomePatternsDF(calls      = chr.nucs,
                                        df         = chr.genes,
                                        col.id     = col.id,
                                        col.pos    = col.pos,
                                        col.strand = col.strand,
                                        ...,
                                        mc.cores   = mc.cores)
            res$chrom <- chrom
            res
        }
    }
    ddply(df, col.chrom, iterFun)
}

###############################################################################
# Lower level wrappers ########################################################

nucleosomePatternsDF <- function (calls, df, col.id="name", col.pos="pos",
                                  col.strand="strand", ..., mc.cores=1)
{   # Works on the calls and genes of a sigle chromosome
    n <- nrow(df)
    iterRows <- function (i, report.every=100) {
        if (i %% report.every == 0) {
            cat(i, "/", n,"\n")
        }
        nucleosomePatterns(calls  = calls,
                           id     = df[i, col.id],
                           pos    = df[i, col.pos],
                           strand = df[i, col.strand],
                           ...)
    }
    do.call(rbind,
            mclapply(1:n,
                     iterRows,
                     report.every=10,
                     mc.cores=mc.cores))
}

emptyChrom <- function (df, chrom, col.id="name", col.pos="pos",
                        col.strand="strand")
    # Return entries filled with NAs for chromosomes with no nucleosomes found
    data.frame(id     = df[, col.id],
               chrom  = chrom,
               strand = df[, col.strand],
               pos    = df[, col.pos],
               p1.pos = NA,
               m1.pos = NA,
               dist   = NA,
               descr  = "NA",
               start  = df[, col.pos],
               end    = df[, col.pos])

###############################################################################
# Function that does the actual work ##########################################

nucleosomePatterns <- function (calls, id, pos, strand="+", window=300,
                                p1.max.merge=3, p1.max.downstream=20,
                                open.thresh=215, max.uncovered=150,
                                position="tss")
{   # The start and end will be the closest nucleosome found nearby the TSS
    # those positions will also be the p1 and m1 positions if they are within
    # a -/+ window
    p1 <- detectNuc(calls, pos, p1.max.downstream, strand, "p1", position)

    no.p1s <- nrow(p1) == 0  # no nucleosome found upstream of the TSS
    if (no.p1s) {
        p1.pos <- pos
    } else {
        p1.pos <- .mid(p1)
    }

    if (!no.p1s && checkPos(p1.pos, pos, window)) {
        # there's a nucleosome upstream and within the window
        p1.nuc <- .mid(p1)
        p1.class <- p1$class

        # look for m1 relative to p1
        m1 <- detectNuc(calls, p1.nuc, 0, strand, "m1", position)

        no.m1s <- nrow(m1) == 0
        if (no.m1s) {
            m1.pos <- pos
        } else {
            m1.pos <- .mid(m1)
        }

        if (!no.m1s && checkPos(m1.pos, pos, window)) {
            # and there's also one downstream
            m1.class <- m1$class
            m1.nuc <- .mid(m1)
        } else {
            m1.nuc <- NA
            m1.class <- NULL
        }
    } else {
        p1.nuc <- NA
        m1.nuc <- NA
        p1.class <- NULL
        m1.class <- NULL

        # no p1... so look for m1 relative to the TSS
        m1.pos <- .mid(detectNuc(calls, pos, 0, strand, "m1", position))
        if (length(m1.pos) == 0) {
            m1.pos <- pos
        }
    }

    dist <- abs(p1.nuc - m1.nuc)
    descr <- getDescr(dist, m1.class, p1.class, open.thresh)

    data.frame(id=id,
               strand=strand,
               pos=pos,
               p1.pos=p1.nuc,
               m1.pos=m1.nuc,
               dist=dist,
               descr=descr,
               start=min(m1.pos, p1.pos),
               end=max(m1.pos, p1.pos))
}
