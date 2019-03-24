#!/usr/bin/Rscript

# Functions to load reads from a BAM file.
# Can process either single-end or paired-end experiments

suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Rsamtools))

source(paste(SOURCE.DIR, "helperfuns.R", sep="/"))

sortBy <- function (xs, a)
    lapply(xs, `[`, sort.list(xs[[a]]))

loadSingleBam <- function (exp)
{
    what <- c("pos", "qwidth", "strand", "rname")
    bam <- scanBam(exp, param=ScanBamParam(what=what))[[1]]

    non.na <- Reduce(`&`, lapply(bam, Negate(is.na)))
    filtered.bam <- lapply(bam, `[`, non.na)

    # IRanges
    GRanges(   filtered.bam$rname,
               IRanges(start = filtered.bam[["pos"]],
                       width = filtered.bam[["qwidth"]]),
               strand = filtered.bam[["strand"]])
}

processStrand <- function (strand, bam)
{
    message(sprintf("processing strand %s", strand))

    p1 <- ifelse(strand == "+", 99, 163)
    p2 <- ifelse(strand == "+", 147, 83)

    unsorted.reads1 <- bam[bam$flag == p1, ]
    unsorted.reads2 <- bam[bam$flag == p2, ]

    rownames(unsorted.reads1) <- as.vector(unsorted.reads1$qname)
    rownames(unsorted.reads2) <- as.vector(unsorted.reads2$qname)

    # Sort by the name of the reads. Assiming the paired reads will have the
    # same name, this will keep the pairs in the same position
    common <- intersect(rownames(unsorted.reads1), rownames(unsorted.reads2))

    reads1 <- unsorted.reads1[common, ]
    reads2 <- unsorted.reads2[common, ]

    # Consistency check
    test <- all(vectorizedAll(reads1$mpos  == reads2$pos,
                              reads2$mpos  == reads1$pos,
                              reads1$rname == reads2$rname))

    if (!test) {
        stop(sprintf("ERROR: Mate selection for %s strand is invalid",
                     strand))
    } else {
        GRanges( as.character(reads1$rname),
                 IRanges(start = reads1$pos,
                         end   = reads2$pos + reads2$qwidth - 1))
    }
}

loadPairedBam <- function (file)
{
    # Read BAM file (only one access to disk, intended for Shared Memory)
    message(sprintf("reading file %s", file))

    what <- c("qname",
              "flag",
              "rname",
              "strand",
              "pos",
              "qwidth",
              "mrnm",
              "mpos")
    bam <- as.data.frame(scanBam(file=file, param=ScanBamParam(what=what))[[1]])

    message("processing flags")
    ## We will process the flags in R
    ## (an alternative is multiple scanBam calls...)
    bam$flag <- bam$flag %% 256

    # Process both strand and return the reads in sorted order
    sortReads(c(processStrand("+", bam),
                processStrand("-", bam)))
}

loadBAM <- function (f, type="single")
{
    if (type == "single") {
        loadSingleBam(f)
    } else if (type == "paired") {
        loadPairedBam(f)
    } else {
        stop("type must be `single` or `paired`")
    }
}
