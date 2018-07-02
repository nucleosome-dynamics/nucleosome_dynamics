#!/usr/bin/Rscript

suppressPackageStartupMessages(library(IRanges))

subsetDyn <- function(dyn, plot.chr, plot.start, plot.end)
{   # This function is intended to be used as a helper to select a
    # region before plotting it.
    # Given:
    #     - a dynamic (a list of lists of reads), where the first level
    #       represents experiments, the second chromosomes and the third types
    #       of reads.
    #     - a region of the genome, specified a chromosome, a stard and
    #       an end.
    # Return:
    #    The part of the dynamic that is in that region of the genome.
    #    Original reads and indels can be directly selected.
    #    In the case of shifts, first we select the ones in the first
    #    experiment that belong to the range, and then we select their pairs
    #    in the second experiment, regardless of their location (if they are
    #    pairs, they won't be very far anyway).

    chr.dyn <- selectChr(dyn, plot.chr)
    shifts <- getShifts(chr.dyn)
    nonshifts <- getNonShifts(chr.dyn)

    sub.shifts <- selectShifts(shifts,
                               plot.start,
                               plot.end)
    sub.nonshifts <- selectNonShifts(nonshifts,
                                     plot.start,
                                     plot.end)

    reJoinReads(sub.shifts, sub.nonshifts)
}

myFilter <- function(x, f, ...)
    # Faster and more `lapply`-friendly version of Filter
    x[f(x, ...)]

selectChr <- function(dyn, chr)
    # Given a dynamic (a list of two lists of reads), subset by chromosome
    lapply(dyn, `[[`, chr)

getShifts <- function(dyn)
    # Given a dynamic already subsetted by chromosome, return the
    # reads belonging to shifts
    lapply(dyn,
           `[`,
           c("left.shifts",
             "right.shifts"))

getNonShifts <- function(dyn)
    # Given a dynamic already subsetted by chromosome, return the
    # reads belonging to the original reads and to the indels
    lapply(dyn,
           `[`,
           c("originals",
             "indels"))

getRanSel <- function(x, start, end)
    # Return what elements of an IRanges are contained in a range specified
    # by `start` and `end`
    start(x) >= start & end(x) <= end

selectShifts <- function(xs, start, end)
{   # Select all the pairs of shifts that belong to a range specified
    # by `start` and `end`.
    # We first find what reads in the first set are contained in the range
    # and then select them and their pairs in the second set.
    sel <- lapply(xs[[1]],
                  getRanSel,
                  start,
                  end)
    lapply(xs,
           function(x) mapply(`[`,
                              x,
                              sel,
                              SIMPLIFY=FALSE))
}

selectNonShifts <- function(xs, start, end)
    # Select al reads that are contained by a range specified by `start`
    # and `end`
    lapply(xs,
           lapply,
           myFilter,
           getRanSel,
           start,
           end)

reJoinReads <- function(xs, ys)
    # Join back two subsets of reads
    mapply(c,
           xs,
           ys,
           SIMPLIFY=FALSE)
