#!/usr/bin/Rscript

## Imports ####################################################################

suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(IRanges))

where <- function () {
    spath <-parent.frame(2)$ofile

    if (is.null(spath)) {
        args <- commandArgs()
        filearg <- args[grep("^--file=", args)]
        fname <- strsplit(filearg, "=")[[1]][2]
    } else {
        fname <- spath
    }

    dirname(normalizePath(fname))
}

SOURCE.DIR <- paste(where(), "../sourced", sep="/")
source(paste(SOURCE.DIR,
             "wig_funs.R",
             sep="/"))

## Parameters and Arguments ###################################################

spec <- matrix(c("input",  "i", 1, "character",
                 "output", "o", 1, "character",
                 "genome", "g", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

## Load RData #################################################################

cover <- get(load(args$input))

## Do it ######################################################################

writeBigWig(lapply(cover,
                   splitAtZeros),
            args$output,
            args$genome)

###############################################################################
