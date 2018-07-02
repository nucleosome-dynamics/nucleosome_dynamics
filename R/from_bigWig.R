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

spec <- matrix(c("input",       "i", 1, "character",
                 "output",      "o", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

## Do it ######################################################################

cover <- readBigWig(args$input)

## Save as an RData ###########################################################

print("string output as RData")
save(cover, file=args[["output"]])

###############################################################################
