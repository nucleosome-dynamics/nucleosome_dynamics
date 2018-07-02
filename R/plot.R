#!/usr/bin/Rscript

## Imports ####################################################################

library(IRanges)
library(getopt)

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
sourced <- c("plot_subset", "make_plot")
for (x in paste0(SOURCE.DIR, "/", sourced, ".R")) {
    source(x)
}

## Parameters and Arguments ###################################################

spec <- matrix(c("input",  "i", 1, "character",
                 "output", "o", 1, "character",
                 "start",  "s", 1, "integer",
                 "end",    "e", 1, "integer",
                 "chr",    "c", 1, "character"),
               byrow=TRUE,
               ncol=4)

args <- getopt(spec)

## Make the plot ##############################################################

dyn <- get(load(args$input))
subdyn <- subsetDyn(dyn, args$chr, args$start, args$end)

png(filename=args[["output"]], width=1000, height=500)
makePlot(subdyn, args$start, args$end)
dev.off()

###############################################################################
