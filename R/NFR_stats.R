#!/usr/bin/Rscript

## Imports ####################################################################

suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(GenomicRanges))

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
sourced <- c("get_genes", "gff_funs")
for (x in sourced) {
    source(paste0(SOURCE.DIR, "/", x, ".R"))
}


## Parameters and Arguments ###################################################

spec <- matrix(c("input",  "a", 1, "character",
                 "genome", "b", 1, "character",
                 "out_gw", "c", 1, "character"),
               byrow=TRUE,
               ncol=4)

params <- getopt(spec)

## Read NFR data  ########################################################

nfr <- readGff(params$input)
nfr_width = nfr$end - nfr$start

## Statistics genome-wide  ####################################################

print("-- computing statistics genome-wide")

stat_nfr = data.frame(NFR=c("Total", "Mean width", "Std. Dev. width"), 
                      Value=c(nrow(nfr),
                              round(mean(nfr_width), 2),
                              round(sd(nfr_width), 2)))

write.csv(stat_nfr, params$out_gw, row.names=F, quote=F)
