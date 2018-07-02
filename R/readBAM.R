#!/usr/bin/Rscript

suppressPackageStartupMessages(library(getopt))

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
             "loadbams.R",
             sep="/"))

spec <- matrix(c("type",   "t", 1, "character",
                 "input",  "i", 1, "character",
                 "output", "o", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

cat("-- loading ", args[["input"]],"\n")
reads <- loadBAM(args[["input"]], args[["type"]])

cat("-- saving ", args[["output"]],"\n")
save(reads, file=args[["output"]])
