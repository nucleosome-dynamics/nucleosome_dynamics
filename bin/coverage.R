#!/usr/bin/Rscript

## Imports ####################################################################

suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(htSeqTools))
suppressPackageStartupMessages(library(nucleR))
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
source(paste(SOURCE.DIR,
             "helperfuns.R",
             sep="/"))

## Parameters and Arguments ###################################################

defaults <- list(type           = "paired",
                 mc.cores       = 1,
                 fdrOverAmp     = 0.05,
                 components     = 1,
                 fragmentLen    = NULL,
                 trim           = 50)

spec <- matrix(c("input",       "a", 1, "character",
                 "output",      "b", 1, "character",
                 "cores",       "c", 1, "integer",
                 "type",        "e", 1, "character",
                 "fdrOverAmp",  "d", 1, "double",
                 "components",  "f", 1, "integer",
                 "fragmentLen", "g", 1, "integer",
                 "trim",        "h", 1, "integer"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

names(args) <- sub("cores", "mc.cores", names(args))

params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

## Pipeline Itself ############################################################

reads <- get(load(params$input))

message("filtering duplicated reads")
f.reads <- filterDuplReads(reads,
                           fdrOverAmp=params$fdrOverAmp,
                           components=params$components,
                           mc.cores=params$mc.cores)

if (is.null(params$fragmentLen)) {
    if (params$type == "single") {
        message("estimating fragment length")
        params$fragmentLen <- fragmentLenDetect(f.reads)
    } else if (params$type == "paired") {
        params$fragmentLen <- 170
    }
}

message("processing reads")
prep <- processReads(f.reads,
                     type=params$type,
                     fragmentLen=params$fragmentLen,
                     trim=params$trim)

message("calculating coverage")
cover <- coverage.rpm(prep)

## Store the Result ###########################################################

message("-- saving ", args[["output"]])
save(cover, file=args[["output"]])

###############################################################################
