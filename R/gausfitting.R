#!/usr/bin/env Rscript

# Score: Stiffness estimation. It's the energy required to move the nucleosome (expressed in cal/bp), it's derived from the gaussian standard deviation.

# nucleR_score: the nucleR score given to that nucleosome
# nucleR_class: the nucleR class given to that nucleosome
# gauss_k: the height of the peak of the gaussian curve
# gauss_m: the position of the peak of the gaussian curve
# gauss_sd: the standard deviation of the gaussian curve

## Imports ####################################################################

library(getopt)
library(IRanges)
library(plyr)

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
sourced <- c("gff_funs", "fp", "gauss_funs", "helperfuns")
for (x in paste0(SOURCE.DIR, "/", sourced, ".R")) {
     source(x)
}

## Command line arguments #####################################################

defaults <- list(t     = 310.15,
                 range = "All")

spec <- matrix(c("calls",  "a", 1, "character",
                 "reads",  "b", 1, "character",
                 "output", "c", 1, "character",
                 "range",  "g", 1, "character",
                 "t",      "t", 1, "double"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

## Subset #####################################################################

range <- parseRange(params$range)

message("loading inputs")
calls <- readGff(params[["calls"]])
reads <- get(load(params[["reads"]]))

reads <- subsetReads(reads, range$chr, range$start, range$end)
calls <- subsetCalls(calls, range$chr, range$start, range$end)

message("performing the fittings")
gauss.df <- fitIt(calls, reads)
gauss.df$score <- sd2stiffness(gauss.df$sd, params$t)

## Save output ################################################################

gauss.df$score_width <- NULL
gauss.df$score_height <- NULL
gauss.df$nmerge <- NULL

names(gauss.df)[names(gauss.df) == "class"] <- "nucleR_class"
names(gauss.df)[names(gauss.df) == "k"] <- "gauss_k"
names(gauss.df)[names(gauss.df) == "m"] <- "gauss_m"
names(gauss.df)[names(gauss.df) == "sd"] <- "gauss_sd"

gauss.df$feature <- "Stiffness estimation"

gff <- df2gff(gauss.df)
writeGff(gff, params[["output"]])

###############################################################################
