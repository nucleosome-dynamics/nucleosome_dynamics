#!/usr/bin/Rscript

# Position: Region between two nucleosomes surrounding the TSS.

# classification: Descriptor of the Transcription Start Site. See the help for possible options.
# distance: Distance in base pairs between the nucleosome +1 and the nucleosome -1.
# id
# nucleosome minus1: Position of the nucleosome -1.
# nucleosome plus1: Position of the nucleosome +1
# TTS_position: Position of the Transcription Start Site.

## Imports ####################################################################

suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(nucleR))

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
sourced <- c("helperfuns", "nucleosome_patterns", "get_genes", "gff_funs","htseqtools_funs")
for (x in paste0(SOURCE.DIR, "/", sourced, ".R")) {
    source(x)
}

## Parameters and Arguments ###################################################

defaults <- list(mc.cores          = 1,
                 window            = 300,
                 #p1.max.merge      = 3,
                 p1.max.downstream = 20,
                 open_thresh       = 215#,
                 #max.uncovered     = 150
                 )

spec <- matrix(c("calls",             "a", 1, "character",
                 "genome",            "c", 1, "character",
                 "output",            "d", 1, "character",
                 "cores",             "e", 1, "integer",
                 "window",            "f", 1, "integer",
                 "p1.max.merge",      "g", 1, "integer",
                 "p1.max.downstream", "h", 1, "integer",
                 "open_thresh",       "i", 1, "integer",
                 "max.uncovered",     "j", 1, "integer"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

names(args) <- sub("cores", "mc.cores", names(args))

params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

## Some function declarations #################################################

message("-- loading inputs")
nucs <- with(readGff(params$calls),
             RangedData(ranges  = IRanges(start = as.numeric(start),
                                          end   = as.numeric(end)),
                        space   = as.character(seqname),
                        score   = as.numeric(score),
                        score_w = as.numeric(score_width),
                        score_h = as.numeric(score_height),
                        #nmerge  = as.numeric(nmerge),
                        class   = as.character(class)))

## Read input #################################################################

message("-- loading used genome")
genes <- getGenes(params$genome)

genes$tss <- as.numeric(genes$tss)
genes$tts <- as.numeric(genes$tts)

## Do it ######################################################################

message("-- checking the classes")
tx.classes <- with(params,
                   patternsByChrDF(calls             = nucs,
                                   df                = genes,
                                   col.id            = "name",
                                   col.chrom         = "chrom",
                                   col.pos           = "tss",
                                   col.strand        = "strand",
                                   window            = window,
                                   p1.max.merge      = p1.max.merge,
                                   p1.max.downstream = p1.max.downstream,
                                   open.thresh       = open_thresh,
                                   max.uncovered     = max.uncovered,
                                   mc.cores          = mc.cores))

tx.classes <- tx.classes[!tx.classes$descr == "NA", ]

## Store output ###############################################################

names(tx.classes)[names(tx.classes) == "descr"] <- "classification"
names(tx.classes)[names(tx.classes) == "dist"] <- "distance"
names(tx.classes)[names(tx.classes) == "m1.pos"] <- "nucleosome_minus1"
names(tx.classes)[names(tx.classes) == "p1.pos"] <- "nucleosome_plus1"
names(tx.classes)[names(tx.classes) == "chrom"] <- "seqname"
names(tx.classes)[names(tx.classes) == "pos"] <- "TSS_position"
#names(tx.classes)[names(tx.classes) == "id"] <- "gene_id"

message("-- saving gff output")
gff <- df2gff(tx.classes,
              source="nucleR",
              feature="TSS classification")
writeGff(gff, params$output)

###############################################################################
