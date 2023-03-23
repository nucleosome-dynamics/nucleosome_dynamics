#!/usr/bin/env Rscript

# nucleosome_first: First nucleosome of the gene.
# nucleosme_last: Last nucleosome of the gene.
# score_phase: Is a measure of the phase between the first and the last nucleosome. A score of 0 means the nucleosome are completely phased and a score of 82 corresponds to totally antiphased nucleosomes.
# score_autocorrelation: It is directly computed from the experimental coverage and is quantitative measure of the periodicity of nucleosomes inside the gene body.

## Imports ####################################################################

suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(nucleR))
suppressPackageStartupMessages(library(htSeqTools))

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
sourced <- c("helperfuns", "wig_funs", "get_genes", "periodicity_funs",
             "gff_funs")
for (x in sourced) {
    source(paste0(SOURCE.DIR, "/", x, ".R"))
}

## Parameters and Arguments ###################################################

defaults <- list(periodicity = 165,
                 mc.cores    = 1)#,
                 #genome      = "R64-1-1") #unused

spec <- matrix(c("calls",       "a", 1, "character",
                 "genes",       "b", 1, "character",
                 "chrom_sizes", "c", 1, "character",
                 "bwOutput",    "d", 1, "character",
                 "gffOutput",   "e", 1, "character",
                 "cores",       "f", 1, "integer",
                 "reads",       "g", 1, "character",
                 "type",        "h", 1, "character",
                 "periodicity", "i", 1, "double"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

names(args) <- subMany("cores", "mc.cores", names(args))

params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

## Some function definitions ##################################################

message("loading genes")
genes <- getGenes(params$genes)
message("reading calls")
calls.df <- readGff(params$calls)
calls.rd <- with(calls.df,
                 RangedData(space=seqname,
                            range=IRanges(start=start,
                                          end=end)))

## Do it ######################################################################

message("calculating coverage")

reads <- get(load(params$reads))
reads <- keepSeqlevels(reads, seqlevelsInUse(reads))
f.reads <- filterDuplReads(reads, fdrOverAmp=0.05, components=1)
if (params$type == "single") {
    fragmentLen <- fragmentLenDetect(f.reads)
} else if (params$type == "paired") {
    fragmentLen <- 170
}

prep <- processReads(f.reads,
                     type=params$type,
                     fragmentLen=fragmentLen,
                     trim=50)
cov <- coverage.rpm(prep)

message("identifying first and last nucleosomes")
genes.nucs <- findGenesNucs(genes, calls.rd, params$mc.cores)
genes.nucs$dfi <- getDfi(genes.nucs$nuc.len, params$periodicity)
genes.nucs$autocor <- autocorFromDf(genes.nucs, cov, params$periodicity)

covPredAll <- getPeriodCov(genes.nucs, params$periodicity, params$mc.cores)

## Store output ###############################################################

names(genes.nucs)[names(genes.nucs) == "chrom"] <- "seqname"
names(genes.nucs)[names(genes.nucs) == "dfi"] <- "score_phase"
names(genes.nucs)[names(genes.nucs) == "autocor"] <- "score_autocorrelation"
names(genes.nucs)[names(genes.nucs) == "first"] <- "nucleosome_first"
names(genes.nucs)[names(genes.nucs) == "last"] <- "nucleosome_last"
genes.nucs$nuc.length <- NULL

message("writing GFF")
gff <- df2gff(genes.nucs,
              source="nucleR",
              feature="nucleosome periodicity")
writeGff(gff, params$gffOutput)

message("writing bigWig output")
# splited <- lapply(covPredAll, splitAtZeros)
# writeBigWig(splited, params$bwOutput, params$chrom_sizes)
writeBigWig_updated(covPredAll,
                    params$bwOutput,
                    params$chrom_sizes)

##############################################################################
