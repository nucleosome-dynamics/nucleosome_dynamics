#!/usr/bin/Rscript

# Position: region where the movement happens
# Type: change in the nucleosome map
# Score: magnitude of the change

# class: type of hotspot (see help for all possible types)

## Imports ####################################################################

suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(NucDyn))
suppressPackageStartupMessages(library(plyr))

where <- function () {
    spath <- parent.frame(2)$ofile

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
sourced <- c("helperfuns", "gff_funs", "nd_funs", "wig_funs")
for (x in sourced) {
    source(paste0(SOURCE.DIR, "/", x, ".R"))
}

### Parameters and Arguments ###################################################

spec <- matrix(c("input1", "a", 1, "character",
                 "input2", "b", 1, "character",
                 "genome", "c", 1, "character",

                 "outputGff",    "d", 1, "character",
                 "plotRData",    "e", 1, "character",
                 "outputBigWig", "f", 1, "character",

                 "cores",      "g", 1, "integer",
                 "maxLen",     "h", 1, "integer",
                 "equal_size", "i", 1, "logical",
                 "readSize",   "k", 1, "integer",
                 "maxDiff",    "l", 1, "integer",

                 "range", "m", 1, "character",

                 "shift_min_nreads", "p", 1, "double",
                 "shift_threshold",  "q", 1, "double",
                 "indel_min_nreads", "r", 1, "double",
                 "indel_threshold",  "s", 1, "double",

                 "roundPow",       "t", 1, "logical",
                 "same_magnitude", "u", 1, "logical"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

defaults <- list(cores      = 1,
                 maxLen     = 170,
                 equal_size = FALSE,
                 plotRData  = NULL,
                 readSize   = 140,
                 maxDiff    = 74,

                 range      = "All",

                 shift_min_nreads = 3,
                 shift_threshold  = 0.1,
                 indel_min_nreads = 3,
                 indel_threshold  = 0.05)

params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

## Pipeline Itself ############################################################

range <- parseRange(params$range)

message("loading and subsetting reads")
rs <- lapply(params[c("input1", "input2")], function(x) get(load(x)))
rs <- lapply(rs, subsetReads, range$chr, range$start, range$end)

message("running NucleosomeDynamics")
dyn <- nucleosomeDynamics(setA      = rs[[1]],
                          setB      = rs[[2]],
                          maxLen    = params$maxLen,
                          equalSize = params$equal_size,
                          readSize  = params$readSize,
                          maxDiff   = params$maxDiff,
                          mc.cores  = params$cores)

if (!is.null(params$plotRData)) {
    plotable <- makePlotable(dyn)
    save(plotable, file=params$plotRData)
}

message("finding hotspots")
hs <- findHotspots(dyn=dyn, mc.cores=params$cores)

## Calculate vector of -log10(p-value)s  ######################################

cov1 <- lapply(coverage(rs[[1]]), as.vector)
cov2 <- lapply(coverage(rs[[2]]), as.vector)

chrs <- intersect(names(cov1), names(cov2))
pvals <- lapply(chrs, function (x) -log10(findPVals(cov1[[x]], cov2[[x]])))
names(pvals) <- chrs

###############################################################################

indel.threshold <- params$indel_threshold
indel.min.nreads <- params$indel_min_nreads
shift.threshold <- params$shift_threshold
shift.min.nreads <- params$shift_min_nreads

is.indel <- hs$type == "EVICTION" | hs$type == "INCLUSION"
is.shift <- hs$type == "SHIFT +"  | hs$type == "SHIFT -"
indel.sel <- hs$score <= indel.threshold & hs$nreads >= indel.min.nreads
shift.sel <- hs$score <= shift.threshold & hs$nreads >= shift.min.nreads
i <- is.indel & indel.sel
j <- is.shift & shift.sel

hs <- hs[i | j, ]

## Store the Result ###########################################################

names(hs)[names(hs) == "type"] <- "class"
names(hs)[names(hs) == "chr"] <- "seqname"
hs$peak <- NULL

message("saving output as gff")
writeGff(df2gff(hs,
                source="NucleosomeDynamics",
                feature="Nucleosome change"),
         params$outputGff)

message("saving bigWig output")
writeBigWig(lapply(pvals, splitAtZeros),
            params$outputBigWig,
            params$genome)

message("done")

###############################################################################
