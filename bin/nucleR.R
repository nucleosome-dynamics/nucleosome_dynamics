#!/usr/bin/Rscript

# Score: Positionning score. It is calculated as the weighted sum of width and height scores.

# score_width: Witdth score. It is a measure of how sharp a peak is. A value of 0 would be an extremely wide peak and a value of 1 a very sharp one.
# score_height: Height score. Tells how large a peak of a nucleosome is. The bigger this number, the higher the peak.
# class: Whether the nucleosome is well-positioned (W) or fuzzy (F) or undetermined. The taken value depends on score_h and score_w. Undetermined means the exact position of the nucleosome cannot be determined due to strong fuzziness.

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
sourced <- c("helperfuns", "gff_funs")
for (x in sourced) {
    source(paste0(SOURCE.DIR, "/", x, ".R"))
}

## Parameters and Arguments ###################################################

defaults <- list(mc.cores            = 1,
                 type                = "paired",
                 fdrOverAmp          = 0.05,
                 components          = 1,
                 fragmentLen         = 170, #NULL,
                 trim                = 50,
                 pcKeepComp          = 0.02,
                 width               = 147,
                 threshold           = TRUE,
                 thresholdPercentage = 35,
                 thresholdValue      = 10,
                 dyad_length         = 50, #NULL
                 min.overlap         = 80, #NULL
                 score_w.thresh      = 0.6,
                 score_h.thresh      = 0.4,
                 start               = NULL,
                 end                 = NULL,
                 chr                 = NULL)

spec <- matrix(c("input",       "a", 1, "character",
                 "output",      "b", 1, "character",
                 "cores",       "c", 1, "integer",
                 "type",        "e", 1, "character",
                 "fdrOverAmp",  "d", 1, "double",
                 "components",  "f", 1, "integer",
                 "fragmentLen", "g", 1, "integer",
                 "trim",        "h", 1, "integer",
                 "pcKeepComp",  "i", 1, "double",
                 "width",       "j", 1, "integer",

                 "threshold",           "k", 1, "logical",
                 "thresholdPercentage", "s", 1, "double",
                 "thresholdValue",      "t", 1, "integer",

                 "dyad_length", "l", 1, "integer",
                 "minoverlap",  "m", 1, "integer",
                 "wthresh",     "n", 1, "double",
                 "hthresh",     "o", 1, "double",
                 "start",       "p", 1, "integer",
                 "end",         "q", 1, "integer",
                 "chr",         "r", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

names(args) <- subMany(c("cores",
                         "dyadlength",
                         "minoverlap",
                         "wthresh",
                         "hthresh"),
                       c("mc.cores",
                         "dyad_length",
                         "min.overlap",
                         "score_w.thresh",
                         "score_h.thresh"),
                       names(args))

params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

if (is.null(params$dyad_length)) {
    params$dyad_length <- params$trim
}
if (is.null(params$min.overlap)) {
    params$min.overlap <- params$trim
}

if (params$threshold) {
    threshold <- paste0(params$thresholdPercentage, "%")
} else {
    threshold <- params$thresholdValue
}

## Pipeline Itself ############################################################

message("loading data")
reads <- get(load(params$input))
#reads <- RangedData(reads$ranges, space  = droplevels(reads$space))
reads <- keepSeqlevels(reads, seqlevelsInUse(reads))

message("filtering duplicated reads")
f.reads <- filterDuplReads(reads,
                           fdrOverAmp = params$fdrOverAmp,
                           components = params$components,
                           mc.cores   = params$mc.cores)

if (is.null(params$fragmentLen)) {
    if (params$type == "single") {
        message("estimating fragment length")
        params$fragmentLen <- fragmentLenDetect(f.reads)
    } else if (params$type == "paired") {
        params$fragmentLen <- 170
    }
}

if (!is.null(params$chr)) {
    f.reads <- f.reads[seqnames(f.reads) == params$chr, ]
    if (!is.null(params$start) && !is.null(params$end)) {
        sel <- isIn(f.reads,
                    IRanges(start=params$start,
                            end=params$end))
        f.reads <- f.reads[sel, ]
    }
}

message("processing reads")
prep <- processReads(f.reads,
                     type        = params$type,
                     fragmentLen = params$fragmentLen,
                     trim        = params$trim)

message("calculating coverage")
cover <- coverage.rpm(as(prep, "GRanges"))


emptyHandler <- function (f)
    function (x, ...)
        if (length(x) > 0) {
            f(x, ...)
        } else {
            numeric()
        }

message("filtering with filterFFT")
fft <- mclapply(cover,
                emptyHandler(filterFFT),
                pcKeepComp     = params$pcKeepComp,
                mc.preschedule = FALSE,
                mc.cores       = params$mc.cores)

message("detecting peaks")
peaks <- peakDetection(fft,
                       width     = params$width,
                       threshold = threshold,
                       score     = FALSE,
                       min.cov   = 0, 
                       mc.cores  = params$mc.cores)

message("scoring peaks")
scores <- peakScoring(peaks,
                      fft,
                      threshold   = threshold,
                      dyad.length = params$dyad_length,
                      mc.cores    = params$mc.cores)

message("merging peaks")
merged <- mergeCalls(scores,
                     min.overlap = params$min.overlap,
                     mc.cores    = params$mc.cores)

merged$class <- getType(merged$score_w,
                        merged$score_h,
                        params$score_w.thresh,
                        params$score_h.thresh,
                        merged$nmerge)

## Store the Result ###########################################################

merged <- rd2df(merged)
names(merged)[names(merged) == "score_h"] <- "score_height"
names(merged)[names(merged) == "score_w"] <- "score_width"
merged$nmerge <- NULL
colnames(merged)[1] = "seqname"

message("saving output as gff")
writeGff(df2gff(merged,
                source  = "nucleR",
                feature = "Nucleosome"),
         params$output)

message("done")

###############################################################################
