#!/usr/bin/Rscript

# Miscelanious helper functions

library(IRanges)
library(GenomicRanges)

source(paste(SOURCE.DIR, "fp.R", sep="/"))

vectorizedAll <- function(...)
    # Helper function that behaves as a vectorized version of the function
    # `all`
    Reduce(`&`, list(...))

getType <- function(score_w, score_h, thresh_w, thresh_h, nmerge=1)
    # Given a vector of score w, a vector of score h, a threshold for w and a
    # threshold for h, return a vector that tells whether those scores are
    # classified as `well-positioned`, or `fuzzy`.
    ifelse(nmerge > 2,
           "uncertain",
           ifelse(`&`(score_w > thresh_w,
                      score_h > thresh_h),
                  "W",
                  "F"))

checkInF <- function(f)
{   # Some checks on the input file specified
    if (is.null(f)) {
        stop("An input file must be specified")
    }
    if (!grepl("\\.bam$", f)) {
        stop("Input file must be in BAM format")
    }
    if (!file.exists(f)) {
        stop("Specified input file doesn't exist")
    }
}

makePlotable <- function(dyn)
{
    message("building structure to be saved for future plotting")

    useful.types <- c("originals", "right.shifts", "left.shifts", "indels")
    chrs <- unique(unlist(lapply(seqnames(set.a(dyn)), levels)))

    lapply(
        list(set.a(dyn), set.b(dyn)),
        function(set) {
            by.chrs <- lapply(
                chrs,
                function(chr)
                    as.list(ranges(set[seqnames(set) == chr, ])[useful.types])
            )
            names(by.chrs) <- chrs
            by.chrs
        }
    )
}

sortDfBy <- function(df, xs)
    # Sort a data.frame by a given value
    do.call(compose,
            lapply(xs,
                   flip2args(partial),
                   orderBy))(df)

orderBy <- function(df, x)
    df[order(df[, x]), ]

getFirstTx <- function(x, df)
{
    entries <- subset(df, GENEID == x)
    f <- `[[`(list("+"=which.min,
                   "-"=which.max),
              unique(entries$TXSTRAND))
    entries[f(entries$TXSTART), ]
}

subMany <- function (patterns, replacements, x)
    # Replace a vector of patterns by a vector of replacements with a 1 to 1
    # equivalence
    do.call(compose,
            mapply(function(p, r) {force(p)
                                   force(r)
                                   function(x) sub(p, r, x)},
                   patterns,
                   replacements
    ))(x)

###############################################################################

.check.mc <- function (mc.cores)
 {
    lib <- "parallel"
    if (mc.cores > 1 && !lib %in% loadedNamespaces()) {
        warning("'",
                lib,
                "' library not available, switching to m.cores=1")
        return(1)
    } else {
        return(mc.cores)
    }
}

xlapply <- function(X, FUN, ..., mc.cores=1)
{   # Choose between lapply pr mclapply accordingly
    actual.cores <- .check.mc(mc.cores)
    if (actual.cores > 1) {
        mclapply(X=X,
                 FUN=FUN,
                 ...=...,
                 mc.cores=actual.cores)
    } else {
        lapply(X=X,
               FUN=FUN,
               ...=...)
    }
}

updateVals <- function (df, vals)
{   # Update some values on a data.frame
    for (i in names(vals)) {
        df[[i]] <- vals[[i]]
    }
    return(df)
}

dyadPos <- function (x)
    # Given an IRanges representing nucleosomes, return a vector of dyad
    # positions
    (start(x) + end(x))/2

iterDf <- function(df, fun, ...)
    # Iterate over the rows of a data.frame.
    # Possibly dlply from the plyr package is more suiting than this.
    lapply(1:nrow(df),
           function(i) do.call(fun,
                               c(unname(as.list(df[i, ])),
                                 list(...))))

###############################################################################

irLs2rd <- function(x)
    # Convert a list of IRanges to a RangedData
    RangedData(ranges=do.call(c, unname(x)),
               space=rep(names(x),
                         sapply(x, length)))

sortReads <- function (reads)
{   # Sort reads RangedData format. Sort them first by chromosome, then by
    # start and then by end
    sortChrs <- function (rans)
        rans[order(names(rans))]
    sortRans <- function (x) {
        tmp <- x[sort.list(end(x))]
        tmp[sort.list(start(tmp))]
    }
    irLs2rd(lapply(sortChrs(ranges(reads)), sortRans))
}

###############################################################################

isIn <- function (x, y)
    # Is x in y?
    start(x) <= end(y) & end(x) >= start(y)

myFilter <- function (x, f, ...)
    x[f(x, ...)]

selectReads <- function (reads, range.df)
    myFilter(reads,
             isIn,
             range(with(range.df,
                        IRanges(start,
                                end))))

makeSubsetter <- function (getChr, getStart, getEnd)
    function (x, chr=NULL, start=NULL, end=NULL) {
        if (!is.null(chr)) {
            x <- x[getChr(x) == chr, ]
            if (!is.null(start) && !is.null(end)) {
                x <- x[getStart(x) >= start & getEnd(x) <= end, ]
            }
        }
        x
    }

subsetCalls <- makeSubsetter(function (x) x$seqname,
                             function (x) x$start,
                             function (x) x$end)
subsetReads <- makeSubsetter(space, start, end)

###############################################################################

parseRange <- function (x)
{   # Parse the range string into a list containing the chromosome, start and
    # end. NULL means everything is to be selected.
    re <- "([[:alnum:]]+):([[:digit:]]+)-([[:digit:]]+)"
    if (x == "All") {
        list(chr=NULL, start=NULL, end=NULL)
    } else if (!grepl(re, x)) {
        list(chr=x, start=NULL, end=NULL)
    } else {
        chr <- gsub(re, '\\1', x)
        start <- as.integer(gsub(re, '\\2', x))
        end   <- as.integer(gsub(re, '\\3', x))
        list(chr=chr, start=start, end=end)
    }
}

###############################################################################
