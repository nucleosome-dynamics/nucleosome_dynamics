#!/usr/bin/Rscript

## Imports ####################################################################

suppressPackageStartupMessages(library(IRanges))
source(paste(SOURCE.DIR, "get_genes.R", sep="/"))

## Binary paths ###############################################################

wig.dir <- paste(where(), "../wig_utils", sep="/")
towig.bin <- paste(wig.dir, "bigWigToWig", sep="/")
tobig.bin <- paste(wig.dir, "wigToBigWig", sep="/")

## Reading functions ##########################################################

readWig <- function(inf)
{
    lines <- readLines(inf)
    sep.idxs <- grep("^fixedStep", lines)

    vals <-  mapply(function(i, j) as.numeric(lines[i:j]),
                    sep.idxs + 1,
                    c(sep.idxs[-1] - 1, length(lines)),
                    SIMPLIFY=FALSE)

    ids <- lines[sep.idxs]

    xs <- unlist(strsplit(ids, " "))

    getVal <- function(xs, a)
        sub(paste0(a, "="), "", grep(a, xs, value=TRUE))

    chrs <- getVal(xs, "chrom")
    pos <- as.numeric(getVal(xs, "start"))

    doChr <- function(chr) {
        chr.loc <- chrs == chr
        chr.pos <- pos[chr.loc]
        chr.vals <- vals[chr.loc]
        max.pos <- max(chr.pos)
        tail.size <- length(chr.vals[[which(chr.pos == max.pos)]]) - 1
        s <- max.pos + tail.size
        x <- rep(0, s)
        for (i in seq_along(chr.vals)) {
            from <- chr.pos[i]
            to <- from + length(chr.vals[[i]]) - 1
            x[from:to] <- chr.vals[[i]]
        }
        Rle(x)
    }

    chroms <- unique(chrs)
    cover <- lapply(chroms, doChr)
    names(cover) <- chroms
    cover
}

readBigWig <- function(inf)
{
    wigf <- sub(".bw$", ".wig", inf)
    system(paste(towig.bin, inf, wigf))
    cover <- readWig(wigf)
    file.remove(wigf)
    cover
}

## Writing functions ##########################################################

splitAtZeros <- function (xs)
{
    if (length(xs) == 0) {
        return(list())
    } else {
        xs <- as.vector(xs)
        counts <- rle(xs != 0)

        jdxs <- cumsum(counts$length)
        idxs <- c(1, jdxs[-length(jdxs)] + 1)

        by.kinds <- mapply(function(i, j) xs[i:j],
                           idxs, jdxs,
                           SIMPLIFY=FALSE)
        names(by.kinds) <- idxs

        return(by.kinds[counts$values])
    }
}

writeWig <- function(x, outf)
{
    x <- x[sapply(x, length) != 0]
    tag.rows <- unlist(lapply(names(x),
                              function(chr)
                                  paste("fixedStep",
                                        paste0("chrom=", chr),
                                        paste0("start=",
                                               format(as.numeric(names(x[[chr]])),
                                                      trim=TRUE,
                                                      scientific=FALSE)),
                                        "step=1",
                                        sep=" ")))
    vals <- unlist(x, recursive=FALSE)
    nonsc.vals <- lapply(vals,
                         format,
                         trim=TRUE,
                         scientific=FALSE)
    idx <- order(c(seq_along(tag.rows), seq_along(nonsc.vals)))
    cat(unlist(c(tag.rows, nonsc.vals)[idx]),
        sep="\n",
        file=outf)
}

writeBigWig <- function (x, outf, chrom.sizes.f)
{
    wigf <- sub(".bw$", ".wig", outf)
    writeWig(x, wigf)
    system2(tobig.bin, shQuote(c(wigf, chrom.sizes.f, outf)))
    file.remove(wigf)
    invisible()
}

###############################################################################
