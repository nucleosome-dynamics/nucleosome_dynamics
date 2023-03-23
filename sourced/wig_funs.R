#!/usr/bin/Rscript

## Imports ####################################################################

suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
source(paste(SOURCE.DIR, "get_genes.R", sep="/"))

## Binary paths ###############################################################

wig.dir <- ""#paste(where(), "../wig_utils", sep="/")
towig.bin <- paste(wig.dir, "bigWigToWig", sep="")
tobig.bin <- paste(wig.dir, "wigToBigWig", sep="")

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
# This function takes a list of vectors (x) representing signal scores across genomic regions and writes the data to a bigWig file format
# Additionally it uses a file path to write the resulting bigWig file (outf),
# and a file path to a file with chromosome sizes in tab seperated form (chrom.sizes.f)
writeBigWig_updated <- function (x, outf, chrom.sizes.f) {
  
  # Map over the input vectors to create a tibble for each vector with four columns:
  # score, start position, width, and strand (* indicates unstranded).
  # Combine all tibbles into a single data frame using map_df().
  # Use .id = "seqnames" to create a new column seqnames, which contains the names of the vectors in x.
  gr <- purrr::map_df(x,
                      function(y) 
                      {dplyr::tibble(score = y,
                                     start = 1:length(y),
                                     width = 1,
                                     strand = "*")},
                      .id = "seqnames") %>% 
    as_granges()
  
  # Read in a file containing chromosome sizes into a data frame.
  # The file is assumed to have two columns separated by tabs: chromosome name and size.
  chrSizes <- readr::read_delim(chrom.sizes.f,
                                col_names = c("seqnames","size"),
                                delim = "\t",
                                show_col_types = FALSE)
  
  # Filter the chromosome sizes data frame to include only the chromosomes that are present in the granges object.
  chrSizes <- filter(chrSizes, seqnames %in% seqnames(gr)) 
  
  # Set the sequence lengths of the granges object to the chromosome sizes in chrSizes.
  seqlengths(gr) <- chrSizes$size
  
  # Write the granges object to a bigWig file at the specified output path.
  write_bigwig(gr, outf)
}

###############################################################################
