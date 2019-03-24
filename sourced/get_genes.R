#!/usr/bin/Rscript

## Imports ####################################################################

sourced <- c("helperfuns", "gff_funs")
for (x in paste0(SOURCE.DIR, "/", sourced, ".R")) {
    source(x)
}

###############################################################################

genomeFile <- function (genome)
    sprintf("/orozco/services/Rdata/Web/refGenomes/%s/genes.gff",
            genome)

chromSizesFile <- function (genome)
    sprintf("/orozco/services/Rdata/Web/refGenomes/%s/%s.fa.chrom.sizes",
            genome,
            genome)

getFirstTx <- function (x, df)
{
    entries <- subset(df, ID == x)
    f <- `[[`(list("+"=which.min,
                   "-"=which.max),
              unique(entries$strand))
    entries[f(entries$start), ]
}

cleanExons <- function (df)
{
    dupls <- myFilter(df$ID, duplicated)
    sortDfBy(rbind(subset(df, !ID %in% dupls),
                   do.call(rbind,
                           lapply(dupls,
                                  getFirstTx,
                                  df))),
             c("seqname", "start"))
}

getGenes <- function (genome)
{
    #f <- genomeFile(genome)
    #gff <- readGff(f)
    gff <- readGff(genome)
    names(gff)[names(gff) == "seqname"] <- "chrom"
    gff
}

getChromSizes <- function (genome)
{
    f <- chromSizesFile(genome)
    df <- read.table(f, sep="\t")
    vals <- df$V2
    names(vals) <- df$V1
    vals
}

###############################################################################
