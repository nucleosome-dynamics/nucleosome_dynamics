#!/usr/bin/Rscript

## Imports ####################################################################

library(getopt)
library(htSeqTools)
library(nucleR)
library(IRanges)
library(GenomicRanges)
library(ggplot2)

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
sourced <- c("get_genes", "gff_funs")
for (x in sourced) {
    source(paste0(SOURCE.DIR, "/", x, ".R"))
}

## Parameters and Arguments ###################################################

spec <- matrix(c("input",     "a", 1, "character",
                 "genome",    "b", 1, "character",
                 "out_genes", "c", 1, "character",
                 "out_gw",    "d", 1, "character"),
               byrow=TRUE,
               ncol=4)

params <- getopt(spec)

## Read genes #################################################################

message("-- loading used genome")
genes <- getGenes(params$genome)

genes$tss <- as.numeric(genes$tss)
genes$tts <- as.numeric(genes$tts)

genes_gr <- GRanges(genes$chrom,
                    IRanges(start=genes$start, end=genes$end),
                    name=genes$name)

## Statistics per gene ########################################################

message("-- computing statistics per gene")
nd <- readGff(params$input)

if (nrow(nd) > 0) {
    nd_gr <- GRanges(nd$seqname,
                     IRanges(start=nd$start, end=nd$end),
                     class=nd$class)
} else {
    nd_gr <- GRanges()
}

incl    <- nd_gr[nd_gr$class == "INCLUSION"]
evic    <- nd_gr[nd_gr$class == "EVICTION"]
shift_p <- nd_gr[nd_gr$class == "SHIFT +"]
shift_m <- nd_gr[nd_gr$class == "SHIFT -"]

genes$nIncl    <- countOverlaps(genes_gr, incl)
genes$nEvic    <- countOverlaps(genes_gr, evic)
genes$nShift_p <- countOverlaps(genes_gr, shift_p)
genes$nShift_m <- countOverlaps(genes_gr, shift_m)

i <- c("name",
       "nIncl",
       "nEvic",
       "nShift_p",
       "nShift_m")
stat_nd <- genes[, i]

stat_nd <- rbind(c("Name",
                   "Inclusions",
                   "Evictions",
                   "Shifts+",
                   "Shifts-"),
                 stat_nd)

write.table(stat_nd,
            params$out_genes,
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE,
            sep=",")

## Statistics genome-wide  ####################################################

#--- Mean and std.dev ---
message("-- computing statistics genome-wide")

nd_tab <- table(nd$class) / sum(table(nd$class))
df <- as.data.frame(nd_tab)

if (nrow(df) > 0) {
    names(df)[names(df) == "Freq"] <- "Proportion"

    levels(df$Var1)[match("EVICTION", levels(df$Var1))] <- "Eviction"
    levels(df$Var1)[match("INCLUSION", levels(df$Var1))] <- "Inclusion"
    levels(df$Var1)[match("SHIFT -", levels(df$Var1))] <- "Shift -"
    levels(df$Var1)[match("SHIFT +", levels(df$Var1))] <- "Shift +"
} else {
    df <- data.frame(Var1=character(), Proportion=numeric())
}

p <- ggplot(df, aes(Var1, Proportion)) +
    geom_bar(stat="identity", aes(fill=Var1)) +
    scale_fill_manual(values=c("Eviction"  = "#04B431",
                               "Inclusion" = "#FF0000",
                               "Shift -"   = "#8000FF",
                               "Shift +"   = "#0040FF")) +
    theme_bw() +
    theme(legend.position="none",
          axis.title.x=element_blank(),
          text=element_text(size=3))

ggsave(filename=params$out_gw,
       plot=p,
       width=41,
       height=41,
       units="mm",
       dpi=300)
