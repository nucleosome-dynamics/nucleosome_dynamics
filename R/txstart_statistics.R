#!/usr/bin/Rscript


## Imports ####################################################################

library(getopt)
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
                 "out_gw",    "d", 1, "character",
                 "out_gw2",   "e", 1, "character"),
               byrow=TRUE,
               ncol=4)
params <- getopt(spec)

## Read genes #################################################################

message("-- loading used genome")
genes <- getGenes(params$genome)

genes$tss <- as.numeric(genes$tss)
genes$tts <- as.numeric(genes$tts)

## Statistics per gene ########################################################

message("-- computing statistics per gene")
tss <- readGff(params$input)

genes_out <- tss[, c("id", "classification", "distance")]
genes_out <- merge(genes[, c("name", "tss")],
                   genes_out,
                   by.x="name",
                   by.y="id",
                   all.x=TRUE)
genes_out <- genes_out[, c(1, 3, 4)]
genes_out <- rbind(c("Name", "TSS class", "Distance from -1 to +1"), genes_out)

write.table(genes_out,
            params$out_genes,
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE,
            sep=",")

## Statistics genome-wide  ####################################################

message("-- computing statistics genome-wide")

#--- Plot 1 ---

tab_tss <- table(tss$classification)
df1 <- as.data.frame(tab_tss)

x <- max(df1$Freq)
lim <- x + (0.05*x)
p1 <- ggplot(df1, aes(x=Var1, y=Freq)) +
    geom_bar(stat="identity", fill="#66A61E") +
    geom_text(aes(label=Freq), hjust=0, size=1) +
    coord_flip() +
    ylim(0, lim) +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          text=element_text(size=4))

ggsave(filename = params$out_gw,
       plot     = p1,
       width    = 41,
       height   = 41,
       units    = "mm",
       dpi      = 300)

#--- Plot 2 ---

x <- as.numeric(tss$distance)
x <- x[!is.na(x)]
df2 <- data.frame(x=x)

p2 <- ggplot(df2, aes(x=x)) +
    geom_line(stat="density", color="#1B9E77", lwd=0.5) +
    ylim(0, 0.01) +
    xlim(0, 450) +
    labs(x="Width", y="Density") +
    ggtitle("Distribution of NFR width around TSS") +
    theme_bw() +
    theme(text=element_text(size=4))

ggsave(filename = params$out_gw2,
       plot     = p2,
       width    = 41,
       height   = 41,
       units    = "mm",
       dpi      = 300)
