#!/usr/bin/Rscript

## Imports ####################################################################

suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(plyranges))

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

genes <- read_gff3(params$genome)

genes$tss <- as.numeric(genes$tss)
genes$tts <- as.numeric(genes$tts)

genes_gr <- GRanges(as.vector(seqnames(genes)),
                    IRanges(start=start(genes), end=end(genes)),
                    name=genes$name)

## Statistics per gene ########################################################

message("-- computing statistics per gene")

nuc <- readGff(params$input)

nuc_gr <- GRanges(nuc$seqname,
                  IRanges(start=nuc$start, end=nuc$end),
                  class=nuc$class)

nuc_well <- nuc_gr[nuc_gr$class == "W"]
nuc_fuzzy <- nuc_gr[nuc_gr$class == "F"]
nuc_uncertain <- nuc_gr[nuc_gr$class == "uncertain"]

genes$TotalNucleosomes <- countOverlaps(genes_gr, nuc_gr)
genes$TotalWellPositioned <- countOverlaps(genes_gr, nuc_well)
genes$TotalFuzzy <- countOverlaps(genes_gr, nuc_fuzzy)
genes$TotalUncertain <- countOverlaps(genes_gr, nuc_uncertain)

cols <- c("name",
          "TotalNucleosomes",
          "TotalWellPositioned",
          "TotalFuzzy",
          "TotalUncertain")
genes_out = genes[, cols]

# cols <- c("Name",
#           "Total Nucleosomes",
#           "Total Well-Positioned",
#           "Total Fuzzy",
#           "Total Uncertain")
# genes_out = rbind(cols, genes_out)

write.table(genes_out,
            params$out_genes,
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE,
            sep=",")

## Statistics genome-wide  ####################################################

message("-- computing statistics genome-wide")
nuc_gr$class = as.character(nuc_gr$class)
class_lab = gsub("F", "Fuzzy", nuc_gr$class)
class_lab = gsub("W", "Well-positioned", class_lab)
class_lab = gsub("uncertain", "Uncertain", class_lab)

gw_stat <- as.data.frame(table(class_lab))
colnames(gw_stat) = c("Class", "Frequency")
gw_stat$Class = as.character(gw_stat$Class)
gw_stat <- rbind(gw_stat, c("Total", length(nuc_gr)))

write.csv(gw_stat, params$out_gw, row.names=FALSE, quote=FALSE)
