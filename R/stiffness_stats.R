#!/usr/bin/Rscript

## Imports ####################################################################

suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ggplot2))

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

cat("-- loading reference genome\n")
genes <- getGenes(params$genome)

genes$tss <- as.numeric(genes$tss)
genes$tts <- as.numeric(genes$tts)

genes_gr <- GRanges(genes$chrom,
                    IRanges(start=genes$start, end=genes$end),
                    name=genes$name)

## Statistics per gene ########################################################

cat("-- computing statistics per gene\n")
stf <- readGff(params$input)

stf_gr <- GRanges(stf$seqname,
                  IRanges(start=stf$start, end=stf$end),
                  score=stf$score)

ovlps <- findOverlaps(genes_gr, stf_gr)
split_ovlps <- split(subjectHits(ovlps), queryHits(ovlps))

res <- lapply(
    names(split_ovlps),
    function(i) {
        tmp <- stf_gr[split_ovlps[[i]]]
        return(c(genes_gr$name[as.numeric(i)],
                 mean(tmp$score, na.rm=TRUE),
                 sd(tmp$score, na.rm=TRUE)))
    }
)

res <- do.call(rbind, res)
colnames(res) <- c("name", "Mean_STF", "StdDev_STF")

stat_stf <- merge(genes, res, by="name", all.x=T)
colnames(stat_stf)[colnames(stat_stf)=="name"] <- "Name"

i <- c("Name","Mean_STF","StdDev_STF")
write.table(stat_stf[, i],
            params$out_genes,
            row.names = FALSE,
            quote     = FALSE,
            sep       = ",")

## Statistics genome-wide  ####################################################

#--- Mean and std.dev ---
cat("-- computing statistics genome-wide\n")

rownames <- c("Mean stiffness", "Std. Dev. stiffness")
a <- mean(stf_gr$score, na.rm=TRUE)
b <- sd(stf_gr$score, na.rm=TRUE)
gw_stat <- data.frame(c("Statistic", rownames),
                      c("Value", round(c(a, b), 4)))

write.table(gw_stat,
            params$out_gw,
            row.names = FALSE,
            col.names = FALSE,
            quote     = FALSE,
            sep       = ",")

#--- Plot distribution ---

stf$class <- NA
stf$class[stf$score <  0.1] <- "0 - 0.1"
stf$class[stf$score >= 0.1 & stf$score < 0.2] <- "0.1 - 0.2"
stf$class[stf$score >= 0.2 & stf$score < 0.3] <- "0.2 - 0.3"
stf$class[stf$score >= 0.3 & stf$score < 0.4] <- "0.3 - 0.4"
stf$class[stf$score >= 0.4] <- "0.4 - 1"

df <- as.data.frame(table(stf$class)/sum(table(stf$class)))

colors <- c("0 - 0.1"   = "#98F5FF",
            "0.1 - 0.2" = "#71DAE2",
            "0.2 - 0.3" = "#4CC0C4",
            "0.3 - 0.4" = "#26A5A8",
            "0.4 - 1"   = "#008B8B")

p <- ggplot(df, aes(Var1, Freq)) +
    geom_bar(stat="identity", aes(fill=Var1)) +
    scale_fill_manual(values=colors) +
    labs(x="Stiffness", y="Proportion of genes") +
    theme_bw() +
    theme(legend.position="none", text=element_text(size=3))

ggsave(filename = params$out_gw2,
       plot     = p,
       width    = 41,
       height   = 41,
       units    = "mm",
       dpi      = 300)
