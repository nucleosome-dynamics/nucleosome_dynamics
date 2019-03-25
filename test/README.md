# Data for running Nucleosome Dynamics

We describe here the files required to run Nucleosome Dynamics and provide sample input data that can be used to test the code. The directory `test/data` includes:

1. Sample MNase-seq reads
2. Reference genome associated data

#### 1. MNase-seq reads

`test/data/cellcycle*` BAM files contain mapped reads to R64-1-1 genome assembly, from an MNaseSeq experiment published in Deniz et al. (2016) for S. cerevisiae chrII in G2 and M cell cycle phases. 

YBD
BAM files should be indexed to be read with nucleR function readBAM function?
RData and BAI files are needed?
Do we want to include output sample data?

#### 2. Reference genome associated data

Two files are requiered for the corresponding genome assembly:

* Chromosome sizes: two-column tab-separated text file containing assembly sequence names and sizes. It can be downloaded from [UCSC](http://hgdownload.soe.ucsc.edu/downloads.html)
* Gene TSS and TTS coordinates: GFF3 file containing information for every gene. For a full description of the fields in GFF3 files see [ensemle](https://www.ensembl.org/info/website/upload/gff3.html) website. The following columns are required:
  * seqid: name of the chromosome in the corresponding genome assembly
  * source: name of the program that generated this feature, or the data source. If not provided, it should be denoted with '.'
  * type: type of feature. If not provided, it should be denoted with '.'
  * start: start point of the feature. Should correspond with TSS for genes in the + strand, TTS for the genes in the + strand
  * end: end point of the feature. Should correspond with TTS for genes in the + strand, TSS for the genes in the + strand
  * score: not required (it should be denoted with '.')
  * strand: defined as + (forward) or - (reverse)
  * phase: not required (it should be denoted with '.')
  * attributes: a semicolon-separated list of tag-value pairs, providing additional information about each feature. At least three values are required: 
    * name: gene identifier
    * tss: TSS of the gene
    * tts: TTS of the gene

We provide the files for S. cerevisiae genome assembly R64-1-1 in `test/data/refGenomes/R64-1-1`. Chromosome sizes were downloaded from [SGD](http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz) and gene TSS and TTS were compiled from Pelechano et al. (2013), Miura et al. (2006), Yassour et al. (2009), Nagalakshmi et al. (2008).


# References 
Deniz, O., Flores, O., Aldea, M., Soler-López, M., and Orozco, M. (2016). Nucleosome architecture throughout the cell cycle. Scientific Reports 6, 19729.

Pelechano, V., Wei, W., Steinmetz, L.M. (2013). "Extensive transcriptional heterogeneity revealed by isoform profiling." Nature, 497, pp. 127–131. doi: 10.1038/nature12121

Miura, F., Kawaguchi, N., Sese, J., Toyoda, A., Hattori, M., Morishita, S., Ito, T. (2006). "A large-scale full-length cDNA analysis to explore the budding yeast transcriptome." Proc. Natl. Acad. Sci. U.S.A., 103, pp. 17846–17851. doi: 10.1073/pnas.0605645103

Yassour, M., Kaplan, T., Fraser, H.B., Levin, J.Z., Pfiffner, J., Adiconis, X., Schroth, G., Luo, S., Khrebtukova, I., Gnirke, A., Nusbaum, C., Thompson, D.-A., Friedman, N., Regev, A. (2009). "Ab initio construction of a eukaryotic transcriptome by massively parallel mRNA sequencing." Proceedings of the National Academy of Sciences, 106, pp. 3264–3269. doi: 10.1073/pnas.0812841106

Nagalakshmi, U., Wang, Z., Waern, K., Shou, C., Raha, D., Gerstein, M., Snyder, M. (2008). "The Transcriptional Landscape of the Yeast Genome Defined by RNA Sequencing." Science, 320, pp. 1344–1349. doi: 10.1126/science.1158441

