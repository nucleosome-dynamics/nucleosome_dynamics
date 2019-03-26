# Nucleosome Dynamics CLI

This repository includes the set of R programs implementing 'Nucleosome Dynamics' analyses. 

<br/>

| Landing page | http://mmb.irbbarcelona.org/NucleosomeDynamics| 
| ------------ |----------------------------------------------|


----

## Table of contents

- [Nucleosome Dynamics](#Nucleosome_Dynamics)
- [Requirements](#Requirements)
- [Installation](#Installation)
- [Input data](#Input_data)
- [Running Nucleosome Dynamics CLI](#Running_Nucleosome_Dynamics_CLI)
    - [Example](#Example)
- [Usage](Usage)
    - [Analyses usage](#Analyses_usage)
        - [ReadBAM](#usage_readBAM)
        - [nucleR](#usage_nucleR)
        - [NucDyn](#usage_nucDyn)
        - [NFR](#usage_NFR)
        - [txstart -- or TSS???](#usage_txstart)
        - [Periodicity](#usage_periodicity)
        - [Stiffness](#usage_stiffness)
    - [Statistics usage](#Statistics_usage)
        - [nucleR](#usage_nucleR_stats)
        - [NucDyn](#usage_nucDyn_stats)
        - [NFR](#usage_NFR_stats)
        - [txstart -- or TSS???](#usage_txstart_stats)
        - [Periodicity](#usage_periodicity_stats)
        - [Stiffness](#usage_stiffness_stats)
- [Results](#Results)
    - [Analyses](#)
        - [NucDyn](#NucDyn)
        - [nucleR](#nucleR)
        - [Stiffness](#Stiffness)
        - [Periodicity](#Periodicity)
        - [TSS classes](#TSS_classes)

----


<a name="Nucleosome_Dynamics"></a>
# Nucleosome Dynamics

Nucleosome positioning plays a major role in transcriptional regulation DNA replication and DNA repair. The Nucleosome Dynamics CLI offers different R scripts to **analyze** nucleosome positioning from **MNase-seq** experimental data and to compare experiments to account for the transient and dynamic nature of nucleosomes under different cellular states.

The main analyses are performed with nucleR and NucDyn R packages:

* [nucleR](http://mmb.pcb.ub.es/gitlab/NuclDynamics/nucleR) performs Fourier transform filtering and peak calling to efficiently and accurately find the location of nucleosomes and classify them according to their fuzziness. 
* [NucDyn](http://mmb.pcb.ub.es/gitlab/NuclDynamics/NucDyn) detects changes in nucleosome architectures between two MNase-seq experiments. It identifies nucleosomes’ insertions, evictions and shifts at the read level.

'Nucleosome Dynamics' also offers other nucleosome-related analyses:
* NFR: location of nucleosome-free regions
* TSS: classification of transcription start sites based on the surrounding nucleosomes
* Periodicity: study of nucleosome phasing at gene level
* Stiffness: computation of stiffness of the nucleosomes derived from fitting a Gaussian function to nucleosome profiles

Additionally, after computing all the analyses, a series of statistics (plots and tables) can be computed to obtain:
* genome-wide values for all calculations
* counts and averages at the gene level

'Nucleosome Dynamics CLI' is not only offered throught the native R interface, but a number of other implementation have been prepared in order to reach also non-R users. These implementations include two web-based platforms ([MuGVRE](http://mmb.pcb.ub.es/gitlab/NuclDynamics/nucleosome_dynamics_MuGVRE) and [Galaxy](http://mmb.pcb.ub.es/gitlab/NuclDynamics/galaxy)) and two containerized installations ([docker](http://mmb.pcb.ub.es/gitlab/NuclDynamics/docker) and [singularity](https://github.com/nucleosome-dynamics/nucleosome_dynamics_singularity))




        ┌──────────────────────────────────────────────────────────┐
        │                   Nucleomose Dynamics                    │
        │                                                          │
        │       ┌───────────────┐         ┌───────────────┐        │
        │       │    nucleR     │         │     NucDyn    │        │
        │       └───────|───────┘         └────────|──────┘        │
        │               │                          │               │
        │           ┌─────────────────────────────────┐            │
        │           │    Nucleosome Dynamics CLI      │            │
        |           └┌──────────────────────┐────────┐┘            │
        │            │                      │        │             │
        │┌───────────└┐ ┌───────────┐ ┌─────└─────┐ ┌└────────────┐│
        ││   MuGVRE   │ │   Galaxy  │-│   Docker  │ │ Singularity ││
        │└────────────┘ └───────────┘ └───────────┘ └─────────────┘│
        └──────────────────────────────────────────────────────────┘


<a name="Requirements"></a>
# Requirements

* R >= 3.5
* R packages
    -  NucDyn
    -  nucleR
* UCSC wig utilities
    - wigToBigWig
    - bigWigToWig 
    - fetchChromSizes

<a name="Installation"></a>
# Installation

The instructions on how to install 'Nucleosome Dynamics CLI' are detailed at [here](INSTALL.md).

<a name="Input_data"></a>
# Input data

The primary data for 'Nucleosome Dynamics CLI' is essencially MNase-seq or ATAC-seq **aligned and sorted reads** in BAM format. Additionally, as detailed in the [Usage](#Usage) section below, some data on the reference genome is also required.
Go to [`test/data`](test) for a detailed description on such data. The folder also contains an example data set for ilustrating the R program usage.

<a name="Running_Nucleosome_Dynamics_CLI"></a>
# Running Nucleosome Dynamics CLI

Simply run each of the analysis like follows:

```sh
Rscript bin/[analysis].R [analysis_arguments]
```

Where `analysis.R` are:

| `[analysis]` | Description |
| -------- | -------- |
| readBAM         | Read Aligned MSase-seq BAM into a RData structure (required for further process) |
| nucleR           | Determine the positions of the nucleosomes across the genome |
| nucDyn          | Comparison of two diferent MNase-seq experiments to nucleosome architecture local changes |
| NFR             | Short regions depleted of nucleosomes |
| txstart         | Classify Transciption start accordint to the properties of surrounding nucleosomes  |
| periodicity  | Periodic properties of nucleosomes inside gene bodies |
| stiffness | Aparent stiffness constant foreach nucleosome obtained by fitting the coverage to a gaussian distribution |


<br/>
<br/>

Additionally, each analysis has an statistics module that creates a report (tabular or graphical) for summarising the calculation. So, after running an analysis:

```sh
Rscript statistics/[analysis_stats].R [analysis_stats_arguments]
```

Available `analysis_stats`.R are:

| `[analysis_stats]` | Description |
| -------- | -------- |
| nucleR_stats|  Nucleosome call statistics|
| NFR_stats|             Nucleosome Free Regions statistics|
| txstart_stats|         TSS and TTS statistics|
| periodicity_stats|  Statistics on Nucleosome periodicity|
| stiffness_stats|      Statistics on stiffness|
| nucDyn_stats|    Statistics on Nucleosome Dynamics analysis|

<br/>

The [Usage](#Usage) section below describes the arguments for each individual analysis and similary, the [Results](#Results) section includes a description of the sequence annotations files (GFF) generated.


<a name="Example"></a>
### Example
For calling nucleosome positions in a pair of aligned MNase-seq BAM files, 'nucleR' is the analysis of interest. As input, it takes the MNase-seq reads in Rdata format, so first, we need load such reads into an RData file using `readBAM`:

```sh
bin/readBAM.R --input test/data/cellcycleG2_chrII.bam --output test/data/cellcycleG2_chrII.RData --type paired
```
Once the RData is generated, run the calling with `nucleR` as follows:

```sh
bin/nucleR.R --input test/data/cellcycleG2_chrII.RData --output test/data/NR_cellcycleG2_chrII.gff --type paired --width 147 --minoverlap 80
```

For extracting an statistical report of the 'nucleR', feed the resulting GFF file into `statistics/nucleR`:

```sh
statistics/nucleR.R --input test/data/NR_cellcycleG2_chrII.gff --genome test/data/refGenomes/R64-1-1/genes.gff --out_genes test/data/NR_cellcycleG2_chrII_genes.csv --out_gw test/data/NR_cellcycleG2_chrII_genes.csv
```

Check the argument descriptions for each particular analysis at the [Usage](#usage_nucleR) section. The 'nucleR' execution will generate the NR_cellcycleG2_chrII.gff track file with the predicted positions an a number of scores and other annotations described in the [Results](#nucleR) section. Most of the sequence browsers can visualized it, for instance, [JBrowse](https://jbrowse.org/).


<a name="Usage"></a>
# Usage

<a name="Analyses_usage"></a>
## Analyses usage

<a name="usage_readBAM"></a>
**readBAM** --input {bam} --output {RData} --type (single|paired)

        --input
                Sequence file in BAM format
        --output
                RData structure
        --type
                Type of sequence data: single|paired

<a name="usage_nucleR"></a>
**nucleR** --input {RData} --output {gff} --type (single|PAIRED) --width {int} --minoverlap {int}
[--dyad_length {int} --thresholdPercentage {double} --hthresh {double} --wthresh {double} --pcKeepComp {double} --fdrOverAmp {double} --components {int} --fragmentLen {int} --trim {int} --threshold {logical} --thresholdValue {int} --chr STR --start {int} --end {int} ]

        --input: 
                Input BAM file (RData format)
        --output:
                Nucleosome calls in GFF format. Annotations: Score_weight (0-1), Score_height (0-1), class (W, F, uncertain)
        --type
                Type of sequence data: single|paired
        --width: 
                Width given to nucleosome calls previous to merging (bp). Default 147 
        --minoverlap: 
                Minimum number of overlapping base pairs in two nucleosome calls for them to be merged into one (bp).
        --dyad_length: 
                Length of the reads that should be used for the nucleosome calling to define the dyad of the nucleosomes, keeping that number of bases around the center of   each fragment (bp). Optional. Default 50 bp.
        --hthresh: 
                Height threshold (between 0 and 1) to classify (in combination with width threshold) a nucleosome as either fuzzy or well-positioned according to the number of reads in the dyad of the nucleosome call. Nucleosomes below this value will be defined as fuzzy. Optional, default 0.4.
        --wthresh: 
                Width threshold (between 0 and 1) to classify (in combination with height threshold) a nucleosome as either fuzzy or well-positioned according to the disperion of the reads around the dyad. Nucleosomes below this value (that is, nucleosome calls not sharp enough) will be defined as fuzzy. Optional, default 0.4
        --pcKeepComp
                Parameter used in the coverage smoothing when Fourier transformation is applied. Number of components to select with respect to the total size of the sample. Allowed values are numeric (in range 0:1) for manual setting, or auto for automatic detection. Optional, default 0.02.
        --fdrOverAmp
                Threshold to filter over-amplified reads , as defined in filterDuplReads function of htSetqTools R package. Optional, default 0.05. 
        --components
                Number of negative binomials that will be used to filter duplicated reads, as defined in filterDuplReads function of htSetqTools R package. Optional, default 1.
        --fragmentLen
                Fragment Length (bp). Optional, default 170.
        --trim
                Optional, default 50
        --threshold
                 Threshold to filter out which nucleosome calls are reported as significant enough. It can be supplied as either an absolute value (--thresholdValue) or as a relative value expressed as a percentage (--thresholdPercentage). If set to TRUE, the percentage value is considered. If set to FALSE, the absolute value does. Optional, default TRUE.
        --thresholdValue
                 Absolute value to filter out nucleosome calls. It is the minimum number of reads (coverage) in a nucleosome call expressed as reads per million of mapped reads. Optional, default 10.
        --thresholdPercentage: 
                Percentile of coverage in the experiment used as threshold to filter out nucleosome calls (i.e., '25%' would mean that only peaks with coverage in the 1st quantile would be considered). Optional, default 35(%)
        --chr
                Chromosome to consider for the analysis in the given input file. By default, all the genomic range is considered. Optional.
        --start
                Start genomic position to consider for the analysis in the given input file. By default, all the genomic range is considered. Optional.
        --end
                End genomic position to consider for the analysis in the given input file. By default, all the genomic range is considered. Optional.


<a name="usage_nucDyn"></a>
**nucDyn**  --input1 {RData} --input2 {RData} --calls1 {gff} --calls2 {gff} --outputGff {gff} --outputBigWig {bw}  --genome {chrom.sizes} --range {str}
[ --plotRData {RData} --maxDiff {int} --maxLen {int} --shift_min_nreads {int} --shift_threshold {double} --indel_min_nreads {int} --indel_threshold {double} --cores {int} --equal_size (logical) --readSize {int} ]

        -- input1, --input2
                Input BAM from MNase-seq in RData format (from readBAM)
        --calls1, --calls2
                Nucleosome calls in GFF format as obtained from NucleR
        -outputGff
                Output of NucleosomeDynamics in GFF format
        --outputBigWig {bw}
                Output of NucleosomeDynamics in BigWig format: -log10 of the p-value of the significance of the differences found.
        --genome {chrom.sizes}
                Chromosome sizes from reference Genome
        --range (All|chr|chr:start-end)
                Genomic range to be analyzed. All: all chromosomes; chr: a Whole chromosome; chr:start-end: a specific region given the coordinates. Default All
        --plotRData {RData}

        --maxDiff
                Maximum distance between the centers of two fragments for them to be paired as shifts. (bp) Optional,default 70bp
        --maxLen
                This value is used in a preliminar filtering. Fragments longer than this will be filtered out, since they are likely the result of MNase under-digestion and represent two or more nucleosomes. (bp) Optional, default 140bp
        --shift_min_nreads
                Minimum number of shifted reads for a shift hostspot to be reported {int}, optional, default 3
        --shift_threshold
                Threshold applied to the shift hostpots. Only hotspots with a score better than the value will be reported. Notice the score has to be lower than the threshold, since these numbers represent p-values. Optional, default 0.1
        --indel_min_nreads
                Minimum number of shifted reads for an indel hostspot to be reported {int}, optional, default 3
        --indel_threshold [float, 0.05]
                Threshold applied to the indel hostpots. Only hotspots with a score better than the value will be reported. Notice the score has to be lower than the threshold, since these numbers represent p-values. Optional, default 0.05
        --cores 
                Number of computer threads. Optional, default 1
        --equal_size
                Set all fragments to the same size. Optional, default FALSE
        --readSize  
                Length to which all reads will be set in case `equalSize` is `TRUE`. It is ignored when `equalSize` is set to `FALSE`. Default 140


<a name="usage_NFR"></a>
**NFR**  --input {gff} --output {gff} [--minwidth {int} --threshold {int}  ]

        --input 
                Nucleosome calls in GFF format as obtained from NucleR
        --output {gff} 
                Nucleosome Free Regions in GFF format
        --minwidth 
                Minimum length (bp). Optional, default 110bp
        --threshold 
                Maximum length (bp). Optional, default 140bp


<a name="usage_txstart"></a>
**txstart**   ---calls {gff} --genome {gff} --output {gff} --window {int} --open_thresh {int} [ --cores {int} --p1.max.downstream {int} ]


        --calls
                Nucleosome calls as obtained form NucleR. GFF format
        --genome
                Gene positions from the reference genome. GFF Format
        --output
                Classification of TSS according to nucleosome types,  in GFF format
        --window 
                Number of nucleotides on each side of the TSS where -1 and +1 nucleosomes are searched for. Default 300
        --open_thresh
                Distance between nucleosomes -1 and +1 to discriminate between 'open' and 'close' classes. Default 215
        --cores
                Number of computer threads. Optional, default 1
        --p1.max.downstream
                Maximum distance to look for +1 nucleosome relative to the TSSOptional. Optional, default 20


<a name="usage_periodicity"></a>
**periodicity** --calls {gff} --reads {RData} --type (single|paired) --gffOutput {gff} --bwOutput {bw} --genes {gff} --chrom_sizes {chrom.sizes}
[--periodicity {int} --cores {int} ]

        --calls
                Nucleosome calls in GFF format as obtained from NucleR
        --reads
                Sequence Reads in RData format as obtained from readBAM
        --type 
                Type of reads (Single|paired)
        --gffOutput
                Periodicity output in GFF format
        --bwOutput
                Periodicity output in BigWig format
        --genes
                Position of genes in reference genome,in GFF format
        --chrom_sizes
                Chromosome sizes file in the reference genome
        --periodicity
                Average distance between two consecutive nucleosomes. It is used as the period of the nucleosome coverage signal. Optional, default 165
        --cores
                Number of computer threads. Optional, default 1


<a name="usage_stiffness"></a>
**stiffness** --calls {gff} --reads {RData} --output {gff} --range {str}  [--t {double}] ]

        --calls
                Nucleosome calls in GFF format as obtained from NucleR
        --reads
                Sequence data in RData format as obtained from readBAM
        --output
                Output stiffness in GFF format
        --range
                Genomic range to consider. Format: [str, All|chr|chr:start-end]  whre 'All' is all chromosomes, 'chr' is a single chomosome, and 'chr:start-end' is the  range indicated by the coordinates. Default All
        --t
                Temperature (K). Optional, default 310.15



<a name="Statistics_usage"></a>
## Statistics usage

<a name="usage_nucleR_stats"></a>
**nucleR_stats** --input {gff} --genome {gff} --out_genes {csv} --out_gw {csv}

        --input
                Nucleosome calls in GFF format as obtained from NucleR
        --genome
                Gene positions from the reference genome. GFF Format
        --out_genes
                Output file containing statistics of all calculations for each gene. CSV format 
        --out_gw
                Output file containing table of genome wide statistics. CSV format 


<a name="usage_nucDyn_stats"></a>
**nucDyn_stats** --input {gff} --genome {gff} --out_genes {csv} --out_gw {png}

        --input 
                Nucleosome calls in GFF format as obtained from stiffness_stats
        --genome
                Gene positions from the reference genome. GFF Format
        --out_genes
                Output file containing statistics of all calculations for each gene. CSV format
        --out_gw
                Output file containing table of genome wide statistics. PNG format


<a name="usage_NFR_stats"></a>
**NFR_stats** --input {gff}  --genome {gff} --out_gw {csv}
        --input
                Nucleosome calls in GFF format as obtained from NFR
        --genome 
                Gene positions from the reference genome. GFF Format
        --out_gw 
                Output file containing table of genome wide statistics. CSV format


<a name="usage_txstart_stats"></a>
**txstart_stats** --input {gff} --genome {gff} --out_genes {csv} --out_gw {png} --out_gw2 {png}

        --input 
                Nucleosome calls in GFF format as obtained from txstart
        --genome
                Gene positions from the reference genome. GFF Format
        --out_genes
                Output file containing statistics of all calculations for each gene. CSV format
        --out_gw
                Output file containing table of genome wide statistics. PNG format
        --out_gw2
                Output file containing table of genome wide statistics. PNG format


<a name="usage_periodicity_stats"></a>
**periodicity_stats** --input {gff} --genome {gff} --out_genes {csv} --out_gw {csv}

        --input 
                Nucleosome calls in GFF format as obtained from periodicity
        --genome
                Gene positions from the reference genome. GFF Format
        --out_genes
                Output file containing statistics of all calculations for each gene. CSV format
        --out_gw
                Output file containing table of genome wide statistics. CSV format


<a name="usage_stiffness_stats"></a>
**stiffness_stats** -input {gff} --genome {gff} --out_genes {csv} --out_gw {csv} --out_gw2 {png}

        --input 
                Nucleosome calls in GFF format as obtained from stiffness_stats
        --genome
                Gene positions from the reference genome. GFF Format
        --out_genes
                Output file containing statistics of all calculations for each gene. CSV format
        --out_gw
                Output file containing table of genome wide statistics. CSV format
        --out_gw2
                Output file containing table of genome wide statistics. PNG format


<a name="Results"></a>
# Results

Most of 'Nucleosome Dynamics CLI' output files are **sequence annotation files**, either in GFF or BW format. They are ready to be visualized as sequence tracks in most modern sequence browsers.
In addition, for each track hit further information is offered to help analysing the rellevance of each individual result. Following, an explanation on what represents each analysis hit, and the list of attributes accompanying then.

<a name="NucDyn"></a>
##### NucDyn

Primary data
* Position: region where a nucleosome movement is detected
* Type: change in the nucleosome map
* Score: magnitude of the change

Attributes
* class: type of hotspot (see help for all possible types)
* nreads: number of reads involved in this movement


<a name="nucleR"></a>
##### nucleR

Primary data
* Score: Positionning score. It is calculated as the weighted sum of width and height scores.

Attributes
* score_width: Witdth score. It is a measure of how sharp a peak is. A value of 0 would be an extremely wide peak and a value of 1 a very sharp one.
* score_height: Height score. Tells how large a peak of a nucleosome is. The bigger this number, the higher the peak.
* class: Whether the nucleosome is well-positioned (W) or fuzzy (F) or undetermined. The taken value depends on score_h and score_w. Undetermined means the exact position of the nucleosome cannot be determined due to strong fuzziness.


<a name="Stiffness"></a>
##### Stiffness

Primary data
* Score: Stiffness estimation. It represents the energy required to move the nucleosome (expressed in kcal/mol/bp), it's derived from the Gaussian standard deviation.

Attributes
* nucleR_score: the nucleR score given to that nucleosome
* nucleR.class: the nucleR class given to that nucleosome
* gauss_k: the height of the peak of the gaussian curve
* gauss_m: the position of the peak of the gaussian curve
* gauss_sd: the standard deviation of the gaussian curve


<a name="Periodicity"></a>
##### Periodicity

Primary data
* Position: gene position (from TSS to TTS) 

Attributes
* nucleosome_first: First nucleosome of the gene.
* nucleosme_last: Last nucleosome of the gene.
* score_phase: Is a measure of the phase between the first and the last nucleosome. A score of 0 means the nucleosome are completely phased and a score of 82 corresponds to totally antiphased nucleosomes.
* score_autocorrelation: It is directly computed from the experimental coverage and is quantitative measure of the periodicity of nucleosomes inside the gene body.


<a name="TSS_classes"></a>
##### TSS classes

Primary data
* Position: Region between the dyads of two nucleosomes surrounding the TSS.

Attributes
* classification: Descriptor of the Transcription Start Site. See the help for possible options.
* distance: Distance in base pairs between the nucleosome +1 and the nucleosome -1.
* nucleosome minus1: Position of the nucleosome -1.
* nucleosome plus1: Position of the nucleosome +1
* TTS_position: Position of the Transcription Start Site.

