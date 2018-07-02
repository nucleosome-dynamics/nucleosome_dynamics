#!/bin/sh

# PATH to Rscript program
RSCRIPT=/usr/bin/Rscript
# PATH to the root of Nucleosome Dynamics package
NDDIR=/home/gelpi/DEVEL/NuclDynamics/source/distpkg
# PATH to reference genome data
PUBLICDIR=/home/gelpi/DEVEL/NuclDynamics/public_dir

echo "running readBAM (bin)"
/usr/bin/Rscript $NDDIR/R/readBAM.R --input input/cellcycleG2_chrII.bam --output cellcycleG2_chrII.RData --type paired 
#==============================================================
echo "running readBAM (bin)"
/usr/bin/Rscript $NDDIR/R/readBAM.R --input input/cellcycleM_chrII.bam --output cellcycleM_chrII.RData --type paired
#==============================================================
echo "running nucleR (bin)"
/usr/bin/Rscript $NDDIR/R/nucleR.R --input cellcycleG2_chrII.RData --output NR_cellcycleG2_chrII.gff --type paired --width 147 --minoverlap 80 --dyad_length 50 --thresholdPercentage 35 --hthresh 0.4 --wthresh 0.6 --pcKeepComp 0.02
#==============================================================
echo "running nucleR (bin)"
/usr/bin/Rscript $NDDIR/R/nucleR.R --input cellcycleM_chrII.RData --output NR_cellcycleM_chrII.gff --type paired --width 147 --minoverlap 80 --dyad_length 50 --thresholdPercentage 35 --hthresh 0.4 --wthresh 0.6 --pcKeepComp 0.02
#==============================================================
echo "running nucleR (statistics)"
/usr/bin/Rscript $NDDIR/R/nucleR_stats.R --input NR_cellcycleG2_chrII.gff --out_genes NR_cellcycleG2_chrII_genes_stats.csv --out_gw NR_cellcycleG2_chrII_stats.csv --genome $PUBLICDIR/refGenomes/R64-1-1/genes.gff
#==============================================================
echo "running nucleR (statistics)"
/usr/bin/Rscript $NDDIR/R/nucleR_stats.R --input NR_cellcycleM_chrII.gff --out_genes NR_cellcycleM_chrII_genes_stats.csv --out_gw NR_cellcycleM_chrII_stats.csv --genome $PUBLICDIR/refGenomes/R64-1-1/genes.gff
#==============================================================
echo "running NFR (bin)"
/usr/bin/Rscript $NDDIR/R/NFR.R --input NR_cellcycleG2_chrII.gff --output NFR_cellcycleG2_chrII.gff --minwidth 110 --threshold 400
#==============================================================
echo "running NFR (bin)"
/usr/bin/Rscript $NDDIR/R/NFR.R --input NR_cellcycleM_chrII.gff --output NFR_cellcycleM_chrII.gff --minwidth 110 --threshold 400
#==============================================================
echo "running NFR (statistics)"
/usr/bin/Rscript $NDDIR/R/NFR_stats.R --input NFR_cellcycleG2_chrII.gff --out_gw NFR_cellcycleG2_chrII_stats.csv --genome $PUBLICDIR/refGenomes/R64-1-1/genes.gff
#==============================================================
echo "running NFR (statistics)"
/usr/bin/Rscript $NDDIR/R/NFR_stats.R --input NFR_cellcycleM_chrII.gff --out_gw NFR_cellcycleM_chrII_stats.csv --genome $PUBLICDIR/refGenomes/R64-1-1/genes.gff
#==============================================================
echo "running txstart (bin)"
/usr/bin/Rscript $NDDIR/R/txstart.R --calls NR_cellcycleG2_chrII.gff --genome $PUBLICDIR/refGenomes/R64-1-1/genes.gff --output TSS_cellcycleG2_chrII.gff --window 300 --open_thresh 215
#==============================================================
echo "running txstart (bin)"
/usr/bin/Rscript $NDDIR/R/txstart.R --calls NR_cellcycleM_chrII.gff --genome $PUBLICDIR/refGenomes/R64-1-1/genes.gff --output TSS_cellcycleM_chrII.gff --window 300 --open_thresh 215
#==============================================================
echo "running txstart (statistics)"
/usr/bin/Rscript $NDDIR/R/txstart_stats.R --input TSS_cellcycleG2_chrII.gff --genome $PUBLICDIR/refGenomes/R64-1-1/genes.gff --out_genes TSS_cellcycleG2_chrII_genes_stats.csv --out_gw TSS_cellcycleG2_chrII_stats1.png --out_gw2 TSS_cellcycleG2_chrII_stats2.png
#==============================================================
echo "running txstart (statistics)"
/usr/bin/Rscript $NDDIR/R/txstart_stats.R --input TSS_cellcycleM_chrII.gff --genome $PUBLICDIR/refGenomes/R64-1-1/genes.gff --out_genes TSS_cellcycleM_chrII_genes_stats.csv --out_gw TSS_cellcycleM_chrII_stats1.png --out_gw2 TSS_cellcycleM_chrII_stats2.png
#==============================================================
echo "running periodicity (bin)"
/usr/bin/Rscript $NDDIR/R/periodicity.R --calls NR_cellcycleG2_chrII.gff --reads cellcycleG2_chrII.RData --type paired --gffOutput P_cellcycleG2_chrII.gff --bwOutput P_cellcycleG2_chrII.bw --genes $PUBLICDIR/refGenomes//R64-1-1/genes.gff --chrom_sizes $PUBLICDIR/refGenomes/R64-1-1/R64-1-1.fa.chrom.sizes --periodicity 165
#==============================================================
echo "running periodicity (bin)"
/usr/bin/Rscript $NDDIR/R/periodicity.R --calls NR_cellcycleM_chrII.gff --reads cellcycleM_chrII.RData --type paired --gffOutput P_cellcycleM_chrII.gff --bwOutput P_cellcycleM_chrII.bw --genes $PUBLICDIR/refGenomes/R64-1-1/genes.gff --chrom_sizes $PUBLICDIR/refGenomes/R64-1-1/R64-1-1.fa.chrom.sizes --periodicity 165
#==============================================================
echo "running periodicity (statistics)"
/usr/bin/Rscript $NDDIR/R/periodicity_stats.R --input P_cellcycleG2_chrII.gff --genome $PUBLICDIR/refGenomes/R64-1-1/genes.gff --out_genes P_cellcycleG2_chrII_genes_stats.csv --out_gw P_cellcycleG2_chrII_stats.csv
#==============================================================
echo "running periodicity (statistics)"
/usr/bin/Rscript $NDDIR/R/periodicity_stats.R --input P_cellcycleM_chrII.gff --genome $PUBLICDIR/refGenomes/R64-1-1/genes.gff --out_genes P_cellcycleM_chrII_genes_stats.csv --out_gw P_cellcycleM_chrII_stats.csv
#==============================================================
echo "running stiffness (bin)"
/usr/bin/Rscript $NDDIR/R/stiffness.R --calls NR_cellcycleG2_chrII.gff --reads cellcycleG2_chrII.RData --output STF_cellcycleG2_chrII.gff --range All
#==============================================================
echo "running stiffness (bin)"
/usr/bin/Rscript $NDDIR/R/stiffness.R --calls NR_cellcycleM_chrII.gff --reads cellcycleM_chrII.RData --output STF_cellcycleM_chrII.gff --range All
#==============================================================
echo "running stiffness (statistics)"
/usr/bin/Rscript $NDDIR/R/stiffness_stats.R --input STF_cellcycleG2_chrII.gff --genome $PUBLICDIR/refGenomes/R64-1-1/genes.gff --out_genes STF_cellcycleG2_chrII_genes_stats.csv --out_gw STF_cellcycleG2_chrII_stats1.csv --out_gw2 STF_cellcycleG2_chrII_stats2.png
#==============================================================
echo "running stiffness (statistics)"
/usr/bin/Rscript $NDDIR/R/stiffness_stats.R --input STF_cellcycleM_chrII.gff --genome $PUBLICDIR/refGenomes/R64-1-1/genes.gff --out_genes STF_cellcycleM_chrII_genes_stats.csv --out_gw STF_cellcycleM_chrII_stats1.csv --out_gw2 STF_cellcycleM_chrII_stats2.png
#==============================================================
echo "running nucDyn (bin)"
/usr/bin/Rscript $NDDIR/R/nucDyn.R --input1 cellcycleG2_chrII.RData --input2 cellcycleM_chrII.RData --outputGff ND_cellcycleG2_chrII_cellcycleM_chrII.gff --outputBigWig ND_cellcycleG2_chrII_cellcycleM_chrII.bw --plotRData ND_cellcycleG2_chrII_cellcycleM_chrII_plot.RData --genome $PUBLICDIR/refGenomes/R64-1-1/R64-1-1.fa.chrom.sizes --range All --maxDiff 70 --maxLen 140 --shift_min_nreads 3 --shift_threshold 0.1 --indel_min_nreads 3 --indel_threshold 0.05
#==============================================================
echo "running nucDyn (statistics)"
/usr/bin/Rscript $NDDIR/R/nucDyn_stats.R --input ND_cellcycleG2_chrII_cellcycleM_chrII.gff --genome $PUBLICDIR/refGenomes/R64-1-1/genes.gff --out_genes ND_cellcycleG2_chrII_cellcycleM_chrII_genes_stats.csv --out_gw ND_cellcycleG2_chrII_cellcycleM_chrII_stats.png
#==============================================================
