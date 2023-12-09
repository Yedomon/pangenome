#!/usr/bin/bash


################################################
#											   #
# Shell script is used for calculating XPCLR,  #
# 											   #
# nucleotide diversity, Fst between C, C1, C2, #
#											   #
# C3 and W population. Top 5% score regions w- #
#											   #
# ere defined as candidate selection sweeps.   #
#											   #
# Note: XP-CLR results were latter smoothed by #
#											   #
# 100-kb windows with 10-kb steps on each chr- #
#											   #
# omosome.                                     #
#											   #
################################################


# Cross-population composite likelihood ratio test(XP-CLRv1.0)
qsub -t 1-18 xpclr.sh
# XPCLR scores were smoothed
Rscript xpclr_smooth_window.R

# vcftools (v0.1.13) was used to calculate the pi, Fst, and tajimaD
qsub -t 1-18 pi_tajimaD_fst.sh

# extract the top 5% regions from pi, and Fst results
Rscript fst_top.R
Rscript pi_top.R

# using XP-CLR results as input
# filter with top 5% fst and pi by BEDtools with -wa
./region_filter.sh

# plotting the selection sweeps
Rscript whole_genome_plot.R
