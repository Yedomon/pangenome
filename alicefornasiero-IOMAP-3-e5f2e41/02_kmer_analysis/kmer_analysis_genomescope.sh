#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=batch
#SBATCH -J genomescope
#SBATCH -o /path/to/logs/genomescope.%J.out
#SBATCH -e /path/to/logs/genomescope.%J.err
#SBATCH --time=48:00:00
#SBATCH --mem=500G

## Data Variables
# outdir: path to output directory
# filelist: path to file with the list of the paths of your fastq files, one per line
# kmer: kmer length
# outprefix: prefix of the output file
# ploidy: ploidy level (diploid=2)

outdir=/path/to/kmer_analysis/genomescope
filelist=/path/to/kmer_analysis/read.list.CLR
kmer=19
nthreads=32
ploidy=4
outprefix=Oryza_species
##################################

mkdir -p ${outdir}
cd ${outdir}

#### 1. Run KMC
# -ci<value> - exclude k-mers occurring less than <value> times (default: 2)
# -cs<value> - maximal value of a counter (default: 255)
# -cx<value> - exclude k-mers occurring more of than <value> times (default: 1e9)
# -k<len> - k-mer length (k from 1 to 256; default: 25)
# -m<size> - max amount of RAM in GB (from 1 to 1024); default: 12

# Software Module
module load kmc/3.1.2

# run command
kmc -k${kmer} -v -t${nthreads} -m100 -cs10000 @${filelist} ${outprefix} ./

#### 2. Export the k-mer count histogram
kmc_tools transform ${outprefix} histogram ${outprefix}.histo

# Software Module unload
module unload kmc/3.1.2
####

#### 3a. Extract genomic kmers using reasonable coverage thresholds. 
# Inspect the kmer spectra and choose the L (lower) and U (upper) coverage thresholds
# You can estimate L and U using command smudgeplot.py cutoff

# Software Module
module load smudgeplot/0.2.4

# run command
L=$(smudgeplot.py cutoff ${outprefix}.histo L)
U=$(smudgeplot.py cutoff ${outprefix}.histo U)
echo $L $U

# Software Module unload
module unload smudgeplot/0.2.4
####

#### 3b. Alternative: you can run genomeScope to estimate L/U via visual inspection of the kmer coverage histogram

# Software Module
module load kmc/3.1.2

# run command
Rscript /ibex/sw/csi/kmc/3.1.2/gnu6.4.0/genomescope2.0/genomescope.R \
-i ${outdir}/${outprefix}.histo \
-n ${outprefix}_p${ploidy}_Gscope \
-o ${outdir} \
-k ${kmer} \
-p ${ploidy} &>${outdir}/${outprefix}_p${ploidy}_Gscope.log

# extract kmers in the coverage range from L (lower) to U (upper) using kmc_tools
kmc_tools transform ${outprefix} -ci"$L" -cx"$U" reduce ${outprefix}_L"$L"_U"$U"

# run kmc_dump on the reduced file to compute the set of kmer pairs.
kmc_tools transform ${outprefix} -ci"$L" -cx"$U" dump -s ${outprefix}_L"$L"_U"$U".dump

# Software Module unload
module unload kmc/3.1.2
####

#### 4. Run smudgeplot.py hetkmers to compute the set of kmer pairs.
# Software Module
module load smudgeplot/0.2.4

# run command
echo "smudgeplot.py hetkmers is running"
smudgeplot.py hetkmers -o ${outprefix}_L"$L"_U"$U" < ${outprefix}_L"$L"_U"$U".dump

# generate the smudgeplot using the coverages of the identified kmer pairs (*_coverages.tsv file).
smudgeplot_plot.R -i ${outprefix}_L"$L"_U"$U"_coverages.tsv \
-o ${outprefix}_smudge \
-k ${kmer} \
--title ${outprefix}

# Software Module unload
module unload smudgeplot/0.2.4
####