##########################################################################################
##########################################################################################
####### Phylogenetic relationship analysis of the chloroplast across the Oryza genus #####
##########################################################################################

####### Supports the findings in the following manuscript ################################
####### *Oryza genome evolution through a tetraploid lens* ###############################
####### by: Fornasiero et al., Nature Genetics, 2025 #####################################

# 0. Set up a Conda environment and install necessary packages
# For example, run these commands once before proceeding
conda create -n iqtree python=3.9
conda activate iqtree
conda install -c bioconda mafft trimal iqtree

### Note that, input chloroplast genomes and output files for each of the following steps 
### are also provided together with this text file
### at https://github.com/nam-hoang/rice_cp_phylo_analysis

##########################################################################################
##########################################################################################
### 1. Rice chloroplast genome assembly, acquisition, and preparation
# The chloroplast genome sequences of 
# 10 Oryza species (e.g., O. malampuzhaensis, O. minuta, O. alta, O. grandiglumis, 
# O. latifolia, O. coarctata, O. schlechteri, O. ridleyi, O. longiglumis, and O. meyeriana) 
# were assembled from whole-genome PacBio sequencing data (by the above study), 
# plus 17 publicly available chloroplast genome sequences 
# from 16 diploid Oryza species and the outgroup Leersia japonica.
# See Supplementary Data for accession numbers.

# The Large Single Copy (LSC) regions of these 27 chloroplast genomes were concatenated
# into a single FASTA file, i.e. "New_27_cp_genomes_correct_A_genomes_LSC.fa", and used 
# for the following analyses

##########################################################################################
##########################################################################################
### 2. Sequence alignment using MAFFT (v7.480)
# --auto: Automatically selects the best alignment strategy
# --thread 96: Uses 96 CPU threads (adjust based on your system resources)

mafft --auto --thread 96 New_27_cp_genomes_correct_A_genomes_LSC.fa \
> New_27_cp_genomes_correct_A_genomes_MAFFT_auto_LSC.fa

##########################################################################################
##########################################################################################
### 3. Removal of poorly aligned regions by TrimAL (v1.4)
# -automated1: Automatically detects and removes poorly aligned regions

trimal -in New_27_cp_genomes_correct_A_genomes_MAFFT_auto_LSC.fa \
-out New_27_cp_genomes_correct_A_genomes_MAFFT_auto_LSC_trimal_auto.fa -automated1

##########################################################################################
##########################################################################################
### 4. Constructing the maximum likelihood (ML) phylogenetic tree by IQ-TREE (v1.6.12)
# -st DNA: Specifies input as DNA sequences
# -m TEST: Automatically selects the best-fit model using ModelFinder
# -mem 120G: Limits memory usage to 120GB (adjust as needed)
# -ntmax 48: Uses up to 48 CPU threads (adjust as needed)
# -bb 1000 -alrt 1000: Bootstrapping and SH-aLRT tests with 1000 replicates

iqtree -s New_27_cp_genomes_correct_A_genomes_MAFFT_auto_LSC_trimal_auto.fa \
-st DNA -m TEST -mem 120G -ntmax 48 -bb 1000 -alrt 1000
