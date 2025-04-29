#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=canu_assembly
#SBATCH --output=/path/to/logs/canu_assembly.%J.out
#SBATCH --error=/path/to/logs/canu_assembly.%J.err
#SBATCH --time=10-00:00:00

#### Run Canu assembler on PacBio CLR reads ####
#    Alice Fornasiero
#    last modified: February 2025
################################################

# Software Modules
module load canu/2.0/gnu8.2.0;

# Data Variables
INDIR=/path/to/Oryza_species
SPECIES=Oryza_species
OUTDIR=${INDIR}/canu
mkdir -p ${OUTDIR}/
mkdir -p ${OUTDIR}/logs

#### Run Canu assembly
time canu -d ${OUTDIR} \
-p ${SPECIES} \
genomeSize=1008m \
-pacbio-raw ${INDIR}/${SPECIES}.fasta \
gnuplot=undef \
usegrid=1 \
gridOptions="--time=10-00:00:00 --partition=batch --mem-per-cpu=16g" \
gridOptionsJobName=Om-using-grid
