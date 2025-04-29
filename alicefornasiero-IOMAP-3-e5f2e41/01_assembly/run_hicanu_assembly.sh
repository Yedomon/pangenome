#!/bin/bash
#SBATCH --job-name=HiCanu_assembly
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --error=/path/to/logs/HiCanu_assembly.err
#SBATCH --output=/path/to/logs/HiCanu_assembly.out 
#SBATCH --time=5-00:00:00
#SBATCH --mem=450gb

#### Run HiCanu assembler on PacBio HiFi reads ####
#    Alice Fornasiero
#    last modified: February 2025
################################################

## Software Modules
module load hicanu/2.2

## Data Variables
CORES=32
WORKDIR=/path/to/Oryza_species
OUTDIR=/path/to/output_dir
PREFIX=Oryza_species
mkdir -p ${OUTDIR}

## Run HiCanu assembly
canu -p ${OUTDIR} \
-d ${PREFIX} \
genomeSize=800m \
-pacbio-hifi ${WORKDIR}/m64068_230209_044215.hifi_reads.fastq \
${WORKDIR}/m64068_230210_134855.hifi_reads.fastq \
${WORKDIR}/m64068_230211_225339.hifi_reads.fastq \
${WORKDIR}/m64068_230213_080152.hifi_reads.fastq \
-useGrid=false
