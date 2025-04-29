#!/bin/bash
#SBATCH --job-name=Hifiasm_assembly
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --error=/path/to/logs/Hifiasm_assembly.err
#SBATCH --output=/path/to/logs/Hifiasm_assembly.out 
#SBATCH --time=5-00:00:00
#SBATCH --mem=450gb

#### Run Hifiasm assembler on PacBio HiFi reads ####
#    Alice Fornasiero
#    last modified: February 2025
################################################

## Software Modules
module load hifiasm/0.19.5

## Data Variables
CORES=32
WORKDIR=/path/to/Oryza_species
OUTDIR=/path/to/output_dir
PREFIX=Oryza_species.asm
mkdir -p ${OUTDIR}

## Run HiCanu assembly
hifiasm -o ${OUTDIR}/${PREFIX} -t ${CORES} \
--primary ${WORKDIR}/m64041_211103_120859.hifi_reads.fasta.gz \
${WORKDIR}/m64313e_211202_171836.hifi_reads.fasta.gz
