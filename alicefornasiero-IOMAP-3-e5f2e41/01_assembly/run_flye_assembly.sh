#!/bin/bash
#### Run Mecat2 assembler on PacBio CLR reads ####
#    Alice Fornasiero
#    last modified: February 2025
##################################################

## Software Modules
module load flye/2.8.1

## Data Variables 
export INDIR=/ibex/scratch/projects/c2079/analysis/O_minuta
export SPECIES=O_minuta
export OUTDIR=${INDIR}/flye
mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}/logs

#### Run Flye on PacBio CLR reads
MEM="450gb"
CORES=32
JOB1_NAME="flye_assemble"
JOB1_TYPE="sbatch --partition=batch --job-name=${JOB1_NAME}.${SPECIES} --time=00:10:00 --output=${OUTDIR}/logs/${JOB1_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB1_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB1_CMD="flye --pacbio-raw ${INDIR}/${SPECIES}.fasta --out-dir ${OUTDIR} --genome-size 1008000000 --threads ${CORES} --asm-coverage 35";
${JOB1_TYPE} --parsable --wrap="${JOB1_CMD}"
