#!/bin/bash
#### Run Mecat2 assembler on PacBio CLR reads ####
#    Alice Fornasiero
#    last modified: February 2025
##################################################

## Software Modules
module load mecat2/20190314

## Data Variables 
export OUTDIR=/ibex/scratch/projects/c2079/analysis/O_minuta/mecat2
export SPECIES=O_minuta
mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}/logs

#### Step 1: Prepare config file
MEM="1gb"
CORES=1
JOB1_NAME="mecat_config"
JOB1_TYPE="sbatch --partition=batch --job-name=${JOB1_NAME}.${SPECIES} --time=00:10:00 --output=${OUTDIR}/logs/${JOB1_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB1_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB1_CMD="mecat.pl config ${OUTDIR}/${SPECIES}_config_file.txt";
${JOB1_TYPE} --parsable --wrap="${JOB1_CMD}"

#### Step 2: Correct raw reads
MEM="450gb"
CORES=32
JOB2_NAME="mecat_correct"
JOB2_TYPE="sbatch --partition=batch --job-name=${JOB2_NAME}.${SPECIES} --time=7-00:00:00 --output=${OUTDIR}/logs/${JOB2_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB2_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB2_CMD="mecat.pl correct ${OUTDIR}/${SPECIES}_config_file.txt";
${JOB2_TYPE} --parsable --wrap="${JOB2_CMD}"

# Step 3: Assemble contigs using the corrected reads
MEM="450gb"
CORES=32
JOB3_NAME="mecat_assemble"
JOB3_TYPE="sbatch --partition=batch --job-name=${JOB3_NAME}.${SPECIES} --time=7-00:00:00 --output=${OUTDIR}/logs/${JOB3_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB3_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB3_CMD="mecat.pl assemble ${OUTDIR}/${SPECIES}_config_file.txt";
${JOB3_TYPE} --parsable --wrap="${JOB3_CMD}"
