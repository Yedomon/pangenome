#!/bin/bash
#### Run Mecat2 assembler on PacBio CLR reads ####
#    Alice Fornasiero
#    last modified: February 2025
##################################################

## Software Modules
module load smrtlink/8.0
module load samtools/1.8
module load bwa/0.7.17/gnu-6.4.0
module load minimap2/2.15
module load pilon/1.23

## Data Variables 
export RAWDATA=/path/to/raw_data/Oryza_species
export INPUT=/path/to/mecat2/assembly/4-fsa/contigs.fasta
export OUTDIR=/path/to/mecat2/polishing
export TMPDIR=${OUTDIR}/tmp
export OUTPREFIX=Oryza_species_mecat2
export PILON=/path/to/pilon
export PICARD=/path/to/picard
export R1_fq=/path/to/illumina/Oryza_species_R1.fastq
export R2_fq=/path/to/illumina/Oryza_species_R2.fastq
mkdir -p ${OUTDIR}
mkdir -p ${TMPDIR}

## Workflow steps 
for SAMPLE in `ls $RAWDATA/*.subreads.bam`;
do
    SUFFIX=`basename $SAMPLE .subreads.bam`
		
    ## Run pmmb2 to align pacbio reads on the contig assembly and sort them
    MEM="450gb"
    CORES=32
    JOB1_NAME="pbmm2_align"
    JOB1_TYPE="sbatch --partition=batch --job-name=${JOB1_NAME}.${SPECIES} --time=3-00:00:00 --output=${OUTDIR}/logs/${JOB1_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB1_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
    JOB1_CMD="pbmm2 align ${INPUT} ${SAMPLE} ${OUTDIR}/${OUTPREFIX}_${SUFFIX}.sort.bam --sort -j 24 -J 8 ${TMPDIR}";
    ${JOB1_TYPE} --parsable --wrap="${JOB1_CMD}"

done

## Merge the sorted .bam files
MEM="450gb"
CORES=1
JOB2_NAME="samtools merge"
JOB2_TYPE="sbatch --partition=batch --job-name=${JOB2_NAME}.${SPECIES} --time=2-00:00:00 --output=${OUTDIR}/logs/${JOB2_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB2_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB2_CMD="samtools merge ${OUTDIR}/${OUTPREFIX}.merged.bam ${OUTDIR}/*.sort.bam";
${JOB2_TYPE} --parsable --wrap="${JOB2_CMD}"

## Create index of merged bam file using pbindex
MEM="450gb"
CORES=1
JOB3_NAME="pbindex_bam"
JOB3_TYPE="sbatch --partition=batch --job-name=${JOB3_NAME}.${SPECIES} --time=20:00:00 --output=${OUTDIR}/logs/${JOB3_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB3_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB3_CMD="pbindex ${OUTDIR}/${OUTPREFIX}.merged.bam";
${JOB3_TYPE} --parsable --wrap="${JOB3_CMD}"

## Create index of merged bam file using samtools
MEM="450gb"
CORES=1
JOB4_NAME="samtools_index_bam"
JOB4_TYPE="sbatch --partition=batch --job-name=${JOB4_NAME}.${SPECIES} --time=20:00:00 --output=${OUTDIR}/logs/${JOB4_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB4_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB4_CMD="samtools index ${OUTDIR}/${OUTPREFIX}.merged.bam ${OUTDIR}/${OUTPREFIX}.merged.bam.bai";
${JOB4_TYPE} --parsable --wrap="${JOB4_CMD}"

## Create index of the contig assembly using samtools
MEM="450gb"
CORES=1
JOB5_NAME="samtools_index_fasta"
JOB5_TYPE="sbatch --partition=batch --job-name=${JOB5_NAME}.${SPECIES} --time=2:00:00 --output=${OUTDIR}/logs/${JOB5_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB5_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB5_CMD="samtools faidx ${INPUT}";
${JOB5_TYPE} --parsable --wrap="${JOB5_CMD}"

## Run Arrow software
MEM="450gb"
CORES=1
JOB6_NAME="run_arrow"
JOB6_TYPE="sbatch --partition=batch --job-name=${JOB6_NAME}.${SPECIES} --time=2-00:00:00 --output=${OUTDIR}/logs/${JOB6_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB6_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB6_CMD="gcpp --algorithm=arrow -j ${CORES} -r ${INPUT} -o ${OUTDIR}/${OUTPREFIX}.arrow.fasta ${OUTDIR}/${OUTPREFIX}.merged.bam";
${JOB6_TYPE} --parsable --wrap="${JOB6_CMD}"

## Create index of the arrowed polished bam files using bwa 
MEM="450gb"
CORES=1
JOB7_NAME="bwa_index"
JOB7_TYPE="sbatch --partition=batch --job-name=${JOB7_NAME}.${SPECIES} --time=20:00:00 --output=${OUTDIR}/logs/${JOB7_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB7_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB7_CMD="bwa index -p ${OUTDIR}/${OUTPREFIX}.arrow.bwa ${OUTDIR}/${OUTPREFIX}.arrow.fasta";
${JOB7_TYPE} --parsable --wrap="${JOB7_CMD}"

## Short reads mapping to the arrow-polished assembled reads using bwa
MEM="450gb"
CORES=32
JOB8_NAME="align_illumina"
JOB8_TYPE="sbatch --partition=batch --job-name=${JOB8_NAME}.${SPECIES} --time=2-00:00:00 --output=${OUTDIR}/logs/${JOB8_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB8_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB8_CMD="bwa mem -t 32 ${OUTDIR}/${OUTPREFIX}.arrow.bwa ${R1_fq} ${R2_fq} | samtools view -bhS - > ${OUTDIR}/${OUTPREFIX}.arrow.bam";
${JOB8_TYPE} --parsable --wrap="${JOB8_CMD}"

## Sort reads using samtools
MEM="450gb"
CORES=32
JOB9_NAME="align_illumina"
JOB9_TYPE="sbatch --partition=batch --job-name=${JOB9_NAME}.${SPECIES} --time=2-00:00:00 --output=${OUTDIR}/logs/${JOB9_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB9_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB9_CMD="samtools sort -@ 32 -o ${OUTDIR}/${OUTPREFIX}.arrow_sorted.bam ${OUTDIR}/${OUTPREFIX}.arrow.bam";
${JOB9_TYPE} --parsable --wrap="${JOB9_CMD}"

## Index sorted reads
MEM="450gb"
CORES=1
JOB10_NAME="pilon"
JOB10_TYPE="sbatch --partition=batch --job-name=${JOB10_NAME}.${SPECIES} --time=20:00:00 --output=${OUTDIR}/logs/${JOB10_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB10_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB10_CMD="samtools index ${OUTDIR}/${OUTPREFIX}.arrow_sorted.bam";
${JOB10_TYPE} --parsable --wrap="${JOB10_CMD}"

## Polish the genome assembly using Pilon
MEM="450gb"
CORES=32
JOB11_NAME="pilon"
JOB11_TYPE="sbatch --partition=batch --job-name=${JOB11_NAME}.${SPECIES} --time=2-00:00:00 --output=${OUTDIR}/logs/${JOB11_NAME}.${SPECIES}.%J.out --error=${OUTDIR}/logs/${JOB11_NAME}.${SPECIES}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM}";
JOB11_CMD="java -Xmx450G -jar $PILON --genome ${OUTDIR}/${OUTPREFIX}.arrow.fasta --frags ${OUTDIR}/${OUTPREFIX}.arrow_sorted.bam --output ${OUTPREFIX}.arrow.pilon --outdir ${OUTDIR} --threads 32 --changes";
${JOB11_TYPE} --parsable --wrap="${JOB11_CMD}"
