#!/bin/bash
#SBATCH --job-name=tophat2
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=40
#SBATCH --error=alignment.err
#SBATCH --output=alignment.out
#SBATCH --mem=500GB


module load fastqc/0.12.0
#leaf
fastqc -o S3_11r  SRR16526864_1.fastq SRR16526864_2.fastq
fastqc -o S3_11r  SRR16526865_1.fastq SRR16526865_2.fastq
fastqc -o S3_11r  SRR16526866_1.fastq SRR16526866_2.fastq
fastqc -o S3_11r  SRR4017523_1.fastq SRR4017523_2.fastq
fastqc -o S3_11r  SRR4017527_1.fastq SRR4017527_2.fastq

module load tophat2/2.1.1

bowtie2-build IRGSP-1.0_genome.fa IRGSP-1.0_genome

tophat2 -p 40 -G transcripts_exon.gff -o leaf_Os IRGSP-1.0_genome SRR16526864_1.fastq,SRR16526865_1.fastq,SRR16526866_1.fastq,SRR4017523_1.fastq,SRR4017527_1.fastq SRR16526864_2.fastq,SRR16526865_2.fastq,SRR16526866_2.fastq,SRR4017523_2.fastq,SRR4017527_2.fastq

mv leaf_Os/accepted_hits.bam leaf_Os.bam

mv leaf_Os/leaf_Os.bam .

#root

tophat2 -p 40 -G transcripts_exon.gff -o root_Os IRGSP-1.0_genome SRR25078452_1.fastq,SRR25078455_1.fastq,SRR25078456_1.fastq SRR25078452_2.fastq,SRR25078455_2.fastq,SRR25078456_2.fastq

mv root_Os/accepted_hits.bam root_Os.bam

mv root_Os/root_Os.bam .


