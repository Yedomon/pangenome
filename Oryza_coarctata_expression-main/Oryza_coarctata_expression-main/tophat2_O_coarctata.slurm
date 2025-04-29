#!/bin/bash
#SBATCH --job-name=tophat2
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=40
#SBATCH --error=alignment.err
#SBATCH --output=alignment.out
#SBATCH --mem=500GB

module load tophat2/2.1.1


bowtie2-build Ocoarctata_ll_genome.fa Ocoarctata_ll_genome
bowtie2-build Ocoarctata_kk_genome.fa Ocoarctata_kk_genome

tophat2 -p 40 -G Ocoa_kk_genome.gff -o leaf_control_data_mapped_kk_subgenome Ocoarctata_kk_genome S0_1l_1.fq.gz,S0_2l_1.fq.gz,S0_3l_1.fq.gz,S0_4l_1.fq.gz,S0_5l_1.fq.gz,S0_6l_1.fq.gz,S0_8l_1.fq.gz,S0_9l_1.fq.gz,S0_10l_1.fq.gz,S0_11l_1.fq.gz,S0_12l_1.fq.gz,S0_13l_1.fq.gz,S0_14l_1.fq.gz,S0_15l_1.fq.gz,S0_16l_1.fq.gz,S0_17l_1.fq.gz,S0_18l_1.fq.gz  S0_1l_2.fq.gz,S0_2l_2.fq.gz,S0_3l_2.fq.gz,S0_4l_2.fq.gz,S0_5l_2.fq.gz,S0_6l_2.fq.gz,S0_8l_2.fq.gz,S0_9l_2.fq.gz,S0_10l_2.fq.gz,S0_11l_2.fq.gz,S0_12l_2.fq.gz,S0_13l_2.fq.gz,S0_14l_2.fq.gz,S0_15l_2.fq.gz,S0_16l_2.fq.gz,S0_17l_2.fq.gz,S0_18l_2.fq.gz

cd leaf_control_data_mapped_kk_subgenome

mv acceptem_hits.bam leaf_control_data_mapped_kk_subgenome.bam

cd ..

tophat2 -p 40 -G Ocoa_ll_genome.gff -o leaf_control_data_mapped_ll_subgenome Ocoarctata_ll_genome S0_1l_1.fq.gz,S0_2l_1.fq.gz,S0_3l_1.fq.gz,S0_4l_1.fq.gz,S0_5l_1.fq.gz,S0_6l_1.fq.gz,S0_8l_1.fq.gz,S0_9l_1.fq.gz,S0_10l_1.fq.gz,S0_11l_1.fq.gz,S0_12l_1.fq.gz,S0_13l_1.fq.gz,S0_14l_1.fq.gz,S0_15l_1.fq.gz,S0_16l_1.fq.gz,S0_17l_1.fq.gz,S0_18l_1.fq.gz  S0_1l_2.fq.gz,S0_2l_2.fq.gz,S0_3l_2.fq.gz,S0_4l_2.fq.gz,S0_5l_2.fq.gz,S0_6l_2.fq.gz,S0_8l_2.fq.gz,S0_9l_2.fq.gz,S0_10l_2.fq.gz,S0_11l_2.fq.gz,S0_12l_2.fq.gz,S0_13l_2.fq.gz,S0_14l_2.fq.gz,S0_15l_2.fq.gz,S0_16l_2.fq.gz,S0_17l_2.fq.gz,S0_18l_2.fq.gz

cd leaf_control_data_mapped_ll_subgenome

mv acceptem_hits.bam leaf_control_data_mapped_ll_subgenome.bam

cd ..

tophat2 -p 40 -G Ocoa_kk_genome.gff -o root_control_data_mapped_kk_subgenome Ocoarctata_kk_genome S0_1r_1.fq.gz,S0_2r_1.fq.gz,S0_3r_1.fq.gz,S0_5r_1.fq.gz,S0_6r_1.fq.gz,S0_7r_1.fq.gz,S0_9r_1.fq.gz,S0_10r_1.fq.gz,S0_11r_1.fq.gz,S0_13r_1.fq.gz,S0_14r_1.fq.gz,S0_17r_1.fq.gz S0_1r_2.fq.gz,S0_2r_2.fq.gz,S0_3r_2.fq.gz,S0_5r_2.fq.gz,S0_6r_2.fq.gz,S0_7r_2.fq.gz,S0_9r_2.fq.gz,S0_10r_2.fq.gz,S0_11r_2.fq.gz,S0_13r_2.fq.gz,S0_14r_2.fq.gz,S0_17r_2.fq.gz

cd root_control_data_mapped_kk_subgenome

mv acceptem_hits.bam root_control_data_mapped_kk_subgenome.bam

cd ..

tophat2 -p 40 -G Ocoa_ll_genome.gff -o root_control_data_mapped_ll_subgenome Ocoarctata_ll_genome S0_1r_1.fq.gz,S0_2r_1.fq.gz,S0_3r_1.fq.gz,S0_5r_1.fq.gz,S0_6r_1.fq.gz,S0_7r_1.fq.gz,S0_9r_1.fq.gz,S0_10r_1.fq.gz,S0_11r_1.fq.gz,S0_13r_1.fq.gz,S0_14r_1.fq.gz,S0_17r_1.fq.gz S0_1r_2.fq.gz,S0_2r_2.fq.gz,S0_3r_2.fq.gz,S0_5r_2.fq.gz,S0_6r_2.fq.gz,S0_7r_2.fq.gz,S0_9r_2.fq.gz,S0_10r_2.fq.gz,S0_11r_2.fq.gz,S0_13r_2.fq.gz,S0_14r_2.fq.gz,S0_17r_2.fq.gz

cd root_control_data_mapped_ll_subgenome

mv acceptem_hits.bam root_control_data_mapped_ll_subgenome.bam

cd ..

