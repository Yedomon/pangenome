
## Identification of homoelog pairs of genes between KK and LL sub genom of O.coarctata

## Blastp between KK and LL subgenomes of Oryza coarctata and between Oryza coarctata and Oryza sativa proteins.  
blastp -db KK_proteins.fa -query LL_proteins.fa -out LL_hits_KK.txt -outfmt 6 -evalue 0.00001 -num_threads 30
blastp -db LL_proteins.fa -query KK_proteins.fa -out KK_hits_LL.txt -outfmt 6 -evalue 0.00001 -num_threads 30


# add extra column to blast output which help in the identification of reciprocal hits

awk '{print $1"_"$2 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" 12}' LL_hits_KK.txt > LL_hits_KK_with_id.txt

awk '{print $2"_"$1 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" 12}' KK_hits_LL.txt > KK_hits_LL_with_id.txt

#identification of reciprocal hits

perl extract_reciprocal_hits.pl LL_hits_KK_with_id.txt KK_hits_LL_with_id.txt > reciprocal_hits_KK_LL.txt

# filter hits with Bits score threshold of 50 and a minimum alignment identity of 50% for aligning regions

awk '$4>=50 && $13>=50 && $17>=50 && $26>=50 {print $2 "\t" $3}' reciprocal_hits_KK_LL.txt > reciprocal_hits_KK_LL_filtered_paired.txt

# create list of genes in pair from KK and LL sub_genomes

awk '{print $2 "\n" $3}' reciprocal_hits_KK_LL_filtered_paired.txt > list_of_genes_in_pair_KK_LL.txt

awk '!x[$0]++' list_of_genes_in_pair.txt > uniq_list_of_genes_in_pair_KK_LL.txt

sed 's/\g//g'  uniq_list_of_genes_in_pair.txt > temp.txt

awk '$1<=27445 {print g$1}'  temp.txt > KK_genes_in_pair.txt

awk '$1>27445 {print g$1}'  temp.txt > LL_genes_in_pair.txt

#genes fractioned in KK and LL

grep '>' KK_proteins.fa | sed 's/>//g' > KK_gene_list.txt

grep '>' LL_proteins.fa | sed 's/>//g' > LL_gene_list.txt

grep -fvxF KK_gene_list.txt KK_genes_in_pair.txt > KK_genes_not_in_pair.txt

grep -fvxF LL_gene_list.txt LL_genes_in_pair.txt > LL_genes_not_in_pair.txt

# create final list for fractioned and unfractioned genes

awk '{print $1 "\tKK(UF)"}' KK_genes_in_pair.txt > KK_genes_in_pair.list
awk '{print $1 "\tLL(UF)"}' LL_genes_in_pair.txt > LL_genes_in_pair.list
awk '{print $1 "\tKK(F)"}'  KK_genes_not_in_pair.txt > KK_genes_not_in_pair.list
awk '{print $1 "\tLL(F)"}'  LL_genes_not_in_pair.txt > LL_genes_not_in_pair.list





##Expression study of coarctata leaf and root

# FastQC for quality checking

perl create_fastqc_run.pl sample.list > run_fastqc.sh

bash run_fastqc.sh

# Align read data to the KK and LL subgenomes

sbatch tophat2.slurm

#run eagle-RC to remove ambigousely aligned reads from alignment


eagle-rc --paired --ngi --splice  --ref1=Ocoarctata_kk_genome.fa --ref2=Ocoarctata_ll_genome.fa --bam1=leaf_control_data_mapped_kk_subgenome.bam  --bam2=leaf_control_data_mapped_ll_subgenome.bam  -o leaf_control_data > leaf_control_data_classified_reads.list

eagle-rc --paired --ngi --splice  --ref1=Ocoarctata_kk_genome.fa --ref2=Ocoarctata_ll_genome.fa --bam1=root_control_data_mapped_kk_subgenome.bam  --bam2=root_control_data_mapped_ll_subgenome.bam  -o root_control_data > root_control_data_classified_reads.list


# TPM calculation

module load tpmcalculator/0.0.3

gffread Ocoa_kk_genome.gff -o Ocoa_kk_genome.gtf -T
gffread Ocoa_ll_genome.gff -o Ocoa_ll_genome.gtf -T

TPMCalculator -a -g Ocoa_kk_genome.gtf -b leaf_control_data1.ref.bam
TPMCalculator -a -g Ocoa_ll_genome.gtf -b root_control_data1.ref.bam
TPMCalculator -a -g Ocoa_kk_genome.gtf -b leaf_control_data2.ref.bam
TPMCalculator -a -g Ocoa_ll_genome.gtf -b root_control_data2.ref.bam


#######Between Oryza sativa and Oryza coarctata


blastp -db Oryza_coarctata_proteins.fa -query Oryza_sativa_proteins.fa -out sative_hits_coarctata_hits.txt -outfmt 6 -evalue 0.00001 -num_threads 30
blastp -db Oryza_sativa_proteins.fa -query Oryza_coarctata_proteins.fa -out coarctata_hits_sativa_hits.txt -outfmt 6 -evalue 0.00001 -num_threads 30


awk '{print $1"_"$2 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" 12}' sative_hits_coarctata_hits.txt > sative_hits_coarctata_hits_id.txt

awk '{print $2"_"$1 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" 12}' coarctata_hits_sativa_hits.txt > coarctata_hits_sativa_hits_id.txt



perl extract_reciprocal_hits.pl sative_hits_coarctata_hits_id.txt coarctata_hits_sativa_hits_id.txt > reciprocal_hits_Os_Oc.txt



awk '$4>=50 && $13>=50 && $17>=50 && $26>=50 {print $2 "\t" $3}' reciprocal_hits_Os_Oc.txt > reciprocal_hits_Os_Oc_filtered_paired.txt



awk '{print $2 "\n" $3}' reciprocal_hits_Os_Oc_filtered_paired.txt > list_of_genes_in_pair_Os_Oc.txt

awk '!x[$0]++' list_of_genes_in_pair_Os_Oc.txt > uniq_list_of_genes_in_pair_Os_Oc.txt

sed 's/g//g' uniq_list_of_genes_in_pair_Os_Oc.txt > temp1.txt

sed 's/\.t1//g' temp1.txt > uniq_list_of_genes_in_pair_Os_Oc.txt

sed 's/g//g' uniq_list_of_genes_in_pair_KK_LL.txt > temp1.txt

sed 's/g//g' temp1.txt > uniq_list_of_genes_in_pair_KK_LL.txt


perl extract_paired_homoeolog_Os_Oc.pl uniq_list_of_genes_in_pair_KK_LL.txt > paired_homologs_Oc_Os.txt

cut -f 3 paired_homologs_Oc_Os.txt > genes_os_homolog_pairs.txt

sed 's/,/\n/g' genes_os_homolog_pairs.txt > genes_os_homolog_pairs1.txt

awk '!x[$0]++' genes_os_homolog_pairs1.txt > genes_os_double_copy_homolog_pairs.txt

awk '{print $1 "\tOs(PH)"}' genes_os_double_copy_homolog_pairs.txt > genes_os_double_copy_homolog_pairs_final.txt

grep -vxFf genes_os_double_copy_homolog_pairs.txt all_Os_genes1.pair > genes_present_in_blast_hit_not_double_copy_homologs.txt

perl extract_single_copy_homolog_Os_Oc.pl genes_present_in_blast_hit_not_double_copy_homologs.txt > single_copy_homolog_Os_Oc.txt

# Manually curate the list of single_copy_homolog_Os_Oc.txt. Some of the Os genes showed hits with both KK and LL, but those KK-LL genes are not in homoeologous pairs. There were around 600 genes, so we removed them from the SCH list.

# further processing of single copy homolog pairs

sed 's/\t/\n/' single_copy_homolog_Os_Oc.txt > single_copy_homolog_Os_Oc1.txt

grep 'Os' single_copy_homolog_Os_Oc1.txt > single_copy_homolog_Os_Oc2.txt

awk '!x[$0]++' single_copy_homolog_Os_Oc2.txt > single_copy_homolog_Os_Oc.txt


awk '{print $1 "\tOs(SCH)"}' single_copy_homolog_Os_Oc.txt > single_copy_homolog_Os_Oc_final.txt


## Map RNAseq raw reads of leaf tissue from NCBI to Rice genome

sbatch tophat2_O_sativa.slurm

cd leaf_Os 
mv accepted_hits.bam leaf_Os.bam
##calculate TPM for leaf
module load tpmcalculator/0.0.3

TPMCalculator -a -g transcripts_exon1.gtf -b leaf_Os.bam

#output leaf_Os_genes.out

## Map RNAseq raw reads of root tissue from NCBI to Rice genome

TPMCalculator -a -g transcripts_exon1.gtf -b root_Os.bam


#Prepare the final table with TPM values

cat KK_genes_in_pair.list LL_genes_in_pair.list KK_genes_not_in_pair.list LL_genes_not_in_pair.list genes_os_double_copy_homolog_pairs_final.txt single_copy_homolog_Os_Oc_final.txt > genes_different_categories.txt 

cat Ocoa_kk_leaf_genes.out Ocoa_ll_leaf_genes.out leaf_Os_genes.out > leaf_tpm_all.txt 
 
cat Ocoa_kk_root_genes.out Ocoa_ll_root_genes.out root_Os_genes.out > root_tpm_all.txt

cut -f 1,7 leaf_tpm_all.txt > leaf_tpm_all_final.txt

cut -f 1,7 root_tpm_all.txt > root_tpm_all_final.txt

perl form_final_table.pl leaf_tpm_all_final.txt genes_different_categories.txt > tpm_table_leaf.txt

perl form_final_table.pl root_tpm_all_final.txt genes_different_categories.txt > tpm_table_root.txt


### boxplot

Rscript boxplot.R














 
