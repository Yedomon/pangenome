#### Generate the Oryza genome type-level pangenome ####
#    Andrea Zuccolo
#    last modified: February 2025
########################################################

# PGGB
# For each genome type (AA, BB, CC, DD), we generated 12 separeted multifasta files containing the sequences for each individual chromosome.

# For each genome, pggb was then run for each chromosome as follows:
# For AA genome types (where ## is the chromosome number):
pggb -i CHR##.fasta -p 90 -s 15000 -n 13 --n-mappings 1 -k 7 -o out_CHR##_genome_type 
# For BB, CC and DD genome types (where ## is the chromosome number):
pggb -i CHR##.fasta -p 80 -s 15000 -n 13 --n-mappings 1 -k 7 -o out_CHR##_genome_type 

# PANACUS
# The CHR##.smooth.final.gfa output files of each pggb run were processed using panacus as follows:
panacus histgrowth -t 40 -c bp -l 1,1,1 -q 0,0.1,1 -S CHR##.smooth.final.gfa &gt; output.tsv

# The output.tsv file was then used to visualize the pangenome result with the following command:
panacus-visualize -e output.tsv &gt; output.pdf
