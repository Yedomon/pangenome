# conda activate busco4
# conda deactivate

prefix=$1
nohup busco -m genome -l embryophyta_odb10 -i ${prefix}.ragtag.fa -o ${prefix}.OUT -c 6 --augustus_species maize &

