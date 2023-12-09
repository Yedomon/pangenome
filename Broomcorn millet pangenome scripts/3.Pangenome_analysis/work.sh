### get nr
for i in `ls *.fa`; do prefix=${i%.nr.fa}; cd-hit -i $i -o $prefix.nr100.fa -c 1 -n 5 -aS 1 -d 0 -M 0 -T 8; done

### clustering
qsub 11.PBS_orthofinder.sh

### calculating Ka, Ks
# get all BC pep
cat ~/Project/81.orthofinder/BCnr100/BC???.genome.maker.pep.nr100.fa |cut -d"." -f1 > ALL.BC.pep.nr100.fa
# get all BC cds
cat ~/Project/81.orthofinder/BCnr100/BC???.genome.maker.cds.nr100.fa |cut -d"." -f1 > ALL.BC.cds.nr100.fa

# get pep of core, dispensable, softcore, private genes
for i in {core,dispensable,softcore,private}; do file=BCnr100.pan.${i}.lst; bioawk -c fastx 'BEGIN{while((getline k <"'${file}'")>0) i[k]=1}{if(i[$name]) print ">"$name"\n"$seq}' ALL.BC.pep.nr100.fa > BCnr100.pan.${i}.pep.fa; done
# get cds
for i in {core,dispensable,softcore,private}; do file=BCnr100.pan.${i}.lst; bioawk -c fastx 'BEGIN{while((getline k <"'${file}'")>0) i[k]=1}{if(i[$name]) print ">"$name"\n"$seq}' ALL.BC.cds.nr100.fa > BCnr100.pan.${i}.cds.fa; done

# format Sit's pep and cds
cat Sitalica_312_v2.2.protein.fa |cut -d" " -f1 |sed -e 's/\.p//' > Sitalica_312_v2.2.protein.fasta
cat Sitalica_312_v2.2.cds.fa |cut -d" " -f1 > Sitalica_312_v2.2.cds.fasta

# get all seqs into one file
cat BCnr100.pan.core.pep.fa BCnr100.pan.dispensable.pep.fa BCnr100.pan.private.pep.fa BCnr100.pan.softcore.pep.fa Sitalica_312_v2.2.protein.fasta > ALL.pep.fa
cat BCnr100.pan.core.cds.fa BCnr100.pan.dispensable.cds.fa BCnr100.pan.private.cds.fa BCnr100.pan.softcore.cds.fa Sitalica_312_v2.2.cds.fasta > ALL.cds.fa

# blastp
makeblastdb -in Sitalica_312_v2.2.protein.fasta -dbtype prot
for i in {core,dispensable,softcore,private}; do blastp -query ../BCnr100.pan.${i}.pep.fa -db Sitalica_312_v2.2.protein.fasta -evalue 1e-5 -max_target_seqs 1 -num_threads 8 -out BCnr100.pan.${i}-Sit.blastp -outfmt 6; done

# get homo pair from blast output
for i in {core,dispensable,softcore,private}; do cat blast/BCnr100.pan.${i}-Sit.blastp |cut -f1-2 |sort |uniq > BCnr100.pan.${i}-Sit.homo; done

# run ParaAT to calculate Ka, Ks
qsub -q fat -t 1-4 21.PBS_paraAT.sh

# summarize the results
for i in {core,dispensable,softcore,private}; do find ./paraAT_result_BCnr100.pan.${i}-Sit/ -name *.axt.kaks |xargs cat |awk -v var="$i" '{if ($0 !~ /Sequence/) print var"\t"$1"\t"$5}' |sed 's/-0/0/' > Ka_Ks.$i; done

### calculating pi
# get all cds
cat ~/Project/00.gene_anno/BC???.genome.maker.cds.fa |cut -d" " -f1 > ALL.BC.cds.fa

# prepare seq
perl prepare_for_align_family.pl -fa ALL.BC.cds.fa -group BCnr100_Orthogroups.num.txt

#
mkdir pan_core pan_dispensable pan_softcore
mv pan_core.O* pan_core
mv pan_softcore.O* pan_softcore
mv pan_dispensable.O* pan_dispensable

# multi align
mkdir pan_core_aligned
for i in `ls pan_core`; do mafft --thread 36 --auto --quiet pan_core/$i > pan_core_aligned/${i}.aligned; done
mkdir pan_softcore_aligned
for i in `ls pan_softcore`; do mafft --thread 36 --auto --quiet pan_softcore/$i > pan_softcore_aligned/${i}.aligned; done
mkdir pan_dispensable_aligned
for i in `ls pan_dispensable`; do mafft --thread 36 --auto --quiet pan_dispensable/$i > pan_dispensable_aligned/${i}.aligned; done

# cal pi
for i in `ls pan_core_aligned`; echo $i; do perl compute_pi_use_aln_fasta_v1.1.pl -i pan_core_aligned/$i |grep "#Nucleotide diversity" |sed 's/#Nucleotide diversity, Pi: //' >> pi.core; done

# summarize
rm pi.all
echo -e "Class\tPi" >> pi.all
for i in `cat pi.pan_core`; do echo -n -e "core\t" >> pi.all; echo -e "$i" >> pi.all; done
for i in `cat pi.pan_softcore`; do echo -n -e "softcore\t" >> pi.all; echo -e "$i" >> pi.all; done
for i in `cat pi.pan_dispensable`; do echo -n -e "dispensable\t" >> pi.all; echo -e "$i" >> pi.all; done
