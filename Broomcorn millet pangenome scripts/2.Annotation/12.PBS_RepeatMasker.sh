#!/bin/bash
#PBS -l nodes=1:ppn=32
##PBS -l mem=100gb
#PBS -l walltime=20:00:00
#PBS -d ./
#PBS -j oe

##qsub -q fat -t 1 RepeatMasker.sh
##ppn=32 ~2h

CPU=$PBS_NP
if [ ! $CPU ]; then CPU=2; fi

N=$PBS_ARRAYID
if [ ! $N ]; then N=1; fi

#----------------
genome=`ls *.ragtag.fa | head -n $N | tail -n 1`
prefix=${genome%%.*}
dbName=${prefix}_db-families.fa

start=`date +%s`
start_date=`date -d @"$start" "+%Y-%m-%d %H:%M:%S"`
echo "<<<<< Start= $start_date"
echo "<<<<< CPU= $CPU"
echo "<<<<< FILE= $genome"
echo "<<<<< FILE_PREFIX= $prefix"

#----------------
echo "<<<<< Loading conda env= RM ..."
source activate
conda deactivate
conda activate RM

# running RepeatMasker
   echo "<<<<<<<<<< Running RepeatMasker ..."
   RepeatMasker -e rmblast -pa $CPU -div 40 -nolow -norna -no_is -gff -lib $dbName -dir ${prefix}_repeatmasker.out $genome

# summary the TE content
   echo "<<<<<<<<<< Summarying TE content ..."
   samtools faidx $genome
   cat ${genome}.fai | cut -f1,2 > ${genome}.tsv
   perl /public1/home/liuyang/bin/buildSummary.pl -maxDiv 40 -genome ${genome}.tsv ${prefix}_repeatmasker.out/${genome}.out > ${genome}.sum
   rm ${genome}.fai ${genome}.tsv

#----------------
end=`date +%s`
end_date=`date -d @"$end" "+%Y-%m-%d %H:%M:%S"`
runtime=$((end-start))
h=$(($runtime/3600))
hh=$(($runtime%3600))
m=$(($hh/60))
s=$(($hh%60))
echo "<<<<< CPU= $CPU"
echo "<<<<< Start= $start_date"
echo "<<<<< End= $end_date"
echo "<<<<< Run time= $h:$m:$s"
echo "<<<<< All Done !!!"

