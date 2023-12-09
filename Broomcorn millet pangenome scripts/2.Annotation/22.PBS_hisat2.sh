#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -d ./
#PBS -j oe

CPU=$PBS_NP
if [ ! $CPU ]; then CPU=2; fi

N=$PBS_ARRAYID
if [ ! $N ]; then N=1; fi

#----------------
genome=`ls *.ragtag.fa | head -n $N | tail -n 1`
prefix=${genome%%.*}
index=${prefix}_index

start=`date +%s`
start_date=`date -d @"$start" "+%Y-%m-%d %H:%M:%S"`
echo "<<<<< Start= $start_date"
echo "<<<<< CPU= $CPU"
echo "<<<<< FILE= $genome"
echo "<<<<< FILE_PREFIX= $prefix"

#----------------
echo "<<<<< Loading conda env= RNAseq ..."
source activate
conda deactivate
conda activate RNAseq

#----------------
if [ ! -e ${index}* ]; then
   echo "<<<<<<<<<< Running hisat2-build ..."
   hisat2-build -p $CPU $genome $index
fi
if [ ! -e ${prefix}.sorted.bam ]; then
   echo "<<<<<<<<<< Running hisat2 ..."
   hisat2 --dta --threads $CPU --min-intronlen 20 --max-intronlen 15000 --fr -x $index \
      -1 ../23.Trimmomatic/out.${prefix}_1.fq.gz -2 ../23.Trimmomatic/out.${prefix}_2.fq.gz -S ${prefix}.sam
   
   samtools view -Su ${prefix}.sam | samtools sort -@ $CPU -o ${prefix}.sorted.bam
   rm ${prefix}.sam
fi

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
echo "<<<<< Done!!!"

