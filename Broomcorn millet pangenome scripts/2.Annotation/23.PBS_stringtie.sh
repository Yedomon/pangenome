#!/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -l walltime=10:00:00
#PBS -d ./
#PBS -j oe

CPU=$PBS_NP
if [ ! $CPU ]; then CPU=2; fi

N=$PBS_ARRAYID
if [ ! $N ]; then N=1; fi

#----------------
bam=`ls *.sorted.bam | head -n $N | tail -n 1`
prefix=${bam%%.*}

start=`date +%s`
start_date=`date -d @"$start" "+%Y-%m-%d %H:%M:%S"`
echo "<<<<< Start= $start_date"
echo "<<<<< CPU= $CPU"
echo "<<<<< FILE= $bam"
echo "<<<<< FILE_PREFIX= $prefix"

#----------------
echo "<<<<< Loading conda env= RNAseq ..."
source activate
conda deactivate
conda activate RNAseq

#----------------
if [ ! -e ${prefix}.stringtie.gtf ]; then
   echo "<<<<<<<<<< Running stringtie ..."
   stringtie -p $CPU -l ${prefix} -o ${prefix}.stringtie.gtf $bam

   ~/program/cufflinks-2.2.1/gffread -g ./genome/${prefix}.ragtag.fa -w ${prefix}.transcript.fa ${prefix}.stringtie.gtf
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

