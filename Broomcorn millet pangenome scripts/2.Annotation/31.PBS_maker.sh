#!/bin/bash
#PBS -N BC
#PBS -l nodes=1:ppn=4
#PBS -l walltime=100:00:00
#PBS -d ./
#PBS -j oe

CPU=$PBS_NP
if [ ! $CPU ]; then CPU=2; fi

N=$PBS_ARRAYID
if [ ! $N ]; then N=1; fi

##----------------
genome=`ls *.fa.masked | head -n $N | tail -n 1`
prefix=${genome%.fa.masked}

start=`date +%s`
start_date=`date -d @"$start" "+%Y-%m-%d %H:%M:%S"`
echo "<<<<< Start= $start_date"
echo "<<<<< CPU= $CPU"
echo "<<<<< FILE= $genome"
echo "<<<<< FILE_PREFIX= $prefix"

##----------------
echo "<<<<< Loading conda env= Maker2 ..."
source activate
conda deactivate
conda activate Maker2
export PATH=/public1/home/liuyang/program/maker-2.31.11/bin:$PATH

## running maker
   echo "<<<<<<<<<< Running maker ..."
   maker -genome $genome -base $prefix -cpus $CPU -nodatastore -quiet
   #mpiexec -n $CPU maker -genome $genome -base $prefix -nodatastore -quiet

## Making gff3 files
   echo "<<<<<<<<<< Making gff3 files ..."
   gff3_merge -d ${prefix}.maker.output/${prefix}_master_datastore_index.log -n
   fasta_merge -d ${prefix}.maker.output/${prefix}_master_datastore_index.log

##----------------
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

