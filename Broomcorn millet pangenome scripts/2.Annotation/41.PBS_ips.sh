#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=100:00:00
#PBS -d ./
#PBS -j oe

CPU=$PBS_NP
if [ ! $CPU ]; then CPU=2; fi

N=$PBS_ARRAYID
if [ ! $N ]; then N=1; fi

##----------------
bindir=/public1/home/liuyang/program/interproscan-5.52-86.0
fa=`ls *.fasta | head -n $N | tail -n 1`

start=`date +%s`
start_date=`date -d @"$start" "+%Y-%m-%d %H:%M:%S"`
echo "<<<<< Start= $start_date"
echo "<<<<< CPU= $CPU"
echo "<<<<< FILE= $fa"

##----------------
echo "<<<<< Loading conda env ..."
source activate
conda deactivate
conda activate python3

## running
   echo "<<<<<<<<<< Running interproscan ..."
   $bindir/interproscan.sh -i $fa -f tsv -cpu $CPU -goterms -dp

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

