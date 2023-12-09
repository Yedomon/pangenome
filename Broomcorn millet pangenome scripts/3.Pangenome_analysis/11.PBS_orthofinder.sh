#!/bin/bash
#PBS -l nodes=1:ppn=80
#PBS -l walltime=100:00:00
#PBS -d ./
#PBS -j oe

CPU=$PBS_NP
if [ ! $CPU ]; then CPU=2; fi

N=$PBS_ARRAYID
if [ ! $N ]; then N=1; fi

##----------------change it, dir contain all pep.fa
dir=/public1/home/liuyang/Project/81.orthofinder/BCnr100/

start=`date +%s`
start_date=`date -d @"$start" "+%Y-%m-%d %H:%M:%S"`
echo "<<<<< Start= $start_date"
echo "<<<<< CPU= $CPU"

##----------------
echo "<<<<< Loading conda env= OrthoFinder ..."
source activate
conda deactivate
conda activate OrthoFinder

## running
   echo "<<<<<<<<<< Running ..."
   orthofinder -f $dir -t $CPU -X -og

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

