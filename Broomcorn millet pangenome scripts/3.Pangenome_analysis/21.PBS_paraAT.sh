#!/bin/bash
#PBS -l nodes=1:ppn=24
##PBS -l mem=100gb
#PBS -l walltime=100:00:00
#PBS -d ./
#PBS -j oe

start=`date +%s`

CPU=$PBS_NP
if [ ! $CPU ]; then CPU=2; fi

N=$PBS_ARRAYID
if [ ! $N ]; then N=1; fi

##----------------
file=`ls BCnr100.pan.*-Sit.homo | head -n $N | tail -n 1`
prefix=${file%.homo}

echo "CPU= $CPU"
echo "FILE= $file"
echo "FILE_PREFIX= $prefix"
  
  #echo $CPU > proc
  ~/program/ParaAT2.0/ParaAT.pl -h $file -n ALL.cds.fa -a ALL.pep.fa -p proc -m mafft -f axt -g -k -o paraAT_result_$prefix

##----------------
end=`date +%s`
runtime=$((end-start))
h=$(($runtime/3600))
hh=$(($runtime%3600))
m=$(($hh/60))
s=$(($hh%60))

echo "Start= $start"
echo "End= $end"
echo "Run time= $h:$m:$s"
echo "Done!"

