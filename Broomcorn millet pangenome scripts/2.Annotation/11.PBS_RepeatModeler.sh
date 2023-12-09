#!/bin/bash
#PBS -l nodes=1:ppn=32
#PBS -l walltime=100:00:00
##PBS -l mem=100gb
#PBS -d ./
#PBS -j oe

#qsub -q fat -t 1 RepeatModeler.sh
#ppn=32	860Mb 26h

CPU=$PBS_NP
if [ ! $CPU ]; then CPU=2; fi

N=$PBS_ARRAYID
if [ ! $N ]; then N=1; fi

#----------------
genome=`ls *.fa | head -n $N | tail -n 1`
prefix=${genome%%.*}
dbName=${prefix}_db
pa=$(($CPU/4))

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

#----------------
if [ ! -e ${dbName}* ]; then
   echo "<<<<<<<<<< Start running BuildDatabase ..."
   BuildDatabase -name $dbName $genome
fi
if [ ! -e *-families.fa ]; then
   echo "<<<<<<<<<< Start running RepeatModeler ..."
   RepeatModeler -database $dbName -pa $pa -LTRStruct
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
echo "<<<<< All Done !!!"

