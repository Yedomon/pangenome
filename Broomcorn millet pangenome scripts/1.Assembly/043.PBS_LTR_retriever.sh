#!/bin/bash
#PBS -l nodes=1:ppn=10
##PBS -l mem=100gb
#PBS -l walltime=100:00:00
#PBS -d ./
#PBS -j oe

start=`date +%s`

CPU=$PBS_NP
if [ ! $CPU ]; then CPU=2; fi

N=$PBS_ARRAYID
if [ ! $N ]; then N=1; fi

#----------------
fa=`ls BC???.genome.fasta | head -n $N | tail -n 1`
prefix=${fa%%.*}

echo "CPU= $CPU"
echo "FILE= $fa"
echo "FILE_PREFIX= $prefix"

#if [ ! -e ${prefix} ]; then
   echo "running LTRharvest..."
   source activate
   conda deactivate   
   conda activate RM

   LTR_retriever -genome $fa -inharvest ../3.combineRawLTR/${fa}.rawLTR.scn -threads $CPU
#fi

#----------------
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

