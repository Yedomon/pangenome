#!/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -l walltime=10:00:00
#PBS -d ./
#PBS -j oe

start=`date +%s`

CPU=$PBS_NP
if [ ! $CPU ]; then CPU=2; fi

N=$PBS_ARRAYID
if [ ! $N ]; then N=1; fi

#----------------
trimmomatic=/public1/home/liuyang/program/Trimmomatic-0.39/trimmomatic-0.39.jar
fq=`ls BC*1.fq.gz | head -n $N | tail -n 1`
prefix=${fq%%_*}

echo "CPU= $CPU"
echo "FILE= $ctg"
echo "FILE_PREFIX= $prefix"

   echo "Trimmomatic"
   java -jar $trimmomatic PE -phred33 ${prefix}_1.fq.gz ${prefix}_2.fq.gz \
      out.${prefix}_1.fq.gz out.trim.${prefix}_1.fq.gz out.${prefix}_2.fq.gz out.trim.${prefix}_2.fq.gz \
      ILLUMINACLIP:/public1/home/liuyang/program/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
      SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50

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

