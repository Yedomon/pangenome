#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l walltime=20:00:00
#PBS -d ./
#PBS -j oe

# cpu=4, 30min

CPU=$PBS_NP
if [ ! $CPU ]; then CPU=2; fi

N=$PBS_ARRAYID
if [ ! $N ]; then N=1; fi

#----------------
minimap2=/public1/home/chenjf/Software/bin/minimap2
reference=Pmlongmi4.chr.fa
ctg=`ls *asm.p_ctg.fa | head -n $N | tail -n 1`
prefix=${ctg%%.*}

start=`date +%s`
start_date=`date -d @"$start" "+%Y-%m-%d %H:%M:%S"`
echo "<<<<< Start= $start_date"
echo "<<<<< CPU= $CPU"
echo "<<<<< FILE= $ctg"
echo "<<<<< FILE_PREFIX= $prefix"

#if [ ! -e ${prefix} ]; then
	source activate
	conda deactivate   
	conda activate RagTag 
	echo "Running ragtag.py ..."
	echo "Command: \
	ragtag.py scaffold --aligner $minimap2 --mm2-params '-x asm5' -t $CPU -u -o ${prefix}.scaf.out $reference $ctg "
	#
	ragtag.py scaffold --aligner $minimap2 --mm2-params '-x asm5' -t $CPU -u -o ${prefix}.scaf.out $reference $ctg
#fi

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

