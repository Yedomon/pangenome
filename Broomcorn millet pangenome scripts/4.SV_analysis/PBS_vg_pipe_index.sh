#!/bin/bash
#PBS -l nodes=comput17:ppn=12
##PBS -l mem=160gb
#PBS -l walltime=100:00:00
#PBS -d ./
#PBS -j oe

echo "NODE= comput17"

start=`date +%s`
start_date=`date -d @"$start" "+%Y-%m-%d %H:%M:%S"`

CPU=$PBS_NP
if [ ! $CPU ]; then CPU=2; fi

N=$PBS_ARRAYID
if [ ! $N ]; then N=1; fi

#----------------
GENOME=Longmi4.genome.Chr.fa

VCF=ALL.PAV.vcf
graph_vg=ALL.PAV.a.vg
graph_xg=ALL.PAV.a.xg

#----------------
echo "CPU= $CPU"
echo "FILE= $fq"
echo "FILE_PREFIX= $prefix"

### Loading environment ----------
#source activate
#conda deactivate   
#conda activate vg
#PATH=/public1/home/chenjf/Software/miniconda3/envs/vg_tools2/bin/:$PATH
PATH=/public1/home/liuyang/program/Miniconda2/envs/vg/bin/:$PATH
#PATH=/public1/home/chenjf/Share/Chenlab_public3/15.pangenome_vg/vg:$PATH

### construct graph
  ### generate *.vg, *.xg files
  echo "<<<<<<<<<< vg construct ..."
  vg construct -a -r $GENOME -v $VCF -t $CPU > $graph_vg
  
  echo "<<<<<<<<<< vg index ..."
  vg index $graph_vg -x $graph_xg -L

  ### generate snarls
  echo "<<<<<<<<<< vg snarls ..."
  if [ ! -e ${graph_xg}.snarls ] && true; then
     # For larger graphs, it is recommended to compute snarls separately
     #vg snarls ${VCF}.a.vg -t $CPU > ${VCF}.a.snarls
     #vg snarls ${AUTO_INDEX}.giraffe.gbz -t $CPU > ${AUTO_INDEX}.giraffe.gbz.snarls
     vg snarls $graph_xg -t $CPU > ${graph_xg}.snarls
  fi

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

