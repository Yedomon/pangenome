#!/bin/bash
#PBS -l nodes=comput19:ppn=10
##PBS -l mem=160gb
#PBS -l walltime=100:00:00
#PBS -d ./
#PBS -j oe

echo "NODE= comput19"

# Run -t 1 very first, for generating "*.giraffe.gbz" !!!

start=`date +%s`
start_date=`date -d @"$start" "+%Y-%m-%d %H:%M:%S"`

CPU=$PBS_NP
if [ ! $CPU ]; then CPU=2; fi

N=$PBS_ARRAYID
if [ ! $N ]; then N=1; fi

#----------------
fq=`ls ./00.BC_516_fq/BC???_1.fq.gz | head -n $N | tail -n 1`
#fq=`ls ./00.BC_32_fq/BC???_1.fq.gz | head -n $N | tail -n 1`
fq_dir=$(dirname $fq)
BCXXX=$(basename $fq _1.fq.gz)
prefix=$BCXXX

GENOME=Longmi4.genome.Chr.fa
VCF=ALL.PAV.vcf
VG_INDEX=ALL.PAV.a

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

### construct graph, use PBS_vg_pipe_index.sh run only once (20min)
#echo "<<<<<<<<<< vg construct ..."
#vg construct -a -r $GENOME -v $VCF -t $CPU > ${VG_INDEX}.vg
#vg construct -a -r $GENOME -v $VCF -t $CPU > $graph_vg
#echo "<<<<<<<<<< vg index ..."
#vg index ${VCF}.a.vg -x ${VCF}.a.xg -L
#vg index ${VG_INDEX}.vg -x ${VG_INDEX}.xg -L

if [ ! -e ./out.vcf/${prefix}.genotypes.vcf ] && true; then
  ### Mapping ---------- take hours
  echo "<<<<<<<<<< Mapping (vg giraffe) ..."
  if [ ! -e ./out.gam/${prefix}.aln.gam ] && true; then
     # map reads against the graph, and get BAM output in linear space
     # -Z, --gbz-name FILE : use this GBZ file (GBWT index + GBWTGraph)
     # -f, --fastq-in FILE : read and align FASTQ-format reads from FILE (two are allowed, one for each mate)
     # -o, --output-format NAME : output the alignments in NAME format (gam / gaf / json / tsv / SAM / BAM / CRAM) [gam]
     #vg giraffe -Z ${AUTO_INDEX}.giraffe.gbz -f $fq_dir/${BCXXX}_1.fq.gz -f $fq_dir/${BCXXX}_2.fq.gz -o gam -t $CPU > ./out.gam/${BCXXX}.aln.gam
     #also generate .giraffe.gbz .min .dist files
     vg giraffe -x ${VG_INDEX}.xg -f $fq_dir/${BCXXX}_1.fq.gz -f $fq_dir/${BCXXX}_2.fq.gz -o gam -t $CPU > ./out.gam/${BCXXX}.aln.gam
  fi

  ###
  echo "<<<<<<<<<< vg pack ..."
  echo "<<<<<<<<<< Compute the reads support (vg pack) -----" # 20min
  if [ ! -e ./out.pack/${prefix}.aln.pack ] && true; then
     # Call only variants that are present in the graph (use -g).
     # Compute the read support from the gam.
     # -Q 5: ignore mapping and base qualitiy < 5
     # -x, --xg FILE : use this basis graph (any format accepted, does not have to be xg)
     # -g, --gam FILE : read alignments from this GAM file (could be '-' for stdin)
     # -o, --packs-out FILE : write compressed coverage packs to this output file
     #vg pack -x ${VCF}.a.vg -g ./out.gam/${BCXXX}.aln.gam -o ./out.pack/${BCXXX}.aln.pack -Q 5 -t $CPU
     #vg pack -x ${AUTO_INDEX}.giraffe.gbz -g ./out.gam/${BCXXX}.aln.gam -o ./out.pack/${BCXXX}.aln.pack -Q 5 -t $CPU
     vg pack -x ${VG_INDEX}.xg -g ./out.gam/${BCXXX}.aln.gam -o ./out.pack/${BCXXX}.aln.pack -Q 5 -t $CPU
  fi

  ###
  echo "<<<<<<<<<< vg snarls ..."
  echo "<<<<<<<<<< Compute the snarls (vg snarls) -----" #
  if [ ! -e ${VG_INDEX}.xg.snarls ] && true; then
     # For larger graphs, it is recommended to compute snarls separately
     #vg snarls ${VCF}.a.vg -t $CPU > ${VCF}.a.snarls
     #vg snarls ${AUTO_INDEX}.giraffe.gbz -t $CPU > ${AUTO_INDEX}.giraffe.gbz.snarls
     vg snarls ${VG_INDEX}.xg -t $CPU > ${VG_INDEX}.xg.snarls
  fi

  ###
  echo "<<<<<<<<<< vg call ..."
  echo "<<<<<<<<<< Genotype the VCF (vg call) -----" # 40min
  if [ ! -e ./out.vcf/${prefix}.genotypes.vcf ] && true; then
     # -k, --pack FILE : Supports created from vg pack for given input graph
     # -v, --vcf FILE : VCF file to genotype (must have been used to construct input graph with -a)
     # -a, --genotype-snarls : Genotype every snarl, including reference calls (use to compare multiple samples), incompatible with -v
     # -s, --sample NAME : Sample name [default=SAMPLE]
     # usage: vg call [options] <graph> > output.vcf
     #vg call ${VCF}.a.vg -r ${VCF}.a.snarls -k ${BCXXX}.aln.pack -v $VCF -s $BCXXX -t $CPU > ${BCXXX}.genotypes.vcf
     #vg call ${AUTO_INDEX}.giraffe.gbz -r ${AUTO_INDEX}.giraffe.gbz.snarls -k ./out.pack/${BCXXX}.aln.pack -v $VCF -s $BCXXX -t $CPU > ./out.vcf/${BCXXX}.genotypes.vcf
     #worked for pri
     vg call ${VG_INDEX}.xg -r ${VG_INDEX}.xg.snarls -k ./out.pack/${BCXXX}.aln.pack -v $VCF -s $BCXXX -t $CPU > ./out.vcf/${BCXXX}.genotypes.vcf
  fi
  ###
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

