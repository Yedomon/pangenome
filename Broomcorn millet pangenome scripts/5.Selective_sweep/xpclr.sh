##/usr/bin/bash
#!/bin/bash
#PBS -N C_W_XPCLR
#PBS -l nodes=1:ppn=1
#PBS -l walltime=9999:00:00
#PBS -d ./
#PBS -j oe
#PBS -V

cd $PBS_O_WORKDIR

start=`date +%s`

CPU=$PBS_NP
if [ ! $CPU ]; then
   CPU=2
fi

N=$PBS_ARRAYID
if [ ! $N ]; then
    N=1
fi

echo "$CPU"
echo "$N"

VCF=/public1/home/guowl/project/XPCLR_calculation/02_input_files/BC_516_merged.SNP.vcf.gz

pop1=W
pop2=C

size=20000
step=2000

xpclr --out ../04_output_files/${pop2}_vs_${pop1}/${pop2}_vs_${pop1}.chr${chr}.xpclr.txt --input $VCF \
      --samplesA ../../PI_tajimaD_Fst_calculation/03_intermediate_files/${pop1}.samples.txt \
      --samplesB ../../PI_tajimaD_Fst_calculation/03_intermediate_files/${pop2}.samples.txt \
      --chr ${chr} --maxsnps 300 --size ${size} --step ${step}


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"