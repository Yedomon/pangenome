##/usr/bin/bash
#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l mem=20b
#PBS -l walltime=40:00:00
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

vcf=BC_516_merged.SNP.vcf.gz
vcftools=/public1/home/guowl/software/vcftools_0.1.13/bin/vcftools


pops=("C1" "C2" "C3" "W")
window_size=(10000 20000 50000 100000)
window_step=(1000 2000 5000 10000)


mkdir ../04_output_files/pi/
# calculate the pi with different levels of pop, window size
for i in ${pops[*]};do
    mkdir ../04_output_files/pi/${i}/
        for j in {0..3};do
            $vcftools --gzvcf ../02_input_files/$vcf --window-pi ${window_size[${j}]} --window-pi-step ${window_step[${j}]} \
            --keep ../03_intermediate_files/${i}.samples.txt --out ../04_output_files/pi/${i}/${i}.win_${window_size[${j}]}
        done;
done;  

# calculate the Fst between the groups for each chromosome
mkdir ../04_output_files/fst
for pair in C1,C2 C1,C3 C1,W C2,C3 C2,W C3,W;do
    str=(${pair//,/ })
    pop1=${str[0]}
    pop2=${str[1]}
    mkdir ../04_output_files/fst/${pop1}_${pop2}
    for i in {0..3};do
        $vcftools --gzvcf ../02_input_files/$vcf --weir-fst-pop ../03_intermediate_files/${pop1}.samples.txt \
        --weir-fst-pop ../03_intermediate_files/${pop2}.samples.txt \
        --out ../04_output_files/fst/${pop1}_${pop2}/${pop1}_vs_${pop2}.win_${window_size[${i}]} \
        --fst-window-size ${window_size[${i}]} --fst-window-step ${window_step[${i}]}
    done
done 

# calculate the tajimaD for each group
mkdir ../04_output_files/tajimaD
for i in ${pops[*]};do
    mkdir ../04_output_files/tajimaD/${i}
    for j in {0..3};do
        $vcftools --gzvcf ../02_input_files/$vcf --TajimaD ${window_size[${j}]} --out ../04_output_files/tajimaD/${i}/${i}.win_${window_size[${j}]}
    done
done


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"