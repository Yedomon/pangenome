#filter
vcftools --vcf BC_SVtyper_paragraph.vcf --maf 0.05 --max-missing 0.5  --recode --recode-INFO-all --out BC_516_SV
#emmax
plink --vcf BC_516_SV.recode.vcf --recode 12 transpose --out BC_516_SV_emmax_in --allow-extra-chr
emmax-kin BC_516_SV_emmax_in -v -d 10  
for i in `ls ./trait` ;do
emmax -v -d 10 -t BC_516_SV_emmax_in \
				-p ./trait/$i  \
				-k BC_516_SV_emmax_in.BN.kinf \
				-c PCA_snp.txt \
				-o result/${i}
done
convert_emax2manhatton.py ./result
mkdir for_manhatton
mv *_for_manhaton.txt for_manhatton
Rscript manhatton.R 