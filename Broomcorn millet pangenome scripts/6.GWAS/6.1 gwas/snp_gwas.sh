vcftools --vcf BC_516_pass_snp.vcf --maf 0.05 --max-missing 0.9  --recode --recode-INFO-all --out BC_516_snp
#beagle
java -Xmx100000m -jar ~/software/beagle.27Jan18.7e1.jar gt=BC_516_snp.recode.vcf out=BC_516_snp_beagle
gunzip BC_516_snp_beagle.vcf.gz
#PCA
plink --vcf BC_516_snp_beagle.vcf --allow-extra-chr  --pca 8 -out PCA
cat PCA.eigenvec | awk '{print $1,$2,"1",$3,$4,$5}' > PCA_snp.txt
#emmax
plink --vcf BC_516_snp_beagle.vcf --recode 12 transpose --out BC_516_snp_emmax_in --allow-extra-chr
emmax-kin BC_516_snp_emmax_in -v -d 10  

for i in `ls ./trait` ;do
emmax -v -d 10 -t BC_516_snp_emmax_in \
				-p ./trait/$i  \
				-k BC_516_snp_emmax_in.BN.kinf \
				-c PCA_snp.txt \
				-o result/${i}
done

convert_emax2manhatton.py ./result
rm -rf for_manhatton
mkdir for_manhatton
mv *_for_manhaton.txt for_manhatton
Rscript manhatton.R 