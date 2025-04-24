gcta64 --bfile SV --make-grm --make-grm-alg 1 --out A1
gcta64 --grm A1 --pheno your_data.phen --reml --out phe
gcta64 --grm A1 --pheno your_data.phen --reml --out phe --reml-pred-rand
