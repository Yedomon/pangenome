import os,re,sys
def convert(gwas_file):
    print(gwas_file)
    outfile=gwas_file+"_for_manhaton.txt"
    OUT=open(outfile,'w')
    with open (gwas_file,'r') as f:
        print("SNP","CHR","POS","PS",sep="\t",file=OUT)
        for line in f:
            line=line.strip()
            SNP,SE,P=line.split("\t")
            chr,pos=SNP.split("_")
            chr=re.sub(r'S',"",chr)
            print(SNP,chr,pos,P,sep="\t",file=OUT)
    OUT.close()

path=sys.argv[1]
for i in os.listdir(path):
    if "result.ps" in i:
        convert("%s/%s"%(path,i))
