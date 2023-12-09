import os,sys,glob
import pandas as pd
datatotal=[]
hash={}
files = glob.glob(sys.argv[1]+'/*.txt')
for fs in files:
    name = os.path.basename(fs).split('.')[0]
    with open (fs,'r' )as f:
        for line in f:
            line=line.strip()
            data=line.split("\t")
            data.append(name)
            datatotal.append(data)
for cell in datatotal:
    peak=cell[0]
    if peak in hash:
        hash[peak].append(cell)
    else:
        hash[peak]=[cell]
        

outfile="%s.txt"%sys.argv[1]
OUT=open(outfile,'w+')
#title="\t".join(["peak","location_num","chr","start","end","count","p","peak","loc"])
for key,value in hash.items():
    if len(value)>1:
        pdvalue=pd.DataFrame(value,columns=["peak","ef","p","loc"])
        minp=pdvalue[pdvalue['p'].astype("float")==pdvalue['p'].astype("float").min()]
        minpvalue="\t".join(minp.iloc[0,].values.tolist())
        print(key,len(value),minpvalue,sep="\t",file=OUT)
OUT.close()