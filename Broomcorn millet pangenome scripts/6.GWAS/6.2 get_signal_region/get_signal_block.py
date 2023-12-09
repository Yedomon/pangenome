import os,sys
import pandas as pd
df = pd.read_csv(sys.argv[1],header=None,sep='\t')
df.columns=["pos1","env_count","pos2","ef","p","loc"]

df['chr'] = df["pos1"].str.split('_',expand=True)[0]
df['pos'] = df["pos1"].str.split('_',expand=True)[1]
df['pos'] = df['pos'].astype('int')
df['chr'] = df['chr'].astype('str')
chroms = list(df['chr'].unique())

for chrs in chroms:
    loci = []
    size = []
    sub1 = df[df['chr']==chrs]
    sub1=sub1.sort_values("pos")
    #print(sub1)
    data=sub1.values.tolist()
    lastchr=""
    laststart=""
    group={}
    start=""
    pvalue={}
    env={}
    ld=int(200000)
    for line in data:
        id=line[0]
        p=line[4]
        chr,pos=id.split("_")
        pvalue[id]=p
        env[id]=line
        if start=="":
            start=pos
            group[start]=[id]
            laststart=pos
            continue
        if int(laststart)+ld>=int(pos):
            group[start].append(id)
            laststart=pos
        else:
            group[start].append(int(laststart)+ld)
            start=pos
            group[start]=[id]
            laststart=pos
    group[start].append(int(laststart)+ld)
    for key,values in group.items():
        p=[float(pvalue[i]) for i in values[:-1]]  
        minp=min(p)
        minindex=p.index(min(p))
        minpos=values[minindex]
        count=int(len(values))-1
        print(int(key)-200000,values[-1],count,minp,minpos,env[minpos][1],env[minpos][5],sep="\t")