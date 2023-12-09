# assembling
qsub -q fat -t 1-32 01.PBS_hifiasm.sh

# anchoring contigs to chromosomes, base Longmi4
qsub -q fat -t 1-32 02.scaffold.sh

# assessing the assembly quality (BUSCO)
./03.run_busco.sh BC027

# assessing the assembly quality (LAI)
qsub -t 1-32 041.PBS_LTRharvest.sh
qsub -t 1-32 042.PBS_LTR_FINDER.sh
for i in {BC027,BC040,BC048,BC100,BC136,BC170,BC188,BC204,BC217,BC235,BC246,BC264,BC292,BC310,BC311,BC315,BC328,BC332,BC350,BC357,BC362,BC382,BC396,BC404,BC407,BC418,BC426,BC434,BC476,BC479,BC494,BC498}; do cat $i.genome.fasta.harvest.scn $i.genome.fasta.finder.combine.scn > $i.genome.fasta.rawLTR.scn; done
qsub -t 1-32 043.PBS_LTR_retriever.sh

# HiFi reads mapping

# HiC

