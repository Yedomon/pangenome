# repeat annotation
qsub 11.PBS_RepeatModeler.sh
qsub 12.PBS_RepeatMasker.sh

# assembling the transcripts
qsub 21.PBS_Trimmomatic.sh
qsub 22.PBS_hisat2.sh
qsub 23.PBS_stringtie.sh

# gene annotation
qsub 31.PBS_maker.sh

# gene function annatation
qsub 41.PBS_ips.sh
