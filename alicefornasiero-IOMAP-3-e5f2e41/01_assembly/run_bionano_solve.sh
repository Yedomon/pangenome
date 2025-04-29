#!/bin/bash
#SBATCH --job-name=bnsolve
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=batch
#SBATCH --error=/path/to/logs/bnsolve.%J.log
#SBATCH --output=/path/to/logs/bnsolve.%J.out
#SBATCH --time=48:00:00
#SBATCH --mem=100GB

## Software Modules
module load bionano/solve3.4

## Data Variables
CORES=8
HYBRID_DIR=/sw/csi/bionano/solve3.4/el7.5_binary/Solve3.4_06042019a/HybridScaffold/06042019
PIPELINE_DIR=/sw/csi/bionano/solve3.4/el7.5_binary/Solve3.4_06042019a/Pipeline/06042019
REF_ALIGNER=/sw/csi/bionano/solve3.4/el7.5_binary/Solve3.4_06042019a/RefAligner/8949.9232rel
CMAP_DIR=/path/to/optical_map/exp_refineFinal1
PREFIX=Oryza_species_contigs
OUTPUT=/path/to/optical_map
INPUT_FASTA=/path/to/Oryza_species_contig.fasta
mkdir -p ${OUTPUT}

#### Step 1: hybridScaffold
# Run the Hybrid Scaffold pipeline

# perl hybridScaffold.pl
# -n <sequence file in FASTA format>
# -b <Bionano CMAP file>
# -c <hybrid scaffold configuration file in XML format>
# -r <RefAligner binary file>
# -o <output directory>
# -B <conflict filter level genome maps: 1, 2, or 3>
# -N <conflict filter level for sequences: 1, 2, or 3>
# -f <a flag to overwrite existing files; optional>

perl ${HYBRID_DIR}/hybridScaffold.pl \
-n ${INPUT_FASTA} \
-b ${CMAP_DIR}/EXP_REFINEFINAL1.cmap \
-c ${HYBRID_DIR}/hybridScaffold_config.xml \
-r ${REF_ALIGNER}/RefAligner \
-o ${OUTPUT}/${PREFIX} \
-B 2 \
-N 2 \
-f

#### Step 2: fa2cmap_multi_color
# A FASTA file can be converted to a CMAP file through in-silico digestion with a specific enzyme motif. Once converted (in-silico digested), the output CMAP file can be imported into Access for further data analysis, such as alignment.

# perl fa2cmap_multi_color.pl 
# -i <input FASTA file> 
# -e <Enzyme Name (s) or Sequence (s), followed by channel #>

cd ${OUTPUT}
perl ${HYBRID_DIR}/scripts/fa2cmap_multi_color.pl \
-i ${INPUT_FASTA} \
-e DLE-1 1

#### Step 3: runCharacterize
# A CMAP (Bionano assembly consensus genome map) file can be aligned to a sequence reference CMAP or another assembly CMAP. The output alignment (trio files, including .xmap, _q.cmap and _r.cmap) can be imported into Access for visualization.

# python3 runCharacterize.py 
# -t <RefAligner binary> 
# -q <query CMAP> 
# -r <sequence reference CMAP> 
# -p <Pipeline directory>
# -a <assembly opt argument used to generate the consensus map> 
# -n <number of threads to use, default 4> 

python ${PIPELINE_DIR}/runCharacterize.py \
-t ${REF_ALIGNER}/RefAligner \
-q ${CMAP_DIR}/EXP_REFINEFINAL1.cmap \
-r $(dirname $INPUT_FASTA)/Omin_GPM.ctg_DLE1_0kb_0labels.cmap \
-p ${PIPELINE_DIR} \
-a ${REF_ALIGNER}/optArguments_nonhaplotype_saphyr.xml \
-n ${CORES}
