### Rice Chloroplast Phylogenetic Analysis

This repository contains a step-by-step script for reconstructing the phylogenetic relationships of chloroplast genomes across the Oryza genus. The analysis supports the findings in:

Fornasiero et al., "Oryza genome evolution through a tetraploid lens," Nature Genetics, 2025.

### Contents

- Rice_Chloroplast_Phylogenetic_Analysis_Script.sh - A script with all necessary commands for sequence alignment, trimming, and phylogenetic tree construction.
- New_27_cp_genomes_correct_A_genomes.fa - Input FASTA file containing concatenated chloroplast genome sequences.
- Output files generated at each step:
  - New_27_cp_genomes_correct_A_genomes_MAFFT_auto_LSC.fa (LSC regions alignment)
  - New_27_cp_genomes_correct_A_genomes_MAFFT_auto_LSC_trimal_auto.fa (trimmed alignment)
  - New_27_cp_genomes_correct_A_genomes_MAFFT_auto_LSC_trimal_auto.treefile (ML tree)
  - New_27_cp_genomes_correct_A_genomes_MAFFT_auto_LSC_trimal_auto.fa.contree (concensus tree)

### Usage

### Clone the repository:
git clone https://github.com/nam-hoang/rice_cp_phylo_analysis.git && cd rice_cp_phylo_analysis
