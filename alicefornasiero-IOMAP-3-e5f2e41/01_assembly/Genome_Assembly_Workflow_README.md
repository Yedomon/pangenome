# Genome Assembly Workflow

Workflow of commands for assembling the genomes of the wild *Oryza* species using different software for PacBio long-read sequencing data. 

The software used are:
- **Mecat**, **Canu**, **Flye** for PacBio CLR (Circular Long Reads) assembly
- **pbmm2**, **Arrow** and **Pilon** for CLR assembly polishing
- **HiFiAsm** and **HiCanu** for PacBio HiFi assembly
- **BioNano Solve** for BioNano hybrid scaffolding using BioNano genome maps

## Table of Contents

- [Requirements](#requirements)

- [1. PacBio CLR Assembly](#1-pacbio-clr-assembly-1)
  - [MECAT](#mecat)
  - [Canu](#canu)
  - [Flye](#flye)

- [2. PacBio CLR Assembly Polishing](#2-pacbio-clr-assembly-polishing-1)

- [3. PacBio HiFi Assembly](#3-pacbio-hifi-assembly)
  - [HiFiasm](#hifiasm)
  - [HiCanu](#hicanu)

- [4. Bionano Genome Mapping](#4-bionano-genome-mapping)

- [References](#references)

- [Author](#author)

---

## Requirements

The following tools are required for this workflow. You can find installation instructions on their respective GitHub pages:

### 1. PacBio CLR Assembly

- **MECAT**: For de novo assembly of CLR reads  
  GitHub: [https://github.com/xiaochuanle/MECAT](https://github.com/xiaochuanle/MECAT)

- **Canu**: For de novo assembly of CLR reads  
  GitHub: [https://github.com/marbl/canu](https://github.com/marbl/canu)

- **Flye**: For de novo assembly of CLR reads  
  GitHub: [https://github.com/fenderglass/Flye](https://github.com/fenderglass/Flye)

### 2. PacBio CLR Assembly Polishing

- **pbmm2**: For aligning PacBio subreads to the reference contig assembly  
  GitHub: [https://github.com/PacificBiosciences/pbmm2](https://github.com/PacificBiosciences/pbmm2)

- **samtools**: For merging and indexing BAM files  
  GitHub: [https://github.com/samtools/samtools](https://github.com/samtools/samtools)

- **pbindex**: For indexing PacBio BAM files  
  GitHub: [https://github.com/PacificBiosciences/pbindex](https://github.com/PacificBiosciences/pbindex)

- **Arrow**: For polishing genome assemblies using PacBio reads  
  GitHub: [https://github.com/PacificBiosciences/gcpp](https://github.com/PacificBiosciences/gcpp)

- **BWA**: To align short Illumina reads  
  GitHub: [https://github.com/lh3/bwa](https://github.com/lh3/bwa)
  
- **Pilon**: For polishing CLR read assembly  
  GitHub: [https://github.com/broadinstitute/pilon](https://github.com/broadinstitute/pilon)

### 3. PacBio HiFi Assembly

- **HiFiAsm**: For de novo assembly of HiFi reads  
  GitHub: [https://github.com/PacificBiosciences/HiFiAsm](https://github.com/PacificBiosciences/HiFiAsm)

- **HiCanu**: For de novo assembly of HiFi reads  
  GitHub: [https://github.com/marbl/canu](https://github.com/marbl/canu)

### 4. Bionano Genome Map

- **BioNano Solve**: For processing Bionano genome mapping data, including hybrid scaffolding and structural variation detection  
  Official Site: [https://bionanogenomics.com/support/software-downloads/](https://bionanogenomics.com/support/software-downloads)

--- 

## 1. PacBio CLR assembly

### MECAT
You can run MECAT with the following commands:

```bash
# Configuring MECAT
mecat.pl config oryza_spp_config_file.txt

# Correcting reads
mecat.pl correct oryza_spp_config_file.txt

# Assembling the genome
mecat.pl assemble oryza_spp_config_file.txt
```

#### MECAT Configuration file example
You need to generate the configuration file (e.g. oryza_spp_config_file.txt) using the "config" command, and then modify it manually:

```txt
PROJECT=Oryza_species
RAWREADS=/path/to/directory/Oryza_species.fasta
GENOME_SIZE=1008000000
THREADS=32
MIN_READ_LENGTH=500
CNS_OVLP_OPTIONS=""
CNS_OPTIONS="-r 0.6 -a 1000 -c 4 -l 2000"
CNS_OUTPUT_COVERAGE=30
TRIM_OVLP_OPTIONS="-B"
ASM_OVLP_OPTIONS="-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400"
FSA_OL_FILTER_OPTIONS="--max_overhang=-1 --min_identity=-1"
FSA_ASSEMBLE_OPTIONS=""
CLEANUP=0
GRID_NODE=0
USE_GRID=false
```

---
### Canu
Use the following command to run Canu genome assembly:

```bash
time canu -d output -p Oryza_species genomeSize=1008m -pacbio-raw Oryza_species.fasta gnuplot=undef usegrid=1 gridOptions="--time=10-00:00:00 --partition=batch --mem-per-cpu=16g" gridOptionsJobName=Oryza-using-grid
```

---
### Flye
Use the following command to run Flye genome assembly:

```bash
flye --pacbio-raw Oryza_species.fasta --out-dir output --genome-size 1008000000 --threads n_threads --asm-coverage 35
```

---

## 2. PacBio CLR Assembly Polishing

The workflow utilizes **pbmm2**, **samtools**, **pbindex**, and **Arrow** for CLR (Circular Long Reads) assembly polishing.

### Workflow Steps

1. Align PacBio subreads to the contig assembly using **pbmm2**:  
   ```bash
   pbmm2 align input subreads_1.bam outdir/subreads_1.sort.bam --sort -j 24 -J 8 tmp_dir
   pbmm2 align input subreads_2.bam outdir/subreads_2.sort.bam --sort -j 24 -J 8 tmp_dir
   ```

2. Merge sorted BAM files using **samtools**  
   ```bash
   samtools merge outdir/output_prefix.merged.bam outdir/*.sort.bam
   ```

3. Index merged BAM file using **pbindex**  
   ```bash
   pbindex outdir/output_prefix.merged.bam
   ```

4. Index merged BAM file using **samtools**  
   ```bash
   samtools index outdir/output_prefix.merged.bam outdir/output_prefix.merged.bam.bai
   ```

5. Index contig assembly using **samtools**
   ```bash
   samtools faidx Oryza_species.fasta
   ```

6. Polish assembly with **Arrow**  
   ```bash
   gcpp --algorithm=arrow -j n_threads -r ${INPUT} -o outdir/output_prefix.arrow.fasta outdir/output_prefix.merged.bam
   ```

7. Create index of the arrowed polished bam files using **BWA** 
   ```bash
   bwa index -p outdir/output_prefix.arrow.bwa outdir/output_prefix.arrow.fasta
   ```

8. Short reads mapping to the arrow-polished assembled reads using **BWA**
   ```bash
   bwa mem -t 32 outdir/output_prefix.arrow.bwa Illumina_R1.fastq Illumina_R2.fastq | samtools view -bhS - > outdir/output_prefix.arrow.bam
   ```  

9. Sort reads using **samtools**
   ```bash
   samtools sort -@ 32 -o outdir/output_prefix.arrow_sorted.bam outdir/output_prefix.arrow.bam
   ```

10. Index sorted using **samtools**
   ```bash
   samtools index outdir/output_prefix.arrow_sorted.bam
   ```

11. Polish the genome assembly using **Pilon**
   ```bash
   java -Xmx450G -jar $PILON --genome outdir/output_prefix.arrow.fasta --frags outdir/output_prefix.arrow_sorted.bam --output output_prefix.arrow.pilon --outdir outdir/ --threads 32 --change
   ```

---

## 3. PacBio HiFi assembly

### HiFiasm 
Use the following command to run HiFiasm genome assembly:

```bash
hifiasm -o output -t 32 --primary m64041_211103_120859.hifi_reads.fasta.gz m64313e_211202_171836.hifi_reads.fasta.gz
```

---

### HiCanu
Use the following command to run HiCanu genome assembly:

```bash
canu -p output_dir -d prefix_name genomeSize=800m -pacbio-hifi m64068_230209_044215.hifi_reads.fastq m64068_230210_134855.hifi_reads.fastq m64068_230211_225339.hifi_reads.fastq m64068_230213_080152.hifi_reads.fastq -useGrid=false
```

---

## 4. Bionano genome mapping
The following commands uses the Bionano assembly consensus genome map (CMAP) to generate hybrid scaffolds and align CMAPs:

```bash
# Generate the hybrid scaffold using the Bionano assembly consensus genome map
perl hybridScaffold.pl -n contig.fasta -b EXP_REFINEFINAL1.cmap -c hybridScaffold_config.xml -r RefAligner -o Oryza_species -B 2 -N 2 -f

# Convert a FASTA file to a CMAP file through in-silico digestion with a specific enzyme motif
perl fa2cmap_multi_color.pl -i contig.fasta -e DLE-1 1

# Align the CMAP to a sequence reference CMAP.
python runCharacterize.py -t RefAligner -q EXP_REFINEFINAL1.cmap -r Oryza_species.cmap -p pipeline_dir -a optArguments_nonhaplotype_saphyr.xml -n 8

```

---

## References

- **MECAT2**: Chuan-Le X. et al. *Nature Methods*, 2017.  
  https://doi.org/10.1038/nmeth.4432

- **Canu**: Koren S. et al. *Genome Res.*, 2017.  
  http://www.genome.org/cgi/doi/10.1101/gr.215087.116.
  
- **Flye**: Kolmogorov M. et al. *Nat. Biotechnol.*, 2019.  
  https://www.nature.com/articles/s41587-019-0072-8

- **pbmm2**: https://github.com/PacificBiosciences/pbbioconda

- **samtools**: Danecek P. et al. *GigaScience*, 2021.  
  https://academic.oup.com/gigascience/article/10/2/giab008/6137722?login=true
  
- **pbindex**: https://github.com/PacificBiosciences/pbbioconda

- **Arrow**: https://github.com/PacificBiosciences/pbbioconda  
  https://github.com/PacificBiosciences/gcpp

- **BWA**: Li H. & Durbin R. *Bioinformatics*, 2010.  
  https://academic.oup.com/bioinformatics/article/26/5/589/211735?login=true

- **Pilon**: Walker B.J. et al. *PLoS ONE*, 2014.  
  https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112963

## Author

Alice Fornasiero - King Abdullah University of Science and Technology

