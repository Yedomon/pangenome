# K-mer Analysis Workflow

Workflow of commands for k-mer analysis using KMC, GenomeScope and SmudgePlot to estimate genome size, heterozygosity and ploidy.

## Table of Contents

- [Requirements](#requirements)  
- [Input Parameters](#input-parameters)  
- [Pipeline Steps](#pipeline-steps)  
  - [1. K-mer Counting with KMC](#1-k-mer-counting-with-kmc)  
  - [2. Estimation of Coverage Thresholds](#2-estimation-of-coverage-thresholds)  
  - [3. Alternative: GenomeScope Analysis](#3-alternative-genomescope-analysis)  
  - [4. Extract and Analyze Genomic K-mers](#4-extract-and-analyze-genomic-k-mers)  
  - [5. Generate Smudge Plot](#5-generate-smudge-plot)  
- [References](#references)  
- [Author](#author)  

## Requirements

The following software modules are required for this workflow. You can find installation instructions on their respective GitHub pages:

- **KMC**:
  GitHub: (https://github.com/refresh-bio/KMC)
  
- **GenomeScope**:
  GitHub: (https://github.com/schatzlab/genomescope)
  
- **Smudgeplot**:
  GitHub: (https://github.com/KamilSJaron/smudgeplot)

## Input Parameters

The following variables should be configured before running the script:

| Parameter   | Description                                                  |
| ----------- | ------------------------------------------------------------ |
| `outdir`    | Path to the output directory                                 |
| `filelist`  | Path to the file containing a list of input FASTQ file paths |
| `kmer`      | K-mer length (default: 19)                                   |
| `nthreads`  | Number of threads for computation (default: 32)              |
| `ploidy`    | Expected ploidy level (e.g., diploid=2, tetraploid=4)        |
| `outprefix` | Prefix for output files                                      |

## Pipeline Steps

### 1. K-mer Counting with KMC

```bash
module load kmc/3.1.2
kmc -k${kmer} -v -t${nthreads} -m100 -cs10000 @${filelist} ${outprefix} ./
kmc_tools transform ${outprefix} histogram ${outprefix}.histo
module unload kmc/3.1.2
```

### 2. Estimation of Coverage Thresholds

- The lower (L) and upper (U) coverage thresholds are estimated using Smudgeplot.

```bash
module load smudgeplot/0.2.4
L=$(smudgeplot.py cutoff ${outprefix}.histo L)
U=$(smudgeplot.py cutoff ${outprefix}.histo U)
echo $L $U
module unload smudgeplot/0.2.4
```

### 3. Alternative: GenomeScope Analysis

- GenomeScope can be used to visually inspect k-mer coverage histograms to estimate L/U values.

```bash
module load kmc/3.1.2
Rscript /ibex/sw/csi/kmc/3.1.2/gnu6.4.0/genomescope2.0/genomescope.R \
  -i ${outdir}/${outprefix}.histo \
  -n ${outprefix}_p${ploidy}_Gscope \
  -o ${outdir} \
  -k ${kmer} \
  -p ${ploidy} &>${outdir}/${outprefix}_p${ploidy}_Gscope.log
module unload kmc/3.1.2
```

### 4. Extract and Analyze Genomic K-mers

- KMC is used to extract k-mers in the range defined by the lower (L) and upper (U) coverage thresholds.

```bash
module load kmc/3.1.2
kmc_tools transform ${outprefix} -ci"$L" -cx"$U" reduce ${outprefix}_L"$L"_U"$U"
kmc_tools transform ${outprefix} -ci"$L" -cx"$U" dump -s ${outprefix}_L"$L"_U"$U".dump
module unload kmc/3.1.2
```

### 5. Generate Smudge Plot

- The extracted k-mers are processed to compute k-mer pairs and generate the smudge plot.

```bash
module load smudgeplot/0.2.4
smudgeplot.py hetkmers -o ${outprefix}_L"$L"_U"$U" < ${outprefix}_L"$L"_U"$U".dump
smudgeplot_plot.R -i ${outprefix}_L"$L"_U"$U"_coverages.tsv \
  -o ${outprefix}_smudge \
  -k ${kmer} \
  --title ${outprefix}
module unload smudgeplot/0.2.4
```

## References

If you use this pipeline in your research, please cite the respective tools:

- **KMC**: Kokot et al., *Bioinformatics*, 2017. [DOI: 10.1093/bioinformatics/btx304](https://doi.org/10.1093/bioinformatics/btx304)
- **GenomeScope**: Ranallo-Benavidez et al., *Nature Communications*, 2020. [DOI: 10.1038/s41467-020-14998-3](https://doi.org/10.1038/s41467-020-14998-3)

## Author

Alice Fornasiero - King Abdullah University of Science and Technology

