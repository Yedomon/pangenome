# Rice Phylogenomics
The code was used in our rice genome evolution [paper](https://doi.org/10.1101/2024.05.29.596369) to generate microsynteny gene clusters and to make phylogeny based on the clusters. The main function is:
1) to process and clean the pangene table from GENESPACE;
2) to identify syntenic gene clusters;
3) to infer gene trees and phylogeny of (sub)genomes

### Installation
Python3 and Pandans are required to run the scripts; In addition, the following tools need to be callable from working envitonment: [GENESPACE ](https://github.com/jtlovell/GENESPACE), [MAFFT](https://mafft.cbrc.jp/alignment/server/index.html), [trimAl](https://trimal.cgenomics.org/), [IQ-TREE](http://www.iqtree.org/), and [Astral-Pro](https://github.com/chaoszhang/A-pro)

### Example
The computed rice pangenome data is used as an example to demonstrate the pipeline
#### 1. run GENESPACE on your genomes, and use the following settings in **query_pangenes()** to output pangenome:
```R
pangenome <- query_pangenes(
  gsParam,
  bed = NULL,
  refGenome = Reference,
  transform = TRUE,
  showArrayMem = FALSE,
  showNSOrtho = FALSE,
  maxMem2Show = Inf,
  showUnPlacedPgs = FALSE)
```
The ouput is a plain text table of all homologous groups identifed by GENESPACE, whcih will be used in next step.
#### 2. With the pangene table generated from GENESPACE, run SynCluster to identify syntenic gene clusters 

```python
pangenome = "pangenome.txt"
all = "pangenome_all.txt"
syn = "pangenome_syn.txt"
synnet = "pangenome_synnet.txt"

clean_pangenome = pangenome_cleaner(pangenome)
parse_pangenome(clean_pangenome,all,syn,synnet)
```
The outputs are three plain text tables: the cleaned homologous groups, the homologous groups with only syntenic genes and the syntenic gene pairs.
#### 3. Explore the syntenic gene clusters in [syntenet](https://github.com/almeidasilvaf/syntenet)
#### 4. Generate FASTA sequences for each syntenic gene cluster
As we have the syntenic clusters generated in step2, now we will extract the protein sequences for each syntenic clusters. 
```python
prepare_fasta_from_synCluster(syn,fasta)
```
syn is the output from step2, fasta is the protein sequences used in GENESPACE run; 
Output is the protein sequence in FASTA format, one fasta file for each syntenic cluster
#### 5. Infer gene tree for each syntenic cluster
we follow the [pipeline](https://bitbucket.org/yangya/phylogenomic_dataset_construction/src/master/) to make gene tree for each syntenic cluster
#### 6. Infer a phylogeny based on the gene trees using coalescent algorithm
```sh
astral-pro -i $geneTrees -o $genomeTree -t 8 -u 1
```
$geneTrees is the tree inferred in last step, and $genomeTree is the coalescent phylogeny of the (sub)genomes

