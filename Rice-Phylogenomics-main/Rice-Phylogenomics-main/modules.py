import os
import pandas as pd
from itertools import combinations
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def prefix_cell(cell, header):
    if pd.notna(cell) and cell.strip():  # Check for non-empty and non-whitespace cells
        items = cell.split("|")
        items = [header + item for item in items]
        return "|".join(items)
    else:
        return "NA"

def pangenome_cleaner(pangenome):
    table = []
    f1 = open(pangenome, "r")
    for line in f1:
        # read text data into a data table
        line = line.strip()
        line = line.split("\t")
        e = [line[0],line[1],line[4]] + line[9:]
        table.append(e)
    
    df = pd.DataFrame(table)
    species_code = df.iloc[0] + "@" # extract species code from the first raw, and prefix it to genes
    df.iloc[1:, 3:] = df.iloc[1:, 3:].apply(lambda col: col.apply(lambda cell: prefix_cell(cell, species_code[col.name])))
    return df

# 2. read in the modified pangenome.txt data and converst it to nonSynteny,tandem and synteny datatable, and also datatable for SynNet package
# this script is different from 2_pangenome_parser.py by categorize genes into syntenic genes and all (syntenic, non-syntenic and tandem)
def parse_pangenome(pangenome,all,syn,synnet):
    with open(all,"w") as f2, open(syn,"w") as f3, open(synnet,"w") as f4:
        list_ = []
        for index, row in pangenome.iterrows():
            res = []
            res_new = []
            list_all = []
            list_syn = []
            list_synnet = []
            character = "+*"
            pgid = str(row.iloc[0]) # extract the pangene ID which is in the first coloum, if pangenome eas read in by pd.read_csv(), we need to convert int to str
            #note: if you don't specify the delimiter in string.split(), the split() will take mutile "\t" as a single "\t"
            list_line = row.iloc[1:].tolist() # extract the member of the given pangene from all species, and store them in list;
            print(list_line)
            for element in list_line:
                if element == "NA": # no gene
                    list_all.append(element)
                    #list_tandem.append(element)
                    list_syn.append(element)
                elif "|" in element: # multiple genes (>=2)
                    member = element.split("|")
                    list1 = []
                    list2 = []
                    for e in member:
                        if not e.endswith("*") and not e.endswith("+"):
                            list1.append(e) # syntenic members, the associated tandems(+) are excluded
                        else:
                            e = e.strip(character)
                            list2.append(e) # other members, including syntenic-associated tandems(+) and non-syntenic members
                    if list1 == []: # all genes are non-syntenic (either * or +), means list2 is not empty
                        list1.append("NA")
                    else: # list1 is not empty, means that there are at least one syntenic genes
                        list2 = list1 + list2 # add syntenic genes to non-synetnic list
                    list_syn.append("|".join(list1))
                    list_all.append("|".join(list2))
                elif element.endswith("*"): # only one non-syntenic gene
                    element = element.strip(character)
                    list_all.append(element)
                    list_syn.append("NA")
#                elif element.endswith("+"): # this is not needed, because if there is "+", "|" must exist
#                    list_tandem.append(element)
#                    list_nonSyn.append("NA")
#                    list_syn.append("NA")
                else: # only one syntenic gene
                    list_syn.append(element)
                    list_all.append(element)

            # write results
            if pgid == "pgID": # write the header line
                f2.write(pgid + "\t" + "\t".join(list_line) + "\n")
            else:
                f2.write(pgid + "\t" + "\t".join(list_all) + "\n")
            f3.write(pgid + "\t" + "\t".join(list_syn) + "\n") # the header line has been included, see line 32     
                                    
            # generate non-redundant syntenic gene pairs of all genes stored in list_syn per chromosome
            chromosome = list_line[0]
            if not chromosome in list_ and not pgid == "pgID" and not chromosome == "":
                list_.append(chromosome)
                line = "\t".join(list_line[2:])
                line = line.replace("|", "\t") # split genes from same species as independent elements
                gene_list = line.split("\t")
                for element in gene_list:
                    if not element.endswith("*") and not element.endswith("+") and element != "NA": # this ensures that non-syntenic genes and syntenic tandem duplicates were removed
                        list_synnet.append(element)
                pairs = list(combinations(list_synnet, 2))
                for element in pairs:
                    string = "\t".join(element)
                    res.append(string)

                # remove gene pairs from same species
                for e in res:
                    gene1 = e.split("\t")[0]
                    #id1 = gene1.split("_")[0:3]
                    name1 = gene1.split("_")[0]
                    gene2 = e.split("\t")[1]
                    #id2 = gene2.split("_")[0:3]
                    name2 = gene2.split("_")[0]
                    if name1 != name2: # this will exclude pairs from the same species, in my case sp21 has 3 varieties and all will be removed if a pair has any of them
                        res_new.append(e)
                        #res.remove(e) # there is probablem with this, the iteration will stop after first match
                if len(res_new) >= 1:
                    res_new = [chromosome + "\t" + element for element in res_new]
                    #f4.write("### the clusters on " + chromosome + "\n")
                    f4.write("\n".join(res_new) + "\n")
            elif chromosome in list_ and not pgid == "pgID" and not chromosome == "":
                #print(list_line)
                line = "\t".join(list_line[2:])
                line = line.replace("|", "\t") # split genes from same species as independent elements
                gene_list = line.split("\t")
                for element in gene_list:
                    if not element.endswith("*") and not element.endswith("+") and element != "NA": # this ensures that non-syntenic genes and syntenic tandem duplicates were removed
                        list_synnet.append(element)
                pairs = list(combinations(list_synnet, 2))
                for element in pairs:
                    string = "\t".join(element)
                    res.append(string)

                # remove gene pairs from same species
                for e in res:
                    gene1 = e.split("\t")[0]
                    #id1 = gene1.split("_")[0:3]
                    name1 = gene1.split("_")[0]
                    gene2 = e.split("\t")[1]
                    #id2 = gene2.split("_")[0:3]
                    name2 = gene2.split("_")[0]
                    if name1 != name2: # this will exclude pairs from the same species, in my case sp21 has 3 varieties and all will be removed if a pair has any of them
                        res_new.append(e)
                        #res.remove(e) # there is probablem with this, the iteration will stop after first match
                if len(res_new) >= 1:
                    res_new = [chromosome + "\t" + element for element in res_new]
                    f4.write("\n".join(res_new) + "\n")
    f2.close()
    f3.close()
    f4.close()

def prepare_fasta_from_synCluster(synnet,fasta):
    "after synnet analysis, we get wide formate synnet data: ClusterID\tsynnet_members;"
    "from here, we want to generate fasta file for each Cluster (synteny constrained HOGs ??), and then use this fasta files to do phylogenomic analysis"
    dic = {}
    dic_ = {}
    for record in SeqIO.parse(fasta, "fasta"):
        dic_[record.id] = record.seq
    for line in open(synnet,'r'):
        line = line.strip()
        if not line.startswith("Cluster"):
            cluster = line.split("\t")[0]
            gene_list = line.split("\t")[1].split(",")
            dic[cluster] = gene_list
    for key in dic.keys():
        sequence = []
        file_name = key + ".fasta"
        genes = dic[key]
        for gene in genes:
            new_record = SeqRecord(dic_[gene],gene)
            sequence.append(new_record)
        SeqIO.write(sequence, file_name, "fasta")