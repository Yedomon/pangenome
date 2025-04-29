# 2. read in the modified pangenome.txt data and converst it to nonSynteny,tandem and synteny datatable, and also datatable for SynNet package

#%%
import os
from itertools import combinations

# add speceis name to gene name
def comBed2taxon(comBed):
    dic = {}
    with open(comBed,"r") as f:
        for line in f:
            geneID = line.split("\t")[3]
            species = line.split("\t")[7]
            dic[geneID] = "rice_" + species + "_" + geneID
    return dic


def parse_pangenome(pangenome,all,syn,synnet,combBed):
    with open(pangenome, "r") as f1, open(all,"w") as f2, open(syn,"w") as f3, open(synnet,"w") as f4:
        dic = comBed2taxon(combBed)
        for line in f1:
            res = []
            res_new = []
            list_all = []
            list_syn = []
            list_synnet = []
            character = "+*"
            line = line.strip()
            pgid = line.split("\t")[0] # extract the pangene ID which is in the first coloum
            #note: if you don't specify the delimiter in string.split(), the split() will take mutile "\t" as a single "\t"
            list_line = line.split("\t")[1:] # extract the member of the given pangene from all species, and store them in list;
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
            if line.startswith("pgID"): # write the header line
                header_line = line.split()[9:]
                f2.write(pgid + "\t" + "\t".join(header_line) + "\n")
            else:
                f2.write(pgid + "\t" + "\t".join(list_all) + "\n")
            f3.write(pgid + "\t" + "\t".join(list_syn) + "\n") # the header line has been included, see line 32     
                                    
            # generate non-redundant syntenic gene pairs of all genes stored in list_syn
            if line.startswith("pgID"):
                print("this is the header line, skip")
            else:
                line = line.replace("|", "\t") # split genes from same species as independent elements
                gene_list = line.split("\t")[9:]
                for element in gene_list:
                    if not element.endswith("*") and not element.endswith("+") and element != "NA": # this ensures that non-syntenic genes and syntenic tandem duplicates were removed
                        gene = dic[element] #append speceis_code to rice gene
                        list_synnet.append(gene)
                pairs = list(combinations(list_synnet, 2))
                for element in pairs:
                    string = "\t".join(element)
                    if string.startswith("rice_"):
                        res.append(string)

            # remove gene pairs from same species
            for e in res:
                #print(e)
                gene1 = e.split("\t")[0]
                id1 = gene1.split("_")[1]
                gene2 = e.split("\t")[1]
                id2 = gene2.split("_")[1]
                if id1 != id2: # this will exclude pairs from the same species
                    res_new.append(e)
                    #res.remove(e) # there is probablem with this, the iteration will stop after first match
            if len(res_new) >= 1:
                f4.write("\n".join(res_new) + "\n")
    f2.close()
    f3.close()
    f4.close()


combBed = "combBed.txt"
pangenome = "pangenome.txt"
all = "og_all.txt"
syn = "og_syn.txt"
synnet = "synnet.txt"


work_dir = "path/to/work_dir"
os.chdir(work_dir)
parse_pangenome(pangenome,all,syn,synnet,combBed)
# %%
