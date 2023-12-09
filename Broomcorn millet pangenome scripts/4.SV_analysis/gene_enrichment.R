library(tidyverse)
library(readr)
library(clusterProfiler)
# R 4.1 or later (need bioconductor to install clusterProfiler)

        # read the standard Gene Ontology and KEGG annotation database
        go_db <- read_tsv("/public1/home/liuyang/Project/00.reference/db.go.tsv")
        ko_db <- read_tsv("/public1/home/liuyang/Project/00.reference/db.kegg.tsv")
        pfam_db <- read_tsv("/public1/home/liuyang/Project/00.reference/db.pfam.tsv")

        # read genes from regions
        set_genes <- read_tsv("PAV.DE.candidate.gwl.lst", col_names = c("gene"))
        set_genes <- as.vector(set_genes$gene)

        # read longmi whole genome Gene Onotology and KEGG annotation
        longmi_go <- read_tsv("/public1/home/liuyang/Project/00.reference/Longmi.go.tsv", col_names = c("gene", "GOID"))
        longmi_ko <- read_tsv("/public1/home/liuyang/Project/00.reference/longmi.kegg.tsv", col_names = c("gene", "KO"))
        longmi_pfam <- read_tsv("/public1/home/liuyang/Project/00.reference/Longmi.pfam.tsv", col_names = c("gene", "PFAMID"))

        # prepare enricher data
	# GO
        go_term2gene <- left_join(longmi_go, go_db, by="GOID") %>% select("GOID", "gene")
        go_term2name <- left_join(longmi_go, go_db, by="GOID") %>% select("GOID", "TERM")

        names(go_term2gene) <- c("go_term", "gene")
        names(go_term2name) <- c("go_term", "name")
	# KEGG
        kegg_term2gene <- left_join(longmi_ko, ko_db, by="KO") %>% select("KO", "gene")
        kegg_term2name <- left_join(longmi_ko, ko_db, by="KO") %>% select("KO", "INFO")

        names(kegg_term2gene) <- c("ko_term", "gene")
        names(kegg_term2name) <- c("ko_term", "name")
	# PFAM
	pfam_term2gene <- left_join(longmi_pfam, pfam_db, by="PFAMID") %>% select("PFAMID", "gene")
	pfam_term2name <- left_join(longmi_pfam, pfam_db, by="PFAMID") %>% select("PFAMID", "INFO")

	names(pfam_term2gene) <- c("pfam_term", "gene")
	names(pfam_term2name) <- c("pfam_term", "name")

        # go enrich
        go_enrich <- enricher(gene = set_genes, pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = go_term2gene, TERM2NAME = go_term2name)
        go_enrich_result <- as_tibble(go_enrich@result) %>% arrange(p.adjust) %>% select("ID", "p.adjust", "Description", "GeneRatio", "BgRatio", "geneID")

        write_tsv(go_enrich_result, file = paste0("enrichment/", "PAV.DE.candidate.gwl.lst.go.enrich.tsv"))
        
        # go DAG plot
        #pdf(file=paste0("enrichment/plot/", j, ".win_", i, ".go.bp.tree.pdf"), width=10, height=15)
        #plotGOgraph(go_enrich)
        #dev.off()
        
        #pdf(file=paste0("enrichment/plot/", j, ".win_", i, ".go.mf.tree.pdf"), width=10, height=15)
        #plotGOgraph(go_enrich.MF)
        #dev.off()
        
        #pdf(file=paste0("enrichment/plot/", j, ".win_", i, ".go.cc.tree.pdf"), width=10, height=15)
        #plotGOgraph(go_enrich.CC)
        #dev.off()
        #barplot(c3_go_enrich, showCategory = 40)

        # kegg enrich
        ko_enrich <- enricher(gene = set_genes, pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = kegg_term2gene, TERM2NAME = kegg_term2name)
        ko_enrich_result <- as_tibble(ko_enrich@result) %>% arrange(p.adjust) %>% select("ID", "p.adjust", "Description", "GeneRatio", "BgRatio", "geneID")

        write_tsv(ko_enrich_result, file = paste0("enrichment/", "PAV.DE.candidate.gwl.lst.ko.enrich.tsv"))

	#barplot(all_ko_enrich, showCategory = 40)

	# pfam enrich
	pfam_enrich <- enricher(gene = set_genes, pvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = pfam_term2gene, TERM2NAME = pfam_term2name)
	pfam_enrich_result <- as_tibble(pfam_enrich@result) %>% arrange(p.adjust) %>% select("ID", "p.adjust", "Description", "GeneRatio", "BgRatio", "geneID")

	write_tsv(pfam_enrich_result, file = paste0("enrichment/", "PAV.DE.candidate.gwl.lst.pfam.enrich.tsv"))
	
	###EOF
