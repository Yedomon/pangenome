library(tidyverse)

if(FALSE) {
    "
    This R script is used to smooth the XP-CLR using 100kb window 10kb step.
    "
}


pops=c("C_vs_W", "C1_vs_W", "C2_vs_W", "C3_vs_W", "C1_vs_C2", "C1_vs_C3")
for(i in pops){
    rolling <- readr::read_delim(paste0("input/xpclr_raw/", i,".win_0.002.txt"),
                                 col_names = c("CHR", "RANK", "SNP_count", "BP", "CM", "P", "MAXS")) %>% filter(P!=Inf)

    genome_length <- readr::read_delim("~/dataset/mizi/Pmlongmi4.genome.fa.fai", col_names = c("CHR", "Length", "A", "B", "C"))

    rolling$MID <- rolling$BP + 10000/2
    rolling$P <- as.double(rolling$P)
    # 100kb window 10kb step 

    res <- tibble(CHR = 1, MID = seq(50000, genome_length$Length[1]-50000, by=10000), AVAE = 0)

    for(j in 2:18){
        res <- rbind(res, tibble(CHR = j, MID = seq(50000, genome_length$Length[j] - 50000, by = 10000), AVAE = 0))
    }

    temp_avae <- c()

    for(j in 1:nrow(res)){
        chr <- res$CHR[j]
        mid <- res$MID[j]
        temp_avae[j] <- mean(filter(rolling, CHR == chr, MID > mid - 50000, MID < mid + 50000)$P)
    }
    res$START <- res$MID - 50000
    res$END <- res$MID + 50000
    res$AVAE <- temp_avae
    write_tsv(res, file=paste0("smooth/", i, ".smooth.100kb.10kb.txt"))
}