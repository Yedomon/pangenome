library(readr)
library(dplyr)


if(FALSE) {
    "
    This R script is used to select top 5% Fst regions.
    "
}

window_size <- c("100000", "10000", "50000", "20000")
pop_pair <- c("C2_vs_C1", "C3_vs_C1")

for(i in pop_pair){
    for(j in window_size){
        fst <- read_tsv(file = paste0("~/project/PI_tajimaD_Fst_calculation/04_output_files/fst/", i, "/", i,".win_",j,".windowed.weir.fst")
                        , col_select = c(1, 2, 3, 6))
        fst_top10_regions <- slice_max(fst, MEAN_FST, prop = 0.1) %>% select(1, 2, 3)
        fst_top5_regions <- slice_max(fst, MEAN_FST, prop = 0.05) %>% select(1, 2, 3)
        fst_top1_regions <- slice_max(fst, MEAN_FST, prop = 0.01) %>% select(1, 2, 3)

        write_tsv(fst_top10_regions, file = paste0("fst_top10/",i,".fst.top10.win_",j,".regions.txt"))
        write_tsv(fst_top5_regions, file = paste0("fst_top5/",i,".fst.top5.win_",j,".regions.txt"))
        write_tsv(fst_top1_regions, file = paste0("fst_top1/",i,".fst.top1.win_",j,".regions.txt"))
    }
}