library(readr)
library(dplyr)


if(FALSE) {
    "
    This R script is used to select top 5% PI regions between populations.
    "
}

window_size <- c("10000", "20000", "50000", "100000")

for(i in window_size){
    C1 <- read_tsv(file = paste0("~/project/PI_tajimaD_Fst_calculation/04_output_files/pi/C1/C1.win_",i,".windowed.pi"), col_types = "iiiid")
    C2 <- read_tsv(file = paste0("~/project/PI_tajimaD_Fst_calculation/04_output_files/pi/C2/C2.win_",i,".windowed.pi"), col_types = "iiiid")
    C3 <- read_tsv(file = paste0("~/project/PI_tajimaD_Fst_calculation/04_output_files/pi/C3/C3.win_",i,".windowed.pi"), col_types = "iiiid")
    W <- read_tsv(file = paste0("~/project/PI_tajimaD_Fst_calculation/04_output_files/pi/W/W.win_",i,".windowed.pi"), col_types = "iiiid")
    C <- read_tsv(file = paste0("~/project/PI_tajimaD_Fst_calculation/04_output_files/pi/C/C.win_",i,".windowed.pi"), col_types = "iiiid")

    C1_vs_C2 <- left_join(C1, C2, by = c("CHROM", "BIN_START", "BIN_END"), na_matches = "never")
    C1_vs_C3 <- left_join(C1, C3, by = c("CHROM", "BIN_START", "BIN_END"), na_matches = "never")
    C1_vs_W <- left_join(C1, W, by = c("CHROM", "BIN_START", "BIN_END"), na_matches = "never")
    C2_vs_W <- left_join(C2, W, by = c("CHROM", "BIN_START", "BIN_END"), na_matches = "never")
    C3_vs_W <- left_join(C3, W, by = c("CHROM", "BIN_START", "BIN_END"), na_matches = "never")
    C_vs_W <- left_join(C, W, by = c("CHROM", "BIN_START", "BIN_END"), na_matches = "never")

    C1_vs_C2$diff <- C1_vs_C2$PI.x / C1_vs_C2$PI.y
    C1_vs_C3$diff <- C1_vs_C3$PI.x / C1_vs_C3$PI.y
    C1_vs_W$diff <-  C1_vs_W$PI.y / C1_vs_W$PI.x
    C2_vs_W$diff <-  C2_vs_W$PI.y / C2_vs_W$PI.x
    C3_vs_W$diff <-  C3_vs_W$PI.y / C3_vs_W$PI.x
    C_vs_W$diff <-  C_vs_W$PI.y / C_vs_W$PI.x
    
    c1_c2_top10_regions <- slice_max(C1_vs_C2, diff, prop=0.5) %>% select(1,2,3)
    c1_c3_top10_regions <- slice_max(C1_vs_C3, diff, prop=0.5) %>% select(1,2,3)
    c1_w_top10_regions <- slice_max(C1_vs_W, diff, prop=0.5) %>% select(1,2,3)
    c2_w_top10_regions <- slice_max(C2_vs_W, diff, prop=0.5) %>% select(1,2,3)
    c3_w_top10_regions <- slice_max(C3_vs_W, diff, prop=0.5) %>% select(1,2,3)
    c_w_top10_regions <- slice_max(C_vs_W, diff, prop=0.5) %>% select(1,2,3)

    write_tsv(c1_c2_top10_regions, file=paste0("pi_top/pi_top50/C1_vs_C2.top50.pi.win_", i,".regions.txt"))
    write_tsv(c1_c3_top10_regions, file=paste0("pi_top/pi_top50/C1_vs_C3.top50.pi.win_", i,".regions.txt"))
    write_tsv(c1_w_top10_regions, file=paste0("pi_top/pi_top50/C1_vs_W.top50.pi.win_", i,".regions.txt"))
    write_tsv(c2_w_top10_regions, file=paste0("pi_top/pi_top50/C2_vs_W.top50.pi.win_", i,".regions.txt"))
    write_tsv(c3_w_top10_regions, file=paste0("pi_top/pi_top50/C3_vs_W.top50.pi.win_", i,".regions.txt"))
    write_tsv(c_w_top10_regions, file=paste0("pi_top/pi_top50/C_vs_W.top50.pi.win_", i,".regions.txt"))
}