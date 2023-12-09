library(ggplot2)
library(ggthemes)
library(readr)
library(dplyr)
library(reshape2)
library(cowplot)
library(patchwork)


if(FALSE) {
    "
    1. Read the XPCLR, pi, and Fst data.
    2. Split the data into subgenome A and B.
    3. Concat mutli chromosome into one subgenome continuous data.
    4. Plot the subgenome continuous XPCLR, pi, Fst data as barplot figure
    "
}


setwd("~/Project/mizi/population_genetics/")
# reads data

sub_genome_A <- c(12, 11, 10, 16, 8, 18, 15, 17, 4)
sub_genome_B <- c(6, 2, 3, 7, 5, 9, 14, 13, 1)

PI_C <- read_tsv("input/plot_input/C.win_20000.windowed.pi")
PI_W <- read_tsv("input/plot_input/W.win_20000.windowed.pi")
PI <- PI_C %>% left_join(PI_W, by = c("CHROM", "BIN_START", "BIN_END")) %>% na.omit()
PI$BIN_MID <- (PI$BIN_START + PI$BIN_END - 1) / 2
PI$ratio <- PI$PI.y / PI$PI.x
PI <- select(PI, CHROM, BIN_MID, ratio)


Fst_C <- read_tsv("input/plot_input/C_vs_W.win_20000.windowed.weir.fst")
Fst_C$BIN_MID <- (Fst_C$BIN_START + Fst_C$BIN_END -11)/2
Fst_C <- select(Fst_C, CHROM, BIN_MID, MEAN_FST)

xpclr <- read_delim(file = paste0("input/plot_input/C_vs_W.win_20000.txt"), delim = " ")
names(xpclr) <- c("CHR", "RANK", "SNP_count", "BP", "CM", "P", "MAXS")
xpclr <- xpclr %>% filter(P!="Inf") %>% select(CHR,BP,P)
xpclr$CHR <- as.numeric(xpclr$CHR)


# subgenome split
PI_sub_genome_A <- tibble(CHROM = sub_genome_A)
PI_sub_genome_A <- inner_join(PI_sub_genome_A, PI, by="CHROM")
PI_sub_genome_B <- tibble(CHROM = sub_genome_B)
PI_sub_genome_B <- inner_join(PI_sub_genome_B, PI, by="CHROM")

Fst_sub_genome_A <- tibble(CHROM = sub_genome_A)
Fst_sub_genome_A <- inner_join(Fst_sub_genome_A, Fst_C, by="CHROM")
Fst_sub_genome_B <- tibble(CHROM = sub_genome_B)
Fst_sub_genome_B <- inner_join(Fst_sub_genome_B, Fst_C, by="CHROM")

xpclr_sub_genome_A <- tibble(CHR = sub_genome_A)
xpclr_sub_genome_A <- inner_join(xpclr_sub_genome_A, xpclr, by="CHR")
xpclr_sub_genome_B <- tibble(CHR = sub_genome_B)
xpclr_sub_genome_B <- inner_join(xpclr_sub_genome_B, xpclr, by="CHR")


# concat mutli chromosome into one single continuous chromosome-like data structure
add_vector <- c()
temp_length <- 0
split_chr_A <- c()
for(i in sub_genome_A){
  temp_mid_vector <- filter(PI_sub_genome_A, CHROM==i)$BIN_MID
  add_vector <- append(add_vector, rep(temp_length, length(temp_mid_vector)))
  temp_length <- temp_mid_vector[length(temp_mid_vector)] + temp_length
  split_chr_A <- append(split_chr_A, temp_length)
}
split_chr_A <- split_chr_A[1:length(split_chr_A)-1]
PI_sub_genome_A$BIN_MID <- PI_sub_genome_A$BIN_MID + add_vector

add_vector <- c()
temp_length <- 0
split_chr_B <- c()
for(i in sub_genome_B){
  temp_mid_vector <- filter(PI_sub_genome_B, CHROM==i)$BIN_MID
  add_vector <- append(add_vector, rep(temp_length, length(temp_mid_vector)))
  temp_length <- temp_mid_vector[length(temp_mid_vector)] + temp_length
  split_chr_B <- append(split_chr_B, temp_length)
}
split_chr_B <- split_chr_B[1:length(split_chr_B)-1]
PI_sub_genome_B$BIN_MID <- PI_sub_genome_B$BIN_MID + add_vector

add_vector_fst <- c()
temp_length <- 0
split_chr_A_fst <- c()
for(i in sub_genome_A){
  temp_mid_vector <- filter(Fst_sub_genome_A, CHROM==i)$BIN_MID
  add_vector_fst <- append(add_vector_fst, rep(temp_length, length(temp_mid_vector)))
  temp_length <- temp_mid_vector[length(temp_mid_vector)] + temp_length
  split_chr_A_fst <- append(split_chr_A_fst, temp_length)
}
split_chr_A_fst <- split_chr_A_fst[1:length(split_chr_A_fst)-1]
Fst_sub_genome_A$BIN_MID <- Fst_sub_genome_A$BIN_MID + add_vector_fst

add_vector_fst <- c()
temp_length <- 0
split_chr_B_fst <- c()
for(i in sub_genome_B){
  temp_mid_vector <- filter(Fst_sub_genome_B, CHROM==i)$BIN_MID
  add_vector_fst <- append(add_vector_fst, rep(temp_length, length(temp_mid_vector)))
  temp_length <- temp_mid_vector[length(temp_mid_vector)] + temp_length
  split_chr_B_fst <- append(split_chr_B_fst, temp_length)
}
split_chr_B_fst <- split_chr_B_fst[1:length(split_chr_B_fst)-1]
Fst_sub_genome_B$BIN_MID <- Fst_sub_genome_B$BIN_MID + add_vector_fst

add_vector_xpclr <- c()
temp_length <- 0
split_chr_A_xpclr <- c()
for(i in sub_genome_A){
  temp_mid_vector <- filter(xpclr_sub_genome_A, CHR==i)$BP
  add_vector_xpclr <- append(add_vector_xpclr, rep(temp_length, length(temp_mid_vector)))
  temp_length <- temp_mid_vector[length(temp_mid_vector)] + temp_length
  split_chr_A_xpclr <- append(split_chr_A_fst, temp_length)
}
split_chr_A_xpclr <- split_chr_A_xpclr[1:length(split_chr_A_xpclr)-1]
xpclr_sub_genome_A$BP <- xpclr_sub_genome_A$BP + add_vector_xpclr

add_vector_xpclr <- c()
temp_length <- 0
split_chr_B_xpclr <- c()
for(i in sub_genome_B){
  temp_mid_vector <- filter(xpclr_sub_genome_B, CHR==i)$BP
  add_vector_xpclr <- append(add_vector_xpclr, rep(temp_length, length(temp_mid_vector)))
  temp_length <- temp_mid_vector[length(temp_mid_vector)] + temp_length
  split_chr_B_xpclr <- append(split_chr_B_fst, temp_length)
}
split_chr_B_xpclr <- split_chr_B_xpclr[1:length(split_chr_B_xpclr)-1]
xpclr_sub_genome_B$BP <- xpclr_sub_genome_B$BP + add_vector_xpclr


# compress data 
PI_A_scaled <- tibble(BIN_MID=apply(matrix(PI_sub_genome_A$BIN_MID, ncol = 130, byrow = TRUE), 1, mean),
                      ratio=apply(matrix(PI_sub_genome_A$ratio, ncol = 130, byrow = TRUE), 1, mean))

PI_B_scaled <- tibble(BIN_MID=apply(matrix(PI_sub_genome_B$BIN_MID, ncol = 135, byrow = TRUE), 1, mean),
                      ratio=apply(matrix(PI_sub_genome_B$ratio, ncol = 135, byrow = TRUE), 1, mean))

Fst_A_scaled <- tibble(BIN_MID=apply(matrix(Fst_sub_genome_A$BIN_MID[1:36660], ncol = 130, byrow = TRUE), 1, mean),
                      MEAN_FST=apply(matrix(Fst_sub_genome_A$MEAN_FST[1:36660], ncol = 130, byrow = TRUE), 1, mean))

Fst_B_scaled <- tibble(BIN_MID=apply(matrix(Fst_sub_genome_B$BIN_MID[1:46845], ncol = 135, byrow = TRUE), 1, mean),
                       MEAN_FST=apply(matrix(Fst_sub_genome_B$MEAN_FST[1:46845], ncol = 135, byrow = TRUE), 1, mean))

xpclr_A_scaled <- tibble(BP=apply(matrix(xpclr_sub_genome_A$BP[1:31950], ncol = 150, byrow = TRUE), 1, mean),
                       P=apply(matrix(xpclr_sub_genome_A$P[1:31950], ncol = 150, byrow = TRUE), 1, mean))

xpclr_B_scaled <- tibble(BP=apply(matrix(xpclr_sub_genome_B$BP[1:41000], ncol = 100, byrow = TRUE), 1, mean),
                         P=apply(matrix(xpclr_sub_genome_B$P[1:41000], ncol = 100, byrow = TRUE), 1, mean))


# add gene annotation
longmi_id <- read_tsv("plot/annotate.ids")
longmi_bed <- read_tsv("longmi.bed")
longmi_bed$mid <- (longmi_bed$start + longmi_bed$end)/2
longmi_bed <- select(longmi_bed, chr, id, mid)
longmi_id <- left_join(longmi_id, longmi_bed, by= c("longmiID"="id"))


gene_sub_genome_A <- tibble(CHROM = sub_genome_A)
gene_sub_genome_A <- inner_join(gene_sub_genome_A, longmi_id, by=c("CHROM" = "chr"))
gene_sub_genome_B <- tibble(CHROM = sub_genome_B)
gene_sub_genome_B <- inner_join(gene_sub_genome_B, longmi_id, by=c("CHROM" = "chr"))

chr_length_A <- tibble(chr=sub_genome_A, length=c(0, split_chr_A_xpclr))
gene_sub_genome_A <- left_join(gene_sub_genome_A, chr_length_A, by = c("CHROM"="chr"))
gene_sub_genome_A$pos <- gene_sub_genome_A$mid + gene_sub_genome_A$length

chr_length_B <- tibble(chr=sub_genome_B, length=c(0, split_chr_B_xpclr))
gene_sub_genome_B <- left_join(gene_sub_genome_B, chr_length_B, by = c("CHROM"="chr"))
gene_sub_genome_B$pos <- gene_sub_genome_B$mid + gene_sub_genome_B$length

# add color for selection sweep
regions <- read_tsv("~/Documents/manuscript/millet/regions/C_vs_W.top5.20kb.regions.txt", col_names = c("chr", "start", "end"))
regions_sub_genome_A <- tibble(CHROM = sub_genome_A)
regions_sub_genome_A <- inner_join(regions_sub_genome_A, regions, by=c("CHROM" = "chr"))
regions_sub_genome_A <- left_join(regions_sub_genome_A, chr_length_A, by = c("CHROM"="chr"))
regions_sub_genome_A$start <- regions_sub_genome_A$start + regions_sub_genome_A$length
regions_sub_genome_A$end <- regions_sub_genome_A$end + regions_sub_genome_A$length

regions_sub_genome_B <- tibble(CHROM = sub_genome_B)
regions_sub_genome_B <- inner_join(regions_sub_genome_B, regions, by=c("CHROM" = "chr"))
regions_sub_genome_B <- left_join(regions_sub_genome_B, chr_length_B, by = c("CHROM"="chr"))
regions_sub_genome_B$start <- regions_sub_genome_B$start + regions_sub_genome_B$length
regions_sub_genome_B$end <- regions_sub_genome_B$end + regions_sub_genome_B$length



# plotting
pi_a_plot <- ggplot(data = PI_A_scaled, aes(x = BIN_MID, y=ratio))+
  geom_bar(stat = "identity", fill="#599B44") +
  geom_vline(xintercept = split_chr_A, linetype = "solid") +
  geom_vline(xintercept = gene_sub_genome_A$pos, linetype = "dashed", size=0.3) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 16)) +
  labs(x = '', y = "Subgenoem A PI ratio", hjust = 0.5) +
  theme_bw() +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), plot.margin = margin(0,10,0,0, "pt"),
        panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8)) +
  annotate("rect", xmin = regions_sub_genome_A$start, xmax = regions_sub_genome_A$end, ymin=-Inf, ymax=Inf,fill='#DF0101', alpha = .3)

pi_b_plot <- ggplot(data = PI_B_scaled, aes(x = BIN_MID, y=ratio))+
  geom_bar(stat = "identity",fill="#599B44") +
  geom_vline(xintercept = split_chr_B, linetype = "solid") +
  geom_vline(xintercept = gene_sub_genome_B$pos, linetype = "dashed", size=0.3) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 16)) +
  labs(x = '', y = "Subgenome B PI ratio", hjust = 0.5) +
  theme_bw() +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), plot.margin = margin(0,10,0,0, "pt"),
        panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8)) +
  annotate("rect", xmin = regions_sub_genome_B$start, xmax = regions_sub_genome_B$end, ymin=-Inf, ymax=Inf,fill='#DF0101', alpha = .3)

fst_a_plot <- ggplot(data = Fst_A_scaled, aes(x = BIN_MID, y=MEAN_FST))+
  geom_bar(stat = "identity",fill="#4E204F") +
  geom_vline(xintercept = split_chr_A_fst, linetype = "solid") +
  geom_vline(xintercept = gene_sub_genome_A$pos, linetype = "dashed", size=0.3) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  labs(x = '', y = "Subgenome A Fst", hjust = 0.5) +
  theme_bw() +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), plot.margin = margin(0,10,0,0, "pt"),
        panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8)) +
  annotate("rect", xmin = regions_sub_genome_A$start, xmax = regions_sub_genome_A$end, ymin=-Inf, ymax=Inf,fill='#DF0101', alpha = .3)

fst_b_plot <- ggplot(data = Fst_B_scaled, aes(x = BIN_MID, y=MEAN_FST))+
  geom_bar(stat = "identity",fill="#4E204F") +
  geom_vline(xintercept = split_chr_B_fst, linetype = "solid") +
  geom_vline(xintercept = gene_sub_genome_B$pos, linetype = "dashed", size=0.3) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  labs(x = '', y = "Subgenome B Fst", hjust = 0.5) +
  theme_bw() +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), plot.margin = margin(0,10,0,0, "pt"),
        panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8))+
  annotate(geom = "text", fontface="bold", x=gene_sub_genome_B$pos, y=0.75, label=gene_sub_genome_B$RiceID, size=2,angle=45) +
  annotate("rect", xmin = regions_sub_genome_B$start, xmax = regions_sub_genome_B$end, ymin=-Inf, ymax=Inf,fill='#DF0101', alpha = .3)

xpclr_a_plot <- ggplot(data = xpclr_A_scaled, aes(x = BP, y=P))+
  geom_bar(stat = "identity",fill="#1973A5") +
  geom_vline(xintercept = split_chr_A_fst, linetype = "solid") +
  geom_vline(xintercept = gene_sub_genome_A$pos, linetype = "dashed", size=0.3) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 40)) +
  labs(x = '', y = "Subgenome B XPCLR", hjust = 0.5) +
  theme_bw() +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), plot.margin = margin(0,10,0,0, "pt"),
        panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8)) +
  annotate(geom = "text", fontface="bold", x=gene_sub_genome_A$pos, y=30, label=gene_sub_genome_A$RiceID, size=2,angle=45) +
  annotate("rect", xmin = regions_sub_genome_A$start, xmax = regions_sub_genome_A$end, ymin=-Inf, ymax=Inf,fill='#DF0101', alpha = .3)

xpclr_b_plot <-  ggplot(data = xpclr_B_scaled, aes(x = BP, y=P))+
  geom_bar(stat = "identity",fill="#1973A5") +
  geom_vline(xintercept = split_chr_B_fst, linetype = "solid") +
  geom_vline(xintercept = gene_sub_genome_B$pos, linetype = "dashed", size=0.3) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 40)) +
  labs(x = '', y = "Subgenome B XPCLR", hjust = 0.5) +
  theme_bw() +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), plot.margin = margin(0,10,0,0, "pt"),
        panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8)) +
  annotate("rect", xmin = regions_sub_genome_B$start, xmax = regions_sub_genome_B$end, ymin=-Inf, ymax=Inf,fill='#DF0101', alpha = .3)


xpclr_a_plot + xpclr_b_plot + 
  pi_a_plot + pi_b_plot + 
  fst_a_plot + fst_b_plot + 
  plot_layout(ncol = 1, heights = c(4, 4, 4, 4, 4, 4))


