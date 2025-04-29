library("svglite")
library("ggsignif")
library("ggplot2")
#svglite(file="boxplot_root_all.svg", width=5, height=5);
pdf(file="boxplot_root_all.pdf", width=5, height=5);
morphData <- read.table("tpm_table_leaf.txt", sep="\t", header=TRUE);
ggplot(morphData, aes(x=Set, y=TPM))+ ylim(0,20)+
geom_boxplot() +
geom_signif(comparisons = list(c("KK(UF)", "KK(F)"), c("LL(UF)", "LL(F)"), c("Os(PH)", "Os(SCH)")),map_signif_level = function(p) sprintf("p = %.2g", p))
dev.off()


library("svglite")
library("ggsignif")
library("ggplot2")
#svglite(file="boxplot_leaf_all.svg", width=5, height=5);
pdf(file="boxplot_leaf_all.pdf", width=5, height=5);
morphData <- read.table("tpm_table_root.txt", sep="\t", header=TRUE);
ggplot(morphData, aes(x=Set, y=TPM))+ ylim(0,20)+
geom_boxplot() + 
geom_signif(comparisons = list(c("KK(UF)", "KK(F)"), c("LL(UF)", "LL(F)"), c("Os(PH)", "Os(SCH)")),map_signif_level = function(p) sprintf("p = %.2g", p))
dev.off()

