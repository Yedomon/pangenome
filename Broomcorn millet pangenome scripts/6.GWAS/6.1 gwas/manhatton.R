#install.packages("CMplot")
library(CMplot)
file_name = list.files("for_manhatton")  
dir = paste("./for_manhatton/",file_name,sep="")
n = length(dir)  
for(i in 1:n){         
        #print(dir[i],i)
data<-read.table(dir[i],header = T)
CMplot(data,plot.type = c("m","q"),
       threshold = 0.05/nrow(data),
       threshold.col=c('grey','black'),
       threshold.lty = c(1,2),threshold.lwd = c(1,1), amplify = F,
       main="",
       cex=0.8,
       cex.axis=1.3,
       cir.chr=TRUE,
       dpi=300,
       cex.lab = 1.8,
       main.font = 2,
       file="pdf",
       memo=file_name[i]) 
}

