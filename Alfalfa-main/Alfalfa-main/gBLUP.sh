source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
myG <- read.table("allSNP.hmp.txt", head = FALSE)
myGAPIT <- GAPIT(G=myG, output.numerical=TRUE,file.output =FALSE)
myGD= myGAPIT$GD 
myGM= myGAPIT$GM

myCV=read.table("PCA.txt",head=T)

data <- as.matrix(read.table(file="plinkg.txt", header=F))

repeat_number=5000

for(i in 1:60){
t1<-data[,1]
t2<-data[,i+1]
myY<- as.matrix(cbind(t1, t2))

acc_testing_all=rep(0,repeat_number)
for (r in 1:repeat_number){
training= as.matrix(sample(1:176,140))
testing<-setdiff(1:176,training)
  
myGAPIT5<- GAPIT(
  Y=myY[training,],
  GD=myGD,
  GM=myGM,
  PCA.total=5,
  CV=myCV,
  model="gBLUP",
  SNP.test=FALSE,
  file.output=F
  )
order=match(myY[,1],myGAPIT5$Pred[,1])
myPred=myGAPIT5$Pred[order,]
ry2.blup=cor(myPred[testing,5],myY[testing,2],use = "complete")

acc_testing_all[r]=ry2.blup} 
write.table (acc_testing_all, paste(i,".txt"))
}