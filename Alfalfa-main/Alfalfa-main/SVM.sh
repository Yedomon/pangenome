Markers=read.table("GAPIT.Genotype.Numerical.txt",head=T)
Markers=data.frame(Markers,row.names=1)
data <-as.matrix(read.table(file ="SVMP.txt", header=TRUE,row.name=1))

library(e1071)

traits=1
cycles=5000
for(i in 1:60){
Pheno=data[,i]
Pheno<-as.matrix(Pheno)

accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles){
train= as.matrix(sample(1:176, 140))
test<-setdiff(1:176,train)
svm_fit <- svm(
    x = Markers[train,], 
    y = Pheno[train,], 
    kernel = "linear", 
    cost = 2^(-9),
    gamma = 1, 
    probability = TRUE
  )

svm_pred <- predict(svm_fit, Markers[test, ])

accuracy[r,1] <-cor(svm_pred,Pheno[test, ], use="complete" )}
write.table (accuracy, paste(i,".txt"))
}