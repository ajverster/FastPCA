
source("../Lib/FastPCA.R")
pr.use <- 50

Data <- read.table("train.csv", sep = ",", header = T)
Data.X <- as.matrix(Data[,-1])
Data.X.pca <- FastPCA(Data.X, top.k = pr.use)
summary(Data.X.pca,Data.X)

