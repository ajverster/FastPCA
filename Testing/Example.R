#!/usr/bin/Rscript

############
#Libraries
############
source("../Lib/FastPCA.R")

############
#Paramters
############
pr.use <- 50
outfile <- "TimeCompare.png"
plots.w <- 2000
plots.h <- 1600

############
#Input
############
Data <- read.table("train.csv", sep = ",", header = T)
Data.X <- as.matrix(Data[,-1])


############
#Computation
############
times <- vector(mode ="numeric",length = 4)
times[1] <- system.time(FastSVD(Data.X,pr.use,iterative = F))[3]
times[2] <- system.time(svd(Data.X))[3]
times[3] <- system.time(FastPCA(Data.X, top.k = pr.use))[3]
times[4] <- system.time(prcomp(Data.X))[3]
names(times) <- c("FastSVD","SVD","FastPCA","PCA")


############
#Output
############
png(outfile, width = plots.w, height = plots.h)
barplot(times,ylab = "Time(sec)")
dev.off()