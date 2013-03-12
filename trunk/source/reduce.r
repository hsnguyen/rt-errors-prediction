#!/usr/bin/env Rscript
library(e1071)
library(FSelector)
setwd("../results")
args <- commandArgs(TRUE)
codes <- as.numeric(args)
metList <- c("MoreauBroto","Moran","Geary","CTDC","CTDT","CTDD","CTriad","SOCN","QSO","PAAC","APAAC")
for (i in 1:length(codes)){
	met <- metList[codes[i]]
	setwd(met)	
	###################################################
	### Data preprocessing ############################
	###################################################
	datastd <- read.table(paste(met, "-full.txt", sep=""), header=TRUE)
	error <- datastd$error
	#Feature selector
	t1 <- Sys.time()
	#weightRF <- random.forest.importance(error~., datastd,importance.type=1)
	weightCS <- chi.squared(error~., datastd)
	subsetCS <- cutoff.k(weightCS,20)
	fCS <- as.simple.formula(subsetCS, "error")
	t2 <- Sys.time()

	csdata <- cbind(error,datastd[,subsetCS])
	write.table(csdata, paste(met, "-cs.txt", sep=""), row.name=FALSE)
	
	sink("./log-reduce.txt")	
	cat(paste("Chi-squared reduction time:",t2-t1))
	cat("\n")	
	subsetCFS <- cfs(error~., datastd)
	fCFS <- as.simple.formula(subsetCFS, "error")
	t3 <- Sys.time()
	
	cfsdata <- cbind(error,datastd[,subsetCFS])
	write.table(cfsdata, paste(met, "-cfs.txt", sep=""),row.name=FALSE)

	cat(paste("CFS reduction time:",t3-t2))
	cat("\n")	
	### Calculate the feature vectors for all amino acid

	len <- nrow(datastd)
	
	index <- 1:len
	testindex <- sample(index, trunc(length(index)*30/100))
	testset <- datastd[testindex,]
	trainset <- datastd[-testindex,]
	
	###################################################
	### Build predictor with SVM ######################
	###################################################
	## After using CFS for feature selection
	#tuning
	bestGamma <- 0.01
	bestC <- 10 
	def.pred <- mean(trainset[,1])
	def.rss <- sum((testset[,1]-def.pred)^2)
	
	cat("Running SVM with CFS-reduced dataset...")
	model <- svm(fCFS, data=trainset, kernel="radial", gamma=bestGamma, cost=bestC)
	summary(model)
	
	cfssvm.pred <- predict(model, testset[,-1])
	cfssvm.rss <- sum((testset[,1]-cfssvm.pred)^2)
	
	cat(paste("Root mean square error:", sqrt(cfssvm.rss/len)))
	cat("\n")	
	cat(paste("Pseudo R2 value for the predictor:",1.0-cfssvm.rss/def.rss))
	cat("\n")	

	cat("Running SVM with CS-reduced dataset...")
	model <- svm(fCS, data=trainset, kernel="radial", gamma=bestGamma, cost=bestC)
	summary(model)
	
	cssvm.pred <- predict(model, testset[,-1])
	cssvm.rss <- sum((testset[,1]-cssvm.pred)^2)
	
	cat(paste("Root mean square error:", sqrt(cssvm.rss/len)))
	cat("\n")	
	cat(paste("Pseudo R2 value for the predictor:",1.0-cssvm.rss/def.rss))
	cat("\n")	
	sink()
	
	pdf("predCS.pdf")
	plot(testset[,1], cssvm.pred, main="Prediction by SVM with CS-reduced data", ylab="prediction",xlab="observation", xlim=c(-2,2), ylim=c(-2,2))
	abline(a=0,b=1,col='red')
	dev.off()
	
	pdf("predCFS.pdf")
	plot(testset[,1], cfssvm.pred, main="Prediction by SVM with CFS-reduced data", ylab="prediction",xlab="observation", xlim=c(-2,2), ylim=c(-2,2))
	abline(a=0,b=1,col='red')
	dev.off()

	setwd("../")
}
