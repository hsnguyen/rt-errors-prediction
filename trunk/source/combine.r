#!/usr/bin/env Rscript
library(e1071)

### Error threshold: 10%
errThres <- 0.1
#minaa <- 6 #min length of aa

args <- commandArgs(trailingOnly=TRUE)
codes <- as.numeric(args)

metList <- c("MoreauBroto","Moran","Geary","CTDC","CTDT","CTDD","CTriad","SOCN","QSO","PAAC","APAAC")
comb <- numeric()
cdir <- ""

setwd("../results")
for(i in 1:length(codes)){
    met <- metList[codes[i]]
    #print(met) 
    setwd(met)
    if (codes[i]==7)
    	data <- read.table(paste(met, "-cfs.txt", sep=""), header=TRUE)	
    else
	data <- read.table(paste(met, "-full.txt", sep=""), header=TRUE)	
    if(i == 1){
	comb <- data
	cdir <- met
    }
    else{
        comb <- cbind(comb, data[,-1])
    	cdir <- paste(cdir, met, sep="-")
    } 
    setwd("../")
}

dir.create(cdir)
setwd(cdir)
###################################################
### Data preprocessing ############################
###################################################

### Calculate the feature vectors for all amino acid

len <- nrow(comb)

fileName <- paste(cdir,"-full.txt",sep="")
write.table(comb, fileName, row.name=FALSE)
### Separate to training test and test set
# (not necessary if using cross-validation)
index <- 1:len
testindex <- sample(index, trunc(length(index)*30/100))
testset <- comb[testindex,]
trainset <- comb[-testindex,]

###################################################
### Build predictor with SVM ######################
###################################################
#write.table(datastd,"2mer-full.txt",row.name=FALSE)

#tuning
#tobj <- tune.svm(error~., data=datastd, gamma=10^(-6:-1), cost= 10^(-1:1))
#bestGamma <- tobj$best.parameters[[1]]
#bestC <- tobj$best.parameters[[2]]
sink("./log.txt")
cat(paste("Running SVM with original dataset of feature ", cdir, sep=""))
bestGamma <- 0.01
bestC <- 10
def.pred <- mean(trainset[,1])
def.rss <- sum((testset[,1]-def.pred)^2)

model <- svm(error~., data=trainset, kernel="radial", gamma=bestGamma, cost=bestC)
summary(model)

esvm.pred <- predict(model, testset[,-1])
esvm.rss <- sum((testset[,1]-esvm.pred)^2)

cat(paste("Root mean square error:", sqrt(esvm.rss/len)))
cat("\n")
cat(paste("Pseudo R2 value for the predictor:",1.0-esvm.rss/def.rss))
sink()

pdf("pred.pdf")
plot(testset[,1], esvm.pred, main=cdir, ylab="prediction",xlab="observation", xlim=c(-2,2), ylim=c(-2,2))
abline(a=0,b=1,col='red')
dev.off()


setwd("../../source")
