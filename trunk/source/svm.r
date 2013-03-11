#!/usr/bin/env Rscript
library(e1071)
library(protr)

### Error threshold: 10%
errThres <- 0.1
#minaa <- 6 #min length of aa

args <- commandArgs(trailingOnly=TRUE)
code <- as.numeric(args[1])

metList <- c("MoreauBroto","Moran","Geary","CTDC","CTDT","CTDD","CTriad","SOCN","QSO","PAAC","APAAC")
funcList <- paste("extract", metList, sep="")

###################################################
### Feature extraction function ###################
###################################################
featureExt <- function(aa)
{
    #Get rid of amino acides at two end
    seq <- strsplit(aa,"\\.")[[1]][2]
    #do.call(funcList[code],list(seq))
    if(code < 4)
      	do.call(funcList[code],list(seq, nlag=5))
    else if(code < 8)
	do.call(funcList[code],list(seq))
    else if (code < 10)
	do.call(funcList[code],list(seq,nlag=5))
    else
	do.call(funcList[code],list(seq,lambda=5))

}

###################################################
### Data preprocessing ############################
###################################################

dataset <- read.table("../data/rt-predictions.txt", header=TRUE)

wdir <- paste("../results/", metList[code], sep="")
dir.create(wdir)
setwd(wdir)
### Assign label to the data: -1 for large error rate and 1 for the others
# used for later classification model
#label <- apply(dataset, 1, function(row){
#                                        errorRate <- (as.numeric(row['Observed_RT'])-as.numeric(row['Predicted_RT']))/as.numeric(row['Observed_RT'])
#                                        if (abs(errorRate) < errThres)
#                                            return(1)
#                                        else
#                                            return(-1)
#                                    })
### used for regression model
error <- apply(dataset, 1, function(row){
                                        (as.numeric(row['Observed_RT'])-as.numeric(row['Predicted_RT']))/as.numeric(row['Observed_RT'])
                                    })

### Calculate the feature vectors for all amino acid

len <- nrow(dataset)
features <- numeric()
for (i in 1:len){
    aa <- as.character(dataset$Peptide[i])
    features <- rbind(features,featureExt(aa))
}
datastd <- data.frame(error=error)
datastd <- cbind(datastd,features)

fileName <- paste(metList[code],"-full.txt",sep="")
write.table(datastd, fileName, row.name=FALSE)
### Separate to training test and test set
# (not necessary if using cross-validation)
index <- 1:len
testindex <- sample(index, trunc(length(index)*30/100))
testset <- datastd[testindex,]
trainset <- datastd[-testindex,]

###################################################
### Build predictor with SVM ######################
###################################################
#write.table(datastd,"2mer-full.txt",row.name=FALSE)

#tuning
#tobj <- tune.svm(error~., data=datastd, gamma=10^(-6:-1), cost= 10^(-1:1))
#bestGamma <- tobj$best.parameters[[1]]
#bestC <- tobj$best.parameters[[2]]
sink("./log.txt")
cat(paste("Running SVM with original dataset of feature ",metList[code],sep=""))
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
plot(testset[,1], esvm.pred, main=metList[code], ylab="prediction",xlab="observation", xlim=c(-2,2), ylim=c(-2,2))
abline(a=0,b=1,col='red')
dev.off()


setwd("../../source")
