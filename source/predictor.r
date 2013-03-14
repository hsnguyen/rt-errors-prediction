#!/usr/bin/env Rscript
library(e1071)
args <- commandArgs(TRUE)
file <- as.character(args)
setwd("../results")
dir <- strsplit(file, "\\.")[[1]][1]
print(dir)
comb <- read.table(file, header=TRUE)
comb$error <- abs(comb$error)
dir.create(dir)
setwd(dir)
###################################################
### Data preprocessing ############################
###################################################

### Calculate the feature vectors for all amino acid
### Separate to training test and test set
# (not necessary if using cross-validation)
len <- nrow(comb)
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
cat(paste("Running SVM with original dataset of feature ", dir, sep=""))
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
plot(testset[,1], esvm.pred, main=dir, ylab="prediction",xlab="observation", xlim=c(-2,2), ylim=c(-2,2))
abline(a=0,b=1,col='red')
dev.off()


