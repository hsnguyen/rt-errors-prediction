library(e1071)
#library(protr)
#library(FSelector)


###################################################
### Data preprocessing ############################
###################################################
datastd <- read.table("../results/2mer/cfs-2mer.txt", header=TRUE)

### Calculate the feature vectors for all amino acid

len <- nrow(datastd)

index <- 1:len
testindex <- sample(index, trunc(length(index)*30/100))
testset <- datastd[testindex,]
trainset <- datastd[-testindex,]

###################################################
### Build predictor with SVM ######################
###################################################
setwd("../results/2mer")
## After using CFS for feature selection
#tuning
bestGamma <- 0.01
bestC <- 10 
def.pred <- mean(trainset[,1])
def.rss <- sum((testset[,1]-def.pred)^2)

print("Running SVM with CFS-reduced dataset...")
model <- svm(error~., data=trainset, kernel="radial", gamma=bestGamma, cost=bestC)
summary(model)

esvm.pred <- predict(model, testset[,-1])
esvm.rss <- sum((testset[,1]-esvm.pred)^2)

print(paste("Root mean square error:", sqrt(esvm.rss/len)))
print(paste("Pseudo R2 value for the predictor:",1.0-esvm.rss/def.rss))
print("=============================================================")

pdf("predCFS.pdf")
plot(testset[,1], esvm.pred, main="Prediction by SVM with CFS-reduced data", ylab="prediction",xlab="observation", xlim=c(-2,2), ylim=c(-2,2))
abline(a=0,b=1,col='red')
dev.off()

