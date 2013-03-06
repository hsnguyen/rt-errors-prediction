library(e1071)
library(protr)
#library(kernlab)
library(FSelector)
### Error threshold: 10%
errThres <- 0.1

###################################################
### Feature extraction function ###################
###################################################
featureExt <- function(aa)
{
    #Get rid of amino acides at two end
    seq <- strsplit(aa,"\\.")[[1]][2]
    return(extractAAC(seq))
}

hydrSeq <- function(aa)
{
    seq <- strsplit(aa,"\\.")[[1]][2]
    base <- "RHKDESTNQAVILMFYWCGP"
    propr <- "CCCCCPPPPHHHHHHHHOOO"
    return(chartr(base,propr,seq))
}
###################################################
### Data preprocessing ############################
###################################################
dataset <- read.table("../data/rt-predictions.txt", header=TRUE)

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
                                        return ((as.numeric(row['Observed_RT'])-as.numeric(row['Predicted_RT']))/as.numeric(row['Observed_RT']))
                                    })

### Calculate the feature vectors for all amino acid

len <- nrow(dataset)
features <- numeric()
#properties <- character()
for (i in 1:len){
    aa <- as.character(dataset$Peptide[i])
    features <- rbind(features,featureExt(aa))
#    properties <- rbind(properties, hydrSeq(aa))
}
#prop <- data.frame(aa=dataset$Peptide, properties=properties, error=error)
#write.table(prop,"../results/aa-properties.txt",row.name=FALSE)

datastd <- data.frame(error=error)
datastd <- cbind(datastd,features)
### Separate to training test and test set
# (not necessary if using cross-validation)
index <- 1:len
testindex <- sample(index, trunc(length(index)*30/100))
testset <- datastd[testindex,]
trainset <- datastd[-testindex,]

###################################################
### Build predictor with SVM ######################
###################################################
# Feature selector
weightRF <- random.forest.importance(error~., datastd,importance.type=1)
subsetRF <- cutoff.k(weightRF,10)
fRF <- as.simple.formula(subsetRF, "error")
subsetCFS <- cfs(error~., datastd)
fCFS <- as.simple.formula(subsetCFS, "error")


#tuning
#tobj <- tune.svm(error~., data=datastd, gamma=10^(-6:-1), cost= 10^(-1:1))
#bestGamma <- tobj$best.parameters[[1]]
#bestC <- tobj$best.parameters[[2]]
print("Running SVM with original dataset...")
bestGamma <- 0.01
bestC <- 10
def.pred <- mean(trainset[,1])
def.rss <- sum((testset[,1]-def.pred)^2)

model <- svm(error~., data=trainset, kernel="radial", gamma=bestGamma, cost=bestC)
summary(model)

esvm.pred <- predict(model, testset[,-1])
esvm.rss <- sum((testset[,1]-esvm.pred)^2)

print(paste("Root mean square error:", sqrt(esvm.rss/len)))
print(paste("Pseudo R2 value for the predictor:",1.0-esvm.rss/def.rss))
print("=============================================================")

pdf("../results/pred.pdf")
plot(testset[,1], esvm.pred, main="Prediction by SVM", ylab="prediction",xlab="observation", xlim=c(-2,2), ylim=c(-2,2))
abline(a=0,b=1,col='red')
dev.off()
## After using RF for feature selection
#tuning
tobj <- tune.svm(fRF, data=datastd, gamma=10^(-6:-1), cost= 10^(-1:1))
bestGamma <- tobj$best.parameters[[1]]
bestC <- tobj$best.parameters[[2]]
print("Running SVM with RF-reduced dataset...")
model <- svm(fRF, data=trainset, kernel="radial", gamma=bestGamma, cost=bestC)
summary(model)

esvm.pred <- predict(model, testset[,-1])
esvm.rss <- sum((testset[,1]-esvm.pred)^2)

print(paste("Root mean square error:", sqrt(esvm.rss/len)))
print(paste("Pseudo R2 value for the predictor:",1.0-esvm.rss/def.rss))
print("=============================================================")

pdf("../results/predRF.pdf")
plot(testset[,1], esvm.pred, main="Prediction by SVM with RF-reduced data", ylab="prediction",xlab="observation", xlim=c(-2,2), ylim=c(-2,2))
abline(a=0,b=1,col='red')
dev.off()
## After using CSF for feature selection
#tuning
tobj <- tune.svm(fCSF, data=datastd, gamma=10^(-6:-1), cost= 10^(-1:1))
bestGamma <- tobj$best.parameters[[1]]
bestC <- tobj$best.parameters[[2]]
print("Running SVM with CSF-reduced dataset...")
model <- svm(fCSF, data=trainset, kernel="radial", gamma=bestGamma, cost=bestC)
summary(model)

esvm.pred <- predict(model, testset[,-1])
esvm.rss <- sum((testset[,1]-esvm.pred)^2)

print(paste("Root mean square error:", sqrt(esvm.rss/len)))
print(paste("Pseudo R2 value for the predictor:",1.0-esvm.rss/def.rss))
print("=============================================================")

pdf("../results/predCSF.pdf")
plot(testset[,1], esvm.pred, main="Prediction by SVM with CSF-reduced data", ylab="prediction",xlab="observation", xlim=c(-2,2), ylim=c(-2,2))
abline(a=0,b=1,col='red')
dev.off()


###########################################################################
#model <- ksvm(label~., data=trainset, type="C-bsvc", kernel="rbfdot", kpar=list(sigma=0.1),C=1)

# built-in 10-fold CV estimate of prediction error
#spread <- rep(0,20)
#for (i in 1:20) {
#    mysvm <- svm(y ~ x,data,cross=10)
#    spread[i] <- mean(mysvm$tot.MSE)
#    }
