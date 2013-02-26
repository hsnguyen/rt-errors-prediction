library(e1071)
library(protr)

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

###################################################
### Data preprocessing ############################
###################################################
dataset <- read.table("../data/rt-predictions.txt", header=TRUE)

### Assign label to the data: -1 for large error rate and 1 for the others
# used for later classification model
label <- apply(dataset, 1, function(row){
                                        errorRate <- (as.numeric(row['Observed_RT'])-as.numeric(row['Predicted_RT']))/as.numeric(row['Observed_RT'])
                                        if (abs(errorRate) < errThres)
                                            return(1)
                                        else
                                            return(-1)
                                    })
### used for regression model
#error <- apply(dataset, 1, function(row){
#                                        return (as.numeric(row['Observed_RT'])-as.numeric(row['Predicted_RT']))/as.numeric(row['Observed_RT'])
#                                    })

### Calculate the feature vectors for all amino acid

len <- nrow(dataset)
features <- numeric()
for (i in 1:len){
    aa <- as.character(dataset$Peptide[i])
    features <- rbind(features,featureExt(aa))
}
datastd <- data.frame(aa=dataset$Peptide,label=label)
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

#tuning
#tobj <- tune.svm(label, y=features, gamma=10^(-6:-1), cost= 10^(-1:1))
#bestGamma <- tobj$best.parameters[[1]]
#bestC <- tobj$best.parameters[[2]]

#model <- svm(label, y=features, kernel="radial")
#summary(model)
model <- svm(label~., data=trainset, kernel="radial")
pred <- predict(model, testset[,-2])
(acc <- table(pred,testset[,2]))
classAgreement(acc)

