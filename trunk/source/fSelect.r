library(FSelector)

#Feature selector
t1 <- Sys.time()
#weightRF <- random.forest.importance(error~., datastd,importance.type=1)
weightCS <- chi.squared(error~., datastd)
subsetCS <- cutoff.k(weightCS,20)
fCS <- as.simple.formula(subsetCS, "error")
t2 <- Sys.time()

csdata <- cbind(error,datastd[,subsetCS])
write.table(csdata,"cs-.txt",row.name=FALSE)

print(paste("Chi-squared reduction time:",t2-t1))

subsetCFS <- cfs(error~., datastd)
fCFS <- as.simple.formula(subsetCFS, "error")
t3 <- Sys.time()

cfsdata <- cbind(error,datastd[,subsetCFS])
write.table(cfsdata,"cfs-.txt",row.name=FALSE)

print(paste("CFS reduction time:",t3-t2))
