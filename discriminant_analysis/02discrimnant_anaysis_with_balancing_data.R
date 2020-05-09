test_training <- read.table("test_training.txt",header=T)
test_training$haps <- as.factor(test_training$haps)
library(MASS)
library(klaR)
library(randomForest)

test_training <- test_training[,-1]

library(smotefamily)

genData <- BLSMOTE(test_training[,-9],test_training[,9],dupSize=30, K=5, C=5)
genData <- BLSMOTE(genData$data[,-9],genData$data[,9],dupSize=30 ,K=5, C=5)

baldata <- genData$data
colnames(baldata)[9] <- "haps"
baldata$haps <- as.factor(baldata$haps)

set.seed(234)

RF <- randomForest(haps~., data=baldata, mtry=6, ntree=600, importance=T, proximity=TRUE, do.trace = T)
save(RF,file="Random_forest_model.Rdata")


