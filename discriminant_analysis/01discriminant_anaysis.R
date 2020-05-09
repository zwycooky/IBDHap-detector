test_training <- read.table("test_training.txt",header=T)
test_training$haps <- as.factor(test_training$haps)
library(MASS)
library(randomForest)
library(smotefamily)
test_training <- test_training[,-1]

reverse_training_3 <- test_training[test_training$haps == 3,]
reverse_training_3 <- reverse_training_3[,c(5:8,1:4,9)]
colnames(reverse_training_3) <- colnames(test_training)

reverse_training_4 <- test_training[test_training$haps == 4,]
reverse_training_4 <- reverse_training_4[,c(5:8,1:4,9)]
colnames(reverse_training_4) <- colnames(test_training)

reverse_training_1 <- test_training[test_training$haps == 1,]
reverse_training_1 <- reverse_training_1[,c(5:8,1:4,9)]
colnames(reverse_training_1) <- colnames(test_training)
reverse_training_1$haps <- 2

reverse_training_2 <- test_training[test_training$haps == 2,]
reverse_training_2 <- reverse_training_2[,c(5:8,1:4,9)]
colnames(reverse_training_2) <- colnames(test_training)
reverse_training_2$haps <- 1

test_training <- rbind(test_training,reverse_training_3,reverse_training_4,reverse_training_1,reverse_training_2)
test_training$haps <- as.factor(test_training$haps)

#test_training[sample(which(test_training$haps ==1),1200),3] <- rnorm(1200,mean=0.9,sd=0.02)
#test_training[sample(which(test_training$haps ==2),1200),3] <- rnorm(1200,mean=0.9,sd=0.02)

genData <- BLSMOTE(test_training[,-ncol(test_training)],test_training[,ncol(test_training)], dupSize=5, K=5, C=5)
genData <- BLSMOTE(genData$data[,-ncol(test_training)],genData$data[,ncol(test_training)], dupSize=5, K=5, C=5)

baldata <- genData$data
colnames(baldata)[ncol(test_training)] <- "haps"
baldata$haps <- as.factor(baldata$haps)

test_training <- baldata

# random forest ##
err_vec <- NULL
for (i in 1:(ncol(test_training)-1)) {
	mtry_fit <- randomForest(haps~.,data=test_training,mtry=i,do.trace=T,ntree=500)
	err <- mean(mtry_fit$err.rate)
	err_vec <- c(err_vec,err)
}

my_mtry <- which(err_vec == min(err_vec))

ntree_fit <- randomForest(haps~.,data=test_training,mtry=my_mtry,ntree=1200, do.trace = T)
pdf ("ntree_fit.pdf")
plot(ntree_fit)
dev.off()

RF <- randomForest(haps~., data=test_training, mtry=my_mtry, ntree=800, importance=T, proximity=TRUE, do.trace = T)

save(RF,file="Random_forest_model.Rdata")

