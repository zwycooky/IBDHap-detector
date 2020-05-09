test_training <- read.table("test_training.txt",header=T)
test_training$haps <- as.factor(test_training$haps)
library(MASS)
library(randomForest)

test_training <- test_training[,-1]

compare_mat <- matrix(0,ncol=2,nrow=20)
for (repeat_test in 1:20) {
	index <- sample(2,nrow(test_training),replace = T,prob=c(0.5,0.5))
	train <- test_training[index==1,]
	test <- test_training[index==2,]

	train[sample(which(train$haps ==1),300),3] <- rnorm(300,mean=0.92,sd=0.01)
	train[sample(which(train$haps ==2),300),3] <- rnorm(300,mean=0.92,sd=0.01)

	# random forest ##
	err_vec <- NULL
	for (i in 1:(ncol(train)-1)) {
		mtry_fit <- randomForest(haps~.,data=train,mtry=i,do.trace=F,ntree=600)
		err <- mean(mtry_fit$err.rate)
		err_vec <- c(err_vec,err)
	}

	my_mtry <- which(err_vec == min(err_vec))

	#ntree_fit <- randomForest(haps~.,data=train,mtry=my_mtry,ntree=1200, do.trace = T)
	#plot(ntree_fit)

	RF <- randomForest(haps~., data=train, mtry=my_mtry, ntree=800, importance=T, proximity=T, do.trace = F)

	pred <- table(predict(RF,newdata=test),test$haps)
	compare_mat[repeat_test,1] <- sum(diag(pred)) / sum(pred)

	error_dat <- cbind(test[predict(RF,newdata=test) != test$haps,], predict(RF,newdata=test)[predict(RF,newdata=test) != test$haps])
	colnames(error_dat)[ncol(error_dat)] <- "predict"
	
##
	library(smotefamily)

	genData <- BLSMOTE(train[,-ncol(train)],train[,ncol(train)],dupSize=20, K=5, C=5)
	genData <- BLSMOTE(genData$data[,-ncol(train)],genData$data[,ncol(train)],dupSize=20 ,K=5, C=5)

	baldata <- genData$data
	colnames(baldata)[ncol(train)] <- "haps"
	baldata$haps <- as.factor(baldata$haps)

	err_vec <- NULL
	for (i in 1:(ncol(baldata)-1)) {
		mtry_fit <- randomForest(haps~.,data=baldata,mtry=i,do.trace=F,ntree=500)
		err <- mean(mtry_fit$err.rate)
		err_vec <- c(err_vec,err)
	}

	my_mtry <- which(err_vec == min(err_vec))

#	ntree_fit <- randomForest(haps~.,data=baldata,mtry=my_mtry,ntree=800, do.trace = T)
#	plot(ntree_fit)

	RF <- randomForest(haps~., data=baldata, mtry=my_mtry, ntree=800, importance=F, proximity=T, do.trace = F)

	pred <- table(predict(RF,newdata=test),test$haps)
	compare_mat[repeat_test,2] <- sum(diag(pred)) / sum(pred)

	error_dat <- cbind(test[predict(RF,newdata=test) != test$haps,], predict(RF,newdata=test)[predict(RF,newdata=test) != test$haps])
	colnames(error_dat)[ncol(error_dat)] <- "predict"
}

write.table(compare_mat,"compare_balance_unbalance.res.txt",row.names=F,col.names=F,sep="\t",quote=F)
