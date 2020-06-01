options(stringsAsFactors=F)
library(HMM)
library(randomForest)

# load Random Forest model #
load("discriminant_analysis/Random_forest_model.Rdata")
# load 'RF' DONE #

haps <- read.table("correct.haps.filtered.formated.txt",header=F)
snps <- read.table("formated_w221_snp.filtered.txt",header=T)
snps <- snps[haps[,2],]

linkage_files <- paste("../linkage_group/","chr",1:15,".map.ordered.txt",sep="")

fuding_phasing <- function(haps, snps, RF) {
	library(HMM)
	library(randomForest)
	
	chromosome_num <- haps[1,5]

	determine_haps <- function(subhap) {
		sum1 <- sum(subhap == 1)
		sum2 <- sum(subhap == 2)
		if (sum1 < sum2) {
			return(2)
		}else{
			return(1)
		}
	}
	
	hmm_correct <- function(fuding_hap,hmm) {
		genoSymbol <- ifelse(fuding_hap==1,hmm$Symbols[1],hmm$Symbols[2])
		correctGeno <- viterbi(hmm, genoSymbol)
		correctGeno <- ifelse(correctGeno==hmm$States[1], 1, 2)
		names(correctGeno) <- names(fuding_hap)
	
		lft <- which( abs(correctGeno[-1] - correctGeno[-length(correctGeno)]) == 1 )
		rht <- as.numeric(names(lft))
		lft <- as.numeric(names(fuding_hap[lft]))
		res <- cbind(lft,rht)
	
		if (nrow(res) == 0) {
			final_hap <- determine_haps(fuding_hap)
			ss <- as.numeric(names(fuding_hap)[1])
			ee <- as.numeric(names(fuding_hap)[length(fuding_hap)])
			#hapreg <- cbind(haps[,1][ss], haps[,1][ee], final_hap)
			
			# het rate #
			bin_select_reg <- haps[,6] >= ss & haps[,6] <= ee
			bin_hetsnp <- sum(nchar(snps[bin_select_reg]) == 2, na.rm=T)
			bin_homsnp <- sum(nchar(snps[bin_select_reg]) == 1, na.rm=T)
			het_rate <- round(bin_hetsnp/(bin_hetsnp+bin_homsnp),3)			
			
			# correct rate #
			correctnum <- sum(correctGeno == fuding_hap)
			allnum <- length(fuding_hap)
			correct_rate <- round(correctnum/allnum,3)
			
			hapreg <- cbind(ss, ee, final_hap, as.numeric(haps[1,5]), het_rate, correct_rate, allnum)
			#return(list(hapreg=hapreg, correcthaps = correctGeno))
			return(hapreg)
		}else{
			bin_vec <- rep(0,nrow(res)*ncol(res))
			bin_vec[c(T,F)] <- res[,1]
			bin_vec[c(F,T)] <- res[,2]
			bin_vec <- as.numeric( c(names(fuding_hap)[1], bin_vec, names(fuding_hap)[length(fuding_hap)]) )
			bin_mat <- matrix(bin_vec,ncol=2,byrow=T)
			hapreg <- matrix(0,nrow=nrow(bin_mat),ncol=7)
			for (b in 1:nrow(bin_mat)) {
				tmp_s <- which(as.numeric(names(fuding_hap)) == bin_mat[b,1])
				tmp_e <- which(as.numeric(names(fuding_hap)) == bin_mat[b,2])
				final_hap <- determine_haps(fuding_hap[tmp_s:tmp_e])
							
				# correct rate #
				correctnum <- sum(correctGeno[tmp_s:tmp_e] == fuding_hap[tmp_s:tmp_e])
				allnum <- length(tmp_s:tmp_e)
				correct_rate <- round(correctnum/allnum,3)
				
				# het rate #
				bin_select_reg <- haps[,6] >= bin_mat[b,1] & haps[,6] <= bin_mat[b,2]
				bin_hetsnp <- sum(nchar(snps[bin_select_reg]) == 2, na.rm=T)
				bin_homsnp <- sum(nchar(snps[bin_select_reg]) == 1, na.rm=T)
				het_rate <- round(bin_hetsnp/(bin_hetsnp+bin_homsnp),3)	
					
				hapreg[b,] <- c(bin_mat[b,1], bin_mat[b,2], final_hap, as.numeric(haps[1,5]), het_rate, correct_rate, allnum)
			}
		
			#return(list(hapreg=hapreg, correcthaps = correctGeno))
			return(hapreg)
		}
	
	}


	match_hap1 <- snps == haps[,3]
	match_hap2 <- snps == haps[,4]

	fuding_hap <- rep(NA,length(match_hap1))

	fuding_hap[which(match_hap1 == TRUE)] <- 1
	fuding_hap[which(match_hap2 == TRUE)] <- 2
	names(fuding_hap) <- haps[,6]
	##
	fuding_hap_phu1 <- fuding_hap
	fuding_hap_phu2 <- fuding_hap
	fuding_hap_phu1[which(nchar(snps) == 2)] <- 1
	fuding_hap_phu2[which(nchar(snps) == 2)] <- 2
	
	#fuding_hap[is.na(match_hap1)] <- NA
	#names(fuding_hap) <- 1:length(fuding_hap)
	
	
	hmm = initHMM(States=c("F","M"), 
        Symbols=c("f","m"), 
        transProbs=matrix(c(0.99999,0.00001,0.00001,0.99999),2),
        emissionProbs=matrix(c(0.99,0.01,0.01,0.99),2), 
        startProbs = c(0.5,0.5))
	
	
	fuding_hap_phu1 <- fuding_hap_phu1[!is.na(fuding_hap_phu1)]
	fuding_hap_phu2 <- fuding_hap_phu2[!is.na(fuding_hap_phu2)]
	
	fuding_hap_phu1_reg <- hmm_correct(fuding_hap_phu1,hmm)
	fuding_hap_phu2_reg <- hmm_correct(fuding_hap_phu2,hmm)

	fir <- 0
	window_res <- NULL
	compare_mat <- matrix(sort(c(fuding_hap_phu1_reg[,1][-1],fuding_hap_phu1_reg[,2][-nrow(fuding_hap_phu1_reg)],fuding_hap_phu2_reg[,1],fuding_hap_phu2_reg[,2])),ncol=2,byrow=T)
	for(i in 1:nrow(compare_mat)) {
		ws <- compare_mat[i,1]
		we <- compare_mat[i,2]
		
		tmpreg1 <- rbind(NULL,fuding_hap_phu1_reg[fuding_hap_phu1_reg[,1] <= ws & fuding_hap_phu1_reg[,2] >= we,])
		tmpreg2 <- rbind(NULL,fuding_hap_phu2_reg[fuding_hap_phu2_reg[,1] <= ws & fuding_hap_phu2_reg[,2] >= we,])
		
		if (fir == 0) {
			
			tmplen1 <- tmpreg1[1,2] - tmpreg1[1,1]
			tmplen2 <- tmpreg2[1,2] - tmpreg2[1,1]
				
			# calculating overlap rate #
			minreg <- min(tmpreg1[1,2], tmpreg1[1,1], tmpreg2[1,2], tmpreg2[1,1])
			maxreg <- max(tmpreg1[1,2], tmpreg1[1,1], tmpreg2[1,2], tmpreg2[1,1])
			overlap_len <- (tmplen1 + tmplen2) - (maxreg - minreg)
			over_rate1 <- round(overlap_len / tmplen1,3)
			over_rate2 <- round(overlap_len / tmplen2,3)
				
			# generate test model for RF #
			testdata <- rbind(NULL,c(tmplen1, tmpreg1[1,5:6], over_rate1, tmplen2, tmpreg2[1,5:6], over_rate2))
			#len1	het_rate1	correct_rate1	overlap_rate1	len2	het_rate2	correct_rate2	overlap_rate2	haps
			colnames(testdata) <- c("len1","het_rate1","correct_rate1","overlap_rate1","len2","het_rate2","correct_rate2","overlap_rate2")
				
				
			# randomForest discriminant analysis #
			pred <- predict(RF,newdata=testdata)
			ibs <- as.numeric(as.character(pred))
				
			if ( (ibs == 3 | ibs == 4) &  tmpreg1[1,3] == tmpreg2[1,3]) {
				ibs <- as.numeric(paste(rep(tmpreg1[1,3],2), collapse=""))
			}
				
			ss <- ws
			
			fir <- 1
		}else{
			
			if (nrow(tmpreg1) > 0 & nrow(tmpreg2) > 0) {
				tmplen1 <- tmpreg1[1,2] - tmpreg1[1,1]
				tmplen2 <- tmpreg2[1,2] - tmpreg2[1,1]
				
				# calculating overlap rate #
				minreg <- min(tmpreg1[1,2], tmpreg1[1,1], tmpreg2[1,2], tmpreg2[1,1])
				maxreg <- max(tmpreg1[1,2], tmpreg1[1,1], tmpreg2[1,2], tmpreg2[1,1])
				overlap_len <- (tmplen1 + tmplen2) - (maxreg - minreg)
				over_rate1 <- round(overlap_len / tmplen1,3)
				over_rate2 <- round(overlap_len / tmplen2,3)
				
				# generate test model for RF #
				testdata <- rbind(NULL,c(tmplen1, tmpreg1[1,5:6], over_rate1, tmplen2, tmpreg2[1,5:6], over_rate2))
				#len1	het_rate1	correct_rate1	overlap_rate1	len2	het_rate2	correct_rate2	overlap_rate2	haps
				colnames(testdata) <- c("len1","het_rate1","correct_rate1","overlap_rate1","len2","het_rate2","correct_rate2","overlap_rate2")
				
				# randomForest discriminant analysis #
				pred <- predict(RF,newdata=testdata)
				nextibs <- as.numeric(as.character(pred))
				
				if ( (nextibs == 3 | nextibs == 4) &  tmpreg1[1,3] == tmpreg2[1,3]) {
					nextibs <- as.numeric(paste(rep(tmpreg1[1,3],2), collapse=""))
				}
				
			}else if (nrow(tmpreg1) == 0 & nrow(tmpreg2) == 0) {
				nextibs <- 0
			}else{
				next
			}
			
			if (nextibs != ibs) {
				ee <- ws
				window_res <- rbind(window_res,c(ss,ee,ibs))
				ss <- ws
				ibs <- nextibs
			}
			
		}
	}
	
	ee <- we
	window_res <- rbind(window_res,c(ss,ee,ibs))
	window_res <- cbind(window_res,rep(chromosome_num,nrow(window_res)),window_res[,2]-window_res[,1]+1)
	colnames(window_res) <- NULL
	rownames(window_res) <- NULL
	return(window_res)
	#start	end	haps	chr	length	accession#
}

library(pbapply)
#library(parallel)
#cl <- makeCluster(1)

fuding_hap_res <- NULL
for (i in 1:length(linkage_files)) {
	tmpgmap <- read.table(linkage_files[i])
	tmphaps <- haps[haps[,1] %in% tmpgmap[,1],]
	tmpsnps <- snps[tmphaps[,2],]
	
#	tmp_res <- pbapply(tmpsnps, 2, fuding_phasing, haps = tmphaps, RF = RF, cl = cl)
	tmp_res <- pbapply(tmpsnps, 2, fuding_phasing, haps = tmphaps, RF = RF)
	
	if (i == 1) {
		fuding_hap_res <- tmp_res
	}else{
		for (f in 1:length(fuding_hap_res)) {
			fuding_hap_res[[f]] <- rbind(fuding_hap_res[[f]],tmp_res[[f]])
		}
	}
	
}

#stopCluster(cl)


fuding_family_pred_haps <- NULL
for (i in 1:length(fuding_hap_res)) {
	tmp <- fuding_hap_res[[i]]
	id <- names(fuding_hap_res)[i]
	tmp <- rbind(NULL,tmp)
	tmp <- cbind(tmp, rep(id,nrow(tmp)))
	fuding_family_pred_haps <- rbind(fuding_family_pred_haps,tmp)
}

colnames(fuding_family_pred_haps) <- c("start","end","haps","chr","length","accession")

save(fuding_family_pred_haps,file="fuding_family_pred_haps.Rdata")


#j <- which(colnames(snps) == "X3W2.19")

