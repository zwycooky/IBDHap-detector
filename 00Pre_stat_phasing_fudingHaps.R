options(stringsAsFactors=F)
library(HMM)

haps <- read.table("correct.haps.filtered.formated.txt",header=F)
snps <- read.table("formated_w221_snp.filtered.txt",header=T)
snps <- snps[haps[,2],]


linkage_files <- paste("linkage_group/","chr",1:15,".map.ordered.txt",sep="")

fuding_phasing <- function(haps, snps) {
	library(HMM)
	
	determine_haps <- function(subhap) {
		sum1 <- sum(subhap == 1)
		sum2 <- sum(subhap == 2)
		if (sum1 < sum2) {
			return(2)
		}else{
			return(1)
		}
	}
	
	classBydisstance <- function(fuding_hap, mydis_cutoff = 1000000) {
		dis <- as.numeric(names(fuding_hap))
		dis_1MB <- as.numeric(dis[-1] - dis[-length(dis)] < mydis_cutoff)
		dis_1MB_vec <- dis_1MB[-1] - dis_1MB[-length(dis_1MB)]
		end_site <- which(dis_1MB_vec == -1) + 1
		start_site <- which(dis_1MB_vec == 1) + 1
		if (dis_1MB[1] == 1) {
			start_site <- c(1,start_site)
		}
		if (dis_1MB[length(dis_1MB)] == 1) {
			end_site <- c(end_site,length(dis))
		}
		res <- cbind(start_site,end_site)
		return(res)
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
			bin_select_reg <- haps[,6] >= ss & haps[,6] <= ee
			
			# het ratio #
			bin_hetsnp <- sum(nchar(snps[bin_select_reg]) == 2, na.rm=T)
			bin_homsnp <- sum(nchar(snps[bin_select_reg]) == 1, na.rm=T)
			ratio_het <- round(bin_hetsnp / (bin_hetsnp + bin_homsnp),3)
			
			# correct ratio #
			correctnum <- sum(correctGeno == fuding_hap)
			allnum <- length(fuding_hap)
			ratio_correct <- round(correctnum/allnum,3)
			
			hapreg <- cbind(ss, ee, final_hap, as.numeric(haps[1,5]), ratio_het, ratio_correct)
			#return(list(hapreg=hapreg, correcthaps = correctGeno))
			rownames(hapreg) <- NULL
			colnames(hapreg) <- NULL
			return(hapreg)
		}else{
			bin_vec <- rep(0,nrow(res)*ncol(res))
			bin_vec[c(T,F)] <- res[,1]
			bin_vec[c(F,T)] <- res[,2]
			bin_vec <- as.numeric( c(names(fuding_hap)[1], bin_vec, names(fuding_hap)[length(fuding_hap)]) )
			bin_mat <- matrix(bin_vec,ncol=2,byrow=T)
			hapreg <- matrix(0,nrow=nrow(bin_mat),ncol=6)
			for (b in 1:nrow(bin_mat)) {
				tmp_s <- which(as.numeric(names(fuding_hap)) == bin_mat[b,1])
				tmp_e <- which(as.numeric(names(fuding_hap)) == bin_mat[b,2])
				final_hap <- determine_haps(fuding_hap[tmp_s:tmp_e])
				bin_select_reg <- haps[,6] >= bin_mat[b,1] & haps[,6] <= bin_mat[b,2]
				
				# het ratio #
				bin_hetsnp <- sum(nchar(snps[bin_select_reg]) == 2, na.rm=T)
				bin_homsnp <- sum(nchar(snps[bin_select_reg]) == 1, na.rm=T)
				ratio_het <- round(bin_hetsnp / (bin_hetsnp + bin_homsnp),3)
				
				# correcet ratio #
				correctnum <- sum(correctGeno[tmp_s:tmp_e] == fuding_hap[tmp_s:tmp_e])
				allnum <- length(tmp_s:tmp_e)
				ratio_correct <- round(correctnum/allnum,3)
				
				hapreg[b,] <- c(bin_mat[b,1], bin_mat[b,2], final_hap, as.numeric(haps[1,5]), ratio_het, ratio_correct)
			}
		
			#return(list(hapreg=hapreg, correcthaps = correctGeno))
			rownames(hapreg) <- NULL
			colnames(hapreg) <- NULL
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
	
	fuding_hap_phu1_reg <- cbind(fuding_hap_phu1_reg,rep(1,nrow(fuding_hap_phu1_reg)))
	fuding_hap_phu2_reg <- cbind(fuding_hap_phu2_reg,rep(2,nrow(fuding_hap_phu2_reg)))
	
	fuding_hap1_2_res <- rbind(fuding_hap_phu1_reg, fuding_hap_phu2_reg)
	rownames(fuding_hap1_2_res) <- NULL
	colnames(fuding_hap1_2_res) <- NULL
	return(fuding_hap1_2_res)
}

library(pbapply)
library(parallel)
cl <- makeCluster(4)

fuding_hap_res <- NULL
for (i in 1:length(linkage_files)) {
	tmpgmap <- read.table(linkage_files[i])
	tmphaps <- haps[haps[,1] %in% tmpgmap[,1],]
	tmpsnps <- snps[tmphaps[,2],]
	
	tmp_res <- pbapply(tmpsnps, 2, fuding_phasing, haps = tmphaps, cl = cl)
	
	if (i == 1) {
		fuding_hap_res <- tmp_res
	}else{
		for (f in 1:length(fuding_hap_res)) {
			fuding_hap_res[[f]] <- rbind(fuding_hap_res[[f]],tmp_res[[f]])
		}
	}
	
}

stopCluster(cl)

fuding_haps_all_for_plot <- NULL
sample_ids <- NULL
for (i in 1:length(fuding_hap_res)) {
	tmp <- fuding_hap_res[[i]]
	id <- names(fuding_hap_res)[i]
	tmp <- rbind(NULL,tmp)
	sample_ids <- c(sample_ids,rep(id,nrow(tmp)))
	fuding_haps_all_for_plot <- rbind(fuding_haps_all_for_plot,tmp)
}
rownames(fuding_haps_all_for_plot) <- sample_ids
colnames(fuding_haps_all_for_plot) <- c("start","end","haps","chr","het_ratio","correcet_ratio","phu_hap")

write.table(fuding_haps_all_for_plot,"00res_fuding_haps_all_for_plot.txt",row.names=T,col.names=T,quote=F,sep="\t")

save(fuding_haps_all_for_plot,file = "00res_fuding_haps_all_for_plot.Rdata")
#save.image(file="20200308.Rdata")


#j <- which(colnames(snps) == "X3W2.19")

