# laad Predict fudingdabai haps # fuding_family_pred_haps #
load("fuding_family_pred_haps.Rdata")

accessions <- unique(fuding_family_pred_haps[,6])

ibs_res <- matrix(0,nrow=length(accessions),ncol=4)
for (i in 1:length(accessions)) {
	tmp <- fuding_family_pred_haps[fuding_family_pred_haps[,6] == accessions[i],]
	#start       end         haps chr length
	tmp1 <- rbind(NULL,tmp[tmp[,3] == 1,])
	tmp2 <- rbind(NULL,tmp[tmp[,3] == 2,])
	tmp3 <- rbind(NULL,tmp[tmp[,3] == 3,])
	tmp0 <- rbind(NULL,tmp[tmp[,3] == 0,])
	tmp4 <- rbind(NULL,tmp[tmp[,3] == 4,])
	tmpother <- rbind(NULL,tmp[nchar(tmp[,3]) == 2,])
	
	if (nrow(tmp1) == 0 & nrow(tmp2) == 0) {
		maxlen <- 0
	}else{
		maxlen <- max(as.numeric(tmp1[,5]),as.numeric(tmp2[,5]))
	}
	
	ibs1_len <- sum(as.numeric(tmp1[,5])) + sum(as.numeric(tmp2[,5])) + sum(as.numeric(tmp4[,5]))
	ibs2_len <- sum(as.numeric(tmp3[,5])) + sum(as.numeric(tmpother[,5]))
	ibs0_len <- sum(as.numeric(tmp0[,5]))
	
	total_len <- sum(ibs1_len+ibs2_len+ibs0_len)
	
	ibs1 <- ibs1_len/total_len
	ibs2 <- ibs2_len/total_len
	ibs0 <- ibs0_len/total_len
	
	ibs_res[i,] <- c(ibs0,ibs1,ibs2,maxlen)
}

rownames(ibs_res) <- accessions
colnames(ibs_res) <- c("IBS0","IBS1","IBS2","maxlen")

write.table(ibs_res,file="03IBS_results.txt",col.names=T,row.names=T,quote=F,sep="\t")

#co_file <- read.table("co_final.filter.txt",header=F)
#chr_vec <- paste("chr",1:15,sep="")
#gmt <- unique(co_file[,2])

#chr_len <- read.table("chromosome_length_re_anchored.txt",header=T)
#co_to_bins <- function(co_file, chr_len) {
#	chr_vec <- paste("chr",1:15,sep="")
#	bin_len <- NULL
#	for (i in 1:length(chr_vec)) {
#		if (chr_vec[i] %in% co_file[,1]) {
#			tmp <- co_file[co_file[,1] == chr_vec[i],]
#			tmp <- rbind(NULL,tmp)
#			bin_vec <- rep(0,nrow(tmp)*2)
#			bin_vec[c(T,F)] <- tmp[,5]
#			bin_vec[c(F,T)] <- tmp[,6]
#			bin_vec <- sort(c(1,bin_vec,chr_len[chr_vec[i],2]))
#			bin_mat <- matrix(bin_vec,ncol=2,byrow=T)
#			len <- bin_mat[,2] - bin_mat[,1]
#			bin_len <- c(bin_len,len)
#		}else{
#			len <- chr_len[chr_vec[i],2]
#			bin_len <- c(bin_len,len)
#		}
#	}
#	return(bin_len)
#}

#bin_len <- NULL
#for (i in 1:500) {
#	tar_gmt <- sample(gmt,size=20,replace=F)
#	tmp_co <- co_file[co_file[,2] %in% tar_gmt,]
#	tmp_bin_len <- co_to_bins(tmp_co,chr_len)
#	bin_len <- c(bin_len, max(tmp_bin_len))
#}
#quantile(bin_len)

