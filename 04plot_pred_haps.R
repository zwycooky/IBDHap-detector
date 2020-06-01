options(stringsAsFactors=F)
load("00res_fuding_haps_all_for_plot.Rdata")
chr_dat <- read.table("chromosome_length_re_anchored.txt",header=T)
chr_len <- chr_dat[,2]

pred_hapdat <- read.table("03IBS_results.txt",header=T)
# laad Predict fudingdabai haps # fuding_family_pred_haps #
load("fuding_family_pred_haps.Rdata")
fuding_family_pred_haps <- fuding_family_pred_haps[fuding_family_pred_haps[,6] != "Fudingdabai",]
# load training samples #
mysamples <- read.table("discriminant_analysis/training_samples.txt",header=F)

# load sample info #
sample_info <- read.table("sample_info.txt",header=F,sep="\t")

# ordered plot samples #
pred_hapdat <- pred_hapdat[order(as.numeric(pred_hapdat[,1])),]
pred_hapdat <- pred_hapdat[rownames(pred_hapdat) != "Fudingdabai",]
#pred_hapdat <- pred_hapdat[!rownames(pred_hapdat) %in% mysamples[,1],]


# xlim = c(- +1,15+1) ylim = c(-1+0,10+1)
# chr_hap_plot(dat_ori, dat_pred, samples_id, tmpibs, chr_len)
format_chr_pos <- function (chr_vec,chr_num,pos) {
	pos + sum(chr_vec[1:chr_num])
}

chr_hap_plot <- function (dat_ori, dat_pred, id, ibs, chr_len, ypos=1) {
	chr_vec <- c(0,chr_len)
	chr_mar <- sum(chr_vec) * 0.01
	
	#cPal1 <- colorRampPalette(c("#CAD5C9","#697767"))
	cPal <- colorRampPalette(c("#F0E3DF","#2E5A71"))
	mycol <- cPal(50)[as.numeric(cut(c(dat_ori[,6],0.5,1),breaks = 50))]
	
	# plot ori #
	for (i in 1:nrow(dat_ori)) {
		start <- dat_ori[i,1]
		end <- dat_ori[i,2]
		haps <- dat_ori[i,3]
		chr <- dat_ori[i,4]
		phu_hap <- dat_ori[i,7]
		
		start <- start + sum(chr_vec[1:chr]) + (chr - 1) * chr_mar
		end <- end + sum(chr_vec[1:chr]) + (chr - 1) * chr_mar
	
		if (phu_hap == 1) {
			ypos_ori <- ypos - 0.18
		}else{
			ypos_ori <- ypos + 0.18
		}
	
		rect(start,ypos_ori - 0.1,end,ypos_ori + 0.1,border="white",col=mycol[i], lwd=0.1)
	}
	
	# plot pred #
	for (i in 1:nrow(dat_pred)) {
		start <- dat_pred[i,1]
		end <- dat_pred[i,2]
		haps <- dat_pred[i,3]
		chr <- dat_pred[i,4]
		
		if (haps == 0) {
			next
		}else{
			start <- start + sum(chr_vec[1:chr]) + (chr - 1) * chr_mar
			end <- end + sum(chr_vec[1:chr]) + (chr - 1) * chr_mar
			if (haps == 1) {
				ypos_pred <- ypos - 0.18
				rect(start,ypos_pred - 0.1,end,ypos_pred + 0.1,border="#E85B51", lwd = 0.2)
			}else if (haps == 2) {
				ypos_pred <- ypos + 0.18
				rect(start,ypos_pred - 0.1,end,ypos_pred + 0.1,border="#FFC31C", lwd = 0.2)
			}else if (haps == 3) {
				ypos_pred <- ypos + 0.18
				rect(start,ypos_pred - 0.1,end,ypos_pred + 0.1,border="#FC8C3D", lwd = 0.2)
				ypos_pred <- ypos - 0.18
				rect(start,ypos_pred - 0.1,end,ypos_pred + 0.1,border="#FC8C3D", lwd = 0.2)
			}else if (haps == 11) {
				ypos_pred <- ypos + 0.18
				rect(start,ypos_pred - 0.1,end,ypos_pred + 0.1,border="#E85B51", lwd = 0.2)
				ypos_pred <- ypos - 0.18
				rect(start,ypos_pred - 0.1,end,ypos_pred + 0.1,border="#E85B51", lwd = 0.2)
			}else if (haps == 22) {
				ypos_pred <- ypos + 0.18
				rect(start,ypos_pred - 0.1,end,ypos_pred + 0.1,border="#FC8C3D", lwd = 0.2)
				ypos_pred <- ypos - 0.18
				rect(start,ypos_pred - 0.1,end,ypos_pred + 0.1,border="#FC8C3D", lwd = 0.2)
			}
		}
	}
	
}

## 3 * 7 ##
pdf("04example_chr_haps_plot.pdf",width=7,height=20)
#par(mfrow=c(66,3))

letters_ext <- letters
for (i in 1:10) {
	tmp <- paste(letters[i],letters,sep="")
	letters_ext <- c(letters_ext,tmp)
}


# xlim = c(-1 +1,15+1) ylim = c(-1+0,10+1)
par(mar=c(3,4,4,2))
plot(0,type="n",ylab="accessions",axes=F,xlim=c(sum(chr_len) * -0.15, sum(chr_len) * 1.15), ylim=c(0,nrow(pred_hapdat)))

# IBS header #
text(sum(chr_len) * -0.1459,221.5,labels="IBD0",cex=0.35,adj=0)
text(sum(chr_len) * -0.1055,221.5,labels="IBD1",cex=0.35,adj=0)
text(sum(chr_len) * -0.065,221.5,labels="IBD2",cex=0.35,adj=0)
lines(x=c(sum(chr_len) * -0.15,sum(chr_len) * -0.035), y=rep(221,2),lwd=0.5)

# chromosome lab #
chr_mar <- sum(chr_len) * 0.01
chr_vec <- c(0,chr_len)
for (i in 1:(length(chr_vec)-1)) {
	pos <- ( sum(chr_vec[1:i]) + sum(chr_vec[1:(i+1)]) ) / 2 + (i - 1) * chr_mar
	text(x=pos,y=222,labels=paste("chr",i),cex=0.4,xpd=NA)
}
median_pos <- (sum(chr_len) + 14*chr_mar) / 2
text(x=median_pos,y=-3,labels="chromosome")


for (i in 1:nrow(pred_hapdat)) {
#for (i in 1:2) {
	dat_ori <- fuding_haps_all_for_plot[rownames(fuding_haps_all_for_plot) == rownames(pred_hapdat)[i],]
	dat_pred <- fuding_family_pred_haps[fuding_family_pred_haps[,6] == rownames(pred_hapdat)[i],]
	dat_pred <- dat_pred[,-6]
	
	dat_pred_numeric <- matrix(0,nrow=nrow(dat_pred),ncol=ncol(dat_pred))
	for (j in 1:5) {
		dat_pred_numeric[,j] <- as.numeric(dat_pred[,j])
	}
	
	tmpibs <- as.numeric(pred_hapdat[rownames(pred_hapdat)[i],1:3])
	tmpibs_for_plot <- NULL
	for(tt in 1:length(tmpibs)) {
		if (tmpibs[tt] <= 0.0005) {
			tmpibs_for_plot <- c(tmpibs_for_plot,"0.000")
		}else if (tmpibs[tt] >= 0.9995) {
			tmpibs_for_plot <- c(tmpibs_for_plot,"1.000")
		}else{
			tmp_ibs <- as.character(round(tmpibs[tt],3))
			sub4 <- 5 - nchar(tmp_ibs)
			if (sub4 > 0) {
				tmp_ibs <- paste(c(tmp_ibs,rep(0,sub4)),collapse="")
			}
			
			tmpibs_for_plot <- c(tmpibs_for_plot, tmp_ibs)
		}
	}
	
	sample_id <- rownames(pred_hapdat)[i]
	# set colors for training data # 
	sample_col <- "black"
	if (sample_id %in% mysamples[,1]) {
		sample_col <- "#E85B51"
	}
	# re-format sample id to 3W1-1 style #
	sample_id <- sub("\\.","-",sub("X","",sample_id))
	
	ypos <- nrow(pred_hapdat) - i + 1
	
	chr_hap_plot(dat_ori, dat_pred_numeric, rownames(pred_hapdat)[i], tmpibs, chr_len, ypos )
	
	# text for sample id #
	text(x=sum(chr_len) * -0.22,y=ypos,labels=sample_id,xpd=NA,font=1,cex=0.35, adj=0, col=sample_col)
	# text for IBS #
	text(x=sum(chr_len) * -0.15,y=ypos,paste(tmpibs_for_plot,collapse="  "),xpd=NA,font=1,cex=0.35, adj=0)
	# text for NO. sample id #
	text(x=sum(chr_len) * -0.24,y=ypos,paste(i,":",sep=""),xpd=NA,font=1,cex=0.35,adj=1)
	
	# define first degree relative #
	if (tmpibs[2] > 0.9) {
		points(x=sum(chr_len) * -0.28,y=ypos,cex=0.4,col="#E85B51",xpd=NA,pch=20)
	}else if (tmpibs[3] > 0.1) {
		#points(x=sum(chr_len) * -0.29,y=ypos,cex=0.4,col="#FC8C3D",xpd=NA,pch=20)
	}
	
	# ancient tree or cultivar #
	sample_type <- sample_info[sample_info[,1] == sample_id,2]
	if (sample_type == "landrace") {
		#points(x=sum(chr_len) * -0.30,y=ypos,cex=0.4,col="#2E5A71",xpd=NA,pch=20)
		#text(x=sum(chr_len) * -0.31,y=ypos,cex=0.35,xpd=NA, labels=sample_info[sample_info[,1] == sample_id,3],adj=1)
	}
	#if (tmpibs[2] > 0.25 & tmpibs[3] <= 0.1) {
	#	points(x=sum(chr_len) * -0.3,y=ypos,cex=0.4,col="#697767",xpd=NA,pch=20)
	#}
}

# plot legend #
rect(sum(chr_len) * -0.24, 228, sum(chr_len) * -0.20, 228.4, border="#E85B51", lwd=0.7, col=NA, xpd=NA)
text(x=sum(chr_len) * -0.18, y=228.2, labels="Predicted Fudingdabai haplotype A",cex=0.35,xpd=NA,adj=0)
rect(sum(chr_len) * -0.24, 227, sum(chr_len) * -0.20, 227.4, border="#FFC31C", lwd=0.7, col=NA, xpd=NA)
text(x=sum(chr_len) * -0.18, y=227.2, labels="Predicted Fudingdabai haplotype B", cex=0.35,xpd=NA,adj=0)
rect(sum(chr_len) * -0.24, 226, sum(chr_len) * -0.20, 226.4, border="#FC8C3D", lwd=0.7, col=NA, xpd=NA)
text(x=sum(chr_len) * -0.18, y=226.2, labels="Predicted Fudingdabai haplotypes (A and B)", cex=0.35,xpd=NA,adj=0)

# haplotype block score # 
cPal <- colorRampPalette(c("#F0E3DF","#2E5A71"))
mycol <- cPal(5)

# ---------------- #
mybase <- sum(chr_len) * 0.14
mylen <- 0.05*sum(chr_len)
# ---------------- #

for (i in 1:length(mycol)) {
	rect(mybase + mylen*(i-1), 227.4 - 0.4, mybase + mylen*i, 227.4 + 0.4, col=mycol[i],border="#2E5A71", xpd=NA, lwd=0.7 )
	
	if (i == 1) {
		text(x=mybase + mylen*(i-1), y=227.4 - 1.2, labels=0.5+0.1*(i-1), cex=0.35, xpd=NA)
		text(x=mybase + mylen*(i), y=227.4 - 1.2, labels=0.5+0.1*i, cex=0.35, xpd=NA)
	}else if (i == 3) {
		median_pos <- (mybase + mylen*(i) + mybase + mylen*(i-1)) / 2
		text(x=median_pos,y=227.4 + 1.2, labels="correct rate for each IBD segment", xpd=NA, cex=0.35)
		text(x=mybase + mylen*(i), y=227.4 - 1.2, labels=0.5+0.1*i, cex=0.35, xpd=NA)
	}else if (i == 5){
		text(x=mybase + mylen*(i), y=227.4 - 1.2, labels="1.0", cex=0.35, xpd=NA)
	}else{
		text(x=mybase + mylen*(i), y=227.4 - 1.2, labels=0.5+0.1*i, cex=0.35, xpd=NA)
	}
	
}

# legend #
points (x=sum(chr_len) * 0.45, y=228.4, pch=20, cex=0.7, col="#E85B51", xpd=NA)
text (x=sum(chr_len) * 0.47, y=228.4, labels="First degree relative with Fudingdabai", cex=0.35, xpd=NA,adj=0)

#points (x=sum(chr_len) * 0.45, y=227.4, pch=20, cex=0.7, col="#FC8C3D", xpd=NA)
#text (x=sum(chr_len) * 0.47, y=227.4, labels="Putative triploid offspring of Fudingdabai", cex=0.35, xpd=NA,adj=0)

#points (x=sum(chr_len) * 0.45, y=226.4, pch=20, cex=0.7, col="#2E5A71", xpd=NA)
#text (x=sum(chr_len) * 0.47, y=226.4, labels="Ancient tree", cex=0.35, xpd=NA,adj=0)

dev.off()



