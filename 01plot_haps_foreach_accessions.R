options(stringsAsFactors=F)
load("00res_fuding_haps_all_for_plot.Rdata")
chr_dat <- read.table("chromosome_length_re_anchored.txt",header=T)
chr_len <- chr_dat[,2]


# xlim = c(- +1,15+1) ylim = c(-1+0,10+1)
chr_hap_plot <- function (dat, id, chr_len) {
	ymax <- max(chr_len)
	yratio <- 10 / ymax
	#cPal1 <- colorRampPalette(c("#CAD5C9","#697767"))
	cPal1 <- colorRampPalette(c("#F0E3DF","#2E5A71"))
	cPal2 <- colorRampPalette(c("#F0E3DF","#BEB9AE"))
	col1 <- cPal2(50)[as.numeric(cut(dat[,5],breaks = 50))]
	col2 <- cPal1(50)[as.numeric(cut(dat[,6],breaks = 50))]
	
	par(mar=c(1,1,3,1))
	plot (0,xlim=c(0,16),ylim=c(-1,11), type="n", axes=F, xlab="", ylab="", main = id)

	for (i in 1:nrow(dat)) {
		start <- dat[i,1]
		end <- dat[i,2]
		haps <- dat[i,3]
		chr <- dat[i,4]
		phu_hap <- dat[i,7]
		
		ytop <- 10 - (start * yratio)
		ybottom <- 10 - (end * yratio)
		
		if (haps == 1) {
			mycol <- "#2E5A71"
		}else{
			mycol <- "#9EB7BB"
		}
	
		if (phu_hap == 1) {
			xmain <- chr - 0.15
			xpara <- chr - 0.3
		}else{
			xmain <- chr + 0.15
			xpara <- chr + 0.3
		}
	
		rect(xmain-0.075,ybottom,xmain+0.075,ytop,border="white",col=col2[i],lwd=0.1)
		#rect(xpara-0.0375,ybottom,xpara+0.0375,ytop,border=NA,col=col2[i])
	}
}

pdf("01example_21_chr_haps_plot.pdf",width=8.36,height=9.86)
mysamples <- read.table("discriminant_analysis/training_samples.txt",header=F)
#samples <- unique(rownames(fuding_haps_all_for_plot))
samples <- mysamples[,1]
samples_rename <- sub("\\.","-",sub("X","",samples))
par(mfrow=c(6,4))

letters_ext <- NULL
for (i in 1:10) {
	tmp <- paste(letters[i],letters,sep="")
	letters_ext <- c(letters_ext,tmp)
}

# xlim = c(-1 +1,15+1) ylim = c(-1+0,10+1)
for (i in 1:length(samples)) {
	dat <- fuding_haps_all_for_plot[rownames(fuding_haps_all_for_plot) == samples[i],]
	chr_hap_plot(dat, samples_rename[i], chr_len)
	text(x=-1,y=12,labels=letters_ext[i],xpd=NA,font=2)
}

# haplotype block score # 
cPal <- colorRampPalette(c("#F0E3DF","#2E5A71"))
mycol <- cPal(5)

plot(0,axes=F,xlab="",ylab="",type="n",xlim=c(1,6),ylim=c(1,10))
# ---- #
mybase <- 1
mylen <- 1
# ---- #

for (i in 1:length(mycol)) {
	rect(mybase + mylen*(i-1), 5 - 0.4, mybase + mylen*i, 5 + 0.4, col=mycol[i],border="#2E5A71", xpd=NA, lwd=0.7 )
	
	if (i == 1) {
		text(x=mybase + mylen*(i-1), y=5 - 1.2, labels=0.5+0.1*(i-1), cex=0.75, xpd=NA)
		text(x=mybase + mylen*(i), y=5 - 1.2, labels=0.5+0.1*i, cex=0.75, xpd=NA)
	}else if (i == 3) {
		median_pos <- (mybase + mylen*(i) + mybase + mylen*(i-1)) / 2
		text(x=median_pos,y=5 + 1.2, labels="correct rate for each haplotype block", xpd=NA, cex=0.75)
		text(x=mybase + mylen*(i), y=5 - 1.2, labels=0.5+0.1*i, cex=0.75, xpd=NA)
	}else if (i == 5){
		text(x=mybase + mylen*(i), y=5 - 1.2, labels="1.0", cex=0.75, xpd=NA)
	}else{
		text(x=mybase + mylen*(i), y=5 - 1.2, labels=0.5+0.1*i, cex=0.75, xpd=NA)
	}
	
}


dev.off()


