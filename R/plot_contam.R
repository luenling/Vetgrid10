setwd("/Volumes/Temp/Lukas/Data/Vienna_2010/Polymorphisms_SimCont/CheckContam")
vienna_2010=read.table("females_pI25_pII25_pIII25_pI18_pII18_pIel_pIIel_pIIIel_unfiltered_q20_simdiv_europe_100_table.txt")
vienna_2010_filt=read.table("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_q20_filtered_simdiv_europe_100_table.txt")
quartz()
col_names <- unlist(strsplit("chrom_pos_ref_pI25_pII25_pIII25_pI18_pII18_pIel_pIIel_pIIIel","_"))
colnames(vienna_2010) <- col_names
samples <- unlist(strsplit("pI25_pI18_pIel_pII25_pII18_pIIel_pIII25_NA_pIIIel","_"))
xmax=0.15
ymax = length(vienna_2010$pI25)
ymax = 5000
par(mfrow=c(3,3))
break_points=seq(0,0.5,0.01)
lablex = expression(paste("freq. of ",italic("D. sim.")," allele"))
labley = "num. of loci"
max_frq = c()
for (sample in samples){
	if(sample=="NA"){plot.new(); next; } 
	max_frq <- c(max_frq, max(vienna_2010[,sample]))
	hpI25 <- hist(vienna_2010[,sample],breaks=break_points,xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley)
	abline(v=median(vienna_2010[,sample]),col="red")
		
}
dev.copy2pdf(file="vienna2010_cont_unfilt.pdf")
quartz()
col_names <- unlist(strsplit("chrom_pos_ref_CCR 1_CCR 2_CCR 3_pI18_pII18_Base 1_Base 2_Base 3","_"))
colnames(vienna_2010) <- col_names
samples <- unlist(strsplit("CCR 1_Base 1_CCR 2_Base 2_CCR 3_Base 3","_"))
xmax=0.15
ymax = length(vienna_2010$pI25)
ymax = 5000
par(mfrow=c(3,2),oma = c(0, 0, 3, 0))
break_points=seq(0,0.5,0.01)
lablex = expression(paste("freq. of ",italic("D. sim.")," allele"))
labley = "num. of loci"
max_frq = c()
for (sample in samples){
	if(sample=="NA"){plot.new(); next; } 
	max_frq <- c(max_frq, max(vienna_2010[,sample]))
	hpI25 <- hist(vienna_2010[,sample],breaks=break_points,xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley)
	abline(v=median(vienna_2010[,sample]),col="red")
		
}
mtext("Vienna, 2010", outer = TRUE, cex = 1,font=2)
dev.copy2pdf(file="vienna2010_cont_unfilt_red.pdf")
quartz()

quartz()
col_names <- unlist(strsplit("chrom_pos_ref_pI25_pII25_pIII25_pIel_pIIel_pIIIel","_"))
colnames(vienna_2010_filt) <- col_names
samples <- unlist(strsplit("pI25_pIel_pII25_pIIel_pIII25_pIIIel","_"))
xmax=0.15
ymax = length(vienna_2010_filt$pI25)
ymax = 5000
par(mfrow=c(3,2))
break_points=seq(0,0.5,0.01)
lablex = expression(paste("freq. of ",italic("D. sim.")," allele"))
labley = "num. of loci"
max_frq = c()
for (sample in samples){
	if(sample=="NA"){plot.new(); next; } 
	max_frq <- c(max_frq, max(vienna_2010_filt[,sample]))
	hpI25 <- hist(vienna_2010_filt[,sample],breaks=break_points,xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley)
	abline(v=median(vienna_2010_filt[,sample]),col="red")
		
}
apply(vienna_2010_filt[samples],2,mean)
apply(vienna_2010_filt[samples],2,median)
# get median for populations
apply(vienna_2010[,-(1:3)],2,median)
dev.copy2pdf(file="vienna2010_cont_filtered.pdf")

setwd("/Volumes/Temp/Lukas/Data/Vienna_2011/Polymorphisms_SimCont/CheckContam")
quartz()
vienna_2011_filt=read.table("fem_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_simdiv_europe_100_table.txt")
vienna_2011=read.table("fem_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_unfiltered_q20_simdiv_europe_100_table.txt")
col_names_11 <- unlist(strsplit("chrom_pos_ref_pI25_pII25_pIII25_pIp_pIIp_pIIIp","_"))
colnames(vienna_2011) <- col_names_11
samples_11 <- unlist(strsplit("pI25_pIp_pII25_pIIp_pIII25_pIIIp","_"))
xmax=0.3
ymax = length(vienna_2011$pI25)
ymax = 5000
par(mfrow=c(3,2))
break_points=seq(0,0.5,0.01)
lablex = expression(paste("freq. of ",italic("D. sim.")," allele"))
labley = "num. of loci"
max_frq = c()
for (sample in samples_11){
	if(sample=="NA"){plot.new(); next; } 
	max_frq <- c(max_frq, max(vienna_2011[,sample]))
	hist(vienna_2011[,sample],breaks=break_points,xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley)
	abline(v=median(vienna_2011[,sample]),col="red")
}
apply(vienna_2011[,-(1:3)],2,median)
dev.copy2pdf(file="vienna2011_cont_unfilt.pdf")
quartz()
colnames(vienna_2011) <- unlist(strsplit("chrom_pos_ref_CCR 1_CCR 2_CCR 3_Base 1_Base 2_Base 3","_"))
samples_11 <- unlist(strsplit("CCR 1_Base 1_CCR 2_Base 2_CCR 3_Base 3","_"))
xmax=0.3
ymax = length(vienna_2011$pI25)
ymax = 5000
par(mfrow=c(3,2),oma = c(0, 0, 3, 0))
break_points=seq(0,0.5,0.01)
lablex = expression(paste("freq. of ",italic("D. sim.")," allele"))
labley = "num. of loci"
max_frq = c()
for (sample in samples_11){
	if(sample=="NA"){plot.new(); next; } 
	max_frq <- c(max_frq, max(vienna_2011[,sample]))
	hist(vienna_2011[,sample],breaks=break_points,xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley)
	abline(v=median(vienna_2011[,sample]),col="red")
}
mtext("Bolzano, 2011", outer = TRUE, cex = 1,font=2)
apply(vienna_2011[,-(1:3)],2,median)
dev.copy2pdf(file="vienna2011_cont_unfilt_red.pdf")



colnames(vienna_2011_filt) <- col_names_11
max=0.3
ymax = length(vienna_2011_filt$pI25)
ymax = 5000
par(mfrow=c(3,2))
break_points=seq(0,0.5,0.01)
lablex = expression(paste("freq. of ",italic("D. sim.")," allele"))
labley = "num. of loci"
max_frq = c()
for (sample in samples_11){
	if(sample=="NA"){plot.new(); next; } 
	max_frq <- c(max_frq, max(vienna_2011_filt[,sample]))
	hist(vienna_2011_filt[,sample],breaks=break_points,xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley)
	abline(v=median(vienna_2011_filt[,sample]),col="red")
}
apply(vienna_2011_filt[,-(1:3)],2,median)
dev.copy2pdf(file="vienna2011_cont_filtered.pdf")
# anova to see whether chrom and prob. simuls. are independent
summary(aov(pIel~chrom,vienna_2010))
aggregate(cbind(pI25,pIel)~chrom,vienna_2010,mean)

setwd("/Volumes/Temp/Lukas/Data/Trident/Contamination")
quartz()
Vie_T_cont=read.table("males_vie_pTL1a_pTL2a_pTD1a_pTD2a_pTL1b_pTL2b_pTD1b_pTD2b_unfiltered_q20.europe_sim_cov100.divergence.txt_table.txt")
col_names_11 <- unlist(strsplit("chrom_pos_ref_pIL_pIIL_pID_pIID_pILb_pIILb_pIDb_pIIDb","_"))
colnames(Vie_T_cont) <- col_names_11
samples_t10 <- unlist(strsplit("pIL_pIIL_pID_pIID_pILb_pIILb_pIDb_pIIDb","_"))
xmax=0.3
ymax = length(Vie_T_cont$pIL)
ymax = 5000
par(mfrow=c(3,2))
break_points=seq(0,0.5,0.01)
lablex = expression(paste("freq. of ",italic("D. sim.")," allele"))
labley = "num. of loci"
max_frq = c()
for (sample in samples_t10){
	if(sample=="NA"){plot.new(); next; } 
	max_frq <- c(max_frq, max(Vie_T_cont[,sample]))
	hist(Vie_T_cont[,sample],breaks=break_points,xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley)
	abline(v=median(Vie_T_cont[,sample]),col="red")
}
apply(vienna_2011[,-(1:3)],2,median)
dev.copy2pdf(file="vie_T_cont_unfilt.pdf")
colnames(vienna_2011) <- col_names_t10
max=0.3
ymax = length(vienna_2011_filt$pI25)
ymax = 5000
par(mfrow=c(3,2))
break_points=seq(0,0.5,0.01)
lablex = expression(paste("freq. of ",italic("D. sim.")," allele"))
labley = "num. of loci"
max_frq = c()
for (sample in samples_11){
	if(sample=="NA"){plot.new(); next; } 
	max_frq <- c(max_frq, max(vienna_2011_filt[,sample]))
	hist(vienna_2011_filt[,sample],breaks=break_points,xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley)
	abline(v=median(vienna_2011_filt[,sample]),col="red")
}
apply(vienna_2011_filt[,-(1:3)],2,median)
dev.copy2pdf(file="vienna2011_cont_filtered.pdf")
# anova to see whether chrom and prob. simuls. are independent
summary(aov(pIel~chrom,vienna_2010))
aggregate(cbind(pI25,pIel)~chrom,vienna_2010,mean)

setwd("/Volumes/Temp/Lukas/Data/Trident/Contamination")
quartz()
vienna_2010_filt=read.table("/Volumes/Temp/Lukas/Data/Trident/Realigned/Polymorphisms/males_vie_pTL1_pTL2_pTD1_pTD2_q20.europe_sim_cov100.divergence.txt_table.txt")
vienna_2010=read.table("males_vie_pTL1a_pTL2a_pTD1a_pTD2a_pTL1b_pTL2b_pTD1b_pTD2b_unfiltered_q20.europe_sim_cov100.divergence.txt_table.txt")
col_names_10 <- unlist(strsplit("chrom_pos_ref_pTL1a_pTL2a_pTD1a_pTD2a_pTL1b_pTL2b_pTD1b_pTD2b","_"))
colnames(vienna_2010) <- col_names_10
samples_10 <- unlist(strsplit("pTL1a_pTL2a_pTD1a_pTD2a_pTL1b_pTL2b_pTD1b_pTD2b","_"))
xmax=0.1
#ymax = length(vienna_2010$pI25)
ymax = 1000
par(mfrow=c(4,2))
break_points=seq(0,0.5,0.005)
lablex = expression(paste("freq. of ",italic("D. sim.")," allele"))
labley = "num. of loci"
max_frq = c()
for (sample in samples_10){
	if(sample=="NA"){plot.new(); next; } 
	max_frq <- c(max_frq, max(vienna_2010[,sample]))
	hist(vienna_2010[,sample],breaks=break_points,xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley)
	abline(v=median(vienna_2010[,sample]),col="red")
}
apply(vienna_2010[,-(1:3)],2,median)
dev.copy2pdf(file="vienna10_trident_cont_unfilt.pdf")
col_names_10 <- unlist(strsplit("chrom_pos_ref_pTL1_pTL2_pTD1_pTD2","_"))
samples_10 <- unlist(strsplit("pTL1_pTL2_pTD1_pTD2","_"))
colnames(vienna_2010_filt) <- col_names_10
max=0.1
#ymax = length(vienna_2010_filt$pI25)
ymax = 1000
par(mfrow=c(2,2))
break_points=seq(0,0.5,0.005)
lablex = expression(paste("freq. of ",italic("D. sim.")," allele"))
labley = "num. of loci"
max_frq = c()
for (sample in samples_10){
	if(sample=="NA"){plot.new(); next; } 
	max_frq <- c(max_frq, max(vienna_2010_filt[,sample]))
	hist(vienna_2010_filt[,sample],breaks=break_points,xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley)
	abline(v=median(vienna_2010_filt[,sample]),col="red")
}
apply(vienna_2010_filt[,-(1:3)],2,median)
dev.copy2pdf(file="vienna10_trident_cont_filtered.pdf")


setwd("/Volumes/Temp/Lukas/Data/SeasonalData/seasonal_raw_data/Polymorphisms")
quartz()
seasonal_unfilt=read.table("/Volumes/Temp/Lukas/Data/SeasonalData/seasonal_raw_data/Polymorphisms/seasonal_July_Sept_Oct_Nov_q20.europe_sim_cov100.divergence.txt_table.txt")
#seasonal_filt=read.table("")
col_names_season <- unlist(strsplit("chrom_pos_ref_July_Sept_Oct_Nov","_"))
colnames(seasonal_unfilt) <- col_names_season
samples_seas <- unlist(strsplit("July_Sept_Oct_Nov","_"))
xmax=0.1
#ymax = length(vienna_2010$pI25)
ymax = 5000
par(mfrow=c(2,2))
break_points=seq(0,0.5,0.005)
lablex = expression(paste("freq. of ",italic("D. sim.")," allele"))
labley = "num. of loci"
max_frq = c()
for (sample in samples_seas){
	if(sample=="NA"){plot.new(); next; } 
	max_frq <- c(max_frq, max(seasonal_unfilt[,sample]))
	hist(seasonal_unfilt[,sample],breaks=break_points,xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley)
	abline(v=median(seasonal_unfilt[,sample]),col="red")
}
apply(seasonal_unfilt[,-(1:3)],2,median)
dev.copy2pdf(file="seasonal_cont_unfilt.pdf")
setwd("/Volumes/Temp/Lukas/Data/SeasonalData/seasonal_raw_data/Poly_filtered/")
quartz()
seasonal_filt=read.table("/Volumes/Temp/Lukas/Data/SeasonalData/seasonal_raw_data/Poly_filtered/seasonal_July_Sept_Oct_Nov_cleaned_q20.europe_sim_cov100.divergence.txt_table.txt")
#seasonal_filt=read.table("")
col_names_season <- unlist(strsplit("chrom_pos_ref_July_Sept_Oct_Nov","_"))
colnames(seasonal_filt) <- col_names_season
samples_seas <- unlist(strsplit("July_Sept_Oct_Nov","_"))
xmax=0.1
#ymax = length(vienna_2010$pI25)
ymax = 5000
par(mfrow=c(2,2))
break_points=seq(0,0.5,0.005)
lablex = expression(paste("freq. of ",italic("D. sim.")," allele"))
labley = "num. of loci"
max_frq = c()
for (sample in samples_seas){
	if(sample=="NA"){plot.new(); next; } 
	max_frq <- c(max_frq, max(seasonal_filt[,sample]))
	hist(seasonal_filt[,sample],breaks=break_points,xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley)
	abline(v=median(seasonal_filt[,sample]),col="red")
}
apply(seasonal_filt[,-(1:3)],2,median)
dev.copy2pdf(file="seasonal_cont_filt.pdf")


setwd("/Volumes/Temp/Lukas/Data/M_Schaeffer_PE_bams/Polymorphisms")
quartz()
lati_unfilt=read.table("m_schaeffer_KEN_ETH_ATH_BAR_DIJ_ESX_BOW_LNZ_q20_filt_mc1.5_chroms_only.europe_sim_cov100.divergence.txt_table.txt")
#lati_filt=read.table("")
col_names_lati <- unlist(strsplit("chrom_pos_ref_KEN_ETH_ATH_BAR_DIJ_ESX_BOW_LNZ","_"))
colnames(lati_unfilt) <- col_names_lati
samples_lati <- unlist(strsplit("KEN_ETH_ATH_BAR_DIJ_ESX_BOW_LNZ","_"))
xmax=0.1
#ymax = length(vienna_2010$pI25)
ymax = 5000
par(mfrow=c(2,4))
break_points=seq(0,0.5,0.005)
lablex = expression(paste("freq. of ",italic("D. sim.")," allele"))
labley = "num. of loci"
max_frq = c()
for (sample in samples_lati){
	if(sample=="NA"){plot.new(); next; } 
	max_frq <- c(max_frq, max(lati_unfilt[,sample]))
	#hist(lati_unfilt[,sample],breaks=break_points,xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley)
	hist(lati_unfilt[,sample],xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley,breaks=50)
	abline(v=median(lati_unfilt[,sample]),col="red")
}
apply(lati_unfilt[,-(1:3)],2,median)
dev.copy2pdf(file="shaeff_cont_unfilt.pdf")

setwd("/Volumes/Temp/Lukas/Data/Florida_Penn_Maine/Polymorphisms/")
quartz(width=7,height=9)
lati_unfilt=read.table("/Volumes/Temp/Lukas/Data/Florida_Penn_Maine/Polymorphisms/lat_FLOR_PENN_MAINE_q20.europe_sim_cov100.divergence.txt_table.txt")
#lati_filt=read.table("")
col_names_lati <- unlist(strsplit("chrom_pos_ref_Florida_Penn_Maine","_"))
colnames(lati_unfilt) <- col_names_lati
samples_lati <- unlist(strsplit("Florida_Penn_Maine","_"))
xmax=0.03
#ymax = length(vienna_2010$pI25)
ymax = 5000
par(mfrow=c(3,1))
break_points=seq(0,1.0,0.001)
lablex = expression(paste("freq. of ",italic("D. sim.")," allele"))
labley = "num. of loci"
max_frq = c()
for (sample in samples_lati){
	if(sample=="NA"){plot.new(); next; } 
	max_frq <- c(max_frq, max(lati_unfilt[,sample]))
	hist(lati_unfilt[,sample],breaks=break_points,xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley)
	#hist(lati_unfilt[,sample],xlim=c(0,xmax),ylim =c(0,ymax), col="gray",main=sample,xlab=lablex, ylab=labley,breaks=seq(0,1,0.001))
	abline(v=median(lati_unfilt[,sample]),col="red")
}
apply(lati_unfilt[,-(1:3)],2,median)
dev.copy2pdf(file="lati_FPM_cont_unfilt.pdf")


for (i in 1:4){
	for (j in 10:12){
		print(j)
		if (j > 10){break; next;}
	}
	print(i)
	
}
