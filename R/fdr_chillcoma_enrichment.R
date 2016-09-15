# get MAF spectra for SNPs
setwd("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/Joined_Analysis/")
vie2010_maf=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25_nowolb.avgMAF.af",header=F)
ita2011_maf=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20_filt_mc8_mcv25_nowolb.avgMAF.af",header=F)
over_maf=read.table("vie2010_ita2011_avgMAF_nowolb.af",header=F)
maf_files=c("CHR","BPS","Alleles","MAF")
colnames(vie2010_maf) = maf_files
colnames(ita2011_maf) = maf_files
colnames(over_maf) = maf_files
over_maf$MAF[over_maf$MAF > 0.5] = 1 - over_maf$MAF[over_maf$MAF > 0.5]

# loading all AF

vie2010_maf=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/females_pI_II_III25_pI_II_IIIel_real_nowolb.avg_frq.af",header=F)
ita2011_maf=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/females11__pI_II_III25_pI_II_IIIhl_real_avg.af",header=F)
#over_maf=read.table("vie2010_ita2011_avgMAF_nowolb.af",header=F)
maf_files=c("CHR","BPS","Alleles","CS_AF","B_AF")
colnames(vie2010_maf) = maf_files
colnames(ita2011_maf) = maf_files
vie2010_maf$CS_MAF = vie2010_maf$CS_AF
vie2010_maf$CS_MAF[vie2010_maf$CS_MAF > 0.5] = 1 - vie2010_maf$CS_MAF[vie2010_maf$CS_MAF > 0.5]
ita2011_maf$CS_MAF = ita2011_maf$CS_AF
ita2011_maf$CS_MAF[ita2011_maf$CS_MAF > 0.5] = 1 - ita2011_maf$CS_MAF[ita2011_maf$CS_MAF > 0.5]
vie2010_maf$B_MAF = vie2010_maf$B_AF
vie2010_maf$B_MAF[vie2010_maf$B_MAF > 0.5] = 1 - vie2010_maf$B_MAF[vie2010_maf$B_MAF > 0.5]
ita2011_maf$B_MAF = ita2011_maf$B_AF
ita2011_maf$B_MAF[ita2011_maf$B_MAF > 0.5] = 1 - ita2011_maf$B_MAF[ita2011_maf$B_MAF > 0.5]
# draw histograms of MAF
quartz()
hist(vie2010_maf$CS_MAF,xlab="MAF",ylab="density",freq=F,breaks=49,col=rgb(1,0,0,1/4), xlim=c(0,0.5),ylim=c(0,15),main="Vienna 2010")
hist(vie2010_maf$B_MAF,xlab="MAF",ylab="density",freq=F,breaks=49,col=rgb(0,1,0,1/4), xlim=c(0,0.5),add=T)
legend("topright", c("CCR", "base"), fill=c(rgb(1,0,0,1/4),rgb(0,1,0,1/4)))
dev.copy2pdf(file="maf_vie2010_cs_base.pdf")
quartz()
hist(ita2011_maf$CS_MAF,xlab="MAF",ylab="density",freq=F,breaks=49,col=rgb(1,0,0,1/4), xlim=c(0,0.5),ylim=c(0,15),main="Bolzano 2011")
hist(ita2011_maf$B_MAF,xlab="MAF",ylab="density",freq=F,breaks=49,col=rgb(0,1,0,1/4), xlim=c(0,0.5),add=T)
legend("topright", c("CCR", "base"), fill=c(rgb(1,0,0,1/4),rgb(0,1,0,1/4)))
dev.copy2pdf(file="maf_ita2011_cs_base.pdf")

# checking the significant snps
vie2010_fdr_af$B_MAF=apply(vie2010_fdr_af[,7:9],1,mean)
vie2010_fdr_af$CS_MAF=apply(vie2010_fdr_af[,4:6],1,mean)
ita2011_fdr_af$B_MAF=apply(ita2011_fdr_af[,7:9],1,mean)
ita2011_fdr_af$CS_MAF=apply(ita2011_fdr_af[,4:6],1,mean)
keeps=c("CHR","BPS","Alleles","CS_MAF","B_MAF")
vie2010_fdr_af=vie2010_fdr_af[,keeps]
ita2011_fdr_af=ita2011_fdr_af[,keeps]
vie2010_fdr_af$CS_MAF[vie2010_fdr_af$B_MAF > 0.5] = 1 - vie2010_fdr_af$CS_MAF[vie2010_fdr_af$B_MAF > 0.5]
vie2010_fdr_af$B_MAF[vie2010_fdr_af$B_MAF > 0.5] = 1 - vie2010_fdr_af$B_MAF[vie2010_fdr_af$B_MAF > 0.5]
ita2011_fdr_af$CS_MAF[ita2011_fdr_af$B_MAF > 0.5] = 1 - ita2011_fdr_af$CS_MAF[ita2011_fdr_af$B_MAF > 0.5] 
ita2011_fdr_af$B_MAF[ita2011_fdr_af$B_MAF > 0.5] = 1 - ita2011_fdr_af$B_MAF[ita2011_fdr_af$B_MAF > 0.5] 
vie2010_fdr_af$AF_change=vie2010_fdr_af$CS_MAF-vie2010_fdr_af$B_MAF
ita2011_fdr_af$AF_change=ita2011_fdr_af$CS_MAF-ita2011_fdr_af$B_MAF
quartz()
hist(vie2010_fdr_af$CS_MAF,xlab="MAF",ylab="density",freq=F,breaks=24,col=rgb(1,0,0,1/4), xlim=c(0,0.5),main="Vienna 2010")
hist(vie2010_fdr_af$B_MAF,xlab="MAF",ylab="density",freq=F,breaks=24,col=rgb(0,1,0,1/4), xlim=c(0,0.5),add=T)
legend("topright", c("CCR", "base"), fill=c(rgb(1,0,0,1/4),rgb(0,1,0,1/4)))
quartz()
hist(ita2011_fdr_af$CS_MAF,xlab="MAF",ylab="density",freq=F,breaks=24,col=rgb(1,0,0,1/4), xlim=c(0,0.5),main="Vienna 2010")
hist(ita2011_fdr_af$B_MAF,xlab="MAF",ylab="density",freq=F,breaks=24,col=rgb(0,1,0,1/4), xlim=c(0,0.5),add=T)
legend("topright", c("CCR", "base"), fill=c(rgb(1,0,0,1/4),rgb(0,1,0,1/4)))

quartz()
hist(vie2010_fdr_af$AF_change,xlab="MAF",ylab="density",freq=F,breaks=24,col=rgb(1,0,0,1/4), xlim=c(-0.5,0.5),main="AF change CCR - base")
hist(ita2011_fdr_af$AF_change,xlab="MAF",ylab="density",freq=F,breaks=24,col=rgb(0,1,0,1/4), xlim=c(-0.5,0.5),add=T)
legend("topright", c("Vie", "Bolz"), fill=c(rgb(1,0,0,1/4),rgb(0,1,0,1/4)))



vie2010_fdr_af = read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25_fdr_0.05.cmhout.af",header=F)
ita2011_fdr_af = read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20_filt_mc8_chroms_only_mct8_mcv25_fdr_0.05.cmhout.af",header=F)
over_fdr_af = read.table("vie2010_ita2011_avgMAF_fdr0.05.af",header=F)
colnames(over_fdr_af) = maf_files
af_files=c("CHR","BPS","Alleles","CS1","CS2","CSI3","B1","B2","B3","P")
colnames(vie2010_fdr_af) = af_files
colnames(ita2011_fdr_af) = af_files
vie2010_fdr_af$MAF=apply(vie2010_fdr_af[,7:9],1,mean)
vie2010_fdr_af$MAF[vie2010_fdr_af$MAF > 0.5] = 1 - vie2010_fdr_af$MAF[vie2010_fdr_af$MAF > 0.5]
ita2011_fdr_af$MAF=apply(ita2011_fdr_af[,7:9],1,mean)
ita2011_fdr_af$MAF[ita2011_fdr_af$MAF > 0.5] = 1 - ita2011_fdr_af$MAF[ita2011_fdr_af$MAF > 0.5] 
keeps=c("CHR","BPS","Alleles","MAF")
vie2010_fdr_af=vie2010_fdr_af[,keeps]
ita2011_fdr_af=ita2011_fdr_af[,keeps]

# draw histograms of MAF
quartz()
hist(vie2010_fdr_af$MAF,xlab="MAF",ylab="density",freq=F,breaks=24,col=rgb(1,0,0,1/4), xlim=c(0,0.5),main="Vienna 2010")
hist(vie2010_maf$MAF,xlab="MAF",ylab="density",freq=F,breaks=24,col=rgb(0,1,0,1/4), xlim=c(0,0.5),add=T)
legend("topright", c("all SNPs", "candidate SNPs"), fill=c(rgb(0,1,0,1/4),rgb(1,0,0,1/4)))
dev.copy2pdf(file="maf_histo_vie2010.pdf")
quartz()
hist(ita2011_fdr_af$MAF,xlab="MAF",ylab="density",freq=F,breaks=24,col=rgb(1,0,0,1/4), xlim=c(0,0.5),main="Bolzano 2011")
hist(ita2011_maf$MAF,xlab="MAF",ylab="density",freq=F,breaks=24,col=rgb(0,1,0,1/4), xlim=c(0,0.5),add=T)
legend("topright", c("all SNPs", "candidate SNPs"), fill=c(rgb(0,1,0,1/4),rgb(1,0,0,1/4)))
dev.copy2pdf(file="maf_histo_ita2011.pdf")
quartz()
hist(over_fdr_af$MAF,xlab="MAF",ylab="density",freq=F,breaks=24,col=rgb(1,0,0,1/4), xlim=c(0,0.5),main="Overlap Vienna 2010 & Bolzano 2011")
hist(over_maf$MAF,xlab="MAF",ylab="density",freq=F,breaks=24,col=rgb(0,1,0,1/4), xlim=c(0,0.5),add=T)
legend("topright", c("all SNPs", "candidate SNPs"), fill=c(rgb(0,1,0,1/4),rgb(1,0,0,1/4)))
dev.copy2pdf(file="maf_histo_over.pdf")

vie2010_cat = read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.snpeff_cat.sync",header=T)
ita2011_cat = read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.snpeff_cats.sync",header=T)
over_cat = read.table("vie2010_bolz2011_all.BDGP5.72.snpeff_cats.sync",header=T)

# enrichment in inversions
# Odds ratio of number of candidate SNPs in inversions and random samples of chr & freq matched general SNPs in inversions : add inversion tag
inversions=read.table("/Volumes/Temp/Lukas/Inversions/inversion_breakpoints.bed",header=F)
colnames(inversions)=c("CHR","BPS1","BPS2","Inversion")
inversions$BPS1 = inversions$BPS1 + 1
inversions$BPS2 = inversions$BPS2 + 1

# cleaning categories
#vie2010
vie2010_cat$INTERGENIC=vie2010_cat$INTERGENIC | vie2010_cat$DOWNSTREAM | vie2010_cat$UPSTREAM
vie2010_cat$NS=vie2010_cat$NON_SYNONYMOUS_CODING | vie2010_cat$NON_SYNONYMOUS_START | vie2010_cat$START_LOST | vie2010_cat$STOP_LOST | vie2010_cat$START_LOST |vie2010_cat$STOP_GAINED | vie2010_cat$CODON_CHANGE
vie2010_cat$S=vie2010_cat$SYNONYMOUS_CODING |  vie2010_cat$SYNONYMOUS_STOP
vie2010_cat$UTR = vie2010_cat$UTR_3_PRIME  | vie2010_cat$UTR_5_PRIME
vie2010_cat$GENE = vie2010_cat$NS | vie2010_cat$S | vie2010_cat$UTR | vie2010_cat$EXON | vie2010_cat$INTRON | vie2010_cat$START_GAINED
vie2010_cat$EXON = vie2010_cat$NS | vie2010_cat$S | vie2010_cat$UTR | vie2010_cat$EXON | vie2010_cat$START_GAINED
vie2010_cat$INTERGENIC = vie2010_cat$INTERGENIC & ! vie2010_cat$GENE
vie2010_cat$CDS = vie2010_cat$EXON & ! vie2010_cat$UTR
colSums(vie2010_cat[,3:length(vie2010_cat)])
#ita2011
ita2011_cat$INTERGENIC=ita2011_cat$INTERGENIC | ita2011_cat$DOWNSTREAM | ita2011_cat$UPSTREAM
ita2011_cat$NS=ita2011_cat$NON_SYNONYMOUS_CODING | ita2011_cat$NON_SYNONYMOUS_START | ita2011_cat$START_LOST | ita2011_cat$STOP_LOST | ita2011_cat$START_LOST |ita2011_cat$STOP_GAINED | ita2011_cat$CODON_CHANGE
ita2011_cat$S=ita2011_cat$SYNONYMOUS_CODING |  ita2011_cat$SYNONYMOUS_STOP
ita2011_cat$UTR = ita2011_cat$UTR_3_PRIME  | ita2011_cat$UTR_5_PRIME
ita2011_cat$GENE = ita2011_cat$NS | ita2011_cat$S | ita2011_cat$UTR | ita2011_cat$EXON | ita2011_cat$INTRON | ita2011_cat$START_GAINED
ita2011_cat$EXON = ita2011_cat$NS | ita2011_cat$S | ita2011_cat$UTR | ita2011_cat$EXON | ita2011_cat$START_GAINED
ita2011_cat$INTERGENIC = ita2011_cat$INTERGENIC & ! ita2011_cat$GENE
ita2011_cat$CDS = ita2011_cat$EXON & ! ita2011_cat$UTR
#over
over_cat$INTERGENIC=over_cat$INTERGENIC | over_cat$DOWNSTREAM | over_cat$UPSTREAM
over_cat$NS=over_cat$NON_SYNONYMOUS_CODING | over_cat$NON_SYNONYMOUS_START | over_cat$START_LOST | over_cat$STOP_LOST | over_cat$START_LOST |over_cat$STOP_GAINED | over_cat$CODON_CHANGE
over_cat$S=over_cat$SYNONYMOUS_CODING |  over_cat$SYNONYMOUS_STOP
over_cat$UTR = over_cat$UTR_3_PRIME  | over_cat$UTR_5_PRIME
over_cat$GENE = over_cat$NS | over_cat$S | over_cat$UTR | over_cat$EXON | over_cat$INTRON | over_cat$START_GAINED
over_cat$EXON = over_cat$NS | over_cat$S | over_cat$UTR | over_cat$EXON | over_cat$START_GAINED
over_cat$INTERGENIC = over_cat$INTERGENIC & ! over_cat$GENE
over_cat$CDS = over_cat$EXON & ! over_cat$UTR
# adding inversion columns to tables
for(i in 1:length(inversions$CHR)){
vie2010_cat[,as.character(inversions$Inversion[i])] = vie2010_cat$CHR == as.character(inversions$CHR[i]) & vie2010_cat$BPS >= inversions$BPS1[i] & vie2010_cat$BPS <= inversions$BPS2[i]
ita2011_cat[,as.character(inversions$Inversion[i])] = ita2011_cat$CHR == as.character(inversions$CHR[i]) & ita2011_cat$BPS >= inversions$BPS1[i] & ita2011_cat$BPS <= inversions$BPS2[i]
over_cat[,as.character(inversions$Inversion[i])] = over_cat$CHR == as.character(inversions$CHR[i]) & over_cat$BPS >= inversions$BPS1[i] & over_cat$BPS <= inversions$BPS2[i]

}
colSums(ita2011_cat[,3:length(ita2011_cat)])
vie2010_joined = merge(vie2010_maf,vie2010_cat,by=c("CHR","BPS"))
ita2011_joined = merge(ita2011_maf,ita2011_cat,by=c("CHR","BPS"))
vie2010_fdr_joined= merge(vie2010_fdr_af,vie2010_cat,by=c("CHR","BPS"))
ita2011_fdr_joined = merge(ita2011_fdr_af,ita2011_cat,by=c("CHR","BPS"))
over_joined = merge(over_maf,over_cat,by=c("CHR","BPS"))
over_fdr_joined= merge(over_fdr_af,over_cat,by=c("CHR","BPS"))

# creating histograms for fdr files, for each chrom indiviually
breaks=seq(0,0.5,0.02)
vie2010_afbins_chrom = list()
ita2011_afbins_chrom = list()
over_afbins_chrom = list()
vie2010_fdr_hist_chrom = list()
ita2011_fdr_hist_chrom = list()
over_fdr_hist_chrom = list()
chroms=levels(vie2010_cat$CHR)
for(j in chroms){
	vie2010_fdr_hist_chrom[[j]]=hist(vie2010_fdr_af$MAF[vie2010_fdr_af$CHR == j],breaks=breaks,plot=F)$counts
	ita2011_fdr_hist_chrom[[j]]=hist(ita2011_fdr_af$MAF[ita2011_fdr_af$CHR == j],breaks=breaks,plot=F)$counts
	over_fdr_hist_chrom[[j]]=hist(over_fdr_af$MAF[over_fdr_af$CHR == j],breaks=breaks,plot=F)$counts
	# initialise lists for chromosomes
	vie2010_afbins_chrom[[j]]=list()
	ita2011_afbins_chrom[[j]]=list()
	over_afbins_chrom[[j]]=list()
	for(i in 1:(length(breaks)-1) ) {
			vie2010_afbins_chrom[[j]][i]=list(which(vie2010_joined$MAF < breaks[i+1] & vie2010_joined$MAF >= breaks[i]& vie2010_joined$CHR == j))
			ita2011_afbins_chrom[[j]][i]=list(which(ita2011_joined$MAF < breaks[i+1] & ita2011_joined$MAF >= breaks[i] & ita2011_joined$CHR == j))	
			over_afbins_chrom[[j]][i]=list(which(over_joined$MAF < breaks[i+1] & over_joined$MAF >= breaks[i] & over_joined$CHR == j))	
	}
}

# creating histograms for fdr files, over all SNPs not chromosome spec
breaks=seq(0,0.5,0.02)
vie2010_afbins = list()
ita2011_afbins = list()
for(i in 1:(length(breaks)-1) ) {
		vie2010_afbins[i]=list(which(vie2010_joined$MAF < breaks[i+1] & vie2010_joined$MAF >= breaks[i]))
		ita2011_afbins[i]=list(which(ita2011_joined$MAF < breaks[i+1] & ita2011_joined$MAF >= breaks[i]))	
}
vie2010_fdr_hist=hist(vie2010_fdr_af$MAF,breaks=breaks,plot=F)$counts
ita2011_fdr_hist=hist(ita2011_fdr_af$MAF,breaks=breaks,plot=F)$counts

# load the foreach library and set it up for parallelisation
library(parallel)
library(doParallel)
registerDoParallel(cores=4)
library(foreach) 
get_samples= function(counts,afbins){
	# takes a vector with counts and a list with indeces in af bins and returns a random sample of indeces with the same size and af spectrum
	vec=foreach(i=1:length(counts), .combine='c') %:% when(counts[i] > 0) %do% sample(afbins[[i]],counts[i])
	return(vec)	
}

get_samples_chrom= function(chrom_counts,afbins){
	# takes a vector with counts and a list with indeces in af bins and returns a random sample of indeces with the same size and af spectrum
	vec0=c()
	for(j in names(chrom_counts)){
		counts=chrom_counts[[j]]
	vec=foreach(i=1:length(counts), .combine='c') %:% when(counts[i] > 0) %do% sample(afbins[[j]][[i]],counts[i])
	vec0=c(vec0,vec)
	}
	return(vec0)	
}

iters=100000
categories=c("INTERGENIC","EXON","INTRON","S","NS","NON_SYNONYMOUS_CODING","SYNONYMOUS_CODING","CDS","UTR","UTR_3_PRIME","UTR_5_PRIME","GENE","IN(2L)t","IN(2R)NS","IN(3L)P","IN(3R)C","IN(3R)P","IN(3R)Mo")
# parallel foreach
vie_samples = foreach(i=1:iters, .combine=rbind) %dopar% {colSums(vie2010_joined[get_samples_chrom(vie2010_fdr_hist_chrom,vie2010_afbins_chrom),categories])}
ita_samples= foreach(i=1:iters, .combine=rbind ) %dopar% colSums(ita2011_joined[get_samples_chrom(ita2011_fdr_hist_chrom,ita2011_afbins_chrom),categories])
over_samples= foreach(i=1:iters, .combine=rbind ) %dopar% colSums(over_joined[get_samples_chrom(over_fdr_hist_chrom,over_afbins_chrom),categories])
#alternative
vie_samples=matrix(nrow=iters,ncol=length(categories),byrow=T)
colnames(vie_samples)=categories
ita_samples=matrix(nrow=iters,ncol=length(categories),byrow=T)
colnames(ita_samples)=categories
over_samples=matrix(nrow=iters,ncol=length(categories),byrow=T)
colnames(over_samples)=categories
for(i=1:iters){	
	vie_samples[i,]=colSums(vie2010_joined[get_samples_chrom(vie2010_fdr_hist_chrom,vie2010_afbins_chrom),categories])
	ita_samples[i,]=colSums(ita2011_joined[get_samples_chrom(ita2011_fdr_hist_chrom,ita2011_afbins_chrom),categories])
	over_samples[i,]=colSums(over_joined[get_samples_chrom(over_fdr_hist_chrom,over_afbins_chrom),categories])
} 

# getting the odds ratios and error bars
vie2010_cat_vals=data.frame(t(sapply(vie2010_fdr_joined[,categories],sum)))
#vie2010_cat_vals$GENES=vie2010_cat_vals$EXON + vie2010_cat_vals$INTRON
ita2011_cat_vals=data.frame(t(sapply(ita2011_fdr_joined[,categories],sum)))
#ita2011_cat_vals$GENES=ita2011_cat_vals$EXON + ita2011_cat_vals$INTRON
over_cat_vals=data.frame(t(sapply(over_fdr_joined[,categories],sum)))

vie2010_samples_quantiles=apply(vie_samples,2,quantile,probs=c(0.025,0.5,0.975))
vie2010_samples_quantiles=data.frame(vie2010_samples_quantiles)
#vie2010_samples_quantiles$GENES=vie2010_samples_quantiles$EXON + vie2010_samples_quantiles$INTRON
ita2011_samples_quantiles=apply(ita_samples,2,quantile,probs=c(0.025,0.5,0.975))
ita2011_samples_quantiles=data.frame(ita2011_samples_quantiles)
#ita2011_samples_quantiles$GENES=ita2011_samples_quantiles$EXON + ita2011_samples_quantiles$INTRON
over_samples_quantiles=apply(over_samples,2,quantile,probs=c(0.025,0.5,0.975))
over_samples_quantiles=data.frame(over_samples_quantiles)
vie2010_fdr_num=length(vie2010_fdr_joined[,1])
ita2011_fdr_num=length(ita2011_fdr_joined[,1])
over_fdr_num=length(over_fdr_joined[,1])
print_cat=c("INTERGENIC","GENE","CDS","INTRON","S","NS")
print_labs=c("Intergenic","Genes","CDS","Introns","Syn.","Nonsyn.")
vie2010_l2OR = foreach(i=1:3,.combine='rbind') %do%  log2( vie2010_cat_vals*(vie2010_fdr_num-vie2010_samples_quantiles[i,])/(vie2010_samples_quantiles[i,]*(vie2010_fdr_num-vie2010_cat_vals)  ) )
ita2011_l2OR = foreach(i=1:3,.combine='rbind') %do%  log2( ita2011_cat_vals*(ita2011_fdr_num-ita2011_samples_quantiles[i,])/(ita2011_samples_quantiles[i,]*(ita2011_fdr_num-ita2011_cat_vals)  ) )
over_l2OR = foreach(i=1:3,.combine='rbind') %do%  log2( over_cat_vals*(over_fdr_num-over_samples_quantiles[i,])/(over_samples_quantiles[i,]*(over_fdr_num-over_cat_vals)  ) )
library(Hmisc)
quartz()
errbar(1:length(print_cat) -0.15,vie2010_l2OR[2,print_cat],sapply(vie2010_l2OR[c(3,1),print_cat],min),sapply(vie2010_l2OR[c(3,1),print_cat],max),ylab=expression(paste("log2(OR)")),main="",xaxt="n",xlab="",pch=".",bg="white",ylim=c(min(min(vie2010_l2OR[,print_cat]),min(ita2011_l2OR[,print_cat])),max(max(vie2010_l2OR[,print_cat]),max(ita2011_l2OR[,print_cat]))),xlim=c(1-0.25,length(print_cat) + 0.25))
points(1:length(print_cat) -0.15,vie2010_l2OR[2,print_cat],pch=22,bg="white",cex=1.5)
errbar(1:length(print_cat) + 0.15,ita2011_l2OR[2,print_cat],sapply(ita2011_l2OR[c(3,1),print_cat],min),sapply(ita2011_l2OR[c(3,1),print_cat],max),ylab="log2 OR",main="",xaxt="n",yaxt="n",xlab="",add=T, pch=16,bg="white",cex=1.5)
axis(1, labels=F,at=1:length(print_cat))
text(1:length(print_cat), par("usr")[3] - 0.15, srt = 45, adj = 1, labels = print_labs, xpd = TRUE)
abline(h=0,lt=3)
legend("topleft",c("Vienna 2010","Bolzano 2011"),pch=c(22,16))
dev.copy2pdf(file="cat_enrichment_vie_ita_p95.pdf")
print_cat=c("INTERGENIC","GENE","CDS","INTRON","UTR_3_PRIME","UTR_5_PRIME","S","NS")
print_labs=c("Intergenic","Genes","CDS","Introns","3\' UTR","5\' UTR","Syn.","Nonsyn.")

#overlap
quartz()
errbar(1:length(print_cat),over_l2OR[2,print_cat],sapply(over_l2OR[c(3,1),print_cat],min),sapply(over_l2OR[c(3,1),print_cat],max),ylab=expression(paste("log2(OR)")),main="",xaxt="n",xlab="",pch=18,cex=1.5)
axis(1, labels=F,at=1:length(print_cat))
text(1:length(print_cat), par("usr")[3] - 0.15, srt = 45, adj = 1, labels = print_labs, xpd = TRUE)
abline(h=0,lt=3)
dev.copy2pdf(file="cat_enrichment_overlap_p95.pdf")


# enrichment in inversions
# Odds ratio of number of candidate SNPs in inversions and random samples of chr & freq matched general SNPs in inversions : add inversion tag
inversions=read.table("/Volumes/Temp/Lukas/Inversions/inversion_breakpoints.bed",header=F)
colnames(inversions)=c("CHR","BPS1","BPS2","Inversion")
inversions$BPS1 = inversions$BPS1 + 1
inversions$BPS2 = inversions$BPS2 + 1
quartz()
print_cat=c("IN.2L.t","IN.2R.NS","IN.3L.P","IN.3R.C","IN.3R.P","IN.3R.Mo")
print_labs=c( "IN(2L)t"  ,"IN(2R)NS" ,"IN(3L)P" , "IN(3R)C" , "IN(3R)Mo" ,"IN(3R)P")
errbar(1:length(print_cat) -0.15,vie2010_l2OR[2,print_cat],sapply(vie2010_l2OR[c(3,1),print_cat],min),sapply(vie2010_l2OR[c(3,1),print_cat],max),ylab=expression(paste("log2(OR)")),main="",xaxt="n",xlab="",pch=".",bg="white",ylim=c(min(min(vie2010_l2OR[,print_cat]),min(ita2011_l2OR[,print_cat])),max(max(vie2010_l2OR[,print_cat]),max(ita2011_l2OR[,print_cat]))),xlim=c(1-0.25,length(print_cat) + 0.25))
points(1:length(print_cat) -0.15,vie2010_l2OR[2,print_cat],pch=22,bg="white",cex=1.5)
errbar(1:length(print_cat) + 0.15,ita2011_l2OR[2,print_cat],sapply(ita2011_l2OR[c(3,1),print_cat],min),sapply(ita2011_l2OR[c(3,1),print_cat],max),ylab="log2 OR",main="",xaxt="n",yaxt="n",xlab="",add=T, pch=16,bg="white",cex=1.5)
axis(1, labels=F,at=1:length(print_cat))
text(1:length(print_cat), par("usr")[3] - 0.05, srt = 45, adj = 1, labels = print_labs, xpd = TRUE)
abline(h=0,lt=3)
legend("topleft",c("Vienna 2010","Bolzano 2011"),pch=c(22,16))
dev.copy2pdf(file="inversion_enrichment_vie_ita_p95.pdf")
quartz()
errbar(1:length(print_cat),over_l2OR[2,print_cat],sapply(over_l2OR[c(3,1),print_cat],min),sapply(over_l2OR[c(3,1),print_cat],max),ylab=expression(paste("log2(OR)")),main="",xaxt="n",xlab="",pch=18,cex=1.5)
axis(1, labels=F,at=1:length(print_cat))
text(1:length(print_cat), par("usr")[3] - 0.15, srt = 45, adj = 1, labels = print_labs, xpd = TRUE)
abline(h=0,lt=3)
dev.copy2pdf(file="inversion_enrichment_overlap_p95.pdf")

quartz()
print_cat=c("IN.2L.t","IN.2R.NS","IN.3L.P","IN.3R.C","IN.3R.P","IN.3R.Mo")
print_labs=c( "IN(2L)t"  ,"IN(2R)NS" ,"IN(3L)P" , "IN(3R)C" , "IN(3R)Mo" ,"IN(3R)P")
errbar(1:length(print_cat) -0.2,vie2010_l2OR[2,print_cat],sapply(vie2010_l2OR[c(3,1),print_cat],min),sapply(vie2010_l2OR[c(3,1),print_cat],max),ylab=expression(paste("log2(OR)")),main="",xaxt="n",xlab="",pch=".",bg="white",ylim=c(min(min(over_l2OR[,print_cat]),min(ita2011_l2OR[,print_cat])),max(max(over_l2OR[,print_cat]),max(ita2011_l2OR[,print_cat]))),xlim=c(1-0.25,length(print_cat) + 0.25))
points(1:length(print_cat) -0.2,vie2010_l2OR[2,print_cat],pch=16,bg="white",cex=1.5)
errbar(1:length(print_cat) + 0.0,ita2011_l2OR[2,print_cat],sapply(ita2011_l2OR[c(3,1),print_cat],min),sapply(ita2011_l2OR[c(3,1),print_cat],max),ylab="log2 OR",main="",xaxt="n",yaxt="n",xlab="",add=T, pch=22,bg="white",cex=1.5)
points(1:length(print_cat),ita2011_l2OR[2,print_cat],pch=22,bg="white",cex=1.5)
errbar(1:length(print_cat) + 0.2,over_l2OR[2,print_cat],sapply(over_l2OR[c(3,1),print_cat],min),sapply(over_l2OR[c(3,1),print_cat],max),ylab="log2 OR",main="",xaxt="n",yaxt="n",xlab="",add=T, pch=18,bg="black",cex=1.5)
axis(1, labels=F,at=1:length(print_cat))
text(1:length(print_cat), par("usr")[3] - 0.1, srt = 45, adj = 1, labels = print_labs, xpd = TRUE)
abline(h=0,lt=3)
legend("topleft",c("Vienna 2010","Bolzano 2011","Overlap"),pch=c(16,22,18))
dev.copy2pdf(file="inversion_enrichment_all_p95.pdf")

# distance to next snp
distances=by(vie2010_fdr_af$BPS,vie2010_fdr_af$CHR,diff)
quartz()
hist(log10(unlist(distances)),xlog=T,breaks=100,xlab="log10 distance",ylab="number of SNPs",main="Vienna sign. SNPs, distance to next")
dev.copy2pdf(file="vie2010_fdr_dist_next_snp.pdf")
distances=by(ita2011_fdr_af$BPS,ita2011_fdr_af$CHR,diff)
quartz()
hist(log10(unlist(distances)),xlog=T,breaks=100,xlab="log10 distance",ylab="number of SNPs",main="Bolzano sign. SNPs, distance to next")
dev.copy2pdf(file="ita2011_fdr_dist_next_snp.pdf")

distances=by(over_fdr_af$BPS,over_fdr_af$CHR,diff)
quartz()
hist(log10(unlist(distances)),xlog=T,breaks=100,xlab="log10 distance",ylab="number of SNPs",main="Overlap sign. SNPs, distance to next")
dev.copy2pdf(file="overlap_fdr_dist_next_snp.pdf")



