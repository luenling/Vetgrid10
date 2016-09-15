library(hexbin)
qqplot_hex <- function(nl_oP,nl_dP1,title,alpha=FALSE,beta=" by coverage"){	
	s_nl_oP <- sort(nl_oP,decreasing=F)
	s_nl_dP <- sort(nl_dP1,decreasing=F)
	q_nl_oP <- quantile(s_nl_oP,probs=seq(0,1,length.out=length(s_nl_oP)))
	q_nl_dP <- quantile(s_nl_dP,probs=seq(0,1,length.out=length(s_nl_oP)))
	quartz()
	bin <- hexbin(q_nl_dP,q_nl_oP,xbins=400,xlab=expression(paste("-log(P",scriptstyle(null),")")),ylab=expression(paste("-log(P",scriptstyle(obs),")")))
	if (alpha != FALSE){
		main_tit=bquote(paste(.(title),alpha,"=",.(alpha)," ,", beta,.(beta)))
		}
	else{
		main_tit=bquote(.(title))
	}
	pp<-plot(bin,legend=FALSE,style="constant.col",main=main_tit)
	
	hvp=hexViewport(bin)
	hexVP.abline(hvp,0,1,col="red",lty=2)
}
qqplot_hex <- function(obs_Pval,dist_Pval,title,alpha,beta=" by coverage"){
	nl_oP = -1 * log10(obs_Pval)
	nl_dP1 = -1 * log10(dist_Pval)
	s_nl_oP <- sort(nl_oP,decreasing=F)
	s_nl_dP <- sort(nl_dP1,decreasing=F)
	bin <- hexbin(s_nl_dP,s_nl_oP,xbins=400,xlab=expression(paste("-log(P",scriptstyle(null),")")),ylab=expression(paste("-log(P",scriptstyle(obs),")")))
	pp<-plot(bin,legend=FALSE,style="constant.col",main=bquote(paste(.(title),alpha,"=",.(alpha)," ", beta,.(beta))))
	hvp=hexViewport(bin)
	hexVP.abline(hvp,0,1,col="red",lty=2)
}

setwd("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_Poly/CMH/Strandbias")
vie2010_pv = scan("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_mct8_mcv25_stranded_q20_chroms_only.fisher_bias.pV")
quartz()
hist.data = hist(log10(vie2010_pv),freq=F,plot=F)
plot(hist.data$mids, hist.data$counts,log="y")
dist_Pval2=scan("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.strandbias_P_gt_0.001.cmhout_a_50.0_beta_r1_71.67_r2_114.74_r3_192.86_nullP.out")
dist_Pval3=scan("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.strandbias_P_gt_0.001.cmhout_a_50.0_beta_r1_71.67_r2_114.74_r3_192.86_nullP.out10x")
obs_Pval=scan("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.strandbias_P_gt_0.001.pV")
dist_Pval=scan("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.strandbias_P_gt_0.001.cmhout_a_45.0_beta_r1_64.50_r2_103.27_r3_173.57_nullP.out")
nl_oP = -1 * log10(obs_Pval)
nl_dP1 = -1 * log10(dist_Pval)
nl_dP2 = -1 * log10(dist_Pval2)
nl_dP3 = -1 * log10(dist_Pval3)

qqplot_hex(nl_oP,nl_dP1,title="Vienna 2010, Strandbias Fisher pV < 0.001",alpha="45.0")
qqplot_hex(nl_oP,nl_dP2,title="Vienna 2010, Strandbias pV < 0.001",alpha="50.0")
qqplot_hex(nl_oP,nl_dP3,title="Vienna 2010, Strandbias pV < 0.001",alpha="50.0")

# add fdr
obs_pval=sort(obs_Pval)
n_ratio=length(obs_pval)/length(dist_Pval3)

ranks = 1:25
for (r in ranks)
	{
	ps20_smaller =  n_ratio*sum(dist_Pval3 < obs_pval[r])
	if (ps20_smaller/r < 0.1){
	cat ("r=", r, obs_pval[r], ps20_smaller/r, "\n")}
# not ps20_smaller/(r + ps20_smaller) as r already contains the FP
	}

allinfo.2L = read.table("2L_snps_allpops.pysamtest",header=F, na.strings = "nan")

colnames(allinfo.2L)=c("CHR","BPS","alleles","SB1","SB2","SB_P","mPos_A","sPos_A","mPos_B","sPos_B","SC","InsA","InsB")
allinfo.2L$nlSB_P = -1 * log10(allinfo.2L$SB_P)
quartz()
boxplot(allinfo.2L[allinfo.2L$SB_P < 1e-3,c("SB1","SB2","nlSB_P","mPos_A","mPos_B")])
boxplot(c(allinfo.2L[allinfo.2L$SB_P < 1e-3,c("mPos_A","mPos_B")],allinfo.2L[allinfo.2L$SB_P >= 1e-3,c("mPos_A","mPos_B")]))

vie_cmh_fsb = read.table("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.cmhout.fisher_bias_pv_only",header=F)
colnames(vie_cmh_fsb)=c("CHR","BPS","alleles","P","FP")
vie_cmh_fsb$nlP=-1*log10(vie_cmh_fsb$P)
vie_cmh_fsb$nlFP=-1*log10(vie_cmh_fsb$FP)
quartz()
plot(vie_cmh_fsb$nlFP,vie_cmh_fsb$nlP)
plot(vie_cmh_fsb$FP,vie_cmh_fsb$P)
cor(vie_cmh_fsb$nlFP,vie_cmh_fsb$nlP)
library(Hmisc)
rcorr(as.matrix(vie_cmh_fsb[,c("nlP","nlFP")]), type="spearman")
rcorr(as.matrix(vie_cmh_fsb[,c("nlP","nlFP")]), type="pearson")
rcorr(as.matrix(vie_cmh_fsb[,c("P","FP")]), type="spearman")
rcorr(as.matrix(vie_cmh_fsb[,c("P","FP")]), type="pearson")

shuffled_sb=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/shuffles/vie2010_shuffles_sb.gwas",header=F, skip = 1)
shuffled_nosb=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/shuffles/vie2010_shuffles_nosb.gwas",header=F, skip = 1)
shuffled=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/shuffles/fem_all_shuf_dgrpandafrorfr_FB_SB.gwas",header=F, skip = 1)
colnames(shuffled)=c("CHR","BPS","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","DGRP","FB","SB")
shuffled$Bad = (shuffled$FB <= 1e-3) & (shuffled$FB <= 1e-2)

hist(shuffled$S6[shuffled$Bad])
hist(shuffled$S6[! shuffled$Bad])

colnames(shuffled_nosb)=c("CHR","BPS","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12")
colnames(shuffled_sb)=c("CHR","BPS","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12")
colSums(shuffled_sb[,3:length(shuffled_sb[1,])] < 1e-5)/length(shuffled_sb[,1])
colSums(shuffled_nosb[,3:length(shuffled_sb[1,])] < 1e-5)/length(shuffled_nosb[,1])
colSums(shuffled_sb[,3:length(shuffled_sb[1,])] < 1e-5)/length(shuffled_sb[,1])/(colSums(shuffled_nosb[,3:length(shuffled_sb[1,])] < 1e-5)/length(shuffled_nosb[,1]))
colSums(shuffled_sb[,3:length(shuffled_sb[1,])] < 1e-8)/length(shuffled_sb[,1])/(colSums(shuffled_nosb[,3:length(shuffled_sb[1,])] < 1e-8)/length(shuffled_nosb[,1]))
quartz()
 hist(shuffled$S2[shuffled$Bad])
quartz()
 hist(shuffled$S2[! shuffled$Bad])
shuffled_sb$SB=1
shuffled_nosb$SB=0
shuffled=rbind(shuffled_sb,shuffled_nosb)
shuffled$SB=factor(shuffled$SB)
wilcox.test(shuffled~shuffled$B, alternative="greater")
apply(shuffled[,3:(length(shuffled)-4)],2,function(x) wilcox.test(x~shuffled$Bad, alternative="greater"))
apply(shuffled[,3:(length(shuffled)-4)],2,function(x) c(length(x[x < 1e-8 & shuffled$Bad]), length(x[x < 1e-8 & ! shuffled$Bad] )))
# load snp file Vie 2010
snp_vie=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/vie2010_gatk_FB_SB_freebayes_afr_dgrpid_union.dat",header=T)
snp_vie$nlFB=-1*log10(snp_vie$FB)
library(ggplot2)

d <- ggplot(data=snp_vie, aes(FS, SB)) 
d+stat_binhex(bins = 75)
p1 = ggplot(data=snp_vie, color="SNPs") +
     stat_binhex(data=snp_vie[snp_vie$DGRP == 1 & snp_vie$AFR_FR_D == 0 ,],aes(x=FS, y=SB,alpha=0.5), fill="blue",bins=75) +
     stat_binhex(data=snp_vie[snp_vie$D_A_F_U == 0,],aes(x=FS, y=SB,alpha=0.5), fill="red", bins=75) +
     stat_binhex(data=snp_vie[snp_vie$AFR_FR_D == 1,],aes(x=FS, y=SB,alpha=0.5), fill="green", bins=75)
p1    
p1 = ggplot(data=snp_vie) +
     stat_binhex(data=snp_vie[snp_vie$DGRP == 1 & snp_vie$AFR_FR_D == 0 ,],aes(x=FS, y=SB,alpha=0.5, color="DGRP"), fill="blue",bins=75) +
     stat_binhex(data=snp_vie[snp_vie$D_A_F_U == 0,],aes(x=FS, y=SB,alpha=0.5, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
     stat_binhex(data=snp_vie[snp_vie$AFR_FR_D == 1,],aes(x=FS, y=SB,alpha=0.5, color="DGRP & Afr or Fr"), fill="green", bins=75) + guides(fill=FALSE,alpha=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
p1
ggsave(filename="fs_sb.pdf")
quartz()
p1 = ggplot(data=snp_vie) +
     stat_binhex(data=snp_vie[snp_vie$DGRP == 1 & snp_vie$AFR_FR_D == 0 ,],aes(x=nlFB, y=SB,alpha=0.5, color="DGRP"), fill="blue",bins=75) +
     stat_binhex(data=snp_vie[snp_vie$D_A_F_U == 0,],aes(x=nlFB, y=SB,alpha=0.5, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
     stat_binhex(data=snp_vie[snp_vie$AFR_FR_D == 1,],aes(x=nlFB, y=SB,alpha=0.5, color="DGRP & Afr or Fr"), fill="green", bins=75) + guides(fill=FALSE,alpha=F,color=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
p1
ggsave(filename="nlfb_sb.pdf")
quartz()
p1 = ggplot(data=snp_vie) +
     stat_binhex(data=snp_vie[snp_vie$DGRP == 1 & snp_vie$AFR_FR_D == 0 ,],aes(x=nlFB, y=ReadPosRankSum,alpha=0.5, color="DGRP"), fill="blue",bins=75) +
     stat_binhex(data=snp_vie[snp_vie$D_A_F_U == 0,],aes(x=nlFB, y=ReadPosRankSum,alpha=0.5, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
     stat_binhex(data=snp_vie[snp_vie$AFR_FR_D == 1,],aes(x=nlFB, y=ReadPosRankSum,alpha=0.5, color="DGRP & Afr or Fr"), fill="green", bins=75) + guides(fill=FALSE,alpha=F,color=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
p1
ggsave(filename="nlfb_readposranksum.pdf")
quartz()
ggplot(data=snp_vie) +
     stat_binhex(data=snp_vie[snp_vie$DGRP == 1 & snp_vie$AFR_FR_D == 0 ,],aes(x=SB, y=ReadPosRankSum,alpha=0.5, color="DGRP"), fill="blue",bins=75) +
     stat_binhex(data=snp_vie[snp_vie$D_A_F_U == 0,],aes(x=SB, y=ReadPosRankSum,alpha=0.5, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
     stat_binhex(data=snp_vie[snp_vie$AFR_FR_D == 1,],aes(x=SB, y=ReadPosRankSum,alpha=0.5, color="DGRP & Afr or Fr"), fill="green", bins=75) + guides(fill=FALSE,alpha=F,color=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
ggsave(filename="sb_readposranksum.pdf")
quartz()
ggplot(data=snp_vie) +
     stat_binhex(data=snp_vie[snp_vie$DGRP == 1 & snp_vie$AFR_FR_D == 0 ,],aes(x=MQRankSum, y=ReadPosRankSum,alpha=0.5, color="DGRP"), fill="blue",bins=75) +
     stat_binhex(data=snp_vie[snp_vie$D_A_F_U == 0,],aes(x=MQRankSum, y=ReadPosRankSum,alpha=0.5, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
     stat_binhex(data=snp_vie[snp_vie$AFR_FR_D == 1,],aes(x=MQRankSum, y=ReadPosRankSum,alpha=0.5, color="DGRP & Afr or Fr"), fill="green", bins=75) + guides(fill=FALSE,alpha=F,color=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
ggsave(filename="MQRankSum_readposranksum.pdf")
quartz()
ggplot(data=snp_vie) +
     stat_binhex(data=snp_vie[snp_vie$DGRP == 1 & snp_vie$AFR_FR_D == 0 ,],aes(x=MQ, y=ReadPosRankSum,alpha=0.5, color="DGRP"), fill="blue",bins=75) +
     stat_binhex(data=snp_vie[snp_vie$D_A_F_U == 0,],aes(x=MQ, y=ReadPosRankSum,alpha=0.5, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
     stat_binhex(data=snp_vie[snp_vie$AFR_FR_D == 1,],aes(x=MQ, y=ReadPosRankSum,alpha=0.5, color="DGRP & Afr or Fr"), fill="green", bins=75) + guides(fill=FALSE,alpha=F,color=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
ggsave(filename="MQ_readposranksum.pdf")
quartz()
ggplot(data=snp_vie) +
     stat_binhex(data=snp_vie[snp_vie$DGRP == 1 & snp_vie$AFR_FR_D == 0 ,],aes(x=EPP, y=ReadPosRankSum,alpha=0.5, color="DGRP"), fill="blue",bins=75) +
     stat_binhex(data=snp_vie[snp_vie$D_A_F_U == 0,],aes(x=EPP, y=ReadPosRankSum,alpha=0.5, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
     stat_binhex(data=snp_vie[snp_vie$AFR_FR_D == 1,],aes(x=EPP, y=ReadPosRankSum,alpha=0.5, color="DGRP & Afr or Fr"), fill="green", bins=75) + guides(fill=FALSE,alpha=F,color=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
ggsave(filename="EPP_readposranksum.pdf")
quartz()
ggplot(data=snp_vie) +
     stat_binhex(data=snp_vie[snp_vie$DGRP == 1 & snp_vie$AFR_FR_D == 0 ,],aes(x=RPP, y=ReadPosRankSum,alpha=0.5, color="DGRP"), fill="blue",bins=75) +
     stat_binhex(data=snp_vie[snp_vie$D_A_F_U == 0,],aes(x=RPP, y=ReadPosRankSum,alpha=0.5, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
     stat_binhex(data=snp_vie[snp_vie$AFR_FR_D == 1,],aes(x=RPP, y=ReadPosRankSum,alpha=0.5, color="DGRP & Afr or Fr"), fill="green", bins=75) + guides(fill=FALSE,alpha=F,color=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
ggsave(filename="RPP_readposranksum.pdf")

quartz()
ggplot(data=snp_vie) +
     stat_binhex(data=snp_vie[snp_vie$DGRP == 1 & snp_vie$AFR_FR_D == 0 ,],aes(x=MQM, y=ReadPosRankSum,alpha=0.5, color="DGRP"), fill="blue",bins=75) +
     stat_binhex(data=snp_vie[snp_vie$D_A_F_U == 0,],aes(x=MQM, y=ReadPosRankSum,alpha=0.5, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
     stat_binhex(data=snp_vie[snp_vie$AFR_FR_D == 1,],aes(x=MQM, y=ReadPosRankSum,alpha=0.5, color="DGRP & Afr or Fr"), fill="green", bins=75) + guides(fill=FALSE,alpha=F,color=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
ggsave(filename="MQM_readposranksum.pdf")
quartz()
ggplot(data=snp_vie) +
     stat_binhex(data=snp_vie[snp_vie$DGRP == 1 & snp_vie$AFR_FR_D == 0 ,],aes(x=FS, y=RPP,alpha=0.5, color="DGRP"), fill="blue",bins=75) +
     stat_binhex(data=snp_vie[snp_vie$D_A_F_U == 0,],aes(x=FS, y=RPP,alpha=0.5, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
     stat_binhex(data=snp_vie[snp_vie$AFR_FR_D == 1,],aes(x=FS, y=RPP,alpha=0.5, color="DGRP & Afr or Fr"), fill="green", bins=75) + guides(fill=FALSE,alpha=F,color=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
ggsave(filename="FS_RPP.pdf")
quartz()
all_snps=1668035
p1=ggplot(data=snp_vie) +
     stat_binhex(data=snp_vie[snp_vie$DGRP == 1 & snp_vie$AFR_FR_D == 0 ,],aes(x=FS, y=RPP,alpha=0.1+0.7*..count../(..count..+1.0), color="DGRP"),fill="blue", bins=75)  +
     stat_binhex(data=snp_vie[snp_vie$D_A_F_U == 0,],aes(x=FS, y=RPP,alpha=0.7*..count../(..count..+1.0)+0.1, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
          stat_binhex(data=snp_vie[snp_vie$AFR_FR_D == 1,],aes(x=FS, y=RPP,alpha=0.1+0.7*..count../(..count..+1.0), color="DGRP & Afr or Frn"), fill="green", bins=75) + 
          guides(fill=FALSE,alpha=F,color=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
p1

library(ggplot2)

setwd("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/")
vie_snps=read.table("all_snps_samtools_biases_self_calc_ident_combined_chroms.tab",na.strings = "NA")
colnames(vie_snps)=c("CHR","BPS","ALL","FS","SB","RTD","ATD","TDB","RPB","DD_i","DD_c")

quartz()
ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 1 & vie_snps$DD_i == 0,],aes(x=FS, y=TDB,alpha=0.5, color="DGRP"), fill="blue",bins=75) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=FS, y=TDB,alpha=0.5, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=FS, y=TDB,alpha=0.5, color="DGRP & Afr or Fr"), fill="green", bins=75) + guides(fill=FALSE,alpha=F,color=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
ggsave(filename="FS_TDB.pdf")

quartz()
ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 1 & vie_snps$DD_i == 0,],aes(x=FS, y=RPB,alpha=0.5, color="DGRP"), fill="blue",bins=75) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=FS, y=RPB,alpha=0.5, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=FS, y=RPB,alpha=0.5, color="DGRP & Afr or Fr"), fill="green", bins=75) + guides(fill=FALSE,alpha=F,color=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
ggsave(filename="FS_RPB.pdf")

quartz()
ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 1 & vie_snps$DD_i == 0,],aes(x=FS, y=ATD,alpha=0.5, color="DGRP"), fill="blue",bins=75) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=FS, y=ATD,alpha=0.5, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=FS, y=ATD,alpha=0.5, color="DGRP & Afr or Fr"), fill="green", bins=75) + guides(fill=FALSE,alpha=F,color=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
ggsave(filename="FS_ATD.pdf")


quartz()
ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 1 & vie_snps$DD_i == 0,],aes(x=TDB, y=ATD,alpha=0.5, color="DGRP"), fill="blue",bins=75) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=TDB, y=ATD,alpha=0.5, color="Not Afr,Fr or DGRP"), fill="red", bins=75) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=TDB, y=ATD,alpha=0.5, color="DGRP & Afr or Fr"), fill="green", bins=75) + guides(fill=FALSE,alpha=F,color=F) + scale_color_manual(values=c("DGRP"="blue", "Not Afr,Fr or DGRP"="red", "DGRP & Afr or Fr"="green"))
ggsave(filename="TDB_ATD.pdf")
library(gridExtra)
mybreaks=c(10,100,1000,10000,100000,1000000)
quartz()
good=ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=TDB, y=ATD),binwidth=c(3,1)) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=TDB, y=ATD),binwidth=c(3,1)) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)     
dev.copy2pdf(file="TDB_ATD.pdf")

quartz()
good=ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=FS, y=RPB),bins=50) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=FS, y=RPB),bins=50) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)     


