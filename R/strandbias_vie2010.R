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
shuffled$Bad = (shuffled$FB <= 1e-3) 

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

library(hexbin)
qqplot_hex <- function(oP,dP1,title){
        nl_oP <- -1*log10(oP)
        nl_dP1 <- -1*log10(dP1)
        q_nl_oP <- oP
        q_nl_dP <- oP
        if (is.data.frame(oP)){
            for(i in 1:length(names(oP))){
                s_nl_oP <- sort(unlist(nl_oP[i]),decreasing=F)
                s_nl_dP <- sort(unlist(nl_dP1[i]),decreasing=F)
                q_nl_oP[i] <- quantile(s_nl_oP,probs=seq(0,1,length.out=length(s_nl_oP)))
                q_nl_dP[i] <- quantile(s_nl_dP,probs=seq(0,1,length.out=length(s_nl_oP)))
            }

        }
        else {
            s_nl_oP <- sort(nl_oP,decreasing=F)
            s_nl_dP <- sort(nl_dP1,decreasing=F)
            q_nl_oP <- quantile(s_nl_oP,probs=seq(0,1,length.out=length(s_nl_oP)))
            q_nl_dP <- quantile(s_nl_dP,probs=seq(0,1,length.out=length(s_nl_oP)))
        }
        bin <- hexbin(as.matrix(q_nl_dP),as.matrix(q_nl_oP),xbins=400,xlab=expression(paste("-log(P",scriptstyle(noFSB),")")),ylab=expression(paste("-log(P",scriptstyle(FSB),")")))

        main_tit=bquote(.(title))
	pp<-plot(bin,legend=FALSE,style="constant.col",main=main_tit)
	hvp=hexViewport(bin)
	hexVP.abline(hvp,0,1,col="red",lty=2)
}
x11()
qqplot_hex(shuffled$S1[shuffled$Bad],shuffled$S1[! shuffled$Bad],"Shuffle 1")
qqplot_hex(shuffled[shuffled$Bad,3:14],shuffled[! shuffled$Bad,3:14],"Shuffles")
dev.copy2pdf(file="12shuffles_qqplot.pdf")
shuffled[,3:14]=-1*log10(shuffled[,3:14])

boxplot(shuffled[shuffled$Bad,3:14], at=1:length(names(shuffled[,3:14]))-0.15,boxwex=0.3, col="red",ylim=c(0,max(shuffled[,3:14])),xlim=c(0,length(names(shuffled[,3:14]))),xaxt="n",ylab=expression(paste("-log(P)")),outpchar=".",outcex=2.0)
boxplot(shuffled[! shuffled$Bad,3:14], at=1:length(names(shuffled[,3:14]))+0.15,boxwex=0.3, col="green", yaxt="n",xaxt="n",outpchar=".",outcex=2.0,add=T)
axis(1,at=1:length(names(shuffled[,3:14])),cex.axis=1,labels=names(shuffled[,3:14]))
legend("topleft",legend=c("strandbias","no strandbias"),fill=c("red","green"),cex=1.0)
dev.print(png,width=600,file="12shuffles_boxplot.png")
dev.copy2pdf(file="12shuffles_boxplot.pdf")

x11()



# load the raw/af/pvalue/bd/wool test file
all.vie=read.table("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.afs.PV.BD.FSB.rawFSB",header=F,na.strings = c("NA","N","NaN","nan","na"))
colnames(all.vie)  =c("CHR","BPS","REF","AF1","AF2","AF3","AF4","AF5","AF6","P","BD","WOOLF","FSB","rFSB")

length(all.vie$P[all.vie$FSB < 1e-3 & all.vie$rFSB < 1e-3]) #33172
length(all.vie$P[all.vie$FSB < 1e-3 & all.vie$rFSB > 1e-3]) #6440
length(all.vie$P[all.vie$FSB < 1e-3 & all.vie$rFSB > 1e-3 & all.vie$P < 1e-3]) #39
length(all.vie$P[all.vie$FSB < 1e-3 & all.vie$rFSB > 1e-3 & all.vie$P < 1e-8]) #39
 
boxplot(-1*log10(all.vie$P[all.vie$FSB < 1e-3 & all.vie$rFSB < 1e-3]),at=1, xlim=c(0.75,2.25),col="green",xaxt="n",ylab=expression(paste("-log(P)")))
boxplot(-1*log10(all.vie$P[all.vie$FSB < 1e-3 & all.vie$rFSB > 1e-3]),at=2,col="red",xaxt="n",yaxt="n",add=T)
legend("topright",legend=c("in raw and processed","only in processed"),title=c("FSB < 1e-3"),fill=c("green","red"),cex=1.0)
dev.copy2pdf(file="raw_processed_boxplot.pdf")

cor.test(all.vie$"FSB",all.vie$"rFSB", method="spearman")
cor.test(all.vie$"FSB",all.vie$"rFSB", method="pearson")

# allele freq variance / mean


boxplot(-log10(all.vie$FSB)/-log10(all.vie$rFSB),ylim=c(0,3))
## all.vie$v1=apply(all.vie[,4:6],1,var)
## all.vie$m1=rowMeans(all.vie[,4:6])
## all.vie$v2=apply(all.vie[,7:9],1,var)
## all.vie$m2=rowMeans(all.vie[,7:9])
## all.vie$maf1=all.vie$m1
## all.vie$maf1[all.vie$maf1 > 0.5]=1-all.vie$maf1[all.vie$maf1 > 0.5]
## all.vie$maf1[all.vie$maf1 < 0.01]=1
## all.vie$maf2=all.vie$m2
## all.vie$maf2[all.vie$maf2 > 0.5]=1-all.vie$maf2[all.vie$maf2 > 0.5]
## all.vie$maf2[all.vie$maf2 < 0.01]=1
## all.vie$sv1=sqrt(all.vie$v1)/all.vie$maf1
## all.vie$sv2=sqrt(all.vie$v2)/all.vie$maf2
str(all.vie)
all.vie[,15:20]=all.vie[,4:9]
names(all.vie)[15:20]=c("AS1","AS2","AS3","AS4","AS5","AS6")
all.vie[,15:20]=asin(sqrt(all.vie[,15:20]))
all.vie=all.vie[,1:20]
all.vie$v1=apply(all.vie[,15:17],1,var)
all.vie$v2=apply(all.vie[,18:20],1,var)

boxplot(all.vie$v2[all.vie$FSB < 1e-3 ],at=1, xlim=c(0.75,2.25),col="red",xaxt="n",ylab=expression(paste("SD/maf")))
boxplot(all.vie$v2[all.vie$FSB > 1e-3],at=2,col="green",xaxt="n",yaxt="n",add=T)
mean(all.vie$v1[all.vie$FSB < 1e-3 ])#0.008902922
mean(all.vie$v1[all.vie$FSB > 1e-3 ])#0.007532555
mean(all.vie$v2[all.vie$FSB < 1e-3 ])#0.004173866
mean(all.vie$v2[all.vie$FSB > 1e-3 ])#0.002946423
all.vie[,4:9]=1-all.vie[,4:9]
all.vie$Sall=apply(all.vie[,4:9],1,function(x) {sum(x > 0.01)})
all.vie$Sc=apply(all.vie[,4:6],1,function(x) {sum(x > 0.01)})
all.vie$Sb=apply(all.vie[,7:9],1,function(x) {sum(x > 0.01)})
all.vie$iSB=rep("N",length(all.vie$FSB))
all.vie$iSB[all.vie$FSB < 1e-3] = "SB"
all.vie$iSB=as.factor(all.vie$iSB)



boxplot(all.vie$Sall[all.vie$FSB < 1e-3 ],at=1, xlim=c(0.75,2.25),col="red",xaxt="n",ylab=expression(paste("SD/maf")))
boxplot(all.vie$Sall[all.vie$FSB > 1e-3],at=2,col="green",xaxt="n",yaxt="n",add=T)
legend("topright",legend=c("in raw and processed","only in processed"),title=c("FSB < 1e-3"),fill=c("green","red"),cex=1.0)
Sall.tab=table(all.vie$Sall,all.vie$iSB)
Sc.tab=table(all.vie$Sc,all.vie$iSB)
Sb.tab=table(all.vie$Sb,all.vie$iSB)
prop.table(Sall.tab,2)
prop.table(Sc.tab,2)
prop.table(Sb.tab,2)
barplot(t(prop.table(Sall.tab,2)),beside=TRUE,col=c("green","red"))
title(main="All Samples")
mtext("samples with MAF > 0.01",side=1,line=2)
mtext("frequency",side=2,line=2.5)
legend("topleft",legend=c("no SB","SB (FSB<1e-3)"),fill=c("green","red"),cex=1.0)
dev.copy2pdf(file="samplecountall.pdf")
dev.off()
x11()
barplot(t(prop.table(Sc.tab,2)),beside=TRUE,col=c("green","red"))
title(main="CCR Samples")
mtext("samples with MAF > 0.01",side=1,line=2)
mtext("frequency",side=2,line=2.5)
legend("topleft",legend=c("no SB","SB (FSB<1e-3)"),fill=c("green","red"),cex=1.0)
dev.copy2pdf(file="samplecountccr.pdf")
barplot(t(prop.table(Sb.tab,2)),beside=TRUE,col=c("green","red"))
title(main="Base Samples")
mtext("samples with MAF > 0.01",side=1,line=2)
mtext("frequency",side=2,line=2.5)
legend("topleft",legend=c("no SB","SB (FSB<1e-3)"),fill=c("green","red"),cex=1.0)
dev.copy2pdf(file="samplecountbase.pdf")
Sallp.tab=table(all.vie$Sall[all.vie$P < 1e-8],all.vie$iSB[all.vie$P < 1e-8])
Scp.tab=table(all.vie$Sc[all.vie$P < 1e-8],all.vie$iSB[all.vie$P < 1e-8])
Sbp.tab=table(all.vie$Sb[all.vie$P < 1e-8],all.vie$iSB[all.vie$P < 1e-8])

barplot(t(prop.table(Sallp.tab,2)),beside=TRUE,col=c("green","red"))
title(main="All Samples")
mtext("samples with MAF > 0.01",side=1,line=2)
mtext("frequency",side=2,line=2.5)
legend("topleft",legend=c("no SB","SB (FSB<1e-3)"),fill=c("green","red"),cex=1.0)
dev.copy2pdf(file="samplecountallp.pdf")
dev.off()
x11()
barplot(t(prop.table(Scp.tab,2)),beside=TRUE,col=c("green","red"))
title(main="CCR Samples")
mtext("samples with MAF > 0.01",side=1,line=2)
mtext("frequency",side=2,line=2.5)
legend("topleft",legend=c("no SB","SB (FSB<1e-3)"),fill=c("green","red"),cex=1.0)
dev.copy2pdf(file="samplecountccrp.pdf")
barplot(t(prop.table(Sbp.tab,2)),beside=TRUE,col=c("green","red"))
title(main="Base Samples")
mtext("samples with MAF > 0.01",side=1,line=2)
mtext("frequency",side=2,line=2.5)
legend("topleft",legend=c("no SB","SB (FSB<1e-3)"),fill=c("green","red"),cex=1.0)
dev.copy2pdf(file="samplecountbasep.pdf")



> prop.table(Sall.tab,2)
   
              N          SB
  1 0.002169674 0.001135532
  2 0.022389883 0.021221832
  3 0.066362315 0.103838098
  4 0.099777384 0.139771380
  5 0.139587283 0.162507255
  6 0.669713461 0.571525903
> prop.table(Sc.tab,2)
   
             N         SB
  0 0.01761525 0.01062353
  1 0.08442636 0.03348558
  2 0.18676757 0.15375104
  3 0.71119083 0.80213985
> prop.table(Sb.tab,2)
   
              N          SB
  0 0.003166103 0.073733882
  1 0.043931297 0.110171844
  2 0.132809967 0.149991168
  3 0.820092632 0.666103106

wilcox.test(all.vie$v2[all.vie$FSB < 1e-3 ],all.vie$sv2[all.vie$FSB > 1e-3], alternative="greater")
# p-value < 2.2e-16
wilcox.test(all.vie$v1[all.vie$FSB < 1e-3 ],all.vie$sv1[all.vie$FSB > 1e-3], alternative="greater")
# p-value < 2.2e-16

length(all.vie$v1[all.vie$FSB < 1e-3 & all.vie$WOOLF < 0.05])#3569
length(all.vie$v1[all.vie$FSB > 1e-3 & all.vie$WOOLF < 0.05])#133696
length(all.vie$v1[all.vie$FSB < 1e-3])#39629
length(all.vie$v1[all.vie$FSB > 1e-3])#1628813
length(all.vie$v1[all.vie$FSB < 1e-3 & all.vie$WOOLF < 0.05])/length(all.vie$v1[all.vie$FSB < 1e-3])#0.06669358
length(all.vie$v1[all.vie$FSB > 1e-3 & all.vie$WOOLF < 0.05])/length(all.vie$v1[all.vie$FSB > 1e-3])#0.06051769







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
colnames(vie_snps)=c("CHR","BPS","ALL","FS","SB","BSB","RTD","ATD","TDB","RPB","DD_i","DD_c")
vie_snps$aRPB=abs(vie_snps$RPB)
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
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=FS, y=ATD),binwidth=c(10,1)) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=FS, y=ATD),binwidth=c(10,1)) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)     
dev.copy2pdf(file="FS_ATD_double.pdf")

quartz()
good=ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=aRPB, y=ATD),binwidth=c(0.5,1)) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=aRPB, y=ATD),binwidth=c(0.5,1)) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)     
dev.copy2pdf(file="aRPB_ATD_double.pdf")

xco="FS" ; yco="TDB"
xlimit=c(min(vie_snps[,xco]),max(vie_snps[,xco]))
ylimit=c(min(vie_snps[,yco]),max(vie_snps[,yco]))
xlimit=c(0,600)
quartz()
good=ggplot(data=vie_snps) + xlim(xlimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=FS, y=TDB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) + xlim(xlimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=FS, y=TDB),bins=40) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)  
dev.copy2pdf(file="FS_TDB_double.pdf")

quartz()
good=ggplot(data=vie_snps) + xlim(xlimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=FS, y=RPB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) + xlim(xlimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=FS, y=RPB),bins=40) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)  
dev.copy2pdf(file="FS_RPB_double.pdf")

quartz()
good=ggplot(data=vie_snps) + xlim(xlimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=FS, y=aRPB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) + xlim(xlimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=FS, y=aRPB),bins=40) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)  
dev.copy2pdf(file="FS_absRPB_double.pdf")

quartz()
good=ggplot(data=vie_snps) + xlim(xlimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=FS, y=SB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) + xlim(xlimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=FS, y=SB),bins=40) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)  
dev.copy2pdf(file="FS_SB_double.pdf")
ylimit=c(0,0.05)
quartz()
good=ggplot(data=vie_snps) + xlim(xlimit) + ylim(ylimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=FS, y=SB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) + xlim(xlimit)  + ylim(ylimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=FS, y=SB),bins=40) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)  
dev.copy2pdf(file="FS_SB_double_blowup.pdf")

ylimit=c(0,0.2)
quartz()
good=ggplot(data=vie_snps) + xlim(xlimit) + ylim(ylimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=FS, y=SB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) + xlim(xlimit)  + ylim(ylimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=FS, y=SB),bins=40) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)  
dev.copy2pdf(file="FS_SB_double_blowup_0.2.pdf")

quartz()
good=ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=TDB, y=aRPB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=TDB, y=aRPB),bins=40) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)  
dev.copy2pdf(file="TDB_aRPB_double.pdf")

quartz()
good=ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=TDB, y=aRPB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=TDB, y=aRPB),bins=40) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)  
dev.copy2pdf(file="TDB_aRPB_double.pdf")

quartz()
good=ggplot(data=vie_snps) + xlim(xlimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=FS, y=BSB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) + xlim(xlimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=FS, y=BSB),bins=40) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)  
dev.copy2pdf(file="FS_BSB_double.pdf")

quartz()
ylimit=c(0,200)
good=ggplot(data=vie_snps) + xlim(xlimit) + ylim(ylimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=FS, y=BSB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) + xlim(xlimit)  + ylim(ylimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=FS, y=BSB),bins=40) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)  
dev.copy2pdf(file="FS_BSB_double_blowup.pdf")

quartz()
ylimit=c(0,10)
xlimit=c(0,200)
good=ggplot(data=vie_snps) + xlim(xlimit) + ylim(ylimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_i == 1,],aes(x=FS, y=ATD),binwidth=c(5,1)) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=vie_snps) + xlim(xlimit) + ylim(ylimit) +
     stat_binhex(data=vie_snps[vie_snps$DD_c == 0,],aes(x=FS, y=ATD),binwidth=c(5,1)) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)     
dev.copy2pdf(file="FS_ATD_blowup_double.pdf")

# for vie_snps$DD_i == 1
good_snps= vie_snps$DD_i == 1
# for vie_snps$DD_c == 0
bad_snps= vie_snps$DD_c == 0
num_good=length(vie_snps$FS[good_snps])
num_bad= length(vie_snps$FS[bad_snps])
num_tot= length(vie_snps$FS)
num_tot #1664671
num_good #1210874
num_bad #225182
length(vie_snps$DD_i[( vie_snps$SB <= 0.01 & vie_snps$FS >= 30 )]) # 13600
length(vie_snps$DD_i[ ( vie_snps$SB <= 0.01 & vie_snps$FS >= 30 ) & bad_snps]) # 11744
length(vie_snps$DD_i[( vie_snps$SB <= 0.01 & vie_snps$FS >= 30 ) & good_snps]) # 493
length(vie_snps$DD_i[(vie_snps$ATD <= 8 & vie_snps$aRPB >= 3.5) & bad_snps]) #45291
length(vie_snps$DD_i[(vie_snps$ATD <= 8 & vie_snps$aRPB >= 3.5) & good_snps]) #624
length(vie_snps$DD_i[(vie_snps$ATD <= 8 & vie_snps$aRPB >= 3.5)]) #52144
length(vie_snps$DD_i[(vie_snps$ATD <= 8 & vie_snps$aRPB >= 3.5) | ( vie_snps$SB <= 0.01 & vie_snps$FS >= 30 )  ]) # 55375
length(vie_snps$DD_i[((vie_snps$ATD <= 8 & vie_snps$aRPB >= 3.5) | ( vie_snps$SB <= 0.01 & vie_snps$FS >= 30 )) & good_snps]) # 1083
length(vie_snps$DD_i[((vie_snps$ATD <= 8 & vie_snps$aRPB >= 3.5) | ( vie_snps$SB <= 0.01 & vie_snps$FS >= 30 ) ) & bad_snps]) # 47568
length(vie_snps$DD_i[! ((vie_snps$ATD <= 8 & vie_snps$aRPB >= 3.5) | ( vie_snps$SB <= 0.01 & vie_snps$FS >= 30 ))]) # 1609296




# create ROC like curves for retrieval SNPs
quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
fs_vals=c(5,10,20,30,60,100,150,200,250)
good_vals = sapply(fs_vals,function (x) length(vie_snps$FS[good_snps & vie_snps$FS >= x] ))
bad_vals =  sapply(fs_vals,function (x) length(vie_snps$FS[bad_snps & vie_snps$FS >= x] ))
plot(fs_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.2),xlab="FS",ylab="fraction filtered")
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=1)
SB_vals=c(0.2,0.1,0.05,0.01)
ltypes=c(2,3,4,5)
for(i in 1:4){
good_vals = sapply(fs_vals,function (x) length(vie_snps$FS[good_snps & vie_snps$FS >= x & vie_snps$SB <= SB_vals[i]] ))
bad_vals =  sapply(fs_vals,function (x) length(vie_snps$FS[bad_snps & vie_snps$FS >= x & vie_snps$SB <=  SB_vals[i]] ))
lines(fs_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","SB <= 1.0","SB <= 0.2","SB <= 0.1","SB <= 0.05","SB <= 0.01"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,1,2,3,4,5))
dev.copy2pdf(file="FS_SB_SNPs_filtered.pdf")
quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
fs_vals=c(5,10,20,30,60,100,150,200,250)
good_vals = sapply(fs_vals,function (x) length(vie_snps$FS[good_snps & vie_snps$FS >= x] ))
bad_vals =  sapply(fs_vals,function (x) length(vie_snps$FS[bad_snps & vie_snps$FS >= x] ))
plot(fs_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.2),xlab="FS",ylab="fraction filtered")
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=1)
bsb_vals=c(25,50,100,150)
ltypes=c(2,3,4,5)
for(i in 1:4){
good_vals = sapply(fs_vals,function (x) length(vie_snps$FS[good_snps & vie_snps$FS >= x & vie_snps$BSB >= bsb_vals[i]] ))
bad_vals =  sapply(fs_vals,function (x) length(vie_snps$FS[bad_snps & vie_snps$FS >= x & vie_snps$BSB >=  bsb_vals[i]] ))
lines(fs_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","BSB >= 0.0","BSB >= 25","BSB >= 50","BSB >= 100","BSB >= 200"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,1,2,3,4,5))
dev.copy2pdf(file="FS_BSB_SNPs_filtered.pdf")

quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
aRPB_vals=c(0,1,2.5,4,5,7.5,10,15)
good_vals = sapply(aRPB_vals,function (x) length(vie_snps$aRPB[good_snps & vie_snps$aRPB >= x] ))
bad_vals =  sapply(aRPB_vals,function (x) length(vie_snps$aRPB[bad_snps & vie_snps$aRPB >= x] ))
plot(aRPB_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.25),xlab="aRPB",ylab="fraction filtered")
lines(aRPB_vals, bad_vals/num_bad,type="l", col="red", lty=1)
ATD_vals=c(15,10,7.5,5,2.5,1)
ltypes=c(2,3,4,5,6,2)
for(i in 1:4){
good_vals = sapply(aRPB_vals,function (x) length(vie_snps$aRPB[good_snps & vie_snps$aRPB >= x & vie_snps$ATD <= ATD_vals[i]] ))
bad_vals =  sapply(aRPB_vals,function (x) length(vie_snps$aRPB[bad_snps & vie_snps$aRPB >= x & vie_snps$ATD <=  ATD_vals[i]] ))
lines(aRPB_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(aRPB_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","ATD <= 15","ATD <= 10","ATD <= 7.5","ATD <= 5","ATD <= 2.5","ATD <= 1.0"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,2,3,4,5,6,2))
dev.copy2pdf(file="aRPB_ATD_SNPs_filtered.pdf")


quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
TDB_vals=c(0,10,25,50,75,100)
good_vals = sapply(TDB_vals,function (x) length(vie_snps$TDB[good_snps & vie_snps$TDB >= x] ))
bad_vals =  sapply(TDB_vals,function (x) length(vie_snps$TDB[bad_snps & vie_snps$TDB >= x] ))
plot(TDB_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.25),xlab="TDB",ylab="fraction filtered")
lines(TDB_vals, bad_vals/num_bad,type="l", col="red", lty=1)
ATD_vals=c(15,10,7.5,5,2.5,1)
ltypes=c(2,3,4,5,6,2)
for(i in 1:4){
good_vals = sapply(TDB_vals,function (x) length(vie_snps$TDB[good_snps & vie_snps$TDB >= x & vie_snps$ATD <= ATD_vals[i]] ))
bad_vals =  sapply(TDB_vals,function (x) length(vie_snps$TDB[bad_snps & vie_snps$TDB >= x & vie_snps$ATD <=  ATD_vals[i]] ))
lines(TDB_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(TDB_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","ATD <= 15","ATD <= 10","ATD <= 7.5","ATD <= 5","ATD <= 2.5","ATD <= 1.0"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,2,3,4,5,6,2))
dev.copy2pdf(file="TDB_ATD_SNPs_filtered.pdf")

