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


setwd("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/FDR")
obs_Pval=scan("females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.obs")
dist_Pval=scan("females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.cmhout_a_19.0_beta_r1_37.62_r2_75.17_r3_144.49_nullP.out10x_sorted")
nl_oP = -1 * log10(obs_Pval)
nl_dP1 = -1 * log10(dist_Pval)
qqplot_hex(nl_oP,nl_dP1,title="Italy 2011, ",alpha="19.0")
dev.copy2pdf(file="ita2011_a_19.0_b_cov_hex.pdf")
dist_Pval=scan("females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.cmhout_a_20.0_beta_r1_39.60_r2_79.13_r3_152.09_nullP.out10x")
nl_oP = -1 * log10(obs_Pval)
nl_dP1 = -1 * log10(dist_Pval)
qqplot_hex(nl_oP,nl_dP1,title="Italy 2011, ",alpha="20.0")
dev.copy2pdf(file="ita2011_a_20.0_b_cov_hex.pdf")

# add fdr
obs_pval=sort(obs_Pval)
n_ratio=length(obs_pval)/length(dist_Pval)

ranks = 1000:1200
for (r in ranks)
	{
	ps20_smaller =  n_ratio*sum(dist_Pval < obs_pval[r])
	if (ps20_smaller/r > 0.04){
	cat ("r=", r, obs_pval[r], ps20_smaller/r, "\n")}
# not ps20_smaller/(r + ps20_smaller) as r already contains the FP
	}
# rank change from 0.049 to 5 % FDR
r= 1023 2.229684e-09 0.04897361 
r= 1024 2.23817e-09 0.04892578 
r= 1025 2.246251e-09 0.04907317 
r= 1026 2.324404e-09 0.05048733 
r= 1027 2.346382e-09 0.05073028 
r= 1028 2.367414e-09 0.05087549 
r= 1029 2.376741e-09 0.05102041 

#does not work 
add_fdr_q <- function(obs_P,s_null_P){	
	rank_obsP = rank(obsP)
	n_ratio=length(obs_P)/length(s_null_P)
	q_value=n_ratio*sapply(obs_P,function(x){sum( s_null_P < x)})/rank_obsP
	return(q_value)
}




setwd("/Volumes/Temp/Lukas/Data/Vienna_2011/Polymorphisms_Filtered/CMH/FDR")
obs_Pval=scan("fem_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_masked_mc8_only_chroms_mct15_mcv25.obsP")
dist_Pval=scan("fem_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_masked_mc8_only_chroms_mct15_mcv25.cmhout_a_21.5_beta_r1_42.57_r2_84.75_r3_163.00_nullP.out")
qqplot_hex(obs_Pval,dist_Pval,title="Vienna 2011, filtered, ",alpha="21.5")
dev.copy2pdf(file="vie2011_filt_a_21.5_b_cov_hex.pdf")

dist_Pval=scan("fem_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_masked_mc8_only_chroms_mct15_mcv25.cmhout_a_10.0_beta__r1_19.8_r2_39.4202898551_r3_75.8139534884_nullP.out")
qqplot_hex(obs_Pval,dist_Pval,title="Vienna 2011, filtered, ",alpha="10")
dev.copy2pdf(file="vie2011_filt_a_10_b_cov_hex.pdf")

dist_Pval=scan("fem_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_masked_mc8_only_chroms_mct15_mcv25.cmhout_a_15.0_beta__r1_29.7_r2_59.1304347826_r3_113.720930233_nullP.out")
qqplot_hex(obs_Pval,dist_Pval,title="Vienna 2011, filtered, ",alpha="15")
dev.copy2pdf(file="vie2011_filt_a_15_b_cov_hex.pdf")

dist_Pval=scan("fem_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_masked_mc8_only_chroms_mct15_mcv25.cmhout_a_15.0_b_54.0_nullP.out")
qqplot_hex(obs_Pval,dist_Pval,title="Vienna 2011, filtered, ",alpha="15", beta="=54.0 ")
dev.copy2pdf(file="vie2011_filt_a_15_b_54_hex.pdf")

dist_Pval=scan("fem_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_masked_mc8_only_chroms_mct15_mcv25.cmhout_a_10.0_b_36.0_nullP.out")
qqplot_hex(obs_Pval,dist_Pval,title="Vienna 2011, filtered, ",alpha="10", beta="=36.0 ")
dev.copy2pdf(file="vie2011_filt_a_15_b_36_hex.pdf")

s_nl_oP <- sort(nl_oP,decreasing=F)
s_nl_dP <- sort(nl_dP1,decreasing=F)
bin <- hexbin(s_nl_dP,s_nl_oP,xbins=400,xlab=expression(paste("-log(P",scriptstyle(null),")")),ylab=expression(paste("-log(P",scriptstyle(obs),")")))
pp<-plot(bin,legend=FALSE,style="constant.col",main=expression(paste("Vienna 2011, filtered, ",alpha,"=20 ", beta," by coverage")))
hvp=hexViewport(bin)
hexVP.abline(hvp,0,1,col="red",lty=2)

# alternative:
pushHexport(pp$plot.vp)
# opens the viewport for all grid commands
grid.abline(0,1,gp = gpar(col=2,lty=2))
#closes the viewport
popViewport()

# adding points to binplot
xy_points=runif(100,0,10)
dim(xy_points) = c(50,2)
grid.points(xy_points[,1],xy_points[,2],pch="x")






# looking at chi2 values
setwd("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH")
chi2_vie2010 = read.table("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.cmhout.chi2",header=T)
chi2_vie2010_avg = read.table("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.cmhout.chi2.avgNull",header=T)
quartz()
par(mfrow=c(2,2))
hist(chi2_vie2010$X1.4,breaks=500,ylim=c(0,10000))
hist(chi2_vie2010$X2.5,breaks=500,ylim=c(0,10000))
hist(chi2_vie2010$X3.6,breaks=500,ylim=c(0,10000))
quartz()
par(mfrow=c(2,2))
hist(chi2_vie2010_avg$X1.4,breaks=500,ylim=c(0,10000))
hist(chi2_vie2010_avg$X2.5,breaks=500,ylim=c(0,10000))
hist(chi2_vie2010_avg$X3.6,breaks=500,ylim=c(0,10000))
uni_dist=runif(length(na.omit(chi2_vie2010$X1.4)),0,1)
qqplot_hex(chi2_vie2010$X1.4,uni_dist,title="",alpha="",beta="")
uni_dist=runif(length(na.omit(chi2_vie2010$X2.5)),0,1)
qqplot_hex(chi2_vie2010$X2.5,uni_dist,title="",alpha="",beta="")
uni_dist=runif(length(na.omit(chi2_vie2010$X3.6)),0,1)
qqplot_hex(chi2_vie2010$X3.6,uni_dist,title="",alpha="",beta="")
# looking at chi2 values
setwd("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH")
chi2_vie2011 = read.table("females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.cmhout.chi2",header=T)
quartz()
par(mfrow=c(2,2))
hist(chi2_vie2011$X1.4,breaks=500,ylim=c(0,10000))
hist(chi2_vie2011$X2.5,breaks=500,ylim=c(0,10000))
hist(chi2_vie2011$X3.6,breaks=500,ylim=c(0,10000))
uni_dist=runif(length(na.omit(chi2_vie2010$X1.4)),0,1)
qqplot_hex(chi2_vie2010$X1.4,uni_dist,title="",alpha="",beta="")
uni_dist=runif(length(na.omit(chi2_vie2010$X2.5)),0,1)
qqplot_hex(chi2_vie2010$X2.5,uni_dist,title="",alpha="",beta="")
uni_dist=runif(length(na.omit(chi2_vie2010$X3.6)),0,1)
qqplot_hex(chi2_vie2010$X3.6,uni_dist,title="",alpha="",beta="")

setwd("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/Strandbias")
vie_2010_pV = read.table("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.strandbias_P_gt_0.001.gwas",header=T)
uni_dist=runif(length(na.omit(vie_2010_pV$P)),0,1)
quartz()
qqplot_hex(-1*log10(vie_2010_pV$P),-1*log10(uni_dist),title="",alpha="",beta="")
length(vie_2010_pV$P)
vie_2010_001_pV = read.table("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.strandbias_P_gt_0.01.gwas",header=T)
uni_dist=runif(length(na.omit(vie_2010_001_pV$P)),0,1)
qqplot_hex(-1*log10(vie_2010_001_pV$P),-1*log10(uni_dist),title="",alpha="",beta="")
vie_2010=read.table("../females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.gwas",header=T)
uni_dist=runif(length(na.omit(vie_2010$P)),0,1)
qqplot_hex(-1*log10(vie_2010$P),-1*log10(uni_dist),title="qqplot Vie2010 against uniform",alpha="",beta="")
colnames(vie_2010_001_pV)=c("CHR","BP","P")
source("/Volumes/Temp/Lukas/Tools/Scripts/R/manhattan_plot_all_chromosomes_function_qqplot.R")
quartz(height=6,width=12)
par(bg="white")
chroms=c("X","2L","2R","3L","3R","4")#,"XHet","2LHet","2RHet","3LHet","3RHet")
manhattan(vie_2010_001_pV,colors=c("black","slategrey"),limitchromosomes=chroms,main="Vienna (2010)", suggestiveline=c(-log10(4.134235e-07)))

setwd("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/Polymorphisms/CMH")
ps =scan("FDR/females_joined_mct20_mcv25_sb_fb_bd.cmhout_a_24.0_beta_r1_34.13_r2_55.17_r3_91.89_r4_48.00_r5_96.70_r6_184.74_nullP.out_10x")
ps = sort(ps)
obs_pval=read.table("females_joined_mct20_mcv25_sb_fb_bd.gwas",header=T)
obs_pval=obs_pval[order(obs_pval$P),]
n_ratio=length(obs_pval$P)/length(ps)
obs_pval$Q = 1.0
obs_pval$Q1 = obs_pval$Q
ranks = 1001:1500
num=1
for (r in ranks)
	{
	while(num <= length(ps) && obs_pval$P[r] >= ps[num] ) {num = num+1 }
	ps_smaller =  n_ratio*(num-1) #
	#ps_smaller2 =  n_ratio* sum(ps < obs_pval$P[r])
	#cat ("r=", r, obs_pval$P[r], ps_smaller/r,ps_smaller2/r, "\n")
	obs_pval$Q[r]=ps_smaller/r
# not ps20_smaller/(r + ps20_smaller) as r already contains the FP
	}
#cutoff at 6.697709e-08
604613   2R 20285559  snp604613 6.615626e-08 0.04958882
943383   3R  3760522  snp943383 6.697709e-08 0.04954807
1272734   4   403451 snp1272734 6.846161e-08 0.05024631
1289929   X  3297596 snp1289929 6.937622e-08 0.05036916
417007   2R  8539881  snp417007 6.948920e-08 0.05032787







