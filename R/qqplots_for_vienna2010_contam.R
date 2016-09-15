library(hexbin)

binned_qq = function(obs_Pval,dist_Pval,main_text){
	nl_oP = -1 * log10(obs_Pval)
	nl_dP1 = -1 * log10(dist_Pval)
	s_nl_oP <- sort(nl_oP,decreasing=F)
	s_nl_dP <- sort(nl_dP1,decreasing=F)
	bin <- hexbin(s_nl_dP,s_nl_oP,xbins=400)
	par(bg="white")
	pp<-plot(bin,legend=FALSE,style="constant.col",main=main_text,xlab=expression(paste("-log(P",scriptstyle(null),")")),ylab=expression(paste("-log(P",scriptstyle(obs),")")))
	# alternative:
	pushHexport(pp$plot.vp)
	# opens the viewport for all grid commands
	grid.abline(0,1,gp = gpar(col="red",lty=2))
	#closes the viewport
	popViewport()
} 


setwd("/Volumes/Temp/Lukas/Data/Vienna_2010/Polymorphisms_SimCont/CMH-Tests/FDR")
obs_Pval_2010=scan("/Volumes/Temp/Lukas/Data/Vienna_2010/Polymorphisms_SimCont/CMH-Tests/FDR/fem_2010_pI25_pII25_pIII25_pIel_pIIel_pIIIel_unfiltered_masked_pC18removed_obsP.cmhout")
dist_Pvalx10=scan("fem_2010_pI25_pII25_pIII25_pIel_pIIel_pIIIel_unfiltered_masked_pC18removed.cmhout_42.0_nullP.out10x")
dist_Pvalx1=scan("/Volumes/Temp/Lukas/Data/Vienna_2010/Polymorphisms_SimCont/CMH-Tests/FDR/fem_2010_pI25_pII25_pIII25_pIel_pIIel_pIIIel_unfiltered_masked_pC18removed.cmhout_42.0_nullP.out")
nl_oP = -1 * log10(obs_Pval_2010)
nl_dP1 = -1 * log10(dist_Pvalx1)
par(bg="white")
qqplot(nl_dP1,nl_oP,xlab=expression(paste("-log(P",scriptstyle(null),")")),ylab=expression(paste("-log(P",scriptstyle(obs),")")))
abline(0,1,col="red",lty=2)
title(main=expression(paste("Vienna 2010, unfiltered, ",alpha,"=42")))
dev.print(png,width=640,file="/Volumes/Temp/Lukas/Data/Vienna_2010/Polymorphisms_SimCont/CMH-Tests/FDR/qqplot_vienna2010_cont_alpha42.png")

#for vienna 2010 filtered
setwd("/Volumes/Temp/Lukas/Data/Vienna_2010/Polymorphisms_Filtered/CMH/FDR")
dist_Pval=scan("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_q20_filtered_masked_mc7_only_chroms.cmhout_4.0_nullP.out")
obs_Pval=scan("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_q20_filtered_masked_mc7_only_chroms.cmhout.obsPval")
nl_oP = -1 * log10(obs_Pval)
nl_dP1 = -1 * log10(dist_Pval)
par(bg="white")
qqplot(nl_dP1,nl_oP,xlab=expression(-log[10](P[null])),ylab=expression(-log[10](P[obs])))
abline(0,1,col="red",lty=2)
title(main=expression(paste("Vienna 2010, filtered, ",alpha,"=4")))
dev.print(png,width=640,file="/Volumes/Temp/Lukas/Data/Vienna_2010/Polymorphisms_SimCont/CMH-Tests/FDR/qqplot_vienna2010_cont_alpha42.png")
hist(obs_Pval)
hist(obs_Pval,main="Vienna 2010, filtered, obs. P values")
qq(obs_Pval,main="qq plot against uniform for Vienna 2010 filtered")
dist_Pval20=scan("above_10/females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_q20_filtered_masked_mc7_only_chroms.cmhout_20.0_nullP.out")
nl_P20 = -log10(dist_Pval20)
qqplot(nl_P20,nl_oP,xlab=expression(-log[10](P[null])),ylab=expression(-log[10](P[obs])))
abline(0,1,col="red",lty=2)
title(main=expression(paste("Vienna 2010 against ",alpha,"=20")))nl_P20 = -log10(dist_Pval20)

qqplot(nl_P16,nl_oP,xlab=expression(-log[10](P[null])),ylab=expression(-log[10](P[obs])))
abline(0,1,col="red",lty=2)
title(main=expression(paste("Vienna 2010 against ",alpha,"=16")))

#for vie2011  filtered
setwd("/Volumes/Temp/Lukas/Data/Vienna_2011/Polymorphisms_Filtered/CMH/FDR")
obs_Pval=scan("fem_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_masked_mc8_only_chroms_mct15_mcv25.obsP")
dist_Pval=scan("fem_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_masked_mc8_only_chroms_mct15_mcv25.cmhout_a_21.0_beta__r1_41.58_r2_82.7826086957_r3_159.209302326_nullP.out")
nl_oP = -1 * log10(obs_Pval)
nl_dP1 = -1 * log10(dist_Pval)
par(bg="white")
qqplot(nl_dP1,nl_oP,xlab=expression(paste("-log(P",scriptstyle(null),")")),ylab=expression(paste("-log(P",scriptstyle(obs),")")))
abline(0,1,col="red",lty=2)
title(main=expression(paste("Vienna 2011, filtered, ",alpha,"=21 ", beta," by cov")))
dev.print(png,width=640,file="qqplot_vienna2011_filt_a21_b_cov.png")
#for fixed beta
dist_Pval=scan("fem_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_masked_mc8_only_chroms_mct15_mcv25.cmhout_a_25.0_b_90.0_nullP.out")
nl_oP = -1 * log10(obs_Pval)
nl_dP1 = -1 * log10(dist_Pval)
par(bg="white")
qqplot(nl_dP1,nl_oP,xlab=expression(paste("-log(P",scriptstyle(null),")")),ylab=expression(paste("-log(P",scriptstyle(obs),")")))
abline(0,1,col="red",lty=2)
title(main=expression(paste("Vienna 2011, filtered, ",alpha,"=25 ", beta,"=90")))
dev.print(png,width=640,file="qqplot_vienna2011_filt_a25_b_90.png")

# plot histograms 
hist(obs_Pval,col=rgb(0,0,1,1/4),breaks=100)
hist(dist_Pval,col=rgb(1,0,0,1/4),breaks=100,add=T)
legend('topright',c('obs. P',expression(paste("null P (",alpha,"=21,",beta,"=90)"))),fill=rgb(0:1,1/4,1:0,1/4))

hist(obs_Pval,col=rgb(0,0,1,1/4),ylim=c(0.8,1.9),freq=FALSE,breaks=250)
hist(dist_Pval,col=rgb(1,0,0,1/4),ylim=c(0.8,1.9),freq=F,breaks=250,add=T)
legend('topright',c('obs. P',expression(paste("null P (",alpha,"=21,",beta,"=90)"))),fill=rgb(0:1,1/4,1:0,1/4))

# Vienna 2010 filtered
setwd("/Volumes/Temp/Lukas/Data/Vienna_2010/Polymorphisms_Filtered/CMH/FDR")
obs_Pval=scan("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_q20_filtered_masked_mc7_only_chroms_noAceto_mct15_mcv25.obsP")
dist_Pval=scan("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_q20_filtered_masked_mc7_only_chroms_noAceto_mct15_mcv25.cmhout_a_40.0_beta__r1_57.3333333333_r2_92.4675324675_r3_154.285714286_nullP.out")
nl_oP = -1 * log10(obs_Pval)
nl_dP1 = -1 * log10(dist_Pval)
par(bg="white")
s_nl_oP <- sort(nl_oP,decreasing=F)
s_nl_dP <- sort(nl_dP1,decreasing=F)
bin <- hexbin(s_nl_dP,s_nl_oP,xbins=400)
par(bg="white")
pp<-plot(bin,legend=FALSE,style="constant.col",main=expression(paste("Vienna 2010, filtered, ",alpha,"=40 ", beta," by coverage")),xlab=expression(paste("-log(P",scriptstyle(null),")")),ylab=expression(paste("-log(P",scriptstyle(obs),")")))
# alternative:
pushHexport(pp$plot.vp)
# opens the viewport for all grid commands
grid.abline(0,1,gp = gpar(col="red",lty=2))
#closes the viewport
popViewport()
dev.copy2pdf(file="vie2010_filt_a_40_b_cov.pdf")




qqplot(nl_dP1,nl_oP,xlab=expression(paste("-log(P",scriptstyle(null),")")),ylab=expression(paste("-log(P",scriptstyle(obs),")")))
abline(0,1,col="red",lty=2)
title(main=expression(paste("Vienna 2010, filtered, ",alpha,"=40 ", beta," by cov")))
dev.print(png,width=640,file="qqplot_vienna2010_filt_a40_b_cov.png")
#for fixed beta
dist_Pval=scan("f")
nl_oP = -1 * log10(obs_Pval)
nl_dP1 = -1 * log10(dist_Pval)
par(bg="white")
qqplot(nl_dP1,nl_oP,xlab=expression(paste("-log(P",scriptstyle(null),")")),ylab=expression(paste("-log(P",scriptstyle(obs),")")))
abline(0,1,col="red",lty=2)
title(main=expression(paste("Vienna 2011, filtered, ",alpha,"=25 ", beta,"=90")))
dev.print(png,width=640,file="qqplot_vienna2011_filt_a25_b_90.png")

dist_Pval=scan("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_q20_filtered_masked_mc7_only_chroms_noAceto_mct15_mcv25.cmhout_a_15.0_beta__r1_21.5_r2_34.6753246753_r3_57.8571428571_nullP.out")
main_text=expression(paste("Vienna 2010, filtered, ",alpha,"=16 ", beta," by coverage"))
binned_qq(obs_Pval,dist_Pval,main_text)


dist_Pval=scan("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_q20_filtered_masked_mc7_only_chroms_noAceto_mct15_mcv25.cmhout_a_44.0_beta_r1_63.07_r2_101.71_r3_169.71_nullP.out")

dev.copy2pdf(file="vie2010_filt_a_44_b_cov.pdf")

