setwd("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/Joined_Analysis")
vie_fdr_ap=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/vie2010_fdr_dsim_afr_zim_flo_port_MAD_bg0s11_bg0s3_cg3s1_cg3s2_cg3s3_noSSA.sync.count",header=T)
vie_ap=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/vie2010_dsim_afr_zim_flo_port_MAD_bg0s11_bg0s3_cg3s1_cg3s2_cg3s3_noSSA.sync.count",header=T)
ita_fdr_ap=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita2011_fdr_dsim_afr_zim_flo_port_MAD_bg0s11_bg0s3_cg3s1_cg3s2_cg3s3_noSSA.count",header=T)
ita_ap=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita2011_dsim_afr_zim_flo_port_MAD_bg0s11_bg0s3_cg3s1_cg3s2_cg3s3_noSSA.count",header=T)

ovr_ap=read.table("common_snps_dsim_afr_zim_flo_port_MAD_port_cages_noSSA.sync.count",header=T)
ovr_fdr_ap=read.table("common_snps_fdr_dsim_afr_zim_flo_port_MAD_port_cages_noSSA.sync.count",header=T)
#ita_all_only=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita2011_simulans_polym_maf_ge_0.05.all_only",header=F)

#recombination rate
vie_fdr_ap_rec=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/vie2010_fdr_dsim_afr_zim_flo_port_MAD_bg0s11_bg0s3_cg3s1_cg3s2_cg3s3_noSSA.sync.comp.rec.count",header=T)
vie_ap_rec=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/vie2010_dsim_afr_zim_flo_port_MAD_bg0s11_bg0s3_cg3s1_cg3s2_cg3s3_noSSA.sync.comp.rec.count",header=T)
ita_fdr_ap_rec=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita2011_fdr_dsim_afr_zim_flo_port_MAD_bg0s11_bg0s3_cg3s1_cg3s2_cg3s3_noSSA.sync.comp.rec.count",header=T)
ita_ap_rec=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita2011_dsim_afr_zim_flo_port_MAD_bg0s11_bg0s3_cg3s1_cg3s2_cg3s3_noSSA.sync.comp.rec.count",header=T)

ovr_ap_rec=read.table("common_snps_dsim_afr_zim_flo_port_MAD_port_cages_noSSA.sync.comp.rec.count",header=T)
ovr_fdr_ap_rec=read.table("common_snps_fdr_dsim_afr_zim_flo_port_MAD_port_cages_noSSA.sync.comp.rec.count",header=T)
#ita_all_only=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita2011_simulans_polym_maf_ge_0.05.all_only",header=F)

get_rec_tot = function(rec_tab){
    tot_rec=by(rec_tab[,-1:-2],rec_tab$RecRate,colSums)
    tot_mat=matrix(nrow=length(tot_rec),ncol=length(tot_rec[[1]])+3)
    rownames(tot_mat)=names(tot_rec)
    for(i in names(tot_rec)){
        nam_temp=names(tot_rec[[i]])
        tot_rec[[i]]=append(tot_rec[[i]],sum(tot_rec[[i]][c("m10","m01")]))
        tot_rec[[i]]=append(tot_rec[[i]],sum(tot_rec[[i]][c("p10","p01")]))
        tot_rec[[i]]=append(tot_rec[[i]],0)
        names(tot_rec[[i]])=c(nam_temp,c("m1","p1","zeros"))        
        tot_mat[i,]=tot_rec[[i]]
    }
    colnames(tot_mat)=names(tot_rec[[1]])
    tot_rows=rownames(tot_mat)
    tot_mat=rbind(tot_mat,tot_mat["zero",]+tot_mat["low",])
    tot_rows=append(tot_rows,"zl")
    tot_mat=rbind(tot_mat,tot_mat["medium",]+tot_mat["high",])
    tot_rows=append(tot_rows,"mh")
    rownames(tot_mat)=tot_rows
    return(as.table(tot_mat))
}

bar_plot_shared_poly = function(rec_tab,rec_tab_fdr,rec_lev,rec_names,title){
    spaces=c(rep(0.1,length(rec_lev)),1,rep(0.1,length(rec_lev)-1))
    labs=c(expression("fixed in "~italic("D. sim.")),expression("polymorph in "~italic("D. sim.")))
    mono_cols=c("m0","m1","zeros")
    poly_cols=c("p00","p1","p12")
    quartz(width=12)
    par(mfrow=c(1,2), mar=c(4.5,3,3,1), oma=c(0,0,2,0))
    barplot(t(as.matrix(rbind(rec_tab[rec_lev,mono_cols],rec_tab[rec_lev,poly_cols]))),col=c("red","green","blue"),names.arg=c(rec_names,rec_names),beside=F,main="all SNPs",space=spaces)
    text(x=c(mean(barplt[1:length(rec_lev)]),mean(barplt[(length(rec_lev)+1):length(barplt)])),y=0, pos=1, offset=2.35, labels=c("rec. rate [cM/Mb]"),xpd=T)
    text(x=c(mean(barplt[1:length(rec_lev)]),mean(barplt[(length(rec_lev)+1):length(barplt)])),y=0, pos=1, offset=3.5, labels=labs,xpd=T)
#    text(x=c(mean(barplt[1:length(rec_lev)]),mean(barplt[(length(rec_lev)+1):length(barplt)])),y=1, pos=1, offset=2.5, labels=labs,xpd=T)
#legend("topright",c("divergent","one allele","two alleles"),fill=c("red","green","blue"),title="Shared alleles")
#dev.copy2pdf(file="vienna_ap_rec_rate.pdf")
#quartz()
    barplt=barplot(t(as.matrix(rbind(rec_tab_fdr[rec_lev,mono_cols],rec_tab_fdr[rec_lev,poly_cols]))),col=c("red","green","blue"),beside=F,names.arg=c(rec_names,rec_names),main="significant SNPs",space=spaces)
    text(x=c(mean(barplt[1:length(rec_lev)]),mean(barplt[(length(rec_lev)+1):length(barplt)])),y=0, pos=1, offset=2.35, labels=c("rec. rate [cM/Mb]"),xpd=T)
    text(x=c(mean(barplt[1:length(rec_lev)]),mean(barplt[(length(rec_lev)+1):length(barplt)])),y=0, pos=1, offset=3.5, labels=labs,xpd=T)
    mtext(title, side=3, line=0, outer=TRUE, cex=1.5,font=2)
    legend("topright",c("divergent","one allele","two alleles"),fill=c("red","green","blue"),title="Shared")

}

rec_lev=c("zero","low","medium","high")
rec_lev=c("zl","mh")
rec_names=c("< 1.7", ">= 1.7")
vie_tot_rec=get_rec_tot(vie_ap_rec)
vie_fdr_tot_rec=get_rec_tot(vie_fdr_ap_rec)
bar_plot_shared_poly(vie_tot_rec,vie_fdr_tot_rec,rec_lev,rec_names,title=c("Vienna, 2010"))
dev.copy2pdf(file="vienna_ap_rec_rate_hl.pdf")
ita_tot_rec=get_rec_tot(ita_ap_rec)
ita_fdr_tot_rec=get_rec_tot(ita_fdr_ap_rec)
bar_plot_shared_poly(ita_tot_rec,ita_fdr_tot_rec,rec_lev,rec_names,title=c("Bolzano, 2011"))
dev.copy2pdf(file="italy_ap_rec_rate_hl.pdf")
ovr_tot_rec=get_rec_tot(ovr_ap_rec)
ovr_fdr_tot_rec=get_rec_tot(ovr_fdr_ap_rec)
bar_plot_shared_poly(ovr_tot_rec,ovr_fdr_tot_rec,rec_lev,rec_names,title=c("Overlap"))
dev.copy2pdf(file="overlap_ap_rec_rate_hl.pdf")



spaces=c(rep(0.1,length(rec_lev)),1,rep(0.1,length(rec_lev)-1))
labs=c(expression("fixed in "~italic("D. sim.")),expression("polymorph in "~italic("D. sim.")))
mono_cols=c("m0","m1","zeros")
poly_cols=c("p00","p1","p12")
quartz(width=12)
par(mfrow=c(1,2), mar=c(4.5,3,3,1), oma=c(0,0,2,0))
barplot(t(as.matrix(rbind(vie_tot_rec[rec_lev,mono_cols],vie_tot_rec[rec_lev,poly_cols]))),col=c("red","green","blue"),names.arg=c(rec_names,rec_names),beside=F,main="all SNPs",space=spaces)
text(x=c(mean(barplt[1:length(rec_lev)]),mean(barplt[(length(rec_lev)+1):length(barplt)])),y=1, pos=1, offset=2.5, labels=labs,xpd=T)
#legend("topright",c("divergent","one allele","two alleles"),fill=c("red","green","blue"),title="Shared alleles")
#dev.copy2pdf(file="vienna_ap_rec_rate.pdf")
#quartz()
barplt=barplot(t(as.matrix(rbind(vie_fdr_tot_rec[rec_lev,mono_cols],vie_fdr_tot_rec[rec_lev,poly_cols]))),col=c("red","green","blue"),beside=F,main="significant SNPs",space=spaces)
text(x=c(mean(barplt[1:length(rec_lev)]),mean(barplt[(length(rec_lev)+1):length(barplt)])),y=0, pos=1, offset=2.35, labels=c("rec. rate [cM/Mb]"),xpd=T)
text(x=c(mean(barplt[1:length(rec_lev)]),mean(barplt[(length(rec_lev)+1):length(barplt)])),y=0, pos=1, offset=3.5, labels=labs,xpd=T)
mtext("Vienna, 2010", side=3, line=0, outer=TRUE, cex=1.5,font=2)
legend("topright",c("divergent","one allele","two alleles"),fill=c("red","green","blue"),title="Shared")
dev.copy2pdf(file="vienna_ap_rec_rate.pdf")
ita_tot_rec=get_rec_tot(ita_ap_rec)
ita_fdr_tot_rec=get_rec_tot(ita_fdr_ap_rec)
quartz(width=14)
par(mfrow=c(1,2))
barplt=barplot(t(as.matrix(rbind(ita_tot_rec[rec_lev,mono_cols],ita_tot_rec[rec_lev,poly_cols]))),col=c("red","green","blue"),beside=F,main="Bolzano, all SNPs",space=spaces)
text(x=c(mean(barplt[1:length(rec_lev)]),mean(barplt[(length(rec_lev)+1):length(barplt)])),y=0, pos=1, offset=2.5, labels=labs,xpd=T)
#legend("topright",c("divergent","one allele","two alleles"),fill=c("red","green","blue"),title="Shared alleles")
#quartz()
barplt=barplot(t(as.matrix(rbind(ita_fdr_tot_rec[rec_lev,mono_cols],ita_fdr_tot_rec[rec_lev,poly_cols]))),col=c("red","green","blue"),beside=F,main="Bolzano, significant SNPs",space=spaces)
text(x=c(mean(barplt[1:length(rec_lev)]),mean(barplt[(length(rec_lev)+1):length(barplt)])),y=0, pos=1, offset=2.5, labels=labs,xpd=T)
legend("topright",c("divergent","one allele","two alleles"),fill=c("red","green","blue"),title="Shared alleles")
dev.copy2pdf(file="italy_ap_rec_rate.pdf")
ovr_tot_rec=get_rec_tot(ovr_ap_rec)
ovr_fdr_tot_rec=get_rec_tot(ovr_fdr_ap_rec)
quartz(width=14)
par(mfrow=c(1,2))
barplt=barplot(t(as.matrix(rbind(ovr_tot_rec[rec_lev,mono_cols],ovr_tot_rec[rec_lev,poly_cols]))),col=c("red","green","blue"),beside=F,main="Overlap,all SNPs",space=spaces)
text(x=c(mean(barplt[1:length(rec_lev)]),mean(barplt[(length(rec_lev)+1):length(barplt)])),y=0, pos=1, offset=2.5, labels=labs,xpd=T)
#legend("topright",c("divergent","one allele","two alleles"),fill=c("red","green","blue"),title="Shared alleles")
#quartz()
barplt=barplot(t(as.matrix(rbind(ovr_fdr_tot_rec[rec_lev,mono_cols],ovr_fdr_tot_rec[rec_lev,poly_cols]))),col=c("red","green","blue"),beside=F,main="Overlap, significant SNPs",space=spaces)
text(x=c(mean(barplt[1:length(rec_lev)]),mean(barplt[(length(rec_lev)+1):length(barplt)])),y=0, pos=1, offset=2.5, labels=labs,xpd=T)
legend("topright",c("divergent","one allele","two alleles"),fill=c("red","green","blue"),title="Shared alleles")
dev.copy2pdf(file="overlap_ap_rec_rate.pdf")


mono_cols=c("m0","m01","m10")
poly_cols=c("p00","p10","p01","p12")
chr_names=c("X","2L","2R","3L","3R","4")
all_cnt=c("None","One","Poly")
pol_cnt=c("pNone","pOne","pPoly")
vie_fdr_ap$Tot=apply(vie_fdr_ap[,c(mono_cols,poly_cols)],1,sum)
vie_fdr_ap$Mt=rowSums(vie_fdr_ap[,mono_cols])  
vie_fdr_ap$Pt=rowSums(vie_fdr_ap[,poly_cols])
vie_fdr_ap$M_fract=vie_fdr_ap$Mt/(vie_fdr_ap$Pt+vie_fdr_ap$Mt)
vie_fdr_ap_ratios=round(vie_fdr_ap[,c(mono_cols,poly_cols)]/vie_fdr_ap$Tot, digits=3)
rownames(vie_fdr_ap_ratios)=as.character(vie_fdr_ap$CHR)
vie_fdr_ap_ratios$None=vie_fdr_ap_ratios$m0+vie_fdr_ap_ratios$p00
vie_fdr_ap_ratios$One=vie_fdr_ap_ratios$m10+vie_fdr_ap_ratios$m01+vie_fdr_ap_ratios$p10 + vie_fdr_ap_ratios$p01
vie_fdr_ap_ratios$Poly=vie_fdr_ap_ratios$p12
vie_fdr_ap_ratios$pNone=(vie_fdr_ap_ratios$p00)*vie_fdr_ap$Tot/vie_fdr_ap$Pt
vie_fdr_ap_ratios$pOne=(vie_fdr_ap_ratios$p10 + vie_fdr_ap_ratios$p01)*vie_fdr_ap$Tot/vie_fdr_ap$Pt
vie_fdr_ap_ratios$pPoly=vie_fdr_ap_ratios$p12*vie_fdr_ap$Tot/vie_fdr_ap$Pt

quartz()
barplot(t(as.matrix(vie_fdr_ap_ratios[chr_names,all_cnt])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main="Significant SNPs Vienna, 2010")
legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="vie_anc_poly_allsim_fdr_barplot.pdf")
quartz()
barplot(t(as.matrix(vie_fdr_ap_ratios[chr_names,pol_cnt])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main=bquote("Significant SNPs Vienna, 2010, polymorph in"~italic("D. sim.") ) )
legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="vie_anc_poly_allsim_fdr_barplot_poly_only.pdf")



vie_ap$Tot=apply(vie_ap[,c(mono_cols,poly_cols)],1,sum)
vie_ap$Mt=rowSums(vie_ap[,mono_cols])  
vie_ap$Pt=rowSums(vie_ap[,poly_cols])
vie_ap$M_fract=vie_ap$Mt/(vie_ap$Pt+vie_ap$Mt)
vie_ap_ratios=round(vie_ap[,c(mono_cols,poly_cols)]/vie_ap$Tot, digits=3)
rownames(vie_ap_ratios)=as.character(vie_ap$CHR)
vie_ap_ratios$None=vie_ap_ratios$m0+vie_ap_ratios$p00
vie_ap_ratios$One=vie_ap_ratios$m10+vie_ap_ratios$m01+vie_ap_ratios$p10 + vie_ap_ratios$p01
vie_ap_ratios$Poly=vie_ap_ratios$p12
vie_ap_ratios$pNone=(vie_ap_ratios$p00)*vie_ap$Tot/vie_ap$Pt
vie_ap_ratios$pOne=(vie_ap_ratios$p10 + vie_ap_ratios$p01)*vie_ap$Tot/vie_ap$Pt
vie_ap_ratios$pPoly=vie_ap_ratios$p12*vie_ap$Tot/vie_ap$Pt
quartz()
barplot(t(as.matrix(vie_ap_ratios[chr_names,all_cnt])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main="All SNPs Vienna, 2010")
#legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="vie_anc_poly_allsim_all_barplot.pdf")
quartz()
barplot(t(as.matrix(vie_ap_ratios[chr_names,pol_cnt])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main=bquote("All SNPs Vienna, 2010, polymorph in"~italic("D. sim.")) )
#legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="vie_anc_poly_allsim_all_barplot_poly_only.pdf")

ita_fdr_ap$Tot=apply(ita_fdr_ap[,c(mono_cols,poly_cols)],1,sum)
ita_fdr_ap$Mt=rowSums(ita_fdr_ap[,mono_cols])  
ita_fdr_ap$Pt=rowSums(ita_fdr_ap[,poly_cols])
ita_fdr_ap$M_fract=ita_fdr_ap$Mt/(ita_fdr_ap$Pt+ita_fdr_ap$Mt)
ita_fdr_ap_ratios=round(ita_fdr_ap[,c(mono_cols,poly_cols)]/ita_fdr_ap$Tot, digits=3)
rownames(ita_fdr_ap_ratios)=as.character(ita_fdr_ap$CHR)
ita_fdr_ap_ratios$None=ita_fdr_ap_ratios$m0+ita_fdr_ap_ratios$p00
ita_fdr_ap_ratios$One=ita_fdr_ap_ratios$m10+ita_fdr_ap_ratios$m01+ita_fdr_ap_ratios$p10 + ita_fdr_ap_ratios$p01
ita_fdr_ap_ratios$Poly=ita_fdr_ap_ratios$p12
ita_fdr_ap_ratios$pNone=(ita_fdr_ap_ratios$p00)*ita_fdr_ap$Tot/ita_fdr_ap$Pt
ita_fdr_ap_ratios$pOne=(ita_fdr_ap_ratios$p10 + ita_fdr_ap_ratios$p01)*ita_fdr_ap$Tot/ita_fdr_ap$Pt
ita_fdr_ap_ratios$pPoly=ita_fdr_ap_ratios$p12*ita_fdr_ap$Tot/ita_fdr_ap$Pt
quartz()
barplot(t(as.matrix(ita_fdr_ap_ratios[chr_names,all_cnt])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main="Significant SNPs Bolzano, 2011")
legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="ita_anc_poly_allsim_fdr_barplot.pdf")
quartz()
barplot(t(as.matrix(ita_fdr_ap_ratios[chr_names,pol_cnt])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main=bquote("Significant SNPs Bolzano, 2010, polymorph in"~italic("D. sim.") ) )
legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="ita_anc_poly_allsim_fdr_barplot_poly_only.pdf")


ita_ap$Tot=apply(ita_ap[,c(mono_cols,poly_cols)],1,sum)
ita_ap$Mt=rowSums(ita_ap[,mono_cols])  
ita_ap$Pt=rowSums(ita_ap[,poly_cols])
ita_ap$M_fract=ita_ap$Mt/(ita_ap$Pt+ita_ap$Mt)
ita_ap_ratios=round(ita_ap[,c(mono_cols,poly_cols)]/ita_ap$Tot, digits=3)
rownames(ita_ap_ratios)=as.character(ita_ap$CHR)
ita_ap_ratios$None=ita_ap_ratios$m0+ita_ap_ratios$p00
ita_ap_ratios$One=ita_ap_ratios$m10+ita_ap_ratios$m01+ita_ap_ratios$p10 + ita_ap_ratios$p01
ita_ap_ratios$Poly=ita_ap_ratios$p12
ita_ap_ratios$pNone=(ita_ap_ratios$p00)*ita_ap$Tot/ita_ap$Pt
ita_ap_ratios$pOne=(ita_ap_ratios$p10 + ita_ap_ratios$p01)*ita_ap$Tot/ita_ap$Pt
ita_ap_ratios$pPoly=ita_ap_ratios$p12*ita_ap$Tot/ita_ap$Pt
quartz()
barplot(t(as.matrix(ita_ap_ratios[chr_names,all_cnt])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main="All SNPs Bolzano, 2011")
#legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="ita_anc_poly_allsim_all_barplot.pdf")
quartz()
barplot(t(as.matrix(ita_ap_ratios[chr_names,pol_cnt])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main=bquote("All SNPs Bolzano, 2010, polymorph in"~italic("D. sim.")) )
#legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="ita_anc_poly_allsim_all_barplot_poly_only.pdf")

ovr_fdr_ap$Tot=apply(ovr_fdr_ap[,c(mono_cols,poly_cols)],1,sum)
ovr_fdr_ap$Mt=rowSums(ovr_fdr_ap[,mono_cols])  
ovr_fdr_ap$Pt=rowSums(ovr_fdr_ap[,poly_cols])
ovr_fdr_ap$M_fract=ovr_fdr_ap$Mt/(ovr_fdr_ap$Pt+ovr_fdr_ap$Mt)
ovr_fdr_ap_ratios=round(ovr_fdr_ap[,c(mono_cols,poly_cols)]/ovr_fdr_ap$Tot, digits=3)
rownames(ovr_fdr_ap_ratios)=as.character(ovr_fdr_ap$CHR)
ovr_fdr_ap_ratios$None=ovr_fdr_ap_ratios$m0+ovr_fdr_ap_ratios$p00
ovr_fdr_ap_ratios$One=ovr_fdr_ap_ratios$m10+ovr_fdr_ap_ratios$m01+ovr_fdr_ap_ratios$p10 + ovr_fdr_ap_ratios$p01
ovr_fdr_ap_ratios$Poly=ovr_fdr_ap_ratios$p12
ovr_fdr_ap_ratios$pNone=(ovr_fdr_ap_ratios$p00)*ovr_fdr_ap$Tot/ovr_fdr_ap$Pt
ovr_fdr_ap_ratios$pOne=(ovr_fdr_ap_ratios$p10 + ovr_fdr_ap_ratios$p01)*ovr_fdr_ap$Tot/ovr_fdr_ap$Pt
ovr_fdr_ap_ratios$pPoly=ovr_fdr_ap_ratios$p12*ovr_fdr_ap$Tot/ovr_fdr_ap$Pt
quartz()
barplot(t(as.matrix(ovr_fdr_ap_ratios[chr_names,all_cnt])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main="Common significant SNPs")
legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="ovr_anc_poly_allsim_fdr_barplot.pdf")
quartz()
barplot(t(as.matrix(ovr_fdr_ap_ratios[chr_names,pol_cnt])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main=bquote("Common significant SNPs, polymorph in"~italic("D. sim.") ) )
legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="ovr_anc_poly_allsim_fdr_barplot_poly_only.pdf")


ovr_ap$Tot=apply(ovr_ap[,c(mono_cols,poly_cols)],1,sum)
ovr_ap$Mt=rowSums(ovr_ap[,mono_cols])  
ovr_ap$Pt=rowSums(ovr_ap[,poly_cols])
ovr_ap$M_fract=ovr_ap$Mt/(ovr_ap$Pt+ovr_ap$Mt)
ovr_ap_ratios=round(ovr_ap[,c(mono_cols,poly_cols)]/ovr_ap$Tot, digits=3)
rownames(ovr_ap_ratios)=as.character(ovr_ap$CHR)
ovr_ap_ratios$None=ovr_ap_ratios$m0+ovr_ap_ratios$p00
ovr_ap_ratios$One=ovr_ap_ratios$m10+ovr_ap_ratios$m01+ovr_ap_ratios$p10 + ovr_ap_ratios$p01
ovr_ap_ratios$Poly=ovr_ap_ratios$p12
ovr_ap_ratios$pNone=(ovr_ap_ratios$p00)*ovr_ap$Tot/ovr_ap$Pt
ovr_ap_ratios$pOne=(ovr_ap_ratios$p10 + ovr_ap_ratios$p01)*ovr_ap$Tot/ovr_ap$Pt
ovr_ap_ratios$pPoly=ovr_ap_ratios$p12*ovr_ap$Tot/ovr_ap$Pt
quartz()
barplot(t(as.matrix(ovr_ap_ratios[chr_names,all_cnt])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main="All comon SNPs")
#legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="ovr_anc_poly_allsim_all_barplot.pdf")
quartz()
barplot(t(as.matrix(ovr_ap_ratios[chr_names,pol_cnt])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main=bquote("All common SNPs, polymorph in"~italic("D. sim.")) )
#legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="ovr_anc_poly_allsim_all_barplot_poly_only.pdf")

# enrichment with lod
# files with results
library(Hmisc)
print_cat=c("X","2R","","INTRON","S","NS")
print_labs=c("Intergenic","Genes","CDS","Introns","Syn.","Nonsyn.")
none_cols=c("m0","p00")
one_cols=c("m01","m10","p10","p01")
poly_cols=c("p12")
chr_names=c("X","2L","2R","3L","3R","4")
all_cnt=c("N","O","T")

clean_smpl_ap <- function(smpl_ap,all_cnt=c("N","O","T")) {
	a=colnames(smpl_ap); a[a=="None"]="N";a[a=="One"]="O";a[a=="Two"]="T"
	a[a=="N_50.0"]="N";a[a=="O_50.0"]="O";a[a=="T_50.0"]="T"
	colnames(smpl_ap)=a
	rownames(smpl_ap)=as.character(smpl_ap$Chrom)
	smpl_ap$Tot = rowSums(smpl_ap[,all_cnt])
	return(smpl_ap)
}
add_totals_ap_tab = function(ap_tab,none_cols=c("m0","p00"),one_cols=c("m01","m10","p10","p01"),poly_cols=c("p12"),chr_names=c("X","2L","2R","3L","3R","4"),all_cnt=c("N","O","T")){
	rownames(ap_tab)=as.character(ap_tab$CHR)
	ap_tab[,"N"] = rowSums(ap_tab[,none_cols])
	ap_tab[,"O"] = rowSums(ap_tab[,one_cols])
	ap_tab[,"T"] = ap_tab[,poly_cols]
	ap_tab["total",] = c(NA,colSums(ap_tab[chr_names,-1]))
	a=as.character(ap_tab$CHR)
	a[length(a)]="total"
	ap_tab$CHR=as.factor(a)
	ap_tab$Tot = rowSums(ap_tab[,all_cnt])
	return(ap_tab)
}
get_lod_tab = function(ap_tab,smp_tab,all_cnt=c("N","O","T"),perc=c("_2.5","","_97.5")){
	p95_cols=NULL; for (i in all_cnt) { p95_cols = rbind(p95_cols,paste(i,perc,sep="")) }
	rnames=rownames(smp_tab)
	lod_ap_tab=data.frame(CHR=smp_tab[,1])
	rownames(lod_ap_tab)=rnames
	for (j in 1:dim(p95_cols)[1]){
		i = p95_cols[j,]
		lod_ap_tab[i] = log2( ap_tab[rnames,i[2]]*(ap_tab[rnames,"Totals"] - smp_tab[rnames,i])/(smp_tab[rnames,i]*(ap_tab[rnames,"Totals"] - ap_tab[rnames,i[2]])  ) )
	}
	return(lod_ap_tab)
}
error_bar_plot_all = function(lod_tab,all_cnt=c("N","O","T"),perc=c("_2.5","","_97.5"),print_chr=c("X","2L","2R","3L","3R","total"),dist=0.2,leg="bottomright",main_title="") {
	p95_cols=NULL; for (i in all_cnt) { p95_cols = rbind(p95_cols,paste(i,perc,sep="")) }
	errbar(1:length(print_chr) - dist,lod_tab[print_chr,p95_cols[1,2]],apply(lod_tab[print_chr,p95_cols[1,]],1,min),apply(lod_tab[print_chr,p95_cols[1,]],1,max),ylab=expression(paste("log2(OR)")),main=main_title,xaxt="n",xlab="",pch=15,bg="white",xlim=c(1-dist - 0.1,length(print_chr) + dist + 0.1), ylim=c(min(lod_tab[print_chr,-1])-0.5,max(lod_tab[print_chr,-1])+0.5),cex=1.5)
	errbar(1:length(print_chr),lod_tab[print_chr,p95_cols[2,2]],apply(lod_tab[print_chr,p95_cols[2,]],1,min),apply(lod_tab[print_chr,p95_cols[2,]],1,max),ylab="log2 OR",main="",xaxt="n",yaxt="n",xlab="",add=T, pch=".",bg="white",cex=1.5)
	points(1:length(print_chr),lod_tab[print_chr,p95_cols[2,2]],pch=22,bg="white",cex=1.5)
	errbar(1:length(print_chr) + dist,lod_tab[print_chr,p95_cols[3,2]],apply(lod_tab[print_chr,p95_cols[3,]],1,min),apply(lod_tab[print_chr,p95_cols[3,]],1,max),ylab="log2 OR",main="",xaxt="n",yaxt="n",xlab="",add=T, pch=16,bg="white",cex=1.5)
	axis(1, labels=F,at=1:length(print_chr))
	text(1:length(print_chr), par("usr")[3] - 0.25, srt = 45, adj = 1, labels = print_chr, xpd = TRUE)
	abline(h=0,lt=3)
	legend(leg,c("none","one","two"),pch=c(15,22,16),pt.cex=1.5,title="shared alleles")
}
error_bar_plot_total = function(lod_tab1,lod_tab2,lod_tab3,all_cnt=c("N","O","T"),xlabs=c(0,1,2),xtext="shared alleles",perc=c("_2.5","","_97.5"),print_chr="total",print_cat=c("Vienna","Bolzano","Common"), dist=0.2,leg="bottomright",main_title="") {
	p95_cols=NULL; for (i in all_cnt) { p95_cols = rbind(p95_cols,paste(i,perc,sep="")) }
	lod_tab=rbind(lod_tab1[print_chr,],lod_tab2[print_chr,],lod_tab3[print_chr,])
	errbar(1:length(all_cnt) - dist,lod_tab[1,p95_cols[,2]],apply(p95_cols,1, function(x) (min(lod_tab[1,x]) )),apply(p95_cols,1, function(x) (max(lod_tab[1,x]) )),ylab=expression(paste("log2(OR)")),main=main_title,xaxt="n",xlab="",pch=15,bg="white",xlim=c(1-dist - 0.1,length(all_cnt) + dist + 0.1), ylim=c(min(lod_tab[,p95_cols]),max(lod_tab[,p95_cols])),cex=1.5)
	
	errbar(1:length(all_cnt),lod_tab[2,p95_cols[,2]],apply(p95_cols,1, function(x) (min(lod_tab[2,x]) )),apply(p95_cols,1, function(x) (max(lod_tab[2,x]) )),ylab="log2 OR",main="",xaxt="n",yaxt="n",xlab="",add=T, pch=".",bg="white",cex=1.5)
	points(1:length(all_cnt),lod_tab[2,p95_cols[,2]],pch=22,bg="white",cex=1.5)
	
	errbar(1:length(all_cnt) + dist,lod_tab[3,p95_cols[,2]],apply(p95_cols,1, function(x) (min(lod_tab[3,x]) )),apply(p95_cols,1, function(x) (max(lod_tab[3,x]) )),ylab="log2 OR",main="",xaxt="n",yaxt="n",xlab="",add=T, pch=16,bg="white",cex=1.5)
	axis(1, labels=F,at=1:length(all_cnt),label=xtext)
	text(1:length(all_cnt), par("usr")[3] - 0.4, srt = 0, adj = 0.5, labels = xlabs, xpd = TRUE)
	text((1+length(all_cnt))/2.0, par("usr")[3] - 0.8, srt = 0, adj = 0.5, labels = xtext, xpd = TRUE)
	abline(h=0,lt=3)
	legend(leg,print_cat,pch=c(15,22,16),pt.cex=1.5)
}

all_cnt=c("N","O","T");perc=c("_2.5","","_97.5")
p95_cols=NULL; for (i in all_cnt) { p95_cols = rbind(p95_cols,paste(i,perc,sep="")) }
print_chr=c("X","2L","2R","3L","3R","total")

# get total number of SNPS on chromosomes
ita_tot_snps=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita2011_fdr_snps_per_chrom.txt",header=F)
rownames(ita_tot_snps)=as.character(ita_tot_snps[,1])
colnames(ita_tot_snps)=c("CHR","Tot")
ita_tot_snps$CHR=NULL
ita_tot_snps["total",]=sum(ita_tot_snps[chr_names,])

vie_tot_snps=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_Poly/CMH/vie2010_fdr_snps_per_chrom.txt",header=F)
rownames(vie_tot_snps)=as.character(vie_tot_snps[,1])
colnames(vie_tot_snps)=c("CHR","Tot")
vie_tot_snps$CHR=NULL
vie_tot_snps["total",]=sum(vie_tot_snps[chr_names,])

ovr_tot_snps=read.table("overlap_fdr_snps_per_chrom.txt",header=F)
rownames(ovr_tot_snps)=as.character(ovr_tot_snps[,1])
colnames(ovr_tot_snps)=c("CHR","Tot")
ovr_tot_snps$CHR=NULL
ovr_tot_snps["total",]=sum(ovr_tot_snps[chr_names,])

ovr_smpl_ap=read.table("common_snps_dsim_afr_flo_port_rnd_sample_1e5_tot.out",header=T,skip=1)
ovr_smpl_ap=clean_smpl_ap(ovr_smpl_ap)
ovr_fdr_ap=add_totals_ap_tab(ovr_fdr_ap)
ovr_fdr_ap$Totals =  ovr_tot_snps[rownames(ovr_fdr_ap),] 
ovr_lod_tab=get_lod_tab(ovr_fdr_ap,ovr_smpl_ap)


vie_smpl_ap=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/vie2010_dsim_afr_flo_port_rnd_sample_1e5_tot.out",header=T,skip=1)
vie_smpl_ap=clean_smpl_ap(vie_smpl_ap)
vie_fdr_ap=add_totals_ap_tab(vie_fdr_ap)
vie_fdr_ap$Totals =  vie_tot_snps[rownames(vie_fdr_ap),] 
vie_lod_tab=get_lod_tab(vie_fdr_ap,vie_smpl_ap)

ita_smpl_ap=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita2011_dsim_afr_flo_port_rnd_sample_1e5_tot.out",header=T,skip=1)
ita_smpl_ap=clean_smpl_ap(ita_smpl_ap)
ita_fdr_ap=add_totals_ap_tab(ita_fdr_ap)
ita_fdr_ap$Totals =  ita_tot_snps[rownames(ita_fdr_ap),] 
ita_lod_tab=get_lod_tab(ita_fdr_ap,ita_smpl_ap)


quartz()
error_bar_plot_all(vie_lod_tab,main_title="Vienna 2010")
dev.copy2pdf(file="vie2010_shared_poly_enrich_all_chrom.pdf")
quartz()
error_bar_plot_all(ita_lod_tab,main_title="Bolzano 2011")
dev.copy2pdf(file="ita2011_shared_poly_enrich_all_chrom.pdf")
quartz()
error_bar_plot_total(vie_lod_tab, ita_lod_tab,ovr_lod_tab,leg="topleft")
dev.copy2pdf(file="vie_ita_overl_shared_poly_enrich_totals.pdf")



# old files
setwd("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/Joined_Analysis")
vie_fdr_ap=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/vie2010_fdr_simulans_polym_maf_ge_0.05.af.count",header=T)
vie_ap=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/vie2010_simulans_polym_maf_ge_0.05.af.count",header=T)
ita_fdr_ap=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita2011_fdr_simulans_polym_maf_ge_0.05.af.count",header=T)
ita_ap=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita2011_simulans_polym_maf_ge_0.05.af.count",header=T)
ovr_fdr_ap=read.table("vie2010_ita2011_fdr_overlap_anc_polym_maf_ge_0.05.af.count",header=T)
ovr_ap=read.table("vie2010_ita2011_overlap_anc_polym_maf_ge_0.05.af.count",header=T)
ita_all_only=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita2011_simulans_polym_maf_ge_0.05.all_only",header=F)


vie_fdr_ap$Tot=apply(vie_fdr_ap[,-1],1,sum)
vie_fdr_ap_ratios=round(vie_fdr_ap[,2:(length(vie_fdr_ap)-1)]/vie_fdr_ap$Tot, digits=3)
rownames(vie_fdr_ap_ratios)=as.character(vie_fdr_ap$CHR)
vie_fdr_ap_ratios$None=vie_fdr_ap_ratios$X0+vie_fdr_ap_ratios$X00
vie_fdr_ap_ratios$One=vie_fdr_ap_ratios$X10+vie_fdr_ap_ratios$X20+vie_fdr_ap_ratios$X1
vie_fdr_ap_ratios$Poly=vie_fdr_ap_ratios$X12+vie_fdr_ap_ratios$X21
vie_fdr_ap_ratios=rbind(vie_fdr_ap_ratios["X",],vie_fdr_ap_ratios[1:5,])
quartz()
barplot(t(as.matrix(vie_fdr_ap_ratios[,8:10])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main="Significant SNPs Vienna, 2010")
legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="vie_anc_poly_fdr_barplot.pdf")

vie_ap$Tot=apply(vie_ap[,-1],1,sum)
vie_ap_ratios=round(vie_ap[,2:(length(vie_ap)-1)]/vie_ap$Tot, digits=3)
rownames(vie_ap_ratios)=as.character(vie_ap$CHR)
vie_ap_ratios$None=vie_ap_ratios$X0+vie_ap_ratios$X00
vie_ap_ratios$One=vie_ap_ratios$X10+vie_ap_ratios$X20+vie_ap_ratios$X1
vie_ap_ratios$Poly=vie_ap_ratios$X12+vie_ap_ratios$X21
vie_ap_ratios=rbind(vie_ap_ratios["X",],vie_ap_ratios[1:5,])
quartz()
barplot(t(as.matrix(vie_ap_ratios[,8:10])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main="All SNPs Vienna, 2010")
#legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="vie_anc_poly_all_barplot.pdf")

ita_fdr_ap$Tot=apply(ita_fdr_ap[,-1],1,sum)
ita_fdr_ap_ratios=round(ita_fdr_ap[,2:(length(ita_fdr_ap)-1)]/ita_fdr_ap$Tot, digits=3)
rownames(ita_fdr_ap_ratios)=as.character(ita_fdr_ap$CHR)
ita_fdr_ap_ratios$None=ita_fdr_ap_ratios$X0+ita_fdr_ap_ratios$X00
ita_fdr_ap_ratios$One=ita_fdr_ap_ratios$X10+ita_fdr_ap_ratios$X20+ita_fdr_ap_ratios$X1
ita_fdr_ap_ratios$Poly=ita_fdr_ap_ratios$X12+ita_fdr_ap_ratios$X21
ita_fdr_ap_ratios=rbind(ita_fdr_ap_ratios["X",],ita_fdr_ap_ratios[1:5,])
quartz()
barplot(t(as.matrix(ita_fdr_ap_ratios[,8:10])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main="Significant SNPs Bolzano, 2011")
legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="ita_anc_poly_fdr_barplot.pdf")

ita_ap$Tot=apply(ita_ap[,-1],1,sum)
ita_ap_ratios=round(ita_ap[,2:(length(ita_ap)-1)]/ita_ap$Tot, digits=3)
rownames(ita_ap_ratios)=as.character(ita_ap$CHR)
ita_ap_ratios$None=ita_ap_ratios$X0+ita_ap_ratios$X00
ita_ap_ratios$One=ita_ap_ratios$X10+ita_ap_ratios$X20+ita_ap_ratios$X1
ita_ap_ratios$Poly=ita_ap_ratios$X12+ita_ap_ratios$X21
ita_ap_ratios=rbind(ita_ap_ratios["X",],ita_ap_ratios[1:5,])
quartz()
barplot(t(as.matrix(ita_ap_ratios[,8:10])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main="All SNPs Bolzano, 2011")
#legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="ita_anc_poly_all_barplot.pdf")

ovr_fdr_ap$Tot=apply(ovr_fdr_ap[,-1],1,sum)
ovr_fdr_ap_ratios=round(ovr_fdr_ap[,2:(length(ovr_fdr_ap)-1)]/ovr_fdr_ap$Tot, digits=3)
rownames(ovr_fdr_ap_ratios)=as.character(ovr_fdr_ap$CHR)
ovr_fdr_ap_ratios$None=ovr_fdr_ap_ratios$X0+ovr_fdr_ap_ratios$X00
ovr_fdr_ap_ratios$One=ovr_fdr_ap_ratios$X10+ovr_fdr_ap_ratios$X20+ovr_fdr_ap_ratios$X1
ovr_fdr_ap_ratios$Poly=ovr_fdr_ap_ratios$X12+ovr_fdr_ap_ratios$X21
ovr_fdr_ap_ratios=rbind(ovr_fdr_ap_ratios["X",],ovr_fdr_ap_ratios[1:5,])
quartz()
barplot(t(as.matrix(ovr_fdr_ap_ratios[,8:10])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main="Significant SNPs Overlap Vienna & Bolzano")
legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="ovr_anc_poly_fdr_barplot.pdf")

ovr_ap$Tot=apply(ovr_ap[,-1],1,sum)
ovr_ap_ratios=round(ovr_ap[,2:(length(ovr_ap)-1)]/ovr_ap$Tot, digits=3)
rownames(ovr_ap_ratios)=as.character(ovr_ap$CHR)
ovr_ap_ratios$None=ovr_ap_ratios$X0+ovr_ap_ratios$X00
ovr_ap_ratios$One=ovr_ap_ratios$X10+ovr_ap_ratios$X20+ovr_ap_ratios$X1
ovr_ap_ratios$Poly=ovr_ap_ratios$X12+ovr_ap_ratios$X21
ovr_ap_ratios=rbind(ovr_ap_ratios["X",],ovr_ap_ratios[1:5,])
quartz()
barplot(t(as.matrix(ovr_ap_ratios[,8:10])),beside=T,col=c("coral3","cornsilk3","aquamarine3"),ylab="frequency",ylim=c(0,1),main="All SNPs Overlap Vienna & Bolzano")
#legend("topleft",c("none","one","both"),fill=c("coral3","cornsilk3","aquamarine3"),title="Shared alleles")
dev.copy2pdf(file="ovr_anc_poly_all_barplot.pdf")

setwd("/Volumes/Temp/Lukas/Sim_Mel_align/mauve_alignments")
mauve_algn=read.table("/Volumes/Temp/Lukas/Sim_Mel_align/mauve_alignments/Dsim_M252_draft_4_3R_align_parse.conserved",header=T)
mauve_ita_snps_algn=read.table("/Volumes/Temp/Lukas/Sim_Mel_align/mauve_alignments/Dmel_Dsim_M252_draft_4_3R_ita2011.conserved",header=T)

quartz()

plot(mauve_algn$win_center/10^6,mauve_algn$mismatch,col="blue",type="l",lty=1,xlab="bps [MB]",ylab="mismatched positions in D. sim./mel. aligment [percent]", main="3R", ylim=c(0,30))
lines(mauve_ita_snps_algn$win_center/10^6,mauve_ita_snps_algn$mismatch,col="red",type="l",lty=2)
legend("topleft",c("all positions","SNPs Italy 2011"),lty=c(1,2), col=c("blue","red"))
dev.copy2pdf(file="3R_mismatches_algn_itaSNPs.pdf")

ita_all_only=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita2011_simulans_polym_maf_ge_0.05.all_only",header=F,colClasses=c("factor",NA,"factor","factor","factor"))
colnames(ita_all_only) = c("CHR","BPS","A_MEL","A_SIM","N_SIM","A_CONS")
get_window_avgs <- function(chrom,winsize=500000) {
	wins=0:floor(max(chrom$BPS)/winsize)
	win_res=c()
	for (win in wins){
		start = win*winsize
		stop = start+winsize
		monomorph=summary(chrom$A_CONS[chrom$N_SIM == 1 & chrom$BPS > start & chrom$BPS <= stop ])[c("00","10","01")]
		polymorph=as.vector(summary(chrom$A_CONS[chrom$N_SIM == 2 & chrom$BPS > start & chrom$BPS <= stop ]))
		win_res=rbind(win_res,c((start+stop)/2,sum(monomorph)+sum(polymorph),sum(monomorph),sum(polymorph), monomorph, polymorph )	)
		
		}
	win_res=data.frame(win_res)
	colnames(win_res)=c("win","tot","mono_tot","poly_tot","mono_none","mono_maj","mono_minor","poly_00","poly_01","poly_02","poly_10","poly_12","poly_20","poly_21")
	return(win_res)
	}
ita_3R_win=get_window_avgs(ita_all_only[ita_all_only$CHR=="3R",],winsize=100000)
quartz()

with(ita_3R_win,plot(win/10^6,mono_none/tot,col="red",type="l",lty=1,ylim=c(0,1.0)))
with(ita_3R_win,lines(win/10^6,mono_maj/tot,col="blue"))
with(ita_3R_win,lines(win/10^6,mono_minor/tot,col="green"))
quartz()
with(ita_3R_win,plot(win/10^6,poly_00/poly_tot,col="red",type="l",lty=1,ylim=c(0,1)))
with(ita_3R_win,lines(win/10^6,(poly_10+poly_20)/poly_tot,col="blue",lty=1))
with(ita_3R_win,lines(win/10^6,(poly_01+poly_02)/poly_tot,col="green",lty=1))
with(ita_3R_win,lines(win/10^6,(poly_21+poly_12)/poly_tot,col="red",lty=2))

