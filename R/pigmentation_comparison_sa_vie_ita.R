setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Joined_SA_Ita")
#setwd("/Volumes/vetgrid10/Data/A7_pigmentation/Joined_SA_Ita")
vie_ita_sa_odds=read.table("vie_LD_ita_LD_eu_LD_sa_VLVD_eusa_VLVD_merged.pV_OR_filt.gz",na.strings = c("NA","nan"),header=F)
vie_ita_sa_odds=read.table("vie_LD_ita_LD_eu_LD_sa_VLVD_eusa_VLVD_merged.pV_OR_filt_3L_1_1.2_X9.1_9.2Mb.gz",na.strings = c("NA","nan"),header=F)
colnames(vie_ita_sa_odds)=c("CHR","BPS","Allele","P_vi","O_vi","Oa_vi","Ob_vi","P_ita","O_ita","Oa_ita","Ob_ita","P_eu","O_eu","Oa_eu","Ob_eu","P_sa","O_sa","Oa_sa","Ob_sa","P_es","O_es","Oa_es","Ob_es")

vie_ita_sa_odds=read.table("vie_LD_ita_LD_eu_LD_sa_VLVD_eusa_VLVD_merged.pV_OR_filt.pV_lt_e_6.gz",na.strings = c("NA","nan"),header=F)
vie_ita_sa_odds=read.table("vie_LD_ita_LD_eu_LD_sa_VLVD_eusa_VLVD_SA_XL_XD_merged.pV_OR_filt.pvlte1_e1.gz",na.strings = c("NA","nan"),header=F)
#xs
colnames(vie_ita_sa_odds)=c("CHR","BPS","Allele","P_vi","O_vi","Oa_vi","Ob_vi","P_ita","O_ita","Oa_ita","Ob_ita","P_eu","O_eu","Oa_eu","Ob_eu","P_sa","O_sa","Oa_sa","Ob_sa","P_es","O_es","Oa_es","Ob_es","P_xs","O_xs","Oa_xs","Ob_xs")
vie_ita_sa_odds$lO_eu=log(vie_ita_sa_odds$O_eu)
vie_ita_sa_odds$lO_sa=log(vie_ita_sa_odds$O_sa)
vie_ita_sa_odds$lO_vi=log(vie_ita_sa_odds$O_vi)
vie_ita_sa_odds$lO_ita=log(vie_ita_sa_odds$O_ita)
vie_ita_sa_odds$rP_eu=(vie_ita_sa_odds$P_eu)
vie_ita_sa_odds$rP_sa=(vie_ita_sa_odds$P_sa)
vie_ita_sa_odds$rP_vi=(vie_ita_sa_odds$P_vi)
vie_ita_sa_odds$rP_ita=(vie_ita_sa_odds$P_ita)
#xs
vie_ita_sa_odds$lO_xs=log(vie_ita_sa_odds$O_xs)
vie_ita_sa_odds$rP_xs=(vie_ita_sa_odds$P_xs)
vie_ita_sa_odds$P_xs=-log10(vie_ita_sa_odds$P_xs)
vie_ita_sa_odds$Oa_xs=log(vie_ita_sa_odds$Oa_xs)
vie_ita_sa_odds$Ob_xs=log(vie_ita_sa_odds$Ob_xs)
vie_ita_sa_odds$sO_xs=sign(vie_ita_sa_odds$Ob_xs) * ( sign(vie_ita_sa_odds$Ob_xs) == sign(vie_ita_sa_odds$Oa_xs) )
# log CI
vie_ita_sa_odds$Oa_eu=log(vie_ita_sa_odds$Oa_eu)
vie_ita_sa_odds$Oa_sa=log(vie_ita_sa_odds$Oa_sa)
vie_ita_sa_odds$Oa_vi=log(vie_ita_sa_odds$Oa_vi)
vie_ita_sa_odds$Oa_ita=log(vie_ita_sa_odds$Oa_ita)
vie_ita_sa_odds$Ob_eu=log(vie_ita_sa_odds$Ob_eu)
vie_ita_sa_odds$Ob_sa=log(vie_ita_sa_odds$Ob_sa)
vie_ita_sa_odds$Ob_vi=log(vie_ita_sa_odds$Ob_vi)
vie_ita_sa_odds$Ob_ita=log(vie_ita_sa_odds$Ob_ita)
# CI sign
vie_ita_sa_odds$sO_eu=sign(vie_ita_sa_odds$Ob_eu) * ( sign(vie_ita_sa_odds$Ob_eu) == sign(vie_ita_sa_odds$Oa_eu) )
vie_ita_sa_odds$sO_sa=sign(vie_ita_sa_odds$Ob_sa) * ( sign(vie_ita_sa_odds$Ob_sa) == sign(vie_ita_sa_odds$Oa_sa) )
vie_ita_sa_odds$sO_vi=sign(vie_ita_sa_odds$Ob_vi) * ( sign(vie_ita_sa_odds$Ob_vi) == sign(vie_ita_sa_odds$Oa_vi) )
vie_ita_sa_odds$sO_ita=sign(vie_ita_sa_odds$Ob_ita) * ( sign(vie_ita_sa_odds$Ob_ita) == sign(vie_ita_sa_odds$Oa_ita) )

length(vie_ita_sa_odds$sO_xs[(vie_ita_sa_odds$sO_xs * vie_ita_sa_odds$sO_eu) > 0 ])
length(vie_ita_sa_odds$sO_xs[abs(vie_ita_sa_odds$sO_xs) == 0 ])

vie_ita_sa_odds$P_eu=-log10(vie_ita_sa_odds$P_eu)
vie_ita_sa_odds$P_sa=-log10(vie_ita_sa_odds$P_sa)
vie_ita_sa_odds$P_vi=-log10(vie_ita_sa_odds$P_vi)
vie_ita_sa_odds$P_ita=-log10(vie_ita_sa_odds$P_ita)

has_Peu = ! (is.na(vie_ita_sa_odds$P_vi) | is.na(vie_ita_sa_odds$P_ita))
has_Pv = ! (is.na(vie_ita_sa_odds$P_sa) | is.na(vie_ita_sa_odds$P_eu) | is.na(vie_ita_sa_odds$P_xs))
has_OR = ! (is.na(vie_ita_sa_odds$O_sa) | is.na(vie_ita_sa_odds$O_eu) | is.na(vie_ita_sa_odds$O_xs))
summary(vie_ita_sa_odds$lO_eu[is.finite(vie_ita_sa_odds$lO_eu)]) # Min: -4.058, Max: 4.125
summary(vie_ita_sa_odds$lO_sa[is.finite(vie_ita_sa_odds$lO_sa)]) # Min: -4.2, Max: 3.827
summary(vie_ita_sa_odds$lO_xs[is.finite(vie_ita_sa_odds$lO_xs)]) # Min: -3.8, Max: 3.6
# set inf: 5.0, -inf: -5.0 (for pearson and such)
## vie_ita_sa_odds$lO_eu_ni=vie_ita_sa_odds$lO_eu
## vie_ita_sa_odds$lO_sa_ni=vie_ita_sa_odds$lO_sa
vie_ita_sa_odds$lO_eu[vie_ita_sa_odds$lO_eu == Inf] = 5
vie_ita_sa_odds$lO_eu[vie_ita_sa_odds$lO_eu == -Inf] = -5
vie_ita_sa_odds$lO_sa[vie_ita_sa_odds$lO_sa == Inf] = 5
vie_ita_sa_odds$lO_sa[vie_ita_sa_odds$lO_sa == -Inf] = -5
vie_ita_sa_odds$lO_vi[vie_ita_sa_odds$lO_vi == Inf] = 5
vie_ita_sa_odds$lO_vi[vie_ita_sa_odds$lO_vi == -Inf] = -5
vie_ita_sa_odds$lO_ita[vie_ita_sa_odds$lO_ita == Inf] = 5
vie_ita_sa_odds$lO_ita[vie_ita_sa_odds$lO_ita == -Inf] = -5
vie_ita_sa_odds$lO_xs[vie_ita_sa_odds$lO_xs == Inf] = 5
vie_ita_sa_odds$lO_xs[vie_ita_sa_odds$lO_xs == -Inf] = -5
# number of same alleles:  2132620
length(vie_ita_sa_odds$lO_sa[has_Pv])
# number of equal sign on odds ratio: 1068935
length(vie_ita_sa_odds$lO_sa[  sign(vie_ita_sa_odds$lO_sa) == sign(vie_ita_sa_odds$lO_eu) & has_Pv ]) # 921815
length(vie_ita_sa_odds$lO_sa[  sign(vie_ita_sa_odds$lO_sa) == sign(vie_ita_sa_odds$lO_eu) & has_Pv &  (is.finite(vie_ita_sa_odds$lO_eu) & is.finite(vie_ita_sa_odds$lO_eu))]) # 1068935
# how many are zero and equal?
length(vie_ita_sa_odds$lO_sa[  sign(vie_ita_sa_odds$lO_sa) == sign(vie_ita_sa_odds$lO_eu) & sign(vie_ita_sa_odds$lO_eu) == 0 & has_Pv ]) # 0

length(vie_ita_sa_odds$lO_sa[  sign(vie_ita_sa_odds$lO_sa) == sign(vie_ita_sa_odds$lO_eu) & vie_ita_sa_odds$P_eu >= 8 & has_Pv] ) # 201
length(vie_ita_sa_odds$lO_sa[  sign(vie_ita_sa_odds$lO_sa) != sign(vie_ita_sa_odds$lO_eu) & vie_ita_sa_odds$P_eu >= 8 & has_Pv] ) # 121

length(vie_ita_sa_odds$lO_sa[  sign(vie_ita_sa_odds$lO_sa) == sign(vie_ita_sa_odds$lO_eu) & vie_ita_sa_odds$P_sa >= 8   & has_Pv] ) # 415
length(vie_ita_sa_odds$lO_sa[  sign(vie_ita_sa_odds$lO_sa) != sign(vie_ita_sa_odds$lO_eu) & vie_ita_sa_odds$P_sa >= 8  & has_Pv ] ) # 308


length(vie_ita_sa_odds$lO_sa[  sign(vie_ita_sa_odds$lO_sa) == sign(vie_ita_sa_odds$lO_eu) & vie_ita_sa_odds$P_eu >= 12  & has_Pv] ) # 56
length(vie_ita_sa_odds$lO_sa[  sign(vie_ita_sa_odds$lO_sa) != sign(vie_ita_sa_odds$lO_eu) & vie_ita_sa_odds$P_eu >= 12  & has_Pv] ) # 25
length(vie_ita_sa_odds$lO_sa[  sign(vie_ita_sa_odds$lO_sa) == sign(vie_ita_sa_odds$lO_eu) & vie_ita_sa_odds$P_sa >= 12  & has_Pv] ) # 70
length(vie_ita_sa_odds$lO_sa[  sign(vie_ita_sa_odds$lO_sa) != sign(vie_ita_sa_odds$lO_eu) & vie_ita_sa_odds$P_sa >= 12  & has_Pv] ) # 44
length(vie_ita_sa_odds$lO_sa[  sign(vie_ita_sa_odds$lO_sa) == sign(vie_ita_sa_odds$lO_eu) & vie_ita_sa_odds$P_sa >= 20  & has_Pv] ) # 11
length(vie_ita_sa_odds$lO_sa[  sign(vie_ita_sa_odds$lO_sa) != sign(vie_ita_sa_odds$lO_eu) & vie_ita_sa_odds$P_sa >= 20  & has_Pv] ) # 11

# spearman between all lods: p<=2e-16, 0.007987311
cor.test(vie_ita_sa_odds$lO_sa[has_Pv],vie_ita_sa_odds$lO_eu[has_Pv],method="spearman")
# spearman between P_eu < 3:  p-value = 7e-8; rho 0.02902994
cor.test(vie_ita_sa_odds$lO_sa[vie_ita_sa_odds$P_eu >= 3 & has_Pv],vie_ita_sa_odds$lO_eu[vie_ita_sa_odds$P_eu >= 3 & has_Pv],method="spearman")
# spearman between P_sa < 3:  p-value  2e-4 rho 0.01540932
cor.test(vie_ita_sa_odds$lO_sa[vie_ita_sa_odds$P_sa >= 3 & has_Pv],vie_ita_sa_odds$lO_eu[vie_ita_sa_odds$P_sa >= 3 & has_Pv],method="spearman")
# spearman between P_eu < 8:  p-value = 1.2e-08 rho=0.3108454
cor.test(vie_ita_sa_odds$lO_sa[vie_ita_sa_odds$P_eu >= 8 & has_Pv],vie_ita_sa_odds$lO_eu[vie_ita_sa_odds$P_eu >= 8 & has_Pv],method="spearman")
# spearman between P_sa < 8:  p-value = 2.8e-11 rho 0.2441134
cor.test(vie_ita_sa_odds$lO_sa[vie_ita_sa_odds$P_sa >= 8 & has_Pv],vie_ita_sa_odds$lO_eu[vie_ita_sa_odds$P_sa >= 8 & has_Pv],method="spearman")
# spearman between P_eu < 12:  p-value = 3.408e-06 rho=0.4960705
cor.test(vie_ita_sa_odds$lO_sa[vie_ita_sa_odds$P_eu >= 12 & has_Pv],vie_ita_sa_odds$lO_eu[vie_ita_sa_odds$P_eu >= 12 & has_Pv],method="spearman")
# spearman between P_sa < 12:  p-value =  1.231e-05 rho 0.4000284
cor.test(vie_ita_sa_odds$lO_sa[vie_ita_sa_odds$P_sa >= 12 & has_Pv],vie_ita_sa_odds$lO_eu[vie_ita_sa_odds$P_sa >= 12 & has_Pv],method="spearman")


#pearson between P_eu < 8:  p-value = 3e-03 rho=0.1613834
cor.test(vie_ita_sa_odds$lO_sa[vie_ita_sa_odds$P_eu >= 8 & has_Pv],vie_ita_sa_odds$lO_eu[vie_ita_sa_odds$P_eu >= 8 & has_Pv],method="pearson")
#pearson between P_sa < 8:  p-value = 3e-10 rho 0.2314178
cor.test(vie_ita_sa_odds$lO_sa[vie_ita_sa_odds$P_sa >= 8 & has_Pv],vie_ita_sa_odds$lO_eu[vie_ita_sa_odds$P_sa >= 8 & has_Pv])

#pearson between P_eu < 12:  p-value = 0.006178 rho=0.3018016
cor.test(vie_ita_sa_odds$lO_sa[vie_ita_sa_odds$P_eu >= 12 & has_Pv],vie_ita_sa_odds$lO_eu[vie_ita_sa_odds$P_eu >= 12 & has_Pv],method="pearson")
#pearson between P_sa < 12:  p-value = 2e-04 rho 0.3403004
cor.test(vie_ita_sa_odds$lO_sa[vie_ita_sa_odds$P_sa >= 12 & has_Pv],vie_ita_sa_odds$lO_eu[vie_ita_sa_odds$P_sa >= 12 & has_Pv])


# manhattan plots

eu_plot=vie_ita_sa_odds[vie_ita_sa_odds$rP_eu <= 0.1 & ! is.na(vie_ita_sa_odds$P_eu),c("CHR","BPS","rP_eu")]
sa_plot=vie_ita_sa_odds[vie_ita_sa_odds$rP_sa <= 0.1 & ! is.na(vie_ita_sa_odds$P_sa),c("CHR","BPS","rP_sa")]
vie_plot=vie_ita_sa_odds[vie_ita_sa_odds$rP_vi <= 0.1 & ! is.na(vie_ita_sa_odds$P_vi),c("CHR","BPS","rP_vi")]
ita_plot=vie_ita_sa_odds[vie_ita_sa_odds$rP_ita <= 0.1 & ! is.na(vie_ita_sa_odds$P_ita),c("CHR","BPS","rP_ita")]
xs_plot=vie_ita_sa_odds[vie_ita_sa_odds$rP_xs <= 0.1 & ! is.na(vie_ita_sa_odds$P_xs),c("CHR","BPS","rP_xs")]
#sa_plot=eu_sa_freqs[eu_sa_freqs$Psa <= 0.01 & ! is.na(eu_sa_freqs$Psa),c("CHR","BPS","Psa")]
colnames(eu_plot)=c("CHR","BP","P")
colnames(sa_plot)=c("CHR","BP","P")
colnames(vie_plot)=c("CHR","BP","P")
colnames(ita_plot)=c("CHR","BP","P")
colnames(xs_plot)=c("CHR","BP","P")
summary(sa_plot)
chroms=c("X","2L","2R","3L","3R","4")
off_eu=get_offset_manh(eu_plot,limitchromosomes=chroms)
off_sa=get_offset_manh(sa_plot,limitchromosomes=chroms)
off_vie=get_offset_manh(vie_plot,limitchromosomes=chroms)
off_ita=get_offset_manh(ita_plot,limitchromosomes=chroms)
off_xs=get_offset_manh(xs_plot,limitchromosomes=chroms)
## eu_plot= eu_plot[! (eu_plot$CHR == "3L" & eu_plot$BP >= 8700000 &  eu_plot$BP <= 8750000 & -log10(eu_plot$P) > 6 ), ]
## ita_plot= ita_plot[!(ita_plot$CHR == "3L" & ita_plot$BP >= 8700000 &  ita_plot$BP <= 8750000 & -log10(ita_plot$P) > 6), ]
## sa_plot= sa_plot[! (sa_plot$CHR == "3L" & sa_plot$BP >= 8700000 &  sa_plot$BP <= 8750000 & -log10(sa_plot$P) > 6 ), ]
## vie_plot= vie_plot[!(vie_plot$CHR == "3L" & vie_plot$BP >= 8700000 &  vie_plot$BP <= 8750000 & -log10(vie_plot$P) > 6), ]
quartz(height=6,width=12)
#x11(height=6,width=12)
#par(bg="white")
#manhattan(eu_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="Vienna (2010) and Bolzano (2011)", suggestiveline=c(22))
par.old <- par()
par(mar=c(1.0,1.0,1.0,0.5))
quartz(height=2.5,width=6)
tiff(filename = "eu_manhattan_pub.tiff",
     width = 6, height = 2.5, units = "in", res=300,
     compression = c("none"), type="cairo" ,
     bg = "white")
par(mar=c(2.0,3.0,1.5,0.5))
manhattan(eu_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="Europe", suggestiveline=c(-log10(3.07403844435e-23)),cex.main=0.75,cex.axis=0.75,cex.lab=1.0,tck=-0.05,ps_max = 0.65, ps_min = 0.40,ylab="")
mtext(expression(-log[10](italic(P))),side=2,line=1.75,cex=0.75)
rect(9111688+off_eu["X"],0,9117290+off_eu["X"],60,lty=3,lwd=1, border="red")
rect(1036369+off_eu["3L"],0,1101089+off_eu["3L"],60,lty=3,lwd=1, border="green")
rect(17062617+off_eu["3R"],0,17062900+off_eu["3R"],60,lty=3,lwd=1, border="blue")
rect(4215005+off_eu["2R"],0,4283785+off_eu["2R"],60,lty=3, border="yellow")
manhattan(eu_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="Europe", suggestiveline=c(-log10(3.07403844435e-23)),cex.main=0.75,cex.axis=0.75,cex.lab=1.0,tck=-0.05,ps_max = 0.65, ps_min = 0.40,ylab="",replot=T)

dev.off()
dev.print(png,width=800,file="eu_manhattan_pub.png")
tiff(filename = "sa_manhattan_pub.tiff",
     width = 6, height = 2.5, units = "in", res=300,
     compression = c("none"), type="cairo" ,
     bg = "white")
par(mar=c(2.0,3.0,1.5,0.5))
manhattan(sa_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="South Africa", suggestiveline=c(-log10(1.916206e-16)),cex.main=0.75,cex.axis=0.75,cex.lab=1.0,tck=-0.05,ps_max = 0.65, ps_min = 0.40,ylab="")
mtext(expression(-log[10](italic(P))),side=2,line=1.75,cex=0.75)
rect(9111688+off_sa["X"],0,9117290+off_sa["X"],60,lty=3,lwd=1, border="red")
rect(1036369+off_sa["3L"],0,1101089+off_sa["3L"],60,lty=3,lwd=1, border="green")
rect(17062617+off_sa["3R"],0,17062900+off_sa["3R"],60,lty=3,lwd=1, border="blue")
rect(4215005+off_sa["2R"],0,4283785+off_sa["2R"],60,lty=3, border="yellow")
manhattan(sa_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="South Africa", suggestiveline=c(-log10(1.916206e-16)),cex.main=0.75,cex.axis=0.75,cex.lab=1.0,tck=-0.05,ps_max = 0.65, ps_min = 0.40,ylab="",replot=T)
dev.off()
dev.print(png,width=800,file="sa_manhattan_pub.png")
quartz(height=6,width=12)
par(bg="white")
chroms=c("X","2L","2R","3L","3R","4")
manhattan(sa_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="South Africa", suggestiveline=c(-log10(1.916206e-16)))
rect(9111688+off_sa["X"],0,9170829+off_sa["X"],60,lty=3,lwd=1.5, border="red")
rect(1071314+off_sa["3L"],0,1120551+off_sa["3L"],60,lty=3,lwd=1.5, border="green")
rect(17055975+off_sa["3R"],0,17069171+off_sa["3R"],60,lty=3,lwd=1.5, border="blue")
#rect(4179149+off_sa["2R"],0,4319221+off_sa["2R"],60,lty=3, border="yellow")
dev.print(png,width=800,file="sa_manhattan_pub.png")
quartz(height=6,width=12)
par(bg="white")
chroms=c("X","2L","2R","3L","3R","4")
manhattan(xs_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="South Africa (2012, extended extremes)", suggestiveline=c(-log10(2.054127e-19)))
rect(9111688+off_xs["X"],0,9170829+off_xs["X"],60,lty=3,lwd=1.5, border="red")
rect(1071314+off_xs["3L"],0,1120551+off_xs["3L"],60,lty=3,lwd=1.5, border="green")
rect(17055975+off_xs["3R"],0,17069171+off_xs["3R"],60,lty=3,lwd=1.5, border="blue")
#rect(4179149+off_xs["2R"],0,4319221+off_xs["2R"],60,lty=3, border="yellow")
dev.print(png,width=800,file="xs_manhattan_pub.png")

quartz(height=6,width=12)
par(bg="white")
manhattan(vie_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="Vienna (2010)", suggestiveline=c(-log10(2.13820160694e-15)))
rect(8867539+off_vie["X"],0,9170829+off_vie["X"],60,lty=3, border="red")
rect(1071314+off_vie["3L"],0,1120551+off_vie["3L"],60,lty=3, border="green")
rect(17055975+off_vie["3R"],0,17069171+off_vie["3R"],60,lty=3, border="blue")
rect(4179149+off_vie["2R"],0,4319221+off_vie["2R"],60,lty=3, border="yellow")
dev.print(png,width=800,file="vie_manhattan.png")
quartz(height=6,width=12)
par(bg="white")
chroms=c("X","2L","2R","3L","3R","4")
manhattan(ita_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="Bolzano (2011)", suggestiveline=c(-log10(2.054127e-19)))
rect(8867539+off_ita["X"],0,9170829+off_ita["X"],60,lty=3, border="red")
rect(1071314+off_ita["3L"],0,1120551+off_ita["3L"],60,lty=3, border="green")
rect(17055975+off_ita["3R"],0,17069171+off_ita["3R"],60,lty=3, border="blue")
rect(4179149+off_ita["2R"],0,4319221+off_ita["2R"],60,lty=3, border="yellow")
dev.print(png,width=800,file="ita_manhattan.png")

a5_vie=read.table("../Vienna_A5/CMH/females_A5_vie_pA5L1_pA5L2_pA5D1_pA5D2_q20_filt_mc2_chroms_only_mct5_mcv15_nlpl1.gwas",header=T)
chroms=c("X","2L","2R","3L","3R","4")
off_a5=get_offset_manh(a5_vie,limitchromosomes=chroms)
quartz(height=6,width=12)
par(bg="white")
manhattan(a5_vie,colors=c("black","slategrey"),limitchromosomes=chroms,main="Vienna (2010) A5 tergite", suggestiveline=c(-log10(2.13820160694e-8)))
rect(8867539+off_a5["X"],0,9170829+off_a5["X"],60,lty=3, border="red")
rect(1071314+off_a5["3L"],0,1120551+off_a5["3L"],60,lty=3, border="green")
rect(17055975+off_a5["3R"],0,17069171+off_a5["3R"],60,lty=3, border="blue")
rect(4179149+off_a5["2R"],0,4319221+off_a5["2R"],60,lty=3, border="yellow")
dev.print(png,width=800,file="a5_vie_manhattan.png")

# correlation plots


## quants=1-10^-(c(seq(0.1,2.6,0.5),seq(3,4.0,0.5),seq(4.15,5.5,0.25)))
## eu_q = quantile( vie_ita_sa_odds$P_eu[has_Pv],probs=quants)
## sa_q = quantile( vie_ita_sa_odds$P_sa[has_Pv],probs=quants)


## sa_cor_sample=read.table("sample_lOR_EU_SA_pV_SA_n_10K.out",na.strings = c("NA","nan"),header=T)
## eu_cor_sample=read.table("sample_lOR_EU_SA_pV_EU_n_10K.out",na.strings = c("NA","nan"),header=T)
## eu.xs_cor_sample=read.table("sample_lOR_EU_xSA_pV_EU_n_10K.out",na.strings = c("NA","nan"),header=T)
## xs.eu_cor_sample=read.table("sample_lOR_EU_xSA_pV_xSA_n_10K.out",na.strings = c("NA","nan"),header=T)
sa_cor_sample=read.table("sample_lOR_EU_SA_pV_SA_n_1K.out",na.strings = c("NA","nan"),header=T)
eu_cor_sample=read.table("sample_lOR_EU_SA_pV_EU_n_1K.out",na.strings = c("NA","nan"),header=T)
eu.xs_cor_sample=read.table("sample_lOR_EU_xSA_pV_EU_n_1K.out",na.strings = c("NA","nan"),header=T)
xs.eu_cor_sample=read.table("sample_lOR_EU_xSA_pV_xSA_n_1K.out",na.strings = c("NA","nan"),header=T)
sa_cor=read.table("sample_lOR_EU_SA_pV_SA_n_1K.out.det",na.strings = c("NA","nan"),header=T)
eu_cor=read.table("sample_lOR_EU_SA_pV_EU_n_1K.out.det",na.strings = c("NA","nan"),header=T)
eu.xs_cor=read.table("sample_lOR_EU_xSA_pV_EU_n_1K.out.det",na.strings = c("NA","nan"),header=T)
xs.eu_cor=read.table("sample_lOR_EU_xSA_pV_xSA_n_1K.out.det",na.strings = c("NA","nan"),header=T)
q.disp=c(10,100,1000,10000,100000,1000000,1500000)
q.disp.eu=eu_cor$pV[eu_cor$rank %in% q.disp]
q.disp.sa=sa_cor$pV[sa_cor$rank %in% q.disp]
q.disp.xs=xs.eu_cor$pV[xs.eu_cor$rank %in% q.disp]
# xs eu
quartz()
plot_cor(eu.xs_cor,xs.eu_cor,"fid_ci",q.disp.eu,q.disp.xs,quants_disp=q.disp,pop.labs=c("Europe","South Africa (ext.)","Rnd. mean","Rnd. 95% CI"),ylims=c(0.0,1.0),ylabs=c("fraction of SNPs with consistent change (95% CI)"))
lines(xs.eu_cor_sample$rank,xs.eu_cor_sample$fr_ci_id_mean,col="blue",lty=2)
lines(xs.eu_cor_sample$rank,xs.eu_cor_sample$fr_ci_id_97.5,col="blue",lty=3)
lines(xs.eu_cor_sample$rank,xs.eu_cor_sample$fr_ci_id_2.5,col="blue",lty=3)
lines(eu.xs_cor_sample$rank,eu.xs_cor_sample$fr_ci_id_mean,col="red",lty=2)
lines(eu.xs_cor_sample$rank,eu.xs_cor_sample$fr_ci_id_97.5,col="red",lty=3)
lines(eu.xs_cor_sample$rank,eu.xs_cor_sample$fr_ci_id_2.5,col="red",lty=3)
dev.copy2pdf(file="xs_eu_lOR_ci_sign_quant_all_n_1K.pdf")
quartz()
plot_cor(eu.xs_cor,xs.eu_cor,"rho",q.disp.eu,q.disp.xs,quants_disp=q.disp,pop.labs=c("Europe","South Africa (ext.)","Rnd. mean","Rnd. 95% CI"),ylims=c(0.0,1.0),ylabs=expression(paste("Spearman\'s ",rho," of ",log(OR))))
lines(xs.eu_cor_sample$rank,xs.eu_cor_sample$rho_id_mean,col="blue",lty=2)
lines(xs.eu_cor_sample$rank,xs.eu_cor_sample$rho_97.5,col="blue",lty=3)
lines(xs.eu_cor_sample$rank,xs.eu_cor_sample$rho_2.5,col="blue",lty=3)
lines(eu.xs_cor_sample$rank,eu.xs_cor_sample$rho_mean,col="red",lty=2)
lines(eu.xs_cor_sample$rank,eu.xs_cor_sample$rho_97.5,col="red",lty=3)
lines(eu.xs_cor_sample$rank,eu.xs_cor_sample$rho_2.5,col="red",lty=3)
dev.copy2pdf(file="xs_eu_lOR_rho_quant_all_n_1K.pdf")
quartz()
plot_cor(eu.xs_cor,xs.eu_cor,"fid",q.disp.eu,q.disp.xs,quants_disp=q.disp,pop.labs=c("Europe","South Africa (ext.)","Rnd. mean","Rnd. 95% CI"),ylims=c(0.4,1.0),ylabs=c("fraction of SNPs with consistent change (mean)"))
lines(xs.eu_cor_sample$rank,xs.eu_cor_sample$fr_id_mean,col="blue",lty=2)
lines(xs.eu_cor_sample$rank,xs.eu_cor_sample$fr_id_97.5,col="blue",lty=3)
lines(xs.eu_cor_sample$rank,xs.eu_cor_sample$fr_id_2.5,col="blue",lty=3)
lines(eu.xs_cor_sample$rank,eu.xs_cor_sample$fr_id_mean,col="red",lty=2)
lines(eu.xs_cor_sample$rank,eu.xs_cor_sample$fr_id_97.5,col="red",lty=3)
lines(eu.xs_cor_sample$rank,eu.xs_cor_sample$fr_id_2.5,col="red",lty=3)
dev.copy2pdf(file="xs_eu_lOR_id_sign_quant_all_n_1K.pdf")
# eu sa
quartz()
plot_cor(eu_cor,sa_cor,"fid_ci",q.disp.eu,q.disp.sa,quants_disp=q.disp,pop.labs=c("Europe","South Africa","Rnd. mean","Rnd. 95% CI"),ylims=c(0.0,1.0),ylabs=c("fraction of SNPs with consistent change (95% CI)"))
lines(sa_cor_sample$rank,sa_cor_sample$fr_ci_id_mean,col="red",lty=2)
lines(sa_cor_sample$rank,sa_cor_sample$fr_ci_id_97.5,col="red",lty=3)
lines(sa_cor_sample$rank,sa_cor_sample$fr_ci_id_2.5,col="red",lty=3)
lines(eu_cor_sample$rank,eu_cor_sample$fr_ci_id_mean,col="blue",lty=2)
lines(eu_cor_sample$rank,eu_cor_sample$fr_ci_id_97.5,col="blue",lty=3)
lines(eu_cor_sample$rank,eu_cor_sample$fr_ci_id_2.5,col="blue",lty=3)
dev.copy2pdf(file="sa_eu_lOR_ci_sign_quant_all_n_1K.pdf")
quartz()
plot_cor(eu_cor,sa_cor,"rho",q.disp.eu,q.disp.xs,quants_disp=q.disp,pop.labs=c("Europe","South Africa","Rnd. mean","Rnd. 95% CI"),ylims=c(0.0,1.0),ylabs=expression(paste("Spearman\'s ",rho," of ",log(OR))))
lines(sa_cor_sample$rank,sa_cor_sample$rho_id_mean,col="red",lty=2)
lines(sa_cor_sample$rank,sa_cor_sample$rho_97.5,col="red",lty=3)
lines(sa_cor_sample$rank,sa_cor_sample$rho_2.5,col="red",lty=3)
lines(eu_cor_sample$rank,eu_cor_sample$rho_mean,col="blue",lty=2)
lines(eu_cor_sample$rank,eu_cor_sample$rho_97.5,col="blue",lty=3)
lines(eu_cor_sample$rank,eu_cor_sample$rho_2.5,col="blue",lty=3)
dev.copy2pdf(file="sa_eu_lOR_rho_quant_all_n_1K.pdf")
quartz()
plot_cor(eu_cor,sa_cor,"fid",q.disp.eu,q.disp.sa,quants_disp=q.disp,pop.labs=c("Europe","South Africa","Rnd. mean","Rnd. 95% CI"),ylims=c(0.4,1.0),ylabs=c("fraction of SNPs with consistent change (mean)"))
lines(sa_cor_sample$rank,sa_cor_sample$fr_id_mean,col="red",lty=2)
lines(sa_cor_sample$rank,sa_cor_sample$fr_id_97.5,col="red",lty=3)
lines(sa_cor_sample$rank,sa_cor_sample$fr_id_2.5,col="red",lty=3)
lines(eu_cor_sample$rank,eu_cor_sample$fr_id_mean,col="blue",lty=2)
lines(eu_cor_sample$rank,eu_cor_sample$fr_id_97.5,col="blue",lty=3)
lines(eu_cor_sample$rank,eu_cor_sample$fr_id_2.5,col="blue",lty=3)
dev.copy2pdf(file="sa_eu_lOR_id_sign_quant_all_n_1K.pdf")

plot_cor(eu_cor,sa_cor,"fid",q.disp.eu,q.disp.sa,quants_disp=q.disp,pop.labs=c("Europe","South Africa","Rnd. mean","Rnd. 95% CI"),ylims=c(0.5,1.0),ylabs=c("fraction of SNPs with consistent change (mean)"))
lines(sa_cor_sample$rank,sa_cor_sample$fr_id_mean,col="red",lty=2)
lines(sa_cor_sample$rank,sa_cor_sample$fr_id_97.5,col="red",lty=3)
lines(sa_cor_sample$rank,sa_cor_sample$fr_id_2.5,col="red",lty=3)
lines(eu_cor_sample$rank,eu_cor_sample$fr_id_mean,col="blue",lty=2)
lines(eu_cor_sample$rank,eu_cor_sample$fr_id_97.5,col="blue",lty=3)
lines(eu_cor_sample$rank,eu_cor_sample$fr_id_2.5,col="blue",lty=3)
dev.copy2pdf(file="sa_eu_lOR_id_sign_quant_all_n_1K_b.pdf")



plot_cor <- function(taba,tabb,typ,a_q_disp,b_q_disp,quants_disp=quants_disp,pop.labs=c("Europe","South Africa (ext.)"),ylabs="fraction of SNPs with consistent change",leg.tit=c("Association"),ylims=c(0.5,1.0),leg.pos="topright"){
    par(mar=c(5, 4, 5, 4) + 0.1)
    plot(taba$rank,taba[,typ],type="l",col="red",lwd=2,xaxt='n',xlab="SNP rank", ylab=ylabs,log="x",ylim=ylims)
    lines(tabb$rank,tabb[,typ], col="blue",lwd=2)
    axis(1,at=quants_disp,labels=quants_disp, col="black", line=0.0)
    legend(leg.pos,pop.labs,lty=c(1,1,2,3),lwd=2,col=c("red","blue","black","black"),title=leg.tit)
    axis(3,col.axis="red",at=quants_disp,labels=round(a_q_disp,digits=1), col="black", line=0.0)
    axis(3,col.axis="blue",at=quants_disp,labels=round(b_q_disp,digits=1), tck=0, lty=0, line=1.0)
    mtext(expression(-log[10](italic(P))),3,line=3.0)
}



library(gplots)
quants=c(10,12,15,17,25,30,75,100,250,500,750,1000,5000,10000,50000,100000,500000,length(vie_ita_sa_odds$P_eu))
quants_disp=c(10,100,1000,10000,100000,1000000)
eu_q = sort( vie_ita_sa_odds$P_eu[has_Pv],decreasing=T)[quants]
sa_q = sort( vie_ita_sa_odds$P_sa[has_Pv],decreasing=T)[quants]
eu_q_disp = sort( vie_ita_sa_odds$P_eu[has_Pv],decreasing=T)[quants_disp]
sa_q_disp = sort( vie_ita_sa_odds$P_sa[has_Pv],decreasing=T)[quants_disp]
eu_x_q <- sapply(eu_q,function(x) {length( vie_ita_sa_odds$P_eu[has_Pv & vie_ita_sa_odds$P_eu >= x])})
sa_x_q <- sapply(sa_q,function(x) {length( vie_ita_sa_odds$P_sa[has_Pv & vie_ita_sa_odds$P_sa >= x])})
eu_y_q <- sapply(eu_q,function(x) {length( vie_ita_sa_odds$P_eu[has_Pv & vie_ita_sa_odds$P_eu >= x &  sign(vie_ita_sa_odds$lO_sa) == sign(vie_ita_sa_odds$lO_eu) ])})
sa_y_q <- sapply(sa_q,function(x) {length( vie_ita_sa_odds$P_sa[has_Pv & vie_ita_sa_odds$P_sa >= x & sign(vie_ita_sa_odds$lO_sa) == sign(vie_ita_sa_odds$lO_eu) ])})
euc_y_q <- sapply(eu_q,function(x) {length( vie_ita_sa_odds$P_eu[has_Pv & vie_ita_sa_odds$P_eu >= x &  ( vie_ita_sa_odds$sO_sa * vie_ita_sa_odds$sO_eu > 0 ) ])})
sac_y_q <- sapply(sa_q,function(x) {length( vie_ita_sa_odds$P_sa[has_Pv & vie_ita_sa_odds$P_sa >= x & (vie_ita_sa_odds$sO_sa * vie_ita_sa_odds$sO_eu > 0 ) ])})

xs_q = sort( vie_ita_sa_odds$P_xs[has_Pv],decreasing=T)[quants]
xs_q_disp = sort( vie_ita_sa_odds$P_xs[has_Pv],decreasing=T)[quants_disp]
xs_x_q <- sapply(xs_q,function(x) {length( vie_ita_sa_odds$P_xs[has_Pv & vie_ita_sa_odds$P_xs >= x])})
xs_euy_q <- sapply(xs_q,function(x) {length( vie_ita_sa_odds$P_xs[has_Pv & vie_ita_sa_odds$P_xs >= x & sign(vie_ita_sa_odds$lO_xs) == sign(vie_ita_sa_odds$lO_eu) ])})
xs_say_q <- sapply(xs_q,function(x) {length( vie_ita_sa_odds$P_xs[has_Pv & vie_ita_sa_odds$P_xs >= x & sign(vie_ita_sa_odds$lO_xs) == sign(vie_ita_sa_odds$lO_sa) ])})
sa_xsy_q <- sapply(sa_q,function(x) {length( vie_ita_sa_odds$P_sa[has_Pv & vie_ita_sa_odds$P_sa >= x & sign(vie_ita_sa_odds$lO_sa) == sign(vie_ita_sa_odds$lO_xs) ])})
eu_xsy_q <- sapply(eu_q,function(x) {length( vie_ita_sa_odds$P_eu[has_Pv & vie_ita_sa_odds$P_eu >= x & sign(vie_ita_sa_odds$lO_eu) == sign(vie_ita_sa_odds$lO_xs) ])})

xsc_euy_q <- sapply(xs_q,function(x) {length( vie_ita_sa_odds$P_xs[has_Pv & vie_ita_sa_odds$P_xs >= x & (vie_ita_sa_odds$sO_xs*vie_ita_sa_odds$sO_eu > 0) ])})
xsc_say_q <- sapply(xs_q,function(x) {length( vie_ita_sa_odds$P_xs[has_Pv & vie_ita_sa_odds$P_xs >= x & (vie_ita_sa_odds$sO_xs*vie_ita_sa_odds$sO_sa > 0) ])})
sac_xsy_q <- sapply(sa_q,function(x) {length( vie_ita_sa_odds$P_sa[has_Pv & vie_ita_sa_odds$P_sa >= x & (vie_ita_sa_odds$sO_sa*vie_ita_sa_odds$sO_xs > 0) ])})
euc_xsy_q <- sapply(eu_q,function(x) {length( vie_ita_sa_odds$P_eu[has_Pv & vie_ita_sa_odds$P_eu >= x & (vie_ita_sa_odds$sO_eu*vie_ita_sa_odds$sO_xs > 0) ])})

quartz()
#x11()
par(mar=c(5, 4, 5, 4) + 0.1)
plot(quants,xs_say_q/xs_x_q,type="l",col="red",lwd=2,xaxt='n',xlab="SNP rank", ylab="fraction of SNPs with consistent mean change",log="x",ylim=c(0.5,1.0))
lines(quants,sa_xsy_q/sa_x_q, col="blue",lwd=2)
axis(1,at=quants_disp,labels=quants_disp, col="black", line=0.0)
legend("bottomleft",c("South Africa (ext.)","South Africa"),lty=c(1,1),lwd=2,col=c("red","blue"),title=c("Association"))
axis(3,col.axis="red",at=quants_disp,labels=round(xs_q_disp,digits=1), col="black", line=0.0)
axis(3,col.axis="blue",at=quants_disp,labels=round(sa_q_disp,digits=1), tck=0, lty=0, line=1.0)
mtext(expression(-log[10](italic(P))),3,line=3.0)
## lines(xs.eu_cor_sample$rank,sa_cor_sample$fr_id_mean,col="blue",lty=2)
## lines(xs.eu_cor_sample$rank,sa_cor_sample$fr_id_97.5,col="blue",lty=3)
## lines(xs.eu_cor_sample$rank,sa_cor_sample$fr_id_2.5,col="blue",lty=3)
## lines(eu.xs_cor_sample$rank,eu_cor_sample$fr_id_mean,col="red",lty=2)
## lines(eu.xs_cor_sample$rank,eu_cor_sample$fr_id_97.5,col="red",lty=3)
## lines(eu.xs_cor_sample$rank,eu_cor_sample$fr_id_2.5,col="red",lty=3)
dev.copy2pdf(file="xs_sa_lOR_sign_quant_all.pdf")
quartz()
#x11()
par(mar=c(5, 4, 5, 4) + 0.1)
plot_correlation(quants,sac_xsy_q,sa_x_q,xsc_say_q,xs_x_q,sa_q_disp,xs_q_disp,quants_disp=quants_disp,pop.labs=c("South Africa","South Africa (ext.)"),ylabs="fraction of SNPs with consistent change (95% CI)",ylims=c(0.0,1.0))
dev.copy2pdf(file="xs_sa_lOR_sign_ci_quant_all.pdf")
quartz()

plot_correlation <- function(quants,a_b_q,a_x_q,b_a_q,b_x_q,a_q_disp,b_q_disp,quants_disp=quants_disp,pop.labs=c("Europe","South Africa (ext.)"),ylabs="fraction of SNPs with consistent change",leg.tit=c("Association"),ylims=c(0.5,1.0),leg.pos="topright"){
    par(mar=c(5, 4, 5, 4) + 0.1)
    plot(quants,a_b_q/a_x_q,type="l",col="red",lwd=2,xaxt='n',xlab="SNP rank", ylab=ylabs,log="x",ylim=ylims)
    lines(quants,b_a_q/b_x_q, col="blue",lwd=2)
    axis(1,at=quants_disp,labels=quants_disp, col="black", line=0.0)
    legend(leg.pos,pop.labs,lty=c(1,1,2,3),lwd=2,col=c("red","blue","black","black"),title=leg.tit)
    axis(3,col.axis="red",at=quants_disp,labels=round(a_q_disp,digits=1), col="black", line=0.0)
    axis(3,col.axis="blue",at=quants_disp,labels=round(b_q_disp,digits=1), tck=0, lty=0, line=1.0)
    mtext(expression(-log[10](italic(P))),3,line=3.0)
}
quartz()
plot_cor(eu.xs_cor,xs.eu_cor,"fid_ci",q.disp.eu,q.disp.xs,quants_disp=q.disp,pop.labs=c("Europe","South Africa (ext.)","Rnd. mean","Rnd. 95% CI"),ylims=c(0.0,1.0))
lines(xs.eu_cor_sample$rank,xs.eu_cor_sample$fr_ci_id_mean,col="blue",lty=2)
lines(xs.eu_cor_sample$rank,xs.eu_cor_sample$fr_ci_id_97.5,col="blue",lty=3)
lines(xs.eu_cor_sample$rank,xs.eu_cor_sample$fr_ci_id_2.5,col="blue",lty=3)
lines(eu.xs_cor_sample$rank,eu.xs_cor_sample$fr_ci_id_mean,col="red",lty=2)
lines(eu.xs_cor_sample$rank,eu.xs_cor_sample$fr_ci_id_97.5,col="red",lty=3)
lines(eu.xs_cor_sample$rank,eu.xs_cor_sample$fr_ci_id_2.5,col="red",lty=3)
dev.copy2pdf(file="xs_eu_lOR_ci_sign_quant_all.pdf")

plot_correlation(quants,eu_xsy_q,eu_x_q,xs_euy_q,xs_x_q,eu_q_disp,xs_q_disp,quants_disp=quants_disp,pop.labs=c("Europe","South Africa (ext.)","Rnd. mean","Rnd. 95% CI"))
lines(xs.eu_cor_sample$rank,sa_cor_sample$fr_id_mean,col="blue",lty=2)
lines(xs.eu_cor_sample$rank,sa_cor_sample$fr_id_97.5,col="blue",lty=3)
lines(xs.eu_cor_sample$rank,sa_cor_sample$fr_id_2.5,col="blue",lty=3)
lines(eu.xs_cor_sample$rank,eu_cor_sample$fr_id_mean,col="red",lty=2)
lines(eu.xs_cor_sample$rank,eu_cor_sample$fr_id_97.5,col="red",lty=3)
lines(eu.xs_cor_sample$rank,eu_cor_sample$fr_id_2.5,col="red",lty=3)
dev.copy2pdf(file="xs_eu_lOR_sign_quant_all.pdf")
# ci
quartz()
par(mar=c(5, 4, 5, 4) + 0.1)
plot(quants,xsc_euy_q/xs_x_q,type="l",col="blue",lwd=2,xaxt='n',xlab="SNP rank", ylab="fraction of SNPs with consistent change",log="x",ylim=c(0.0,1.0))
lines(quants,euc_xsy_q/eu_x_q, col="red",lwd=2)
axis(1,at=quants_disp,labels=quants_disp, col="black", line=0.0)
legend("topright",c("Europe","South Africa (ext.)"),lty=c(1,1,3),lwd=2,col=c("red","blue"),title=c("Association"))
lines(xs.eu_cor_sample$rank,sa_cor_sample$fr_ci_id_mean,col="blue",lty=2)
lines(xs.eu_cor_sample$rank,sa_cor_sample$fr_ci_id_97.5,col="blue",lty=3)
lines(xs.eu_cor_sample$rank,sa_cor_sample$fr_ci_id_2.5,col="blue",lty=3)
lines(eu.xs_cor_sample$rank,eu_cor_sample$fr_ci_id_mean,col="red",lty=2)
lines(eu.xs_cor_sample$rank,eu_cor_sample$fr_ci_id_97.5,col="red",lty=3)
lines(eu.xs_cor_sample$rank,eu_cor_sample$fr_ci_id_2.5,col="red",lty=3)
axis(3,col.axis="red",at=quants_disp,labels=round(eu_q_disp,digits=1), col="black", line=0.0)
axis(3,col.axis="blue",at=quants_disp,labels=round(xs_q_disp,digits=1), tck=0, lty=0, line=1.0)
mtext(expression(-log[10](italic(P))),3,line=3.0)
dev.copy2pdf(file="xs_eu_lOR_sign_ci_quant_all.pdf")

eu_xscor_q <- sapply(eu_q,function(x) {cor.test(vie_ita_sa_odds$O_xs[vie_ita_sa_odds$P_eu >= x & has_Pv],vie_ita_sa_odds$O_eu[vie_ita_sa_odds$P_eu >= x & has_Pv],method="spearman")$estimate})
sa_xscor_q <- sapply(sa_q,function(x) {cor.test(vie_ita_sa_odds$O_sa[vie_ita_sa_odds$P_sa >= x & has_Pv & has_OR],vie_ita_sa_odds$O_xs[vie_ita_sa_odds$P_sa >= x & has_Pv & has_OR],method="spearman")$estimate})
xs_eucor_q <- sapply(xs_q,function(x) {cor.test(vie_ita_sa_odds$O_xs[vie_ita_sa_odds$P_xs >= x & has_Pv],vie_ita_sa_odds$O_eu[vie_ita_sa_odds$P_xs >= x & has_Pv],method="spearman")$estimate})
xs_sacor_q <- sapply(xs_q,function(x) {cor.test(vie_ita_sa_odds$O_sa[vie_ita_sa_odds$P_xs >= x & has_Pv & has_OR],vie_ita_sa_odds$O_xs[vie_ita_sa_odds$P_xs >= x & has_Pv & has_OR],method="spearman")$estimate})

par(mar=c(5, 4, 5, 4) + 0.1)
plot(quants,xs_sacor_q,type="l",col="red",lwd=2,xaxt='n',xlab="SNP rank", ylab=expression(paste("Spearman\'s ",rho," of ",log(OR))),log="x",ylim=c(0,1))
lines(quants,sa_xscor_q, col="blue",lwd=2)
axis(1,at=quants_disp,labels=quants_disp, col="black", line=0.0)
legend("bottomleft",c("South Africa (ext.)","South Africa"),lty=c(1,1),lwd=2,col=c("red","blue"),title=c("Association"))
axis(3,col.axis="red",at=quants_disp,labels=round(eu_q_disp,digits=1), col="black", line=0.0)
axis(3,col.axis="blue",at=quants_disp,labels=round(xs_q_disp,digits=1), tck=0, lty=0, line=1.0)
mtext(expression(-log[10](italic(P))),3,line=3.0)
dev.copy2pdf(file="xs_sa_lORspearman_sign_quant_all.pdf")

par(mar=c(5, 4, 5, 4) + 0.1)
plot(quants,eu_xscor_q,type="l",col="red",lwd=2,xaxt='n',xlab="SNP rank", ylab=expression(paste("Spearman\'s ",rho," of ",log(OR))),log="x",ylim=c(0,1))
lines(quants,xs_eucor_q, col="blue",lwd=2)
axis(1,at=quants_disp,labels=quants_disp, col="black", line=0.0)
legend("bottomleft",c("Europe","South Africa (ext.)"),lty=c(1,1),lwd=2,col=c("red","blue"),title=c("Association"))
axis(3,col.axis="red",at=quants_disp,labels=round(eu_q_disp,digits=1), col="black", line=0.0)
axis(3,col.axis="blue",at=quants_disp,labels=round(xs_q_disp,digits=1), tck=0, lty=0, line=1.0)
mtext(expression(-log[10](italic(P))),3,line=3.0)
dev.copy2pdf(file="xs_eu_lORspearman_sign_quant_all.pdf")



quartz()
par(mar=c(5, 4, 5, 4) + 0.1)
plot(quants,eu_y_q/eu_x_q,type="l",col="red",lwd=2,xaxt='n',xlab="SNP rank", ylab="fraction of SNPs with consistent change",log="x",ylim=c(0.5,1))
lines(quants,sa_y_q/sa_x_q, col="blue",lwd=2)
axis(1,at=quants_disp,labels=quants_disp, col="black", line=0.0)
legend("topright",c("Europe","South Africa","Rnd. mean","Rnd. 95% CI"),lty=c(1,1,2,3),lwd=2,col=c("red","blue","black","black"),title=c("Association"))
axis(3,col.axis="red",at=quants_disp,labels=round(eu_q_disp,digits=1), col="black", line=0.0)
axis(3,col.axis="blue",at=quants_disp,labels=round(sa_q_disp,digits=1), tck=0, lty=0, line=1.0)
mtext(expression(-log[10](italic(P))),3,line=3.0)
lines(sa_cor_sample$rank,sa_cor_sample$fr_id_mean,col="blue",lty=2)
lines(sa_cor_sample$rank,sa_cor_sample$fr_id_97.5,col="blue",lty=3)
lines(sa_cor_sample$rank,sa_cor_sample$fr_id_2.5,col="blue",lty=3)
lines(eu_cor_sample$rank,eu_cor_sample$fr_id_mean,col="red",lty=2)
lines(eu_cor_sample$rank,eu_cor_sample$fr_id_97.5,col="red",lty=3)
lines(eu_cor_sample$rank,eu_cor_sample$fr_id_2.5,col="red",lty=3)
dev.copy2pdf(file="eu_sa_lOR_sign_quant_all.pdf")

quartz()
par(mar=c(5, 4, 5, 4) + 0.1)
plot_correlation(quants,euc_y_q,eu_x_q,sac_y_q,sa_x_q,eu_q_disp,sa_q_disp,quants_disp=quants_disp,pop.labs=c("Europe","South Africa","Rnd. mean","Rnd. 95% CI"),ylims=c(0.0,1.0),ylabs="fraction of SNPs with consistent change (95% CI)")
lines(sa_cor_sample$rank,sa_cor_sample$fr_ci_id_mean,col="blue",lty=2)
lines(sa_cor_sample$rank,sa_cor_sample$fr_ci_id_97.5,col="blue",lty=3)
lines(sa_cor_sample$rank,sa_cor_sample$fr_ci_id_2.5,col="blue",lty=3)
lines(eu_cor_sample$rank,eu_cor_sample$fr_ci_id_mean,col="red",lty=2)
lines(eu_cor_sample$rank,eu_cor_sample$fr_ci_id_97.5,col="red",lty=3)
lines(eu_cor_sample$rank,eu_cor_sample$fr_ci_id_2.5,col="red",lty=3)
dev.copy2pdf(file="eu_sa_lOR_ci_sign_quant_all.pdf")


eu_cor_q <- sapply(eu_q,function(x) {cor.test(vie_ita_sa_odds$O_sa[vie_ita_sa_odds$P_eu >= x & has_Pv],vie_ita_sa_odds$O_eu[vie_ita_sa_odds$P_eu >= x & has_Pv],method="spearman")$estimate})
sa_cor_q <- sapply(sa_q,function(x) {cor.test(vie_ita_sa_odds$O_sa[vie_ita_sa_odds$P_sa >= x & has_Pv & has_OR],vie_ita_sa_odds$O_eu[vie_ita_sa_odds$P_sa >= x & has_Pv & has_OR],method="spearman")$estimate})
quartz()
x11()
par(mar=c(5, 4, 5, 4) + 0.1)
plot(quants,eu_cor_q,type="l",col="red",lwd=2,xaxt='n',xlab="SNP rank", ylab=expression(paste("Spearman\'s ",rho," of ",log(OR))),log="x",ylim=c(0,1))
lines(quants,sa_cor_q, col="blue",lwd=2)
axis(1,at=quants_disp,labels=quants_disp, col="black", line=0.0)
legend("topright",c("Europe","South Africa","Rnd. mean","Rnd. 95% CI"),lty=c(1,1,2,3),lwd=2,col=c("red","blue","black","black"),title=c("Association"))
axis(3,col.axis="red",at=quants_disp,labels=round(eu_q_disp,digits=1), col="black", line=0.0)
axis(3,col.axis="blue",at=quants_disp,labels=round(sa_q_disp,digits=1), tck=0, lty=0, line=1.0)
mtext(expression(-log[10](italic(P))),3,line=3.0)
lines(sa_cor_sample$rank,sa_cor_sample$rho_id_mean,col="blue",lty=2)
lines(sa_cor_sample$rank,sa_cor_sample$rho_97.5,col="blue",lty=3)
lines(sa_cor_sample$rank,sa_cor_sample$rho_2.5,col="blue",lty=3)
lines(eu_cor_sample$rank,eu_cor_sample$rho_id_mean,col="red",lty=2)
lines(eu_cor_sample$rank,eu_cor_sample$rho_97.5,col="red",lty=3)
lines(eu_cor_sample$rank,eu_cor_sample$rho_2.5,col="red",lty=3)
dev.copy2pdf(file="eu_sa_lOR_spear_quant_all.pdf")

quants=c(10,12,15,17,25,30,75,100,250,500,750,1000,5000,10000,50000,100000,500000)
quants_disp=c(10,100,1000,10000,100000)
vie_q = sort( vie_ita_sa_odds$P_vi[has_Peu],decreasing=T)[quants]
ita_q = sort( vie_ita_sa_odds$P_ita[has_Peu],decreasing=T)[quants]
vie_q_disp = sort( vie_ita_sa_odds$P_vi[has_Peu],decreasing=T)[quants_disp]
ita_q_disp = sort( vie_ita_sa_odds$P_ita[has_Peu],decreasing=T)[quants_disp]
vie_x_q <- sapply(vie_q,function(x) {length( vie_ita_sa_odds$P_vi[has_Peu & vie_ita_sa_odds$P_vi >= x])})
ita_x_q <- sapply(ita_q,function(x) {length( vie_ita_sa_odds$P_ita[has_Peu & vie_ita_sa_odds$P_ita >= x])})
vie_y_q <- sapply(vie_q,function(x) {length( vie_ita_sa_odds$P_vi[has_Peu & vie_ita_sa_odds$P_vi >= x &  sign(vie_ita_sa_odds$lO_ita) == sign(vie_ita_sa_odds$lO_vi) ])})
ita_y_q <- sapply(ita_q,function(x) {length( vie_ita_sa_odds$P_ita[has_Peu & vie_ita_sa_odds$P_ita >= x & sign(vie_ita_sa_odds$lO_ita) == sign(vie_ita_sa_odds$lO_vi) ])})




quartz()
par(mar=c(5, 4, 5, 4) + 0.1)
ylimit=c(min(c(vie_y_q/vie_x_q,ita_y_q/ita_x_q)),max(c(vie_y_q/vie_x_q,ita_y_q/ita_x_q)))
plot(quants,vie_y_q/vie_x_q,type="l",col="green",lwd=2,xaxt='n',ylim=ylimit,xlab="SNP rank", ylab="fraction of SNPs with consistent change",log="x")
lines(quants,ita_y_q/ita_x_q, col="orange",lwd=2)
axis(1,at=quants_disp,labels=quants_disp, col="black", line=0.0)
legend("topright",c("Vienna, Austria","Bolzano, Italy"),lty=1,lwd=2,col=c("green","orange"),title=c("Association"))
axis(3,col.axis="green",at=quants_disp,labels=round(vie_q_disp,digits=1), col="black", line=0.0)
axis(3,col.axis="orange",at=quants_disp,labels=round(ita_q_disp,digits=1), tck=0, lty=0, line=1.0)
mtext(expression(-log[10](italic(P))),3,line=3.0)
dev.copy2pdf(file="vie_ita_lOR_sign_quant_all.pdf")

vie_cor_q <- sapply(vie_q,function(x) {cor.test(vie_ita_sa_odds$O_ita[vie_ita_sa_odds$P_vi >= x & has_Peu],vie_ita_sa_odds$O_vi[vie_ita_sa_odds$P_vi >= x & has_Peu],method="spearman")$estimate})
ita_cor_q <- sapply(ita_q,function(x) {cor.test(vie_ita_sa_odds$O_ita[vie_ita_sa_odds$P_ita >= x & has_Peu & has_OR],vie_ita_sa_odds$O_vi[vie_ita_sa_odds$P_ita >= x & has_Peu & has_OR],method="spearman")$estimate})
quartz()
par(mar=c(5, 4, 5, 4) + 0.1)
ylimit=c(min(c(vie_cor_q,ita_cor_q)),max(c(vie_cor_q,ita_cor_q)))
plot(quants,vie_cor_q,type="l",col="green",lwd=2,xaxt='n',xlab="SNP rank",ylim=ylimit, ylab=expression(paste("Spearman\'s ",rho," of ",log(OR))),log="x")
lines(quants,ita_cor_q, col="orange",lwd=2)
axis(1,at=quants_disp,labels=quants_disp, col="black", line=0.0)
legend("topright",c("Vienna, Austria","Bolzano, Italy"),lty=1,lwd=2,col=c("green","orange"),title=c("Association"))
axis(3,col.axis="green",at=quants_disp,labels=round(vie_q_disp,digits=1), col="black", line=0.0)
axis(3,col.axis="orange",at=quants_disp,labels=round(ita_q_disp,digits=1), tck=0, lty=0, line=1.0)
mtext(expression(-log[10](italic(P))),3,line=3.0)
dev.copy2pdf(file="vie_ita_lOR_spear_quant_all.pdf")


q.loc=c(5,10,15,25,30,50,75,100,250,500) #,250,500,750,1000,5000)
q.loc.disp=c(5,10,50,100,500) # ,1000,5000)
tan=c("X",8867539,9122860)
tan=c("X",9117290,9122860) # tan upstream until GR8 end
tan=c("X",9119470,9122860) # just CG.., MSE and GR8
ebony=c("3R",17055975,17069171)
bab1=c("3L",1036369,1101089) # whole of bab1
bab1=c("3L",1040354,1099101) # intron1 of bab1
tan_cor=with(vie_ita_sa_odds, vie_ita_sa_odds[which(CHR == tan[1] & (tan[2] <= BPS & BPS <= tan[3]) & has_Pv),] )
bab1_cor=with(vie_ita_sa_odds, vie_ita_sa_odds[which(CHR == bab1[1] & (bab1[2] <= BPS & BPS <= bab1[3]) & has_Pv),] )
ebony_cor=with(vie_ita_sa_odds, vie_ita_sa_odds[which(CHR == ebony[1] & (ebony[2] <= BPS & BPS <= ebony[3]) & has_Pv),] )

eu_qt = sort( tan_cor$P_eu,decreasing=T)[q.loc]
sa_qt = sort( tan_cor$P_sa,decreasing=T)[q.loc]
eu_qt_disp = sort( tan_cor$P_eu,decreasing=T)[q.loc.disp]
sa_qt_disp = sort( tan_cor$P_sa,decreasing=T)[q.loc.disp]
eu_x_qt <- sapply(eu_qt,function(x) {length( tan_cor$P_eu[ tan_cor$P_eu >= x])})
sa_x_qt <- sapply(sa_qt,function(x) {length( tan_cor$P_sa[ tan_cor$P_sa >= x])})
eu_y_qt <- sapply(eu_qt,function(x) {length( tan_cor$P_eu[tan_cor$P_eu >= x &  sign(tan_cor$lO_sa) == sign(tan_cor$lO_eu) ])})
sa_y_qt <- sapply(sa_qt,function(x) {length( tan_cor$P_sa[tan_cor$P_sa >= x & sign(tan_cor$lO_sa) == sign(tan_cor$lO_eu) ])})
euc_y_qt <- sapply(eu_qt,function(x) {length( tan_cor$P_eu[ tan_cor$P_eu >= x &  ( tan_cor$sO_sa * tan_cor$sO_eu > 0 ) ])})
sac_y_qt <- sapply(sa_qt,function(x) {length( tan_cor$P_sa[ tan_cor$P_sa >= x & (tan_cor$sO_sa * tan_cor$sO_eu > 0 ) ])})

eu_cor_qt <- sapply(eu_qt,function(x) {cor.test(tan_cor$O_xs[tan_cor$P_eu >= x ],tan_cor$O_eu[tan_cor$P_eu >= x ],method="spearman")$estimate})
sa_cor_qt <- sapply(sa_qt,function(x) {cor.test(tan_cor$O_sa[tan_cor$P_sa >= x ],tan_cor$O_xs[tan_cor$P_sa >= x ],method="spearman")$estimate})


plot_correlation(q.loc,euc_y_qt,eu_x_qt,sac_y_qt,sa_x_qt,eu_qt_disp,sa_qt_disp,quants_disp=q.loc.disp,pop.labs=c("Europe","South Africa"),ylabs="fraction of SNPs with consistent change (95% CI)",ylims=c(0.0,1.0))
dev.copy2pdf(file="eu_sa_lOR_sign_CI_tan.pdf")
plot_correlation(q.loc,eu_y_qt,eu_x_qt,sa_y_qt,sa_x_qt,eu_qt_disp,sa_qt_disp,quants_disp=q.loc.disp,pop.labs=c("Europe","South Africa"),ylabs="fraction of SNPs with consistent change",ylims=c(0.5,1.0))
dev.copy2pdf(file="eu_sa_lOR_sign_tan.pdf")






quants=10^-(c(seq(0,2.5,0.5),seq(3,4.0,0.25),seq(4.15,5.5,0.125)))

eu_q = quantile( eu_sa_odds$P_eu[has_Pv],probs=quants)
sa_q = quantile( eu_sa_odds$P_sa[has_Pv],probs=quants)
eu_x_q <- sapply(eu_q,function(x) {length( eu_sa_odds$P_eu[has_Pv & eu_sa_odds$P_eu <= x])})
sa_x_q <- sapply(sa_q,function(x) {length( eu_sa_odds$P_sa[has_Pv & eu_sa_odds$P_sa <= x])})
eu_y_q <- sapply(eu_q,function(x) {length( eu_sa_odds$P_eu[has_Pv & eu_sa_odds$P_eu <= x &  sign(eu_sa_odds$O_sa) == sign(eu_sa_odds$O_eu) ])})
sa_y_q <- sapply(sa_q,function(x) {length( eu_sa_odds$P_sa[has_Pv & eu_sa_odds$P_sa <= x & sign(eu_sa_odds$O_sa) == sign(eu_sa_odds$O_eu) ])})

quartz()
par(mar=c(5, 4, 5, 4) + 0.1)
plot(quants,eu_y_q/eu_x_q,type="l",col="red",lwd=2,xlab="quantile of P values", ylab="fraction of SNPs with consistent log OR sign",log="x")
lines(quants,sa_y_q/sa_x_q, col="blue",lwd=2)
legend("topright",c("Europe","South Africa"),lty=1,lwd=2,col=c("red","blue"),title=c("association"))
axis(3,col.axis="red",at=quants,labels=round(-log10(eu_q),digits=1), col="black", line=0.0)
axis(3,col.axis="blue",at=quants,labels=round(-log10(sa_q),digits=1), tck=0, lty=0, line=1.0)
mtext(expression(-log[10](italic(P))),3,line=3.0)
dev.copy2pdf(file="eu_sa_lOR_sign_quant_pv_lt_1e3.pdf")


# correlation (OR)

cor.test(eu_sa_odds$O_sa_ni[eu_sa_odds$P_eu <= 1e-3 & has_Pv & has_OR],eu_sa_odds$O_eu_ni[eu_sa_odds$P_eu <= 1e-3 & has_Pv & has_OR],method="spearman")$estimate
eu_q = quantile( eu_sa_odds$P_eu[has_Pv],probs=quants)
sa_q = quantile( eu_sa_odds$P_sa[has_Pv],probs=quants)
eu_cor_q <- sapply(eu_q,function(x) {cor.test(eu_sa_odds$O_sa_ni[eu_sa_odds$P_eu <= x & has_Pv & has_OR],eu_sa_odds$O_eu_ni[eu_sa_odds$P_eu <= x & has_Pv & has_OR],method="spearman")$estimate})
sa_cor_q <- sapply(sa_q,function(x) {cor.test(eu_sa_odds$O_sa_ni[eu_sa_odds$P_sa <= x & has_Pv & has_OR],eu_sa_odds$O_eu_ni[eu_sa_odds$P_sa <= x & has_Pv & has_OR],method="spearman")$estimate})

quartz()
par(mar=c(5, 4, 5, 4) + 0.1)
plot(quants,eu_cor_q,type="l",col="red",lwd=2,xlab="quantile of P values", ylab=expression(paste(rho," Spearman ranked correlation of ",log(OR))),log="x")
lines(quants,sa_cor_q, col="blue",lwd=2)
legend("topright",c("Europe","South Africa"),lty=1,lwd=2,col=c("red","blue"),title=c("association"))
axis(3,col.axis="red",at=quants,labels=round(-log10(eu_q),digits=1), col="black", line=0.0)
axis(3,col.axis="blue",at=quants,labels=round(-log10(sa_q),digits=1), tck=0, lty=0, line=1.0)
mtext(expression(-log[10](italic(P))),3,line=3.0)
dev.copy2pdf(file="eu_sa_lOR_cor_OR_quant_pv_lt_1e3.pdf")

eu_sa_odds$lPeu=-log10(eu_sa_odds$P_eu)
eu_sa_odds$lPsa=-log10(eu_sa_odds$P_sa)


plot_OR_all_pV <- function(region,coords,lpVals,ylab_pv=expression(-log[10](italic(P))) ){
    region_idx = which(eu_sa_odds$CHR == coords[1] & (coords[2] <= eu_sa_odds$BPS & eu_sa_odds$BPS <= coords[3] ) & (eu_sa_odds$lPeu >= lpVals[1] | eu_sa_odds$lPsa >= lpVals[2]))
    length(region_idx)
    layout(matrix(1:2, ncol = 1), widths = 1, heights = c(0.75,1.25), respect = FALSE)
    par(mar = c(0, 4.1, 4.1, 2.1))
    main_tit=bquote("SNPs around" ~.(region) ~"with" ~ Peu <= 10^-.(lpVals[1]) ~ " or " ~  Psa <= 10 ^-.(lpVals[2]))
    xlimit=c(1-0.25,length(region_idx)+0.15)  
    plot(1:length(region_idx),eu_sa_odds$lPsa[region_idx],col="blue",cex=1,pch=16,ylim=c(0,max(eu_sa_odds[region_idx,c("lPsa","lPeu")])),ylab=ylab_pv, xaxt="n",main=main_tit,xlim=xlimit )
    abline(v=(1:length(region_idx))-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+0.25,col="grey",lty="dotted",lw=2)
    points(1:length(region_idx),eu_sa_odds$lPeu[region_idx],col="red",cex=1,pch=16)
    legend("topright",legend=c("SA","Europe"),pch=c(16,20),col=c("blue","red"))
    par(mar = c(4.1, 4.1, 0, 2.1))
    ors=c("O_el","O_ed","O_svl","O_svd","O_sd","O_sl")
    vals=unlist(eu_sa_odds[region_idx,ors])
    ymin=min(vals[is.finite(vals)],na.rm=T)
    ymax=max(vals[is.finite(vals)],na.rm=T)
    plot((1:length(region_idx))-0.15,eu_sa_odds[region_idx,"O_svl"],col="blue",cex=1,pch=6,ylim=c(ymin,ymax),ylab="log(OR)",xlab="SNPs",xlim=xlimit )
    abline(v=(1:length(region_idx))-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+0.25,col="grey",lty="dotted",lw=2)
    abline(h=seq(-4,4,by=1),col="grey",lty="dotted",lw=2)
    #points((1:length(region_idx))-0.15,eu_sa_odds[region_idx,O_svl],col="blue",cex=1,pch=6)
    points((1:length(region_idx))-0.15,eu_sa_odds[region_idx,"O_svd"],col="blue",cex=1,pch=17)
     points((1:length(region_idx))-0.15,eu_sa_odds[region_idx,"O_sl"],col="blue",cex=0.75,pch=23)
    points((1:length(region_idx))-0.15,eu_sa_odds[region_idx,"O_sd"],col="blue",cex=0.75,pch=23,bg="cyan")

    points((1:length(region_idx))+0.15,eu_sa_odds[region_idx,"O_el"],col="red",cex=1,pch=6)
    points((1:length(region_idx))+0.15,eu_sa_odds[region_idx,"O_ed"],col="red",cex=1,pch=17)
    legend("bottomright",legend=c("very light","light","dark","very dark"),title="log(OR) to base",pch=c(6,23,23,17),col=c("black"),pt.bg=c("white","white","grey16","black"),pt.cex =c(1,0.75,0.75,1))
}

tan=c("X",8867539,9170829)
tan_p=c(20,12)
quartz(width=18,height=8)
plot_OR_all_pV("tan",tan,tan_p)
dev.copy2pdf(file="eu_sa_OR_comp_tan_20_12.pdf")


tan=c("X",8867539,9170829)
tan_p=c(20,12)
x11(width=18,height=8)
plot_OR_xs_pV("tan",tan,tan_p)
dev.copy2pdf(file="eu_sa_xs_OR_comp_tan_20_12.pdf")

tan=c("X",8867539,9170829)
tan_p=c(12,12)
quartz(width=18,height=8)
plot_OR_all_pV("tan",tan,tan_p)
dev.copy2pdf(file="eu_sa_OR_comp_tan.pdf")

bab1=c("3L",1071314,1120551 )
bab1_p=c(10,16)
quartz(width=18,height=8)
plot_OR_all_pV("bab1",bab1,bab1_p)
dev.copy2pdf(file="eu_all_OR_comp_bab1_lpeu10_psa16.pdf")
bab1_p=c(10,18)
quartz(width=18,height=8)
plot_OR_all_pV("bab1",bab1,bab1_p)
dev.copy2pdf(file="eu_all_OR_comp_bab1_lpeu10_psa18.pdf")
bab1_p=c(12,22)
quartz(width=18,height=8)
plot_OR_all_pV("bab1",bab1,bab1_p)
dev.copy2pdf(file="eu_all_OR_comp_bab1_lpeu12_psa22.pdf")

ebony=c("3R",17055975,17069171)
ebony_p=c(6,6)
quartz(width=18,height=8)
plot_OR_all_pV("ebony",ebony,ebony_p)
dev.copy2pdf(file="eu_all_OR_comp_ebony_lpeu6_lpsa6.pdf")

bab2=c("3L",1128600,1181056)
bab2_p=c(8,8)
quartz(width=18,height=8)
plot_OR_all_pV("bab2",bab2,bab2_p)
dev.copy2pdf(file="eu_all_OR_comp_bab2_lpeu8_lpsa8.pdf")
bab2_p=c(6,6)
quartz(width=18,height=8)
plot_OR_all_pV("bab2",bab2,bab2_p)
dev.copy2pdf(file="eu_all_OR_comp_bab2_lpeu6_lpsa6.pdf")


pdm3=c("2R",4179149,4319221)
pdm3_p=c(6,6)
quartz(width=18,height=8)
plot_OR_all_pV("pdm3",pdm3,pdm3_p)
dev.copy2pdf(file="eu_all_OR_comp_pdm3_lpeu6_lpsa6.pdf")

pdm3=c("2R",4179149,4319221)
pdm3_p=c(6,8)
quartz(width=18,height=8)
plot_OR_all_pV("pdm3",pdm3,pdm3_p)
dev.copy2pdf(file="eu_all_OR_comp_pdm3_lpeu6_lpsa8.pdf")

# manhattan plots

eu_plot=eu_sa_odds[eu_sa_odds$P_eu <= 0.1 & ! is.na(eu_sa_odds$P_eu),c("CHR","BPS","P_eu")]
sa_plot=eu_sa_odds[eu_sa_odds$P_sa <= 0.1 & ! is.na(eu_sa_odds$P_sa),c("CHR","BPS","P_sa")]
#sa_plot=eu_sa_freqs[eu_sa_freqs$Psa <= 0.01 & ! is.na(eu_sa_freqs$Psa),c("CHR","BPS","Psa")]
colnames(eu_plot)=c("CHR","BP","P")
colnames(sa_plot)=c("CHR","BP","P")
summary(sa_plot)
chroms=c("X","2L","2R","3L","3R","4")
off_eu=get_offset_manh(eu_plot,limitchromosomes=chroms)
off_sa=get_offset_manh(sa_plot,limitchromosomes=chroms)

quartz(height=6,width=12)

x11(height=6,width=12)
par(bg="white")
manhattan(eu_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="Vienna (2011) and Bolzano (2012)", suggestiveline=c(22))
rect(8867539+off_eu["X"],0,9170829+off_eu["X"],60,lty=3, border="red")
rect(1071314+off_eu["3L"],0,1120551+off_eu["3L"],60,lty=3, border="green")
rect(17055975+off_eu["3R"],0,17069171+off_eu["3R"],60,lty=3, border="blue")
rect(4179149+off_eu["2R"],0,4319221+off_eu["2R"],60,lty=3, border="yellow")
dev.print(png,width=800,file="eu_manhattan.png")
quartz(height=6,width=12)
par(bg="white")
chroms=c("X","2L","2R","3L","3R","4")
manhattan(sa_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="South Africa (2012)", suggestiveline=c(-log10(2.054127e-19)))
rect(8867539+off_eu["X"],0,9170829+off_eu["X"],60,lty=3, border="red")
rect(1071314+off_eu["3L"],0,1120551+off_eu["3L"],60,lty=3, border="green")
rect(17055975+off_eu["3R"],0,17069171+off_eu["3R"],60,lty=3, border="blue")
rect(4179149+off_eu["2R"],0,4319221+off_eu["2R"],60,lty=3, border="yellow")
dev.print(png,width=800,file="sa_manhattan.png")

# closeups
#tan

tan
X:9111688-9117290
CG15370:
X:9119470-9120721
Gr8a
X:9121170-9122860
CG12121
X:9122829-9126411

plot_region <- function(chrom,a,b, title,ylab=expression(-log[10](italic(P))),plus=10,cex=0.75,exact=F){
    region_sa = which(sa_plot$CHR == chrom & sa_plot$BP >= a  & sa_plot$BP <= b )
    region_eu = which(eu_plot$CHR == chrom & eu_plot$BP >= a  & eu_plot$BP <= b )
    ymax = max(-log10(sa_plot$P[region_sa]),-log10(eu_plot$P[region_eu])) + plus
    if(exact != F) {ymax=exact}
    tit=bquote(italic(.(title)))
    plot(sa_plot$BP[region_sa],-log10(sa_plot$P[region_sa]), col="blue", pch=20, main=tit, ylim=c(0,ymax),xlab=chrom,ylab=ylab,cex=cex )
    points(eu_plot$BP[region_eu],-log10(eu_plot$P[region_eu]), col="red", pch=20,cex=cex)
}

plot_region2 = function(chrom,a,b, title,ylab=expression(-log[10](italic(P)))){
  region_sa = which(sa_plot$CHR == chrom & sa_plot$BP >= a  & sa_plot$BP <= b )
  region_eu = which(eu_plot$CHR == chrom & eu_plot$BP >= a  & eu_plot$BP <= b )
  par(mfrow=c(2,1))
  xlims=c(min(c(eu_plot$BP[region_eu],sa_plot$BP[region_sa])),max(c(eu_plot$BP[region_eu],sa_plot$BP[region_sa])))
  tit=bquote(italic(.(title)))
  ymax = max(-log10(eu_plot$P[region_eu]))+10
  par(mar=c(1.0,3.5,3.1,1))
  plot(eu_plot$BP[region_eu],-log10(eu_plot$P[region_eu]), col="red", pch=20, main=tit, xlim=xlims,ylim=c(0,ymax),xlab="",ylab="",xaxt="n",cex=0.75 )
  legend("topleft",legend=c("Europe","South Africa"),pch=c(20,20),col=c("red","blue"))
  par(mar=c(4.1,3.5,0,1))  
  ymax = max(-log10(sa_plot$P[region_sa]))
  plot(sa_plot$BP[region_sa],-log10(sa_plot$P[region_sa]), col="blue", pch=20, main="", xlim=xlims,ylim=c(0,ymax),xlab=chrom,ylab="",cex=0.75 )
  mtext
}

plot_region_indel = function(chrom,a,b, title,indels_sa=indels_sa7,ylab=expression(-log[10](italic(P)))){
    region_sa = which(indels_sa$CHR == chrom & indels_sa$BP >= a  & indels_sa$BP <= b )
    #region_eu = which(eu_plot$CHR == chrom & eu_plot$BP >= a  & eu_plot$BP <= b )
    ymax = max(-log10(indels_sa$P[region_sa]))#,-log10(eu_plot$P[region_eu]))
    plot(indels_sa$BP[region_sa],-log10(indels_sa$P[region_sa]), col="blue", pch=20, main=title, ylim=c(0,ymax),xlab=chrom,xlim=c(a,b),ylab=ylab,cex=0.75 )
    
    #points(eu_plot$BP[region_eu],-log10(eu_plot$P[region_eu]), col="red", pch=20,cex=0.75)
}

require(grDevices)
quartz()
#x11()
fdr.eu=-log10(3.07403844435e-23)
fdr.sa=-log10(1.916206e-16)
par(mfrow=c(1,1))
par(mar=c(4.1,4.1,4.1,1))
plot_region("X",9116995,9121764,"tan",plus=10,cex=1.5,exact=85)
rect(9111688,60+5,9117290,62+5,col="grey")
text(9117100,63.5+5,labels="tan")
text(9117290-50,61+5,labels="> > >",srt=180,adj = c(0.0, NA))
rect(9121170,60+5,9122860,62+5,col="grey")
text(9121170+50,63.5+5,labels="Gr8a",adj = c(0.0, NA))
text(9121170+50,61+5,labels="> > >",adj = c(0.0, NA))
rect(9119470,60+5,9120721,62+5,col="grey")
text((9119470+9120721)/2,63.5+5,labels="CG15370")
text(9119470+50,61+5,labels="> > >",adj = c(0.0, NA))
rect(9120721,60+5,9121170,62+5,col="white")
text((9120721+9121170)/2,63.5+5,labels="MSE")
abline(h=fdr.eu,col="red",lty=2)
abline(h=fdr.sa,col="blue",lty=2)
#legend("topleft",legend=c("Europe","South Africa"),pch=c(20,20),col=c("red","blue"))
#dev2bitmap("eu_sa_tan_region.tiff_nl",type="tiff24nc",res=300)
dev.copy2pdf(file="eu_sa_tan_region_85bg_nl.pdf")
legend("topleft",legend=c("Europe","South Africa"),pch=c(20,20),pt.cex=1.5,col=c("red","blue"),bty = "n")
dev.copy2pdf(file="eu_sa_tan_region_85bg_lg.pdf")

quartz()
x11()

plot_region("3L",1039508,1111138,"bab1",plus=4,cex=1)
segments(1036369,36+9,1101089,36+9)
rect(c(1099101,1039796,1036369,1036369),60+9,c(1101089,1040223,1039215,1039215),62+9,col="grey")
#rect(1079062,35.5,1079672,36.5,col="grey")
rect(1084814,35.5+9,1085477,36.5+9,col="white")
text((1096651+1098007)/2,34.5+9,labels="AE")
text((1084814+1085477)/2,34.5+9,labels="DME")
rect(1096651,35.5+9,1098007,36.5+9,col="white")
text(1060000,38+9,labels="bab1")
text(1060000,36+9,labels="< < <")
#dev2bitmap("eu_sa_bab1_big_region.tiff",type="tiff24nc",res=300)
dev.copy2pdf(file="eu_sa_bab1_big_region.pdf")
dev.copy2pdf(file="eu_sa_bab1_big_region_bg.pdf")
quartz()
x11()
plot_region("3L",1036369,1177276,"bab1/2",cex=1.5,plus=15,exact=85)
segments(1036369,61+5,1101089,61+5)
segments(1140429,61+5,1177276,61+5)
rect(c(1099101,1039796,1036369,1036369),60+5,c(1101089,1040223,1039215,1039215),62+5,col="grey")
rect(c(1177276,1174607,1143929,1143780,1142921),60+5, c(1176981,1172758,1143851,1143350,1140429),62+5,col="grey")
#rect(1079062,35.5,1079672,36.5,col="grey") removed from flybase and such
rect(1084814,60+5,1085477,62+5,col="white") # DME
#rect(1169506,35.5,1173256,36.5,col="white") # Bab2 reg region
text((1084814+1085477)/2,63.5+5,labels="DME")
rect(1096651,60+5,1098007,62+5,col="white") # DAE
text((1096651+1098007)/2,63.5+5,labels="AE")
text(1060000,63.5+5,labels="bab1")
text((1140429+1177276)/2,63.5+5,labels="bab2")
text((1140429+1177276)/2,61+5,labels="< < <")
text(1060000,61+5,labels="< < <")
abline(h=fdr.eu,col="red",lty=2)
abline(h=fdr.sa,col="blue",lty=2)

#dev2bitmap("eu_sa_bab12_region.tiff",type="tiff24nc",res=300)
dev.copy2pdf(file="eu_sa_bab12_region_85.pdf")
legend("topleft",legend=c("Europe","South Africa"),pch=c(20,20),pt.cex=1.5,col=c("red","blue"),bty = "n")
dev.copy2pdf(file="eu_sa_bab12_region_85lg.pdf")
quartz()
x11()
plot_region("3L",1070000,1111138,"bab1",plus=5,cex=1.5,exact=80)
segments(1036369,36+9,1101089,36+9)
rect(c(1099101,1039796,1036369,1036369),35+9,c(1101089,1040223,1039215,1039215),37+9,col="grey")
#rect(1079062,35.5,1079672,36.5,col="grey")
rect(1084814,35.5+9,1085477,36.5+9,col="white")
text((1084814+1085477)/2,34.5+9,labels="DME")
rect(1096651,35.5+9,1098007,36.5+9,col="white") # DAE
text((1096651+1098007)/2,34.5+9,labels="AE")
text(1075000,38+9,labels="bab1")
text(1075000,36+9,labels="< < <")
#legend("topright",legend=c("Europe","South Africa"),pch=c(20,20),col=c("red","blue"))
dev2bitmap("eu_sa_bab1_small_region_nl.tiff",type="tiff24nc",res=300)
dev.copy2pdf(file="eu_sa_bab1_small_region_bg_nl.pdf")
legend("topright",legend=c("Europe","South Africa"),pch=c(20,20),col=c("red","blue"))
dev.copy2pdf(file="eu_sa_bab1_small_region_bg.pdf")

#pdm3
quartz()
x11()
plot_region("2R",4210633,4235314,"pdm3",plus=1.0,cex=1)
segments(4214995,0.15,42837859,0.15)
rect(c(4214995,4235006),0,c(4215512,4235799),0.3,col="grey")
text(4225000,0.4,labels="pdm3")
text(4225000,0.14,labels="> > >")
legend("topright",legend=c("Europe","South Africa"),pch=c(20,20),col=c("red","blue"))
dev2bitmap("eu_sa_pdm3_region.tiff",type="tiff24nc",res=300)
dev.copy2pdf(file="eu_sa_pdm3_region_bg.pdf")

# ebony
quartz()
x11()
plot_region("3R",17062142,17069000,"ebony",plus=2.0,cex=1)
yc=15
rect(17055561,yc+0,17062899,yc+0.5,col="grey")
text(17062700,yc+0.7,labels="e")
text(17062650,yc+0.25,labels="< < <")
text((17066518+17067319)/2,yc+0.7,labels="aCRE")
rect(17066518,yc+0,17067319,yc+0.5,col="white")
rect(17068351,yc+0,17071140,yc+0.5,col="grey")
text(17068351,yc+0.7,labels="CG5892",adj = c(0.0, NA))
text(17068351+50,yc+0.25,labels="> > >",adj = c(0.0, NA))
legend("right",legend=c("Europe","South Africa"),pch=c(20,20),col=c("red","blue"))
dev2bitmap("eu_sa_ebony_region_nl.tiff",type="tiff24nc",res=300)
dev.copy2pdf(file="eu_sa_ebony_region_bg.pdf")


#dimorphic element: approx 1084336-1085897
bab1:
3L:1036369-1101089
region: 3L:1039508-1111138   
Exon number: 1
3L:1099101-1101089
<b>Dmel:r5:3L:1040354:1099101:-</b>
3L:1040354-1099101
Type = exon_junction
id = FBsf0000145493
Exon number: 1
3L:1099101-1101089
bab1:7
3L:1099101-1101089
bab1:6
3L:1079062-1079672
Exon number: 3
3L:1039796-1040223
bab1:4
3L:1039796-1040223
Exon number: 4
3L:1036369-1039215
bab1:1
3L:1036369-1039215

bab2
3L:1140429-1177276

pdm3
2R:4214995-4283785
Exon number: 1
2R:4214995-4215512
pdm3:1
Exon number: 2
2R:4235006-4235799
pdm3:2


e
3R:17055561-17062899
CG5892
3R:17068340-17071140

indels.eu=read.table("/Volumes/Temp/Lukas/Data//A7_pigmentation/Polymorphisms/Indels/vie_ita_indels_all_filter_95_15_maxAF_gt_0.05_minDP_gt_10_allele_depths_L_D.cmh.gz",na.strings = c("NA","nan"),header=F)
indels.sa=read.table("/Volumes/Temp/Lukas/Data/SA_A7/Polymorphisms/Indels/sa_a7_indels_all_filter_95_15_maxAF_gt_0.05_minDP_gt_15_allele_depths_VL_VVD.cmh.gz",na.strings = c("NA","nan"),header=F)
colnames(indels.eu)=c("CHR","BPS","A1","A2","P","O","Oa","Ob")
colnames(indels.sa)=c("CHR","BPS","A1","A2","P","O","Oa","Ob")
indels.eu$lO=log2(indels.eu$O)
indels.sa$lO=log2(indels.sa$O)
median(abs(indels.sa[order(indels.sa$P)[1:10],]$lO))
median(abs(indels.eu[order(indels.eu$P)[1:10],]$lO))


indels_eu7=read.table("/Volumes/vetgrid10/Data//A7_pigmentation/Polymorphisms/Indels/vie_ita_indels_all_filter_95_15_maxAF_gt_0.05_minDP_gt_10_allele_depths_L_D.gwas",header=T)

indels_sa7=read.table("/Volumes/vetgrid10/Data/SA_A7/Polymorphisms/Indels/sa_a7_indels_all_filter_95_15_maxAF_gt_0.05_minDP_gt_15_allele_depths_VL_VVD.gwas",header=T)
off_sa_indels=get_offset_manh(indels_sa7,limitchromosomes=chroms)
off_eu_indels=get_offset_manh(indels_eu7,limitchromosomes=chroms)


quartz(height=6,width=12)
par(bg="white")
chroms=c("X","2L","2R","3L","3R","4")
manhattan(indels_sa7,colors=c("black","slategrey"),limitchromosomes=chroms,main="South Africa, InDels",suggestiveline=c(-log10(7.89290778441e-16)))
rect(9100000+off_sa_indels["X"],0,9170829+off_sa_indels["X"],60,lty=3, border="red")
rect(1071314+off_sa_indels["3L"],0,1120551+off_sa_indels["3L"],60,lty=3, border="green")
rect(17055975+off_sa_indels["3R"],0,17069171+off_sa_indels["3R"],60,lty=3, border="blue")
rect(4179149+off_sa_indels["2R"],0,4319221+off_sa_indels["2R"],60,lty=3, border="yellow")

#legend("topright",legend=c("tan","bab","ebony","pdm3"),lty=1,col=c("red","green","blue","yellow"))
dev.print(png,width=800,file="sa_indel_manhattan.png")

quartz(height=6,width=12)
par(bg="white")
chroms=c("X","2L","2R","3L","3R","4")
manhattan(indels_eu7,colors=c("black","slategrey"),limitchromosomes=chroms,main="Europe, InDels", suggestiveline=c(-log10(5.92726940102e-12)))
rect(9100000+off_eu_indels["X"],0,9170829+off_eu_indels["X"],60,lty=3, border="red")
rect(1071314+off_eu_indels["3L"],0,1120551+off_eu_indels["3L"],60,lty=3, border="green")
rect(17055975+off_eu_indels["3R"],0,17069171+off_eu_indels["3R"],60,lty=3, border="blue")
rect(4179149+off_eu_indels["2R"],0,4319221+off_eu_indels["2R"],60,lty=3, border="yellow")

#legend("topright",legend=c("tan","bab","ebony","pdm3"),lty=1,col=c("red","green","blue","yellow"))
dev.print(png,width=800,file="eu_indel_manhattan.png")

quartz()
plot_region_indel("3L",1070000,1111138,"bab1, SA (2012) InDels")
segments(1036369,18,1101089,18)
rect(c(1099101,1039796,1036369,1036369),17.5,c(1101089,1040223,1039215,1039215),18.5,col="grey")
rect(1079062,17.75,1079672,18.25,col="grey")
rect(1084814,17.5,1085477,18.5,col="white")
text((1084814+1085477)/2,19,labels="DME")
text(1075000,19,labels="bab1")
text(1075000,18,labels="< < <")
dev.copy2pdf(file="sa_indels_bab1_small_region.pdf")
quartz()
plot_region_indel("3L",1039508,1111138,"bab1, SA (2012) InDels")
segments(1036369,18,1101089,18)
rect(c(1099101,1039796,1036369,1036369),17.5,c(1101089,1040223,1039215,1039215),18.5,col="grey")
rect(1079062,17.75,1079672,18.25,col="grey")
rect(1084814,17.5,1085477,18.5,col="white")
text((1084814+1085477)/2,19,labels="DME")
text(1075000,19,labels="bab1")
text(1075000,18,labels="< < <")
dev.copy2pdf(file="sa_indels_bab1_big_region.pdf")

quartz()
plot_region_indel("X",9116995,9121764,"tan, SA (2012) InDels")
rect(9111688,4.5,9117290,4.7,col="grey")
text(9117100,4.8,labels="tan")
text(9117050,4.6,labels="< < <")
rect(9121170,4.5,9122860,4.7,col="grey")
text(9121300,4.8,labels="Gr8a")
text(9121350,4.6,labels="> > >")
rect(9119470,4.5,9120721,4.7,col="grey")
text((9119470+9120721)/2,4.8,labels="CG15370")
text((9119470+9120721)/2,4.6,labels="> > >")
rect(9120721,4.5,9121170,4.7,col="white")
text((9120721+9121170)/2,4.8,labels="MSE")
dev.copy2pdf(file="sa_indel_tan_region.pdf")

#pdm3
quartz()
plot_region_indel("2R",4210633,4235314,"pdm3, SA (2012) InDels")
segments(4214995,6.5,42837859,6.5)
rect(c(4214995,4235006),6.3,c(4215512,4235799),6.7,col="grey")
text(4225000,6.75,labels="pdm3")
text(4225000,6.5,labels="> > >")
dev.copy2pdf(file="sa_indels_pdm3_region.pdf")

# ebony
quartz()
plot_region_indel("3R",17062142,17065683,"ebony, SA (2012) InDels")
rect(17055561,0,17062899,0.1,col="grey")
text(17062700,0.15,labels="e")
text(17062650,0.05,labels="< < <")
dev.copy2pdf(file="sa_indels_ebony_region.pdf")

# plot fst values
ms_fst_5K=read.table("/Volumes/Temp/Lukas/Data/M_Schaeffer_PE_bams/Polymorphisms/m_schaeffer_KEN_ETH_ATH_BAR_DIJ_ESX_BOW_LNZ_mc_6_mcov10_win_5K.fst.tab", header=T)

quartz(width=10,height=10)
chroms=c("X","2L","2R","3L","3R")
ylimit=c(0,0.1)
layout(matrix(c(1,6,2,3,4,5),byrow=T, ncol = 2), widths = 1, heights = 1, respect = FALSE)

for (chrom in chroms){
par(mar = c(4.1, 4.1, 4.1, 2.1))
aa=subset(ms_fst_5K,CHR==chrom)
plot(aa$BP, predict(loess(aa$X3.4~aa$BP,span=0.01),newdata=aa$BP), lwd=1,lty=2, col="olivedrab1", ylim=ylimit, xlab="base position", main=chrom, ylab=bquote(F[ST]))
lines(aa$BP, predict(loess(aa$X3.5~aa$BP,span=0.01),newdata=aa$BP), lwd=1,lty=2, col="rosybrown1")
lines(aa$BP,predict(loess(aa$X3.6~aa$BP,span=0.01),newdata=aa$BP), col='green', lwd=1)
lines(aa$BP,predict(loess(aa$X3.7~aa$BP,span=0.01),newdata=aa$BP), col='red', lwd=1)

lines(aa$BP,predict(loess(aa$X3.8~aa$BP,span=0.01),newdata=aa$BP), col='cyan', lwd=1,lty=2)
lines(aa$BP,predict(loess(aa$X4.5~aa$BP,span=0.01),newdata=aa$BP), col='blue', lwd=1, lty=2)
lines(aa$BP,predict(loess(aa$X4.6~aa$BP,span=0.01),newdata=aa$BP), col="darkgoldenrod1", lwd=1,lty=2)
lines(aa$BP,predict(loess(aa$X4.7~aa$BP,span=0.01),newdata=aa$BP), col="darkorange2", lwd=1, lty=2)

## points(aa$BPS, aa$fst_vie, pch=".", col=c("darkgoldenrod1"), cex=0.25)
## points(aa$BPS, aa$fst_ita, pch=".", col=c("darkorange2"), cex=0.25)
## points(aa$BPS, aa$fst_vs, pch=".", col=c("cyan"), cex=0.25)
## points(aa$BPS, aa$fst_is, pch=".", col=c("darkcyan"), cex=0.25)

}
dev.copy2pdf(file="fst_some_ms_athens_barc_against_rest_5K.pdf")

# plot fst values
setwd("/Volumes/vetgrid10/Data/A7_pigmentation/Joined_SA_Ita")

eu.sa.fst.200K=read.table("vie_el_ita_hl_sa_base_2012_joined_pos.win200K.fst.tab", header=T)
eu.sa.fst.5K=read.table("vie_el_ita_hl_sa_base_2012_joined_pos.win5K.fst.tab", header=T)
eu.sa.fst.200K$eu.sa=(eu.sa.fst.200K$X1.3+eu.sa.fst.200K$X2.3)/2

quartz(width=10,height=10)
chroms=c("X","2L","2R","3L","3R")
ylimit=c(0,0.15)
layout(matrix(c(1,6,2,3,4,5),byrow=T, ncol = 2), widths = 1, heights = 1, respect = FALSE)

for (chrom in chroms){
par(mar = c(4.1, 4.1, 4.1, 2.1))
aa=subset(eu.sa.fst.200K,CHR==chrom)
plot(aa$BP, aa$X1.3, lwd=1.5,lty=1, type="l", col="red", ylim=ylimit, xlab="base position", main=chrom, ylab=bquote(F[ST]))
lines(aa$BP, aa$X2.3,lwd=1.5, col="blue")
lines(aa$BP, aa$X1.2,lwd=1.5, col="green")
abline(h=median(aa$X1.3), lwd=1.5,lty=3, col="red")
abline(h=median(aa$X2.3),lwd=1.5,lty=3, col="blue")
abline(h=median(aa$X1.2),lwd=1.5,lty=3, col="green")
}
plot(0:2, 0:2, pch = 1, lty = 1, ylim=c(0,2), type = "n", axes = FALSE, ann = FALSE)
legend(0.75,1.15, legend=c(expression(F[ST] ~ Vienna ~ and ~ Bolzano),expression(F[ST] ~ Vienna ~ and ~ South ~ Africa),expression(F[ST] ~ Bolzano ~ and ~ South ~ Africa),"median"), lty=c(1,1,1,3), lwd=c(1.5,1.5,1.5,1.5),col=c("green","red","blue","black"),xpd=TRUE)
dev.copy2pdf(file="fst_vie_ita_sa_200K.pdf")

summary(eu.sa.fst.200K[eu.sa.fst.200K$CHR %in% chroms,])

eu.sa.fst.5K$eu.sa=(eu.sa.fst.5K$X1.3+eu.sa.fst.5K$X2.3)/2

# get conf intervals:
library("boot")
boot.median <- function(data,indices){
  d <- data[indices]   
  return(median(d, na.rm = TRUE))
}
bfit.aut.ita.200K=boot(eu.sa.fst.200K[eu.sa.fst.200K$CHR %in% chroms,c("X1.2")],statistic=boot.median ,R=10000)
boot.ci(bfit.aut.ita.200K)
# Intervals : 
#   Level      Normal              Basic         
# 95%   ( 0.0104,  0.0110 )   ( 0.0104,  0.0109 )  
# Level     Percentile            BCa          
# 95%   ( 0.0105,  0.0110 )   ( 0.0105,  0.0110 )  
median(eu.sa.fst.200K[eu.sa.fst.200K$CHR %in% chroms,c("X1.2")])
#0.01070555
median(eu.sa.fst.200K[eu.sa.fst.200K$CHR %in% chroms,c("X1.3")])
#0.04393359
bfit.aut.sa.200K=boot(eu.sa.fst.200K[eu.sa.fst.200K$CHR %in% chroms,c("X1.3")],statistic=boot.median ,R=10000)
boot.ci(bfit.aut.sa.200K)
#   Level      Normal              Basic         
# 95%   ( 0.0426,  0.0454 )   ( 0.0429,  0.0455 )  
# Level     Percentile            BCa          
# 95%   ( 0.0424,  0.0450 )   ( 0.0423,  0.0449 )  

median(eu.sa.fst.200K[eu.sa.fst.200K$CHR %in% chroms,c("X2.3")])
#0.04017762
bfit.ita.sa.200K=boot(eu.sa.fst.200K[eu.sa.fst.200K$CHR %in% chroms,c("X2.3")],statistic=boot.median ,R=10000)
boot.ci(bfit.ita.sa.200K)
# medianLevel      Normal              Basic         
# 95%   ( 0.0386,  0.0414 )   ( 0.0386,  0.0413 )  
# Level     Percentile            BCa          
# 95%   ( 0.0391,  0.0418 )   ( 0.0390,  0.0418 )  
(eu.sa.fst.200K[eu.sa.fst.200K$CHR %in% chroms,c("eu.sa")])
#0.0425312
bfit.eu.sa.200K=boot(eu.sa.fst.200K[eu.sa.fst.200K$CHR %in% chroms,c("eu.sa")],statistic=boot.median ,R=10000)
boot.ci(bfit.eu.sa.200K)
# Level      Normal              Basic         
# 95%   ( 0.0411,  0.0441 )   ( 0.0412,  0.0439 )  
# Level     Percentile            BCa          
# 95%   ( 0.0412,  0.0439 )   ( 0.0409,  0.0438 )  



a=list()
for(chrom in chroms){
  bfit.sa=boot(fit.win100K$Fit[fit.win100K$CHR == chrom],statistic=boot.median ,R=1000)
  a[[chrom]]=boot.ci(bfit.sa) 
}



quartz(width=10,height=10)
chroms=c("X","2L","2R","3L","3R")
ylimit=c(0,0.15)
layout(matrix(c(1,6,2,3,4,5),byrow=T, ncol = 2), widths = 1, heights = 1, respect = FALSE)

for (chrom in chroms){
par(mar = c(4.1, 4.1, 4.1, 2.1))
aa=subset(eu.sa.fst.5K,CHR==chrom)
plot(aa$BP, aa$X1.3, lwd=1.5,lty=1, type="l", col="red", ylim=ylimit, xlab="base position", main=chrom, ylab=bquote(F[ST]))
lines(aa$BP, aa$X2.3,lwd=1.5, col="blue")
lines(aa$BP, aa$X1.2,lwd=1.5, col="green")
abline(h=median(aa$X1.3), lwd=1.5,lty=3, col="red")
abline(h=median(aa$X2.3),lwd=1.5,lty=3, col="blue")
abline(h=median(aa$X1.2),lwd=1.5,lty=3, col="green")
}
plot(0:2, 0:2, pch = 1, lty = 1, ylim=c(0,2), type = "n", axes = FALSE, ann = FALSE)
legend(0.75,1.15, legend=c(expression(F[ST] ~ Vienna ~ and ~ Italy),expression(F[ST] ~ Vienna ~ and ~ South ~ Africa),expression(F[ST] ~ Italy ~ and ~ South ~ Africa),"median"), lty=c(1,1,1,3), lwd=c(1.5,1.5,1.5,1.5),col=c("green","red","blue","black"),xpd=TRUE)
dev.copy2pdf(file="fst_vie_ita_sa_5K.pdf")

summary(eu.sa.fst.5K[eu.sa.fst.5K$CHR %in% chroms,])


setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Joined_SA_Ita")
setwd("/Volumes/vetgrid10/Data/A7_pigmentation/Joined_SA_Ita")
a7_inv_data =  read.table("vie_pL1_pL2_pL3_pD1_pD2_pD3_ita_IL1_IL2_IL3_ID1_ID2_ID3_vie_pIel_pIIel_pIIIel_ita_pIhl_pIIhl_pIIIhl_sa_VL1_VL4_VL5_L1_L4_L5_DVD1_DVD4_DVD5_VVD1_VVD4_VVD5_sa_base_merged.syncinversions",header=TRUE)
library(doBy)
head(a7_inv_data)
group_names = c( "IN(2L)t"  ,"IN(2R)NS" ,"IN(3L)P" , "IN(3R)C" , "IN(3R)Mo" ,"IN(3R)P")
colnames(a7_inv_data)=c("CHR","INV","BPS","REF","ALT","ViL1","ViL1c","ViL2","ViL2c","ViL3","ViL3c","ViD1","ViD1c","ViD2","ViD2c","ViD3","ViD3c","IL1","IL1C","IL2","IL2C","IL3","IL3C","ID1","ID1C","ID2","ID2C","ID3","ID3c","VipIel","VipIelc","VipIIel","VipIIelc","VipIIIel","VipIIIelc","IpIhl","IpIhlc","IpIIhl","IpIIhlc","IpIIIhl","IpIIIhlc","SaVL1","SaVL1c","SaVL4","SaVL4c","SaVL5","SaVL5C","SaL1","SaL1c","SaL4","SaL4c","SaL5","SaL5c","SaDVD1","SaDVD1c","SaDVD4","SaDVD4c","SaDVD5","SaDVD5c","SaVVD1","SaVVD1c","SaVVD4","SaVVD4c","SaVVD5","SaVVD5c","SAb","SAbc")


a7_base_vie=c("VipIel","VipIIel","VipIIIel")
a7_base_ita=c("IpIhl","IpIIhl","IpIIIhl")
a7l_vie=c("ViL1","ViL2","ViL3")
a7d_vie=c("ViD1","ViD2","ViD3")
a7l_ita=c("IL1","IL2","IL3")
a7d_ita=c("ID1","ID2","ID3")
a7l_sa=c("SaVL1","SaVL4","SaVL5")
a7d_sa=c("SaVVD1","SaVVD4","SaVVD5")
a7_base_sa=c("SAb")




a7_inv_list <- splitBy(formula= ~ CHR + INV, data = a7_inv_data )


plot.invs <- function(idx){
    mtit=substitute(expr = expression(paste(italic(i))), env = list(i=group_names[idx]) ) 
    boxplot(unlist(a7_inv_list[[idx]][,a7l_vie]),unlist(a7_inv_list[[idx]][,a7_base_vie]),unlist(a7_inv_list[[idx]][a7d_vie]),unlist(a7_inv_list[[idx]][,a7l_ita]),unlist(a7_inv_list[[idx]][,a7_base_ita]),unlist(a7_inv_list[[idx]][a7d_ita]),unlist(a7_inv_list[[idx]][,a7l_sa]),unlist(a7_inv_list[[idx]][,a7_base_sa]),unlist(a7_inv_list[[idx]][a7d_sa]),data=a7_inv_list[[idx]],notch=T,col=c("cornsilk","yellow","brown4"),at=c(0.5,1,1.5,2.5,3,3.5,4.5,5,5.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark", "light", "base","dark"), main=mtit)
                                        #text(c(1,3),-0.05,labels=c("Vienna 2010", "Italy 2011"))
    mtext("Vienna", at=1,line=2,side=1)
    mtext("Italy", at=3,line=2,side=1)
    mtext("South Africa", at=5,line=2,side=1)
}

x11(width=10,height=10)
idx=1
#ylimit=c(0,0.1)
layout(matrix(c(1,2,3,4,5,6),byrow=T, ncol = 2), widths = 1, heights = 1, respect = FALSE)

for (idx in 1:length(group_names)){
par(mar = c(4.1, 4.1, 4.1, 2.1))
plot.invs(idx)
}
dev.copy2pdf(file="inv_all_vie_ita_sa_A7.pdf")

fill_data = function(raw_df,reps=3,pops=c("x")){
    # takes a list of inversions with replicates (each with chr, bps, allele ..., freq1,count1, ,,,) with different populations and returns a list with count,cov,rep,conditions
    # the sequence of the populations needs to be: pop1|r1|cond1 pop1|r2|cond1 pop1|r1|cond2 pop1|r2|cond2 pop2|r1|cond1 pop2|r2|cond1 ...
    raw_list=splitBy(formula=~ CHR + INV, data = raw_df)
    data_list = list()
    for(inversion in names(raw_list)) {
        int_tab=raw_list[[inversion]]
        for(i in seq(from=6,to=length(int_tab[1,]),by=2)){int_tab[,i]=int_tab[,i]*int_tab[,i+1]}
        marl_mat=matrix(nrow=0,ncol=5)
	colnames(marl_mat)=c("count","cov","rep","cs","pop")
        for (pop in 1:length(pops)){
            for(rep in 1:reps){
                pst=6+(pop-1)*(reps*4)
		j=(rep-1)*2
		loc_mat=cbind(int_tab[,(pst+j):(pst+j+1)],rep,0,pops[pop])
		colnames(loc_mat)=c("count","cov","rep","cs","pop")
		marl_mat=rbind(marl_mat,loc_mat)
                j = j + reps*2
		loc_mat=cbind(int_tab[,(pst+j):(pst+j+1)],rep,1,pops[pop])
		colnames(loc_mat)=c("count","cov","rep","cs","pop")
		marl_mat=rbind(marl_mat,loc_mat)
            }
        }
	rownames(marl_mat)=NULL
	data=as.data.frame(marl_mat)	
	data$freq<-round((data$count/data$cov)*100,digits=0)
	data$inv_count<-round(data$count,digits=0)
	data$not_inv_count<- data$cov-data$inv_count
	data$experiment<-as.factor(as.numeric(data$pop)*100+ data$rep*10+data$cs)
	data$rep_fac<-as.factor(data$rep+as.numeric(data$pop)*100)
	data$rep<-as.factor(data$rep)
	data$cs_fac<-as.factor(data$cs)
	data_list[[inversion]] = data		
	}		
	return(data_list)	
}	


library(doBy)
a7.vie=c("CHR","INV","BPS","REF","ALT","ViL1","ViL1c","ViL2","ViL2c","ViL3","ViL3c","ViD1","ViD1c","ViD2","ViD2c","ViD3","ViD3c")
a7.ita=c("CHR","INV","BPS","REF","ALT","IL1","IL1C","IL2","IL2C","IL3","IL3C","ID1","ID1C","ID2","ID2C","ID3","ID3c")
a7.sa=c("CHR","INV","BPS","REF","ALT","SaVL1","SaVL1c","SaVL4","SaVL4c","SaVL5","SaVL5C","SaVVD1","SaVVD1c","SaVVD4","SaVVD4c","SaVVD5","SaVVD5c")
a7.eu=c("CHR","INV","BPS","REF","ALT","ViL1","ViL1c","ViL2","ViL2c","ViL3","ViL3c","ViD1","ViD1c","ViD2","ViD2c","ViD3","ViD3c","IL1","IL1C","IL2","IL2C","IL3","IL3C","ID1","ID1C","ID2","ID2C","ID3","ID3c")
a7.all=c("CHR","INV","BPS","REF","ALT","ViL1","ViL1c","ViL2","ViL2c","ViL3","ViL3c","ViD1","ViD1c","ViD2","ViD2c","ViD3","ViD3c","IL1","IL1C","IL2","IL2C","IL3","IL3C","ID1","ID1C","ID2","ID2C","ID3","ID3c","SaVL1","SaVL1c","SaVL4","SaVL4c","SaVL5","SaVL5C","SaVVD1","SaVVD1c","SaVVD4","SaVVD4c","SaVVD5","SaVVD5c")

a7.inv.vie=a7_inv_data[,a7.vie]
a7.inv.ita=a7_inv_data[,a7.ita]
a7.inv.sa=a7_inv_data[,a7.sa]
a7.inv.eu=a7_inv_data[,a7.eu]
a7.inv.all=a7_inv_data[,a7.all]
a7.ic.vie <- splitBy(formula=~ CHR + INV, data = a7.inv.vie )
data.ic.vie = fill_data(a7.inv.vie)
data.ic.ita = fill_data(a7.inv.ita)
data.ic.sa = fill_data(a7.inv.sa)
data.ic.eu = fill_data(a7.inv.eu,pops=c("vie","ita"))
data.ic.all = fill_data(a7.inv.all,pops=c("vie","ita","sa"))


library("lme4")
lmer.res.vie=list()
lmer.res.ita=list()
lmer.res.sa=list()
lmer.res.eu=list()
lmer.res.all=list()
## for (inversion in names(data.ic.vie)){
## 	lmer.res.vie[[inversion]]<-lmer(freq ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data.ic.vie[[inversion]], family=poisson)
## 	lmer.res.ita[[inversion]]<-lmer(freq ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data.ic.ita[[inversion]], family=poisson)
## 	lmer.res.sa[[inversion]]<-lmer(freq ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data.ic.sa[[inversion]], family=poisson)
## 	lmer.res.eu[[inversion]]<-lmer(freq ~ 1+  (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data.ic.eu[[inversion]], family=poisson)
## 	lmer.res.all[[inversion]]<-lmer(freq ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data.ic.all[[inversion]], family=poisson)
## }
for (inversion in names(data.ic.vie)){
	lmer.res.vie[[inversion]]<-glmer(cbind(inv_count, not_inv_count) ~ 1+ (1|rep_fac) + (1|rep_fac/experiment) + cs_fac, data=data.ic.vie[[inversion]], family=binomial,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e9)))
	lmer.res.ita[[inversion]]<-glmer(cbind(inv_count, not_inv_count) ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data.ic.ita[[inversion]], family=binomial,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e9)))
	lmer.res.sa[[inversion]]<-glmer(cbind(inv_count, not_inv_count) ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data.ic.sa[[inversion]], family=binomial,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e9)))
#	lmer.res.eu[[inversion]]<-glmer(cbind(inv_count, not_inv_count) ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data.ic.eu[[inversion]], family=binomial,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e9)))
#	lmer.res.all[[inversion]]<-glmer(cbind(inv_count, not_inv_count) ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data.ic.all[[inversion]], family=binomial,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e9)))
}
warnings()
lmer.res.eu1=list()
lmer.res.all2=list()

for (inversion in names(data.ic.vie)){
	lmer.res.eu1[[inversion]]<-glmer(cbind(inv_count, not_inv_count) ~ 1 +  (1|pop) + (1|rep_fac/experiment) + cs_fac, data=data.ic.eu[[inversion]], family=binomial,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e9)))
	lmer.res.all2[[inversion]]<-glmer(cbind(inv_count, not_inv_count) ~ 1 +  (1|pop) + (1|rep_fac/experiment) + cs_fac, data=data.ic.all[[inversion]], family=binomial,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e9)))
}

summary(lmer.res.eu1[[inversion]])
summary(lmer.res.all2[[inversion]])
for (inversion in names(lmer.res.vie)){ print(inversion); print(summary(lmer.res.vie[[inversion]])$coefficients) }
for (inversion in names(lmer.res.vie)){ print(inversion); print(summary(lmer.res.ita[[inversion]])$coefficients) }
for (inversion in names(lmer.res.vie)){ print(inversion); print(summary(lmer.res.sa[[inversion]])$coefficients) }
#for (inversion in names(lmer.res.vie)){ print(inversion); print(summary(lmer.res.eu[[inversion]])$coefficients) }
#for (inversion in names(lmer.res.vie)){ print(inversion); print(summary(lmer.res.all[[inversion]])$coefficients) }
for (inversion in names(lmer.res.vie)){ print(inversion); print(summary(lmer.res.eu1[[inversion]])$coefficients) }
for (inversion in names(lmer.res.vie)){ print(inversion); print(summary(lmer.res.all2[[inversion]])$coefficients) }


x11()
boxplot(unlist(a7_inv_list[[6]][,a7l_vie]),unlist(a7_inv_list[[6]][,a7_base_vie]),unlist(a7_inv_list[[6]][a7d_vie]),unlist(a7_inv_list[[6]][,a7l_ita]),unlist(a7_inv_list[[6]][,a7_base_ita]),unlist(a7_inv_list[[6]][a7d_ita]),unlist(a7_inv_list[[6]][,a7l_sa]),unlist(a7_inv_list[[6]][,a7_base_sa]),unlist(a7_inv_list[[6]][a7d_sa]),data=a7_inv_list[[6]],notch=T,col=c("cornsilk","yellow","brown4"),at=c(0.5,1,1.5,2.5,3,3.5,4.5,5,5.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark", "light", "base","dark"), main="Inversion 3R(p)")
boxplot(c(a7_inv_list[[6]][,c(a7l_vie,a7_base_vie,a7d_vie,a7l_ita,a7_base_ita,a7d_ita,a7l_sa,a7_base_sa,a7d_sa)]),data=a7_inv_list[[6]],notch=T,col=c(rep("cornsilk",3),rep("yellow",3),rep("brown4",3)))
        at=c(0.5,1,1.5,2.5,3,3.5,4.5,5,5.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark", "light", "base","dark"), main="Inversion 3R(p)")
#text(c(1,3),-0.05,labels=c("Vienna 2010", "Italy 2011"))
mtext("Vienna 2010", at=1,line=2,side=1)
mtext("Italy 2011", at=3,line=2,side=1)
mtext("South Africa 2012", at=5,line=2,side=1)
dev.copy2pdf(file="inv_3Rp_vie_ita_sa_A7.pdf")

quartz()
boxplot(unlist(a7_inv_list[[1]][,a7l_vie]),unlist(a7_inv_list[[1]][,a7_base_vie]),unlist(a7_inv_list[[1]][a7d_vie]),unlist(a7_inv_list[[1]][,a7l_ita]),unlist(a7_inv_list[[1]][,a7_base_ita]),unlist(a7_inv_list[[1]][a7d_ita]),unlist(a7_inv_list[[1]][,a7l_sa]),unlist(a7_inv_list[[1]][,a7_base_sa]),unlist(a7_inv_list[[1]][a7d_sa]),data=a7_inv_list[[1]],notch=T,col=c("cornsilk","yellow","brown4"),at=c(0.5,1,1.5,2.5,3,3.5,4.5,5,5.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark", "light", "base","dark"), main="Inversion 2L(t)")
#text(c(1,3),-0.05,labels=c("Vienna 2010", "Italy 2011"))
mtext("Vienna 2010", at=1,line=2,side=1)
mtext("Italy 2011", at=3,line=2,side=1)
mtext("South Africa 2012", at=5,line=2,side=1)
dev.copy2pdf(file="inv_2Lt_vie_ita_sa_A7.pdf")

quartz()
boxplot(unlist(a7_inv_list[[2]][,a7l_vie]),unlist(a7_inv_list[[2]][,a7_base_vie]),unlist(a7_inv_list[[2]][a7d_vie]),unlist(a7_inv_list[[2]][,a7l_ita]),unlist(a7_inv_list[[2]][,a7_base_ita]),unlist(a7_inv_list[[2]][a7d_ita]),unlist(a7_inv_list[[2]][,a7l_sa]),unlist(a7_inv_list[[2]][,a7_base_sa]),unlist(a7_inv_list[[2]][a7d_sa]),data=a7_inv_list[[2]],notch=T,col=c("cornsilk","yellow","brown4"),at=c(0.5,1,1.5,2.5,3,3.5,4.5,5,5.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark", "light", "base","dark"), main="Inversion 2R(ns)")
#text(c(1,3),-0.05,labels=c("Vienna 2010", "Italy 2011"))
mtext("Vienna 2010", at=1,line=2,side=1)
mtext("Italy 2011", at=3,line=2,side=1)
mtext("South Africa 2012", at=5,line=2,side=1)
dev.copy2pdf(file="inv_2Rns_vie_ita_sa_A7.pdf")

quartz()
boxplot(unlist(a7_inv_list[[3]][,a7l_vie]),unlist(a7_inv_list[[3]][,a7_base_vie]),unlist(a7_inv_list[[3]][a7d_vie]),unlist(a7_inv_list[[3]][,a7l_ita]),unlist(a7_inv_list[[3]][,a7_base_ita]),unlist(a7_inv_list[[3]][a7d_ita]),unlist(a7_inv_list[[3]][,a7l_sa]),unlist(a7_inv_list[[3]][,a7_base_sa]),unlist(a7_inv_list[[3]][a7d_sa]),data=a7_inv_list[[3]],notch=T,col=c("cornsilk","yellow","brown4"),at=c(0.5,1,1.5,2.5,3,3.5,4.5,5,5.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark", "light", "base","dark"), main="Inversion 3L(p)")
#text(c(1,3),-0.05,labels=c("Vienna 2010", "Italy 2011"))
mtext("Vienna 2010", at=1,line=2,side=1)
mtext("Italy 2011", at=3,line=2,side=1)
mtext("South Africa 2012", at=5,line=2,side=1)
dev.copy2pdf(file="inv_3Lp_vie_ita_sa_A7.pdf")

quartz()
boxplot(unlist(a7_inv_list[[4]][,a7l_vie]),unlist(a7_inv_list[[4]][,a7_base_vie]),unlist(a7_inv_list[[4]][a7d_vie]),unlist(a7_inv_list[[4]][,a7l_ita]),unlist(a7_inv_list[[4]][,a7_base_ita]),unlist(a7_inv_list[[4]][a7d_ita]),unlist(a7_inv_list[[4]][,a7l_sa]),unlist(a7_inv_list[[4]][,a7_base_sa]),unlist(a7_inv_list[[4]][a7d_sa]),data=a7_inv_list[[4]],notch=T,col=c("cornsilk","yellow","brown4"),at=c(0.5,1,1.5,2.5,3,3.5,4.5,5,5.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark", "light", "base","dark"), main="Inversion 3R(c)")
#text(c(1,3),-0.05,labels=c("Vienna 2010", "Italy 2011"))
mtext("Vienna 2010", at=1,line=2,side=1)
mtext("Italy 2011", at=3,line=2,side=1)
mtext("South Africa 2012", at=5,line=2,side=1)
dev.copy2pdf(file="inv_3Rc_vie_ita_sa_A7.pdf")

quartz()
boxplot(unlist(a7_inv_list[[5]][,a7l_vie]),unlist(a7_inv_list[[5]][,a7_base_vie]),unlist(a7_inv_list[[5]][a7d_vie]),unlist(a7_inv_list[[5]][,a7l_ita]),unlist(a7_inv_list[[5]][,a7_base_ita]),unlist(a7_inv_list[[5]][a7d_ita]),unlist(a7_inv_list[[5]][,a7l_sa]),unlist(a7_inv_list[[5]][,a7_base_sa]),unlist(a7_inv_list[[5]][a7d_sa]),data=a7_inv_list[[5]],notch=T,col=c("cornsilk","yellow","brown4"),at=c(0.5,1,1.5,2.5,3,3.5,4.5,5,5.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark", "light", "base","dark"), main="Inversion 3R(mo)")
#text(c(1,3),-0.05,labels=c("Vienna 2010", "Italy 2011"))
mtext("Vienna 2010", at=1,line=2,side=1)
mtext("Italy 2011", at=3,line=2,side=1)
mtext("South Africa 2012", at=5,line=2,side=1)
dev.copy2pdf(file="inv_3Rmo_vie_ita_sa_A7.pdf")

# load pV and OR files for vie, ita, eu, sa and eusa
setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Joined_SA_Ita")
eu_sa_odds=read.table("vie_LD_ita_LD_eu_LD_sa_VLVD_eusa_VLVD_merged.one_pV_gt_1e_3.pV_OR.gz",na.strings = c("NA","nan"),header=F)
colnames(eu_sa_odds)=c("CHR","BPS","Allele","P_vi","O_vi","Oa_vi","Ob_vi","P_ita","O_ita","Oa_ita","Ob_ita","P_eu","O_eu","Oa_eu","Ob_eu","P_sa","O_sa","Oa_sa","Ob_sa","P_esi","O_esi","Oa_esi","Ob_esi")
P=c("P_vi","P_ita","P_eu","P_sa")
lP=c("lP_vi","lP_ita","lP_eu","lP_sa")
O=c("O_vi","O_ita","O_eu","O_sa")
lO=c("lO_vi","lO_ita","lO_eu","lO_sa")
eu_sa_odds[,lP]=-log10(eu_sa_odds[,P])
eu_sa_odds[,lO]=log(eu_sa_odds[,O])
summary(eu_sa_odds$lO_eu[is.finite(eu_sa_odds$lO_eu)]) # Min: -4.67500, Max: 4.11500
summary(eu_sa_odds$lO_sa[is.finite(eu_sa_odds$lO_sa)]) # Min: -4.20200, Max: 3.94600
for(i in  lO){
    #eu_sa_odds[is.na(eu_sa_odds[,i]),i] = 0
    eu_sa_odds[(eu_sa_odds[,i] == -Inf) & !(is.na(eu_sa_odds[,i])),i] = -5
    eu_sa_odds[(eu_sa_odds[,i] == Inf) & !(is.na(eu_sa_odds[,i])),i] = 5
}
eu_sa_odds_base=read.table("vie_LB_DB_ita_LB_DB_eu_LB_DB_sa_VLB_LB_DB_VDB_merged.one_pV_gt_1e_3.OR.gz",na.strings = c("NA","nan"),header=F)
colnames(eu_sa_odds_base)=c("CHR","BPS","Allele","Ol_vi","Od_vi","Ol_ita","Od_ita","Ol_eu","Od_eu","Ol_sa","Od_sa")
eu_sa_odds_base=merge(eu_sa_odds,eu_sa_odds_base)

Ob=c("Ol_vi","Od_vi","Ol_ita","Od_ita","Ol_eu","Od_eu","Ol_sa","Od_sa")
lOb=c("Ol_vi","Od_vi","Ol_ita","Od_ita","Ol_eu","Od_eu","Ol_sa","Od_sa")
eu_sa_odds_base[,lOb]=log(eu_sa_odds_base[,Ob])
for(i in  lOb){
    #eu_sa_odds_base[is.na(eu_sa_odds_base[,i]),i] = 0
    eu_sa_odds_base[(eu_sa_odds_base[,i] == -Inf) & !(is.na(eu_sa_odds_base[,i])),i] = -5
    eu_sa_odds_base[(eu_sa_odds_base[,i] == Inf) & !(is.na(eu_sa_odds_base[,i])),i] = 5
}

eu_sa_freqs=read.table('vie_pL1_pL2_pL3_pD1_pD2_pD3_ita_IL1_IL2_IL3_ID1_ID2_ID3_vie_pIel_pIIel_pIIIel_ita_pIhl_pIIhl_pIIIhl_sa_VL1_VL4_VL5_L1_L4_L5_DVD1_DVD4_DVD5_VVD1_VVD4_VVD5_sa_base_merged.one_pV_lt_1e3.af.gz',header=F,na.strings=c("NA","NaN","nan","NAN","N"))
colnames(eu_sa_freqs)=c("CHR","BPS","Allele","VL1","VL2","VL3","VD1","VD2","VD3","IL1","IL2","IL3","ID1","ID2","ID3","Vb1","Vb2","Vb3","Ib1","Ib2","Ib3","SVL1","SVL4","SVL5","SL1","SL4","SL5","SD1","SD4","SD5","SVD1","SVD4","SVD5","Sb")
# adding extended tails
xs.ind=c(61,82,62,78,60,65,93,131,100,67,55,70)
xs.pops=c("SVL1","SVL4","SVL5","SL1","SL4","SL5","SD1","SD4","SD5","SVD1","SVD4","SVD5")
xs.ind.l=sapply(c(1,2,3,7,8,9),function(x) xs.ind[x] + xs.ind[x+3] )
xs.freqs=sapply(c(1,2,3,7,8,9),function(x) (eu_sa_freqs[,xs.pops[x]]*xs.ind[x] + eu_sa_freqs[,xs.pops[x+3]]*xs.ind[x+3])/(xs.ind[x] + xs.ind[x+3]) )
xs.popname=c("XSL1","XSL4","XSL5","XSD1","XSD4","XSD5")
eu_sa_freqs[,xs.popname]=xs.freqs
pop.means=list(c("VL1","VL2","VL3"),c("VD1","VD2","VD3"),c("IL1","IL2","IL3"),c("ID1","ID2","ID3"),c("Vb1","Vb2","Vb3"),c("Ib1","Ib2","Ib3"),c("SVL1","SVL4","SVL5"),c("SL1","SL4","SL5"),c("SD1","SD4","SD5"),c("SVD1","SVD4","SVD5"),c("VL1","VL2","VL3","IL1","IL2","IL3"),c("VD1","VD2","VD3","ID1","ID2","ID3"),c("Vb1","Vb2","Vb3","Ib1","Ib2","Ib3"),c("XSL1","XSL4","XSL5"),c("XSD1","XSD4","XSD5"))
pop.means.names=c("VL","VD","IL","ID","Vb","Ib","SVL","SL","SD","SVD","EUL","EUD","EUb","XSL","XSD")
for(i in 1:length(pop.means.names)){
    eu_sa_freqs[,pop.means.names[i]]=rowMeans(eu_sa_freqs[,unlist(pop.means[i])])
}

col_freqs=c("CHR","BPS","Allele","VL","VD","IL","ID","Vb","Ib","SVL","SL","SD","SVD","XSL","XSD","EUL","EUD","EUb","Sb")
col_freqs=c("CHR","BPS","Allele","VL","VD","IL","ID","Vb","Ib","SVL","SL","SD","SVD","EUL","EUD","EUb","Sb")

col_odds=c("CHR","BPS","Allele","P_vi","P_ita","P_eu","P_sa","P_xs","lO_vi","lO_ita","lO_sa","lO_eu","lO_xs")
col_odds=c("CHR","BPS","Allele","P_vi","P_ita","P_eu","P_sa","lO_vi","lO_ita","lO_sa","lO_eu")
eu.sa.freqs=merge(eu_sa_freqs[,col_freqs],vie_ita_sa_odds[,col_odds])

# including haplotype information for vie, ita, eu and sa
hap_vie=read.table('Haplotypes/vie_haplo_pV_lt_1_e_8_LL_DD_base_light_dark_total.pairs_mc_75_maxP0.05_mr2_0.75.hap_sync',header=F)
colnames(hap_vie)=c("CHR","BPS","HAP")
hap_ita=read.table('Haplotypes/ita_haplo_pV_lt_1_e_8_LL_DD_base_light_dark_total.pairs_mc_75_maxP0.05_mr2_0.75.hap_sync',header=F)
colnames(hap_ita)=c("CHR","BPS","HAP")
hap_eu=read.table('Haplotypes/eu_haplo_pV_lt_1_e_8_LL_DD_base_light_dark_total.pairs_mc_75_maxP0.05_mr2_0.75.hap_sync',header=F)
colnames(hap_eu)=c("CHR","BPS","HAP")
hap_sa=read.table('Haplotypes/sa_haplo_pV_lt_1_e_8_LL_DD_base_light_dark_total.pairs_mc_100_maxP0.05_mr2_0.75.hap_sync',header=F)
colnames(hap_sa)=c("CHR","BPS","HAP")

# add factors to data frame
haps=c("hap_vie","hap_ita","hap_eu","hap_sa")
eu_sa_odds_base[,haps]=as.factor(1:length(eu_sa_odds_base[,1]))
levels(eu_sa_odds_base$hap_vie) <- c(levels(eu_sa_odds_base$hap_vie),levels(hap_vie$HAP))
for(i in 1:length(hap_vie$CHR)){
eu_sa_odds_base$hap_vie[eu_sa_odds_base$CHR == as.character(hap_vie$CHR[i]) & eu_sa_odds_base$BPS == hap_vie$BPS[i] ] = as.character(hap_vie$HAP[i])
}
levels(eu_sa_odds_base$hap_ita) <- c(levels(eu_sa_odds_base$hap_ita),levels(hap_ita$HAP))
for(i in 1:length(hap_ita$CHR)){
eu_sa_odds_base$hap_ita[eu_sa_odds_base$CHR == as.character(hap_ita$CHR[i]) & eu_sa_odds_base$BPS == hap_ita$BPS[i] ] = as.character(hap_ita$HAP[i])
}
levels(eu_sa_odds_base$hap_eu) <- c(levels(eu_sa_odds_base$hap_eu),levels(hap_eu$HAP))
for(i in 1:length(hap_eu$CHR)){
eu_sa_odds_base$hap_eu[eu_sa_odds_base$CHR == as.character(hap_eu$CHR[i]) & eu_sa_odds_base$BPS == hap_eu$BPS[i] ] = as.character(hap_eu$HAP[i])
}
levels(eu_sa_odds_base$hap_sa) <- c(levels(eu_sa_odds_base$hap_sa),levels(hap_sa$HAP))
for(i in 1:length(hap_sa$CHR)){
eu_sa_odds_base$hap_sa[eu_sa_odds_base$CHR == as.character(hap_sa$CHR[i]) & eu_sa_odds_base$BPS == hap_sa$BPS[i] ] = as.character(hap_sa$HAP[i])
}

# for eu.sa.odds
#hap_eu=read.table('Haplotypes/eu_haplo_pV_lt_1_e_6_LL_DD_base_light_dark_total.pairs_mc_25_maxP0.015_mr2_0.85.hap_sync',header=F)
#colnames(hap_eu)=c("CHR","BPS","HAP")
#hap_sa=read.table('Haplotypes/sa_haplo_pV_lt_1_e_6_LL_DD_base_light_dark_total.pairs_mc_12_maxP0.015_mr2_0.85.hap_sync',header=F)
#colnames(hap_sa)=c("CHR","BPS","HAP")


eu.sa.odds.cols=c("CHR","BPS","Allele","P_eu","O_eu","Oa_eu","Ob_eu","P_sa","O_sa","Oa_sa","Ob_sa","P_xs","O_xs","Oa_xs","Ob_xs")
eu.sa.odds.cols=c("CHR","BPS","Allele","P_eu","O_eu","Oa_eu","Ob_eu","P_sa","O_sa","Oa_sa","Ob_sa")
eu.sa.odds=vie_ita_sa_odds[,eu.sa.odds.cols]
eu.sa.odds=eu_sa_odds
eu.sa.odds$P_eu=-log10(eu.sa.odds$P_eu)
eu.sa.odds$P_sa=-log10(eu.sa.odds$P_sa)
eu.sa.odds$O_eu=log( eu.sa.odds$O_eu)
eu.sa.odds$O_sa=log( eu.sa.odds$O_sa)
eu.sa.odds$Ob_eu=log( eu.sa.odds$Ob_eu)
eu.sa.odds$Ob_sa=log( eu.sa.odds$Ob_sa)
eu.sa.odds$Oa_eu=log( eu.sa.odds$Oa_eu)
eu.sa.odds$Oa_sa=log( eu.sa.odds$Oa_sa)
pV.high=( ( ! is.na(eu.sa.odds$P_eu) & eu.sa.odds$P_eu >= 6) | ( ! is.na(eu.sa.odds$P_sa) & eu.sa.odds$P_sa >= 6)  )
eu.sa.odds=eu.sa.odds[pV.high ,]
haps=c("hap_eu","hap_sa")
eu.sa.odds[,haps]=as.factor(1:length(eu.sa.odds[,1]))
levels(eu.sa.odds$hap_eu) <- c(levels(eu.sa.odds$hap_eu),levels(hap_eu$HAP))
for(i in 1:length(hap_eu$CHR)){
eu.sa.odds$hap_eu[eu.sa.odds$CHR == as.character(hap_eu$CHR[i]) & eu.sa.odds$BPS == hap_eu$BPS[i] ] = as.character(hap_eu$HAP[i])
}
levels(eu.sa.odds$hap_sa) <- c(levels(eu.sa.odds$hap_sa),levels(hap_sa$HAP))
for(i in 1:length(hap_sa$CHR)){
eu.sa.odds$hap_sa[eu.sa.odds$CHR == as.character(hap_sa$CHR[i]) & eu.sa.odds$BPS == hap_sa$BPS[i] ] = as.character(hap_sa$HAP[i])
}

head(eu.sa.odds)
eu.sa.odds$LF=0
maxPs=apply(eu.sa.odds[,c("P_eu","P_sa")],1,function(x) max(x,na.rm = T))
eu.sa.odds$LF[! is.na(eu.sa.odds$P_eu) & eu.sa.odds$P_eu >= maxPs ]=sign(eu.sa.odds$O_eu[! is.na(eu.sa.odds$P_eu) & eu.sa.odds$P_eu >= maxPs])
#eu.sa.odds$LF[! is.na(eu.sa.odds$O_xs) & ! is.na(eu.sa.odds$P_xs) & eu.sa.odds$P_xs  >= maxPs]=sign(eu.sa.odds$O_xs[! is.na(eu.sa.odds$O_xs) & ! is.na(eu.sa.odds$P_xs) & eu.sa.odds$P_xs  >= maxPs])
eu.sa.odds$LF[! is.na(eu.sa.odds$O_sa) & ! is.na(eu.sa.odds$P_sa) & eu.sa.odds$P_sa  >= maxPs]=sign(eu.sa.odds$O_sa[! is.na(eu.sa.odds$O_sa) & ! is.na(eu.sa.odds$P_sa) & eu.sa.odds$P_sa  >= maxPs])
eu.sa.odds[,c("lO_sa","lOa_sa","lOb_sa","lO_eu","lOa_eu","lOb_eu")] = eu.sa.odds[,c("O_sa","Oa_sa","Ob_sa","O_eu","Oa_eu","Ob_eu")]*eu.sa.odds$LF

# for frequencies
eu_sa_freqs[,haps]=as.factor(1:length(eu_sa_freqs[,1]))
levels(eu_sa_freqs$hap_vie) <- c(levels(eu_sa_freqs$hap_vie),levels(hap_vie$HAP))
for(i in 1:length(hap_vie$CHR)){
eu_sa_freqs$hap_vie[eu_sa_freqs$CHR == as.character(hap_vie$CHR[i]) & eu_sa_freqs$BPS == hap_vie$BPS[i] ] = as.character(hap_vie$HAP[i])
}
levels(eu_sa_freqs$hap_ita) <- c(levels(eu_sa_freqs$hap_ita),levels(hap_ita$HAP))
for(i in 1:length(hap_ita$CHR)){
eu_sa_freqs$hap_ita[eu_sa_freqs$CHR == as.character(hap_ita$CHR[i]) & eu_sa_freqs$BPS == hap_ita$BPS[i] ] = as.character(hap_ita$HAP[i])
}
levels(eu_sa_freqs$hap_eu) <- c(levels(eu_sa_freqs$hap_eu),levels(hap_eu$HAP))
for(i in 1:length(hap_eu$CHR)){
eu_sa_freqs$hap_eu[eu_sa_freqs$CHR == as.character(hap_eu$CHR[i]) & eu_sa_freqs$BPS == hap_eu$BPS[i] ] = as.character(hap_eu$HAP[i])
}
levels(eu_sa_freqs$hap_sa) <- c(levels(eu_sa_freqs$hap_sa),levels(hap_sa$HAP))
for(i in 1:length(hap_sa$CHR)){
eu_sa_freqs$hap_sa[eu_sa_freqs$CHR == as.character(hap_sa$CHR[i]) & eu_sa_freqs$BPS == hap_sa$BPS[i] ] = as.character(hap_sa$HAP[i])
}

eu.sa.freqs$P_eu=-log10(eu.sa.freqs$P_eu)
eu.sa.freqs$P_sa=-log10(eu.sa.freqs$P_sa)

pV.high=( ( ! is.na(eu.sa.freqs$P_eu) & eu.sa.freqs$P_eu >= 6) | ( ! is.na(eu.sa.freqs$P_sa) & eu.sa.freqs$P_sa >= 6)  )
eu.sa.freqs.p6=eu.sa.freqs[pV.high,]
hap_eu=read.table('Haplotypes/eu_haplo_pV_lt_1_e_6_LL_DD_base_light_dark_total.pairs_mc_25_maxP0.015_mr2_0.75.hap_sync',header=F)
colnames(hap_eu)=c("CHR","BPS","HAP")
hap_sa=read.table('Haplotypes/sa_haplo_pV_lt_1_e_6_LL_DD_base_light_dark_total.pairs_mc_12_maxP0.015_mr2_0.75.hap_sync',header=F)
colnames(hap_sa)=c("CHR","BPS","HAP")
haps=c("hap_eu","hap_sa")
eu.sa.freqs.p6[,haps]=as.factor(1:length(eu.sa.freqs.p6[,1]))
levels(eu.sa.freqs.p6$hap_eu) <- c(levels(eu.sa.freqs.p6$hap_eu),levels(hap_eu$HAP))
for(i in 1:length(hap_eu$CHR)){
eu.sa.freqs.p6$hap_eu[eu.sa.freqs.p6$CHR == as.character(hap_eu$CHR[i]) & eu.sa.freqs.p6$BPS == hap_eu$BPS[i] ] = as.character(hap_eu$HAP[i])
}
levels(eu.sa.freqs.p6$hap_sa) <- c(levels(eu.sa.freqs.p6$hap_sa),levels(hap_sa$HAP))
for(i in 1:length(hap_sa$CHR)){
eu.sa.freqs.p6$hap_sa[eu.sa.freqs.p6$CHR == as.character(hap_sa$CHR[i]) & eu.sa.freqs.p6$BPS == hap_sa$BPS[i] ] = as.character(hap_sa$HAP[i])
}



# adapt graphics formats
    
plot_af_all_pV <- function(region,coords,lpVals,ylab_pv=expression(-log[10](italic(P))) ) {
    region_idx = which(eu_sa_freqs$CHR == coords[1] & (coords[2] <= eu_sa_freqs$BPS & eu_sa_freqs$BPS <= coords[3] ) & (eu_sa_freqs$lP_eu >= lpVals[1] | eu_sa_freqs$lP_sa >= lpVals[2]))
    length(region_idx)
    layout(matrix(1:2, ncol = 1), widths = 1, heights = c(0.75,1.25), respect = FALSE)
    par(mar = c(0, 4.1, 4.1, 2.1))
    main_tit=bquote("SNPs around" ~.(region) ~"with" ~ Peu <= 10^-.(lpVals[1]) ~ " or " ~  Psa <= 10 ^-.(lpVals[2]))
    xlimit=c(1-0.25,length(region_idx)+0.15)
    plot(1:length(region_idx),eu_sa_freqs$lP_sa[region_idx],col="blue",cex=1,pch=16,ylim=c(0,max(eu_sa_freqs[region_idx,c("lP_sa","lP_eu")],na.rm=T)),ylab=ylab_pv, xaxt="n",main=main_tit,xlim=xlimit )
    abline(v=(1:length(region_idx))-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+0.25,col="grey",lty="dotted",lw=2)
    points(1:length(region_idx),eu_sa_freqs$lP_eu[region_idx],col="red",cex=1,pch=16)
    legend("topright",legend=c("SA","Europe"),pch=c(16,20),col=c("blue","red"))
    par(mar = c(4.1, 4.1, 0, 2.1))
    ymin=min(eu_sa_freqs[region_idx,c("SVL","SL","SD","SVD","EUL","EUD","EUb","Sb")],na.rm=T)
    plot((1:length(region_idx))-0.15,eu_sa_freqs[region_idx,c("Sb")],col="blue",cex=1,pch=".",ylim=c(ymin,1),ylab="AF",xlab="SNPs",xlim=xlimit )
    abline(v=(1:length(region_idx))-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+0.25,col="grey",lty="dotted",lw=2)
    abline(h=seq(0,1,by=0.1),col="grey",lty="dotted",lw=2)
    points((1:length(region_idx))-0.15,eu_sa_freqs[region_idx,c("Sb")],col="blue",cex=1,pch=16)
    points((1:length(region_idx))-0.15,eu_sa_freqs[region_idx,c("SVL")],col="blue",cex=1,pch=6)
    points((1:length(region_idx))-0.15,eu_sa_freqs[region_idx,c("SVD")],col="blue",cex=1,pch=17)
     points((1:length(region_idx))-0.15,eu_sa_freqs[region_idx,c("SL")],col="blue",cex=0.75,pch=23)
    points((1:length(region_idx))-0.15,eu_sa_freqs[region_idx,c("SD")],col="blue",cex=0.75,pch=23,bg="cyan")

    points((1:length(region_idx))+0.15,eu_sa_freqs[region_idx,c("EUb")],col="red",cex=1,pch=16)
    points((1:length(region_idx))+0.15,eu_sa_freqs[region_idx,c("EUL")],col="red",cex=1,pch=6)
    points((1:length(region_idx))+0.15,eu_sa_freqs[region_idx,c("EUD")],col="red",cex=1,pch=17)
    legend("bottomright",legend=c("very light","light","base","dark","very dark"),title="mean AF",pch=c(6,23,16,23,17),col=c("black"),pt.bg=c("white","white","white","grey16","black"),pt.cex =c(1,0.75,1,0.75,1))
}


plot_af_all_pV_vie_ita_sa <- function(region,coords,lpVals,ylab_pv=expression(-log[10](italic(P))) ) {
    region_idx = which(eu_sa_freqs$CHR == coords[1] & (coords[2] <= eu_sa_freqs$BPS & eu_sa_freqs$BPS <= coords[3] ) & (eu_sa_freqs$lP_ita >= lpVals[1] |  eu_sa_freqs$lP_vi >= lpVals[1] | eu_sa_freqs$lP_sa >= lpVals[2]))
    length(region_idx)
    layout(matrix(1:2, ncol = 1), widths = 1, heights = c(0.75,1.25), respect = FALSE)
    par(mar = c(0, 4.1, 4.1, 2.1))
    main_tit=bquote("SNPs around" ~.(region) ~"with" ~ Pvie ~ and ~ Pbolz <= 10^-.(lpVals[1]) ~ " or " ~  Psa <= 10 ^-.(lpVals[2]))
    xlimit=c(1-0.25,length(region_idx)+0.15)
    plot(1:length(region_idx),eu_sa_freqs$lP_sa[region_idx],col="blue",cex=1,pch=16,ylim=c(0,max(eu_sa_freqs[region_idx,c("lP_sa","lP_vi","lP_ita")],na.rm=T)),ylab=ylab_pv, xaxt="n",main=main_tit,xlim=xlimit )
    tot_dist=0.35
    plot_dist=0.25
    abline(v=(1:length(region_idx))-tot_dist,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+tot_dist,col="grey",lty="dotted",lw=2)
    points(1:length(region_idx),eu_sa_freqs$lP_vi[region_idx],col="orange",cex=1,pch=16)
    points(1:length(region_idx),eu_sa_freqs$lP_ita[region_idx],col="green",cex=1,pch=16)
    legend("topright",legend=c("SA","Vienna","Bolzano"),pch=c(16,20),col=c("blue","orange","green"))
    par(mar = c(4.1, 4.1, 0, 2.1))
    ymin=min(eu_sa_freqs[region_idx,c("SVL","SL","SD","SVD","VL","VD","IL","ID","Vb","Ib","Sb")],na.rm=T)
    plot((1:length(region_idx))-plot_dist,eu_sa_freqs[region_idx,c("Sb")],col="blue",cex=1,pch=".",ylim=c(ymin,1),ylab="AF",xlab="SNPs",xlim=xlimit )
    abline(v=(1:length(region_idx))-tot_dist,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+tot_dist,col="grey",lty="dotted",lw=2)
    abline(h=seq(0,1,by=plot_dist),col="grey",lty="dotted",lw=2)
    points((1:length(region_idx))-plot_dist,eu_sa_freqs[region_idx,c("Sb")],col="blue",cex=1,pch=16)
    points((1:length(region_idx))-plot_dist,eu_sa_freqs[region_idx,c("SVL")],col="blue",cex=1,pch=6)
    points((1:length(region_idx))-plot_dist,eu_sa_freqs[region_idx,c("SVD")],col="blue",cex=1,pch=17)
     points((1:length(region_idx))-plot_dist,eu_sa_freqs[region_idx,c("SL")],col="blue",cex=0.75,pch=23)
    points((1:length(region_idx))-plot_dist,eu_sa_freqs[region_idx,c("SD")],col="blue",cex=0.75,pch=23,bg="cyan")

    points((1:length(region_idx))+0,eu_sa_freqs[region_idx,c("Vb")],col="orange",cex=1,pch=16)
    points((1:length(region_idx))+0,eu_sa_freqs[region_idx,c("VL")],col="orange",cex=1,pch=6)
    points((1:length(region_idx))+0,eu_sa_freqs[region_idx,c("VD")],col="orange",cex=1,pch=17)
    
    points((1:length(region_idx))+plot_dist,eu_sa_freqs[region_idx,c("Ib")],col="green",cex=1,pch=16)
    points((1:length(region_idx))+plot_dist,eu_sa_freqs[region_idx,c("IL")],col="green",cex=1,pch=6)
    points((1:length(region_idx))+plot_dist,eu_sa_freqs[region_idx,c("ID")],col="green",cex=1,pch=17)
    legend("bottomright",legend=c("very light","light","base","dark","very dark"),title="mean AF",pch=c(6,23,16,23,17),col=c("black"),pt.bg=c("white","white","white","grey16","black"),pt.cex =c(1,0.75,1,0.75,1))
}

tan=c("X",8867539,9170829)
tan_p=c(14,12)
quartz(width=18,height=8)
plot_af_all_pV_vie_ita_sa("tan",tan,tan_p)
dev.copy2pdf(file="vie_ita_sa_freqs_comp_tan.pdf")
bab1_p=c(10,22)
quartz(width=18,height=8)
plot_af_all_pV_vie_ita_sa("bab1",bab1,bab1_p)
dev.copy2pdf(file="vie_ita_sa_freq_comp_bab1_lpeu10_psa22.pdf")
ebony=c("3R",17055975,17069171)
ebony_p=c(6,6)
quartz(width=18,height=8)
plot_af_all_pV_vie_ita_sa("ebony",ebony,ebony_p)
dev.copy2pdf(file="vie_ita_sa_freq_comp_ebony_lpeu6_lpsa6.pdf")

tan=c("X",8867539,9170829)
tan_p=c(12,12)
quartz(width=18,height=8)
plot_af_all_pV("tan",tan,tan_p)
dev.copy2pdf(file="eu_sa_all_freqs_comp_tan.pdf")
tan_p=c(20,12)
quartz(width=18,height=12)
plot_af_all_pV("tan",tan,tan_p)
tan=c("X",8867539,9170829)
tan_p=c(12,12)
quartz(width=18,height=8)
plot_af_all_pV("tan",tan,tan_p)
dev.copy2pdf(file="eu_sa_all_freqs_comp_tan.pdf")
tan_p=c(20,12)
quartz(width=18,height=12)
plot_af_all_pV("tan",tan,tan_p)
dev.copy2pdf(file="eu_sa_all_freqs_comp_tan_p20_12.pdf")

tan=c("X",8867539,9170829)
tan_p=c(12,12)
quartz(width=18,height=8)
plot_af_all_pV_vie_ita_sa("tan",tan,tan_p)
dev.copy2pdf(file="vie_ita_sa_all_freqs_comp_tan.pdf")
tan_p=c(20,12)
quartz(width=18,height=12)
plot_af_all_pV_vie_ita_sa("tan",tan,tan_p)
tan=c("X",8867539,9170829)
tan_p=c(12,12)
quartz(width=18,height=8)
plot_af_all_pV_vie_ita_sa("tan",tan,tan_p)
dev.copy2pdf(file="vie_ita_sa_all_freqs_comp_tan.pdf")
tan_p=c(20,12)
quartz(width=18,height=12)
plot_af_all_pV_vie_ita_sa("tan",tan,tan_p)
dev.copy2pdf(file="vie_ita_sa_all_freqs_comp_tan_p20_12.pdf")


bab1=c("3L",1071314,1120551 )
bab1_p=c(10,16)
quartz(width=18,height=8)
plot_af_all_pV("bab1",bab1,bab1_p)
dev.copy2pdf(file="eu_sa_all_freq_comp_bab1_lpeu10_psa16.pdf")
bab1_p=c(12,22)
quartz(width=18,height=12)
plot_af_all_pV("bab1",bab1,bab1_p)
dev.copy2pdf(file="eu_sa_all_freq_comp_bab1_lpeu12_psa22.pdf")

ebony=c("3R",17055975,17069171)
ebony_p=c(6,6)
quartz(width=18,height=8)
plot_af_all_pV("ebony",ebony,ebony_p)
dev.copy2pdf(file="eu_sa_all_freq_comp_ebony_lpeu6_lpsa6.pdf")

bab2=c("3L",1128600,1181056)
bab2_p=c(8,10)
quartz(width=18,height=8)
plot_af_all_pV("bab2",bab2,bab2_p)
dev.copy2pdf(file="eu_sa_all_freq_comp_bab2_lpeu8_lpsa10.pdf")
bab2_p=c(6,6)

quartz(width=18,height=8)
plot_af_all_pV("bab2",bab2,bab2_p)
dev.copy2pdf(file="eu_sa_all_freq_comp_bab2_lpeu6_lpsa6.pdf")


pdm3=c("2R",4179149,4319221)
pdm3_p=c(6,6)
quartz(width=18,height=8)
plot_af_all_pV("pdm3",pdm3,pdm3_p)
dev.copy2pdf(file="eu_sa_all_freq_comp_pdm3_lpeu6_lpsa6.pdf")

pdm3_p=c(6,8)
quartz(width=18,height=8)
plot_af_all_pV("pdm3",pdm3,pdm3_p)
dev.copy2pdf(file="eu_sa_all_freq_comp_pdm3_lpeu6_lpsa8.pdf")


bab1=c("3L",1071314,1120551 )
bab1_p=c(10,16)
quartz(width=18,height=8)
plot_af_all_pV_vie_ita_sa("bab1",bab1,bab1_p)
dev.copy2pdf(file="vie_ita_sa_all_freq_comp_bab1_lpeu10_psa16.pdf")
bab1_p=c(12,22)
quartz(width=18,height=12)
plot_af_all_pV_vie_ita_sa("bab1",bab1,bab1_p)
dev.copy2pdf(file="vie_ita_sa_all_freq_comp_bab1_lpeu12_psa22.pdf")

ebony=c("3R",17055975,17069171)
ebony_p=c(6,6)
quartz(width=18,height=8)
plot_af_all_pV_vie_ita_sa("ebony",ebony,ebony_p)
dev.copy2pdf(file="vie_ita_sa_all_freq_comp_ebony_lpeu6_lpsa6.pdf")

bab2=c("3L",1128600,1181056)
bab2_p=c(8,10)
quartz(width=18,height=8)
plot_af_all_pV_vie_ita_sa("bab2",bab2,bab2_p)
dev.copy2pdf(file="vie_ita_sa_all_freq_comp_bab2_lpeu8_lpsa10.pdf")
bab2_p=c(6,6)

quartz(width=18,height=8)
plot_af_all_pV_vie_ita_sa("bab2",bab2,bab2_p)
dev.copy2pdf(file="vie_ita_sa_all_freq_comp_bab2_lpeu6_lpsa6.pdf")


pdm3=c("2R",4179149,4319221)
pdm3_p=c(6,6)
quartz(width=18,height=8)
plot_af_all_pV_vie_ita_sa("pdm3",pdm3,pdm3_p)
dev.copy2pdf(file="vie_ita_sa_all_freq_comp_pdm3_lpeu6_lpsa6.pdf")

pdm3_p=c(6,8)
quartz(width=18,height=8)
plot_af_all_pV_vie_ita_sa("pdm3",pdm3,pdm3_p)
dev.copy2pdf(file="vie_ita_sa_all_freq_comp_pdm3_lpeu6_lpsa8.pdf")


# lOR Plots
tan=c("X",8867539,9170829)
tan_p=c(22,12)
quartz(width=18,height=8)
plot_OR_xs_pV("tan",tan,tan_p)
dev.copy2pdf(file="tan_OR_eu_sa_xs_p22_12.pdf")
bab1=c("3L",1071314,1120551 )
bab1_p=c(10,24)
plot_OR_xs_pV("bab1",bab1,bab1_p)
dev.copy2pdf(file="bab_OR_eu_sa_xs_p10_24.pdf")
ebony=c("3R",17055975,17069171)
ebony_p=c(6,6)
plot_OR_xs_pV("ebony",ebony,ebony_p)
dev.copy2pdf(file="tan_OR_eu_sa_xs_p6_6.pdf")


plot_OR_xs_pV <- function(region,coords,lpVals,ylab_pv=expression(-log[10](italic(P))) ){
    region_idx = which(vie_ita_sa_odds$CHR == coords[1] & (coords[2] <= vie_ita_sa_odds$BPS & vie_ita_sa_odds$BPS <= coords[3] ) & (vie_ita_sa_odds$P_eu >= lpVals[1] | vie_ita_sa_odds$P_sa >= lpVals[2] | vie_ita_sa_odds$P_xs >= lpVals[2]))
    length(region_idx)
    layout(matrix(1:2, ncol = 1), widths = 1, heights = c(0.75,1.25), respect = FALSE)
    par(mar = c(0, 4.1, 4.1, 2.1))
    main_tit=bquote("SNPs around" ~.(region) ~"with either" ~ Peu <= 10^-.(lpVals[1]) ~ "or" ~  Psa <= 10 ^-.(lpVals[2]))
    xlimit=c(1-0.25,length(region_idx)+0.15)  
    plot((1:length(region_idx))+0.2,vie_ita_sa_odds$P_sa[region_idx],col="blue",cex=1,pch=16,ylim=c(0,max(vie_ita_sa_odds[region_idx,c("P_sa","P_eu","P_xs")])),ylab=ylab_pv, xaxt="n",main=main_tit,xlim=xlimit )
    abline(v=(1:length(region_idx))-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+0.25,col="grey",lty="dotted",lw=2)
    points((1:length(region_idx))-0.2,vie_ita_sa_odds$P_eu[region_idx],col="red",cex=1,pch=16)
    points(1:length(region_idx),vie_ita_sa_odds$P_xs[region_idx],col="green",cex=1,pch=16)
    legend("topright",legend=c("Europe","SA xt","SA"),pch=c(16,20),col=c("red","green","blue"))
    par(mar = c(4.5, 4.1, 0, 2.1))
    ors=c("lO_eu","lO_xs","lO_sa")
    vals=unlist(vie_ita_sa_odds[region_idx,ors])
    ymin=min(vals[is.finite(vals)],na.rm=T)
    ymax=max(vals[is.finite(vals)],na.rm=T)
    plot((1:length(region_idx))-0.2,vie_ita_sa_odds[region_idx,"lO_eu"],col="red",cex=1,pch=16,ylim=c(ymin,ymax),ylab="log(OR)",xlab="",xlim=xlimit,xaxt="n")
    abline(h=0,col="darkgrey",lty="dashed",lw=2)
    abline(v=(1:length(region_idx))-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+0.25,col="grey",lty="dotted",lw=2)
    abline(h=seq(-4,4,by=1),col="grey",lty="dotted",lw=2)
    #abline(h=0,col="grey",lty=0,lw=2)
    points((1:length(region_idx))-0.2,vie_ita_sa_odds[region_idx,"lO_eu"],col="red",cex=1,pch=16)
    points((1:length(region_idx))+0.2,vie_ita_sa_odds[region_idx,"lO_sa"],col="blue",cex=1,pch=16)
    points((1:length(region_idx))+0.0,vie_ita_sa_odds[region_idx,"lO_xs"],col="green",cex=1,pch=16)
    # axes
    axis(1, at=1:length(region_idx), labels = FALSE)
    axis(2,at=seq(-4,4,by=1))
    # plot haplotypes
    text(x = 1:length(region_idx), par("usr")[3] - (0.35)*(ymax-ymin)/6, labels = vie_ita_sa_odds$BPS[region_idx], srt = 35, pos = 1, xpd = TRUE,cex=0.75)
    mtext(coords[1],side=1,line=3, xpd = TRUE,font=2)
    return(region_idx)

}

plot_OR_pV <- function(region,coords,lpVals,ylab_pv=expression(-log[10](italic(P))) ){
    region_idx = which(eu_sa_odds$CHR == coords[1] & (coords[2] <= eu_sa_odds$BPS & eu_sa_odds$BPS <= coords[3] ) & (eu_sa_odds$lP_eu >= lpVals[1] | eu_sa_odds$lP_sa >= lpVals[2]))
    length(region_idx)
    layout(matrix(1:2, ncol = 1), widths = 1, heights = c(0.75,1.25), respect = FALSE)
    par(mar = c(0, 4.1, 4.1, 2.1))
    main_tit=bquote("SNPs around" ~.(region) ~"with" ~ Peu <= 10^-.(lpVals[1]) ~ " or " ~  Psa <= 10 ^-.(lpVals[2]))
    xlimit=c(1-0.25,length(region_idx)+0.15)  
    plot(1:length(region_idx),eu_sa_odds$lP_sa[region_idx],col="blue",cex=1,pch=16,ylim=c(0,max(eu_sa_odds[region_idx,c("lP_sa","lP_eu")])),ylab=ylab_pv, xaxt="n",main=main_tit,xlim=xlimit )
    abline(v=(1:length(region_idx))-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+0.25,col="grey",lty="dotted",lw=2)
    points(1:length(region_idx),eu_sa_odds$lP_eu[region_idx],col="red",cex=1,pch=16)
    legend("topright",legend=c("SA","Europe"),pch=c(16,20),col=c("blue","red"))
    par(mar = c(4.1, 4.1, 0, 2.1))
    ors=c("lO_eu","lO_sa")
    vals=unlist(eu_sa_odds[region_idx,ors])
    ymin=min(vals[is.finite(vals)],na.rm=T)
    ymax=max(vals[is.finite(vals)],na.rm=T)
    plot((1:length(region_idx))-0.15,eu_sa_odds[region_idx,"lO_sa"],col="blue",cex=1,pch=16,ylim=c(ymin,ymax),ylab="log(OR)",xlab="SNPs",xlim=xlimit )
    abline(v=(1:length(region_idx))-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+0.25,col="grey",lty="dotted",lw=2)
    abline(h=seq(-4,4,by=1),col="grey",lty="dotted",lw=2)
    points((1:length(region_idx))-0.15,eu_sa_odds[region_idx,"lO_sa"],col="blue",cex=1,pch=16)
    ## points((1:length(region_idx))-0.15,eu_sa_odds[region_idx,"O_svd"],col="blue",cex=1,pch=17)
    ##  points((1:length(region_idx))-0.15,eu_sa_odds[region_idx,"O_sl"],col="blue",cex=0.75,pch=23)
    ## points((1:length(region_idx))-0.15,eu_sa_odds[region_idx,"O_sd"],col="blue",cex=0.75,pch=23,bg="cyan")

    points((1:length(region_idx))+0.15,eu_sa_odds[region_idx,"lO_eu"],col="red",cex=1,pch=16)
    ## points((1:length(region_idx))+0.15,eu_sa_odds[region_idx,"O_el"],col="red",cex=1,pch=6)
    ## points((1:length(region_idx))+0.15,eu_sa_odds[region_idx,"O_ed"],col="red",cex=1,pch=17)
    #legend("bottomright",legend=c("very light","light","dark","very dark"),title="log(OR) to base",pch=c(6,23,23,17),col=c("black"),pt.bg=c("white","white","grey16","black"),pt.cex =c(1,0.75,0.75,1))
}

plot_OR_pV_haps <- function(region,coords,lpVals,ylab_pv=expression(-log[10](italic(P))),df=eu_sa_odds_base,haplos=c("hap_eu","hap_sa"),pVs=c("lP_eu","lP_sa"),ors=c("lO_eu","lO_sa"),pVname=c("P_eu","P_sa"),ors.ci=F,ycorr=0.35){
    extended=F
    if(length(pVs)==3){
        extended=T
    }
    if (extended) {
      region_idx = which(df$CHR == coords[1] & (coords[2] <= df$BPS & df$BPS <= coords[3] ) & 
                           (df[,pVs[1]] >= lpVals[1] | df[,pVs[2]] >= lpVals[2] |  df[,pVs[3]] >= lpVals[2] ) )
    }
    else {
      region_idx = which(df$CHR == coords[1] & (coords[2] <= df$BPS & df$BPS <= coords[3] ) & (df[,pVs[1]] >= lpVals[1] | df[,pVs[2]] >= lpVals[2] ) ) 
    }
    
    #length(region_idx)

    df_region=df[region_idx,]
                                        #a=table(df_region[,c("hap_eu")])
    df_region$xpos=seq(1,length(df_region[,1]))
                                        #df_region$hap1=0
                                        #df_region$hap2=0
    start_hap1=c()
    stop_hap1=c()
    start_hap2=c()
    stop_hap2=c()
    hap1=F
    hap2=F
    for(i in seq(1,length(df_region[,1])-1)){
        if(df_region[i,haplos[1]] == df_region[i+1,haplos[1]]){
            if(! hap1){
                start_hap1=c(start_hap1,i)
                hap1=T
            }
        }
        else if (hap1 ){
            stop_hap1=c(stop_hap1,i)
            hap1=F
        }
        if(df_region[i,haplos[2]] == df_region[i+1,haplos[2]]){
            if(! hap2){
                start_hap2=c(start_hap2,i)
                hap2=T
            }
        }
        else if (hap2){
            stop_hap2=c(stop_hap2,i)
            hap2=F    
        }
        if(df_region[i,haplos[2]] == df_region[i+1,haplos[2]] | df_region[i,haplos[1]] == df_region[i+1,haplos[1]]){
            df_region$xpos[i+1]=df_region$xpos[i]+0.5
        }
        else{
            df_region$xpos[i+1]=df_region$xpos[i]+1
        }
    }
    i=length(length(df_region[,1]))
    if(hap1){
        stop_hap1=c(stop_hap1,i)
    }
    if(hap2){
        stop_hap2=c(stop_hap2,i)
    }    
    main_tit=bquote("SNPs around" ~italic(.(region)) ~"with" ~ P[.(pVname[1])] <= 10^-.(lpVals[1]) ~ " or " ~  P[.(pVname[2])] <= 10 ^-.(lpVals[2]))
    xlimit=c(df_region$xpos[1]-0.25,df_region$xpos[length(df_region$xpos)]+0.15)
    layout(matrix(1:2, ncol = 1), widths = 1, heights = c(0.75,1.25), respect = FALSE)
    par(mar = c(0, 4.1, 4.1, 2.1))
    plot(df_region$xpos,df_region[,pVs[2]],col="white",cex=1,pch=NULL,ylim=c(0,max(df_region[,pVs])),ylab=ylab_pv, xaxt="n",main=main_tit,xlim=xlimit )
    print(start_hap1)
    print(start_hap2)
    print(stop_hap1)
    print(stop_hap2)
    abline(v=(df_region$xpos)-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(df_region$xpos)+0.25,col="grey",lty="dotted",lw=2)
    # plot haplotypes
    if(length(start_hap1) > 0){
        plot_mean(df_region,start_hap1,stop_hap1,var=pVs[1],col="red")
    }
    if(length(start_hap2) > 0){
        plot_mean(df_region,start_hap2,stop_hap2,var=pVs[2],col="blue")
        if (extended) {plot_mean(df_region,start_hap2,stop_hap2,var=pVs[3],col="green")}
    }
    # plot pV
    
    if (extended) {
        points(df_region$xpos-0.15,df_region[,pVs[2]],col="blue",cex=1,pch=16)
        points(df_region$xpos,df_region[,pVs[3]],col="green",cex=1,pch=16)
        points(df_region$xpos+0.15,df_region[,pVs[1]],col="red",cex=1,pch=16)
    }
    else{
        points(df_region$xpos,df_region[,pVs[2]],col="blue",cex=1,pch=16)
        points(df_region$xpos,df_region[,pVs[1]],col="red",cex=1,pch=16)
    }
    if (! extended) {
        legend("topright",legend=c("SA","Europe"),pch=c(16,16),col=c("blue","red"))
    }
    else{
        legend("topright",legend=c("SA","SA ext.","Europe"),pch=c(16,16),col=c("blue","green","red"))
    }
    par(mar = c(5.1, 4.1, 0, 2.1))
    if (length(ors.ci) > 1){
        vals=c(unlist(df_region[,c(ors,unlist(ors.ci))]))
    }
    else{
        vals=c(unlist(df_region[,ors])) 
    }
    ymin=min(vals[is.finite(vals)],na.rm=T)
    ymax=max(vals[is.finite(vals)],na.rm=T)    
    plot((df_region$xpos)-0.15,df_region[,ors[2]],col="blue",cex=1,pch=16,ylim=c(ymin,ymax+0.15),ylab="log(OR)",yaxt="n",xaxt="n",xlab="",xlim=xlimit )
    # plot errorbars
    if (length(ors.ci) > 1){
        for (i in 1:length(df_region$xpos)) {
            #draw sa
            lines(c(df_region$xpos[i]+0.15,df_region$xpos[i]+0.15),df_region[i,ors.ci[[1]]],lty="dotted",col="red",lw=2)
            lines(c(df_region$xpos[i]-0.15,df_region$xpos[i]-0.15),df_region[i,ors.ci[[2]]],lty="dotted",col="blue",lw=2)
            points(c(df_region$xpos[i]+0.15,df_region$xpos[i]+0.15),df_region[i,ors.ci[[1]]],pch=4,col="red")
            points(c(df_region$xpos[i]-0.15,df_region$xpos[i]-0.15),df_region[i,ors.ci[[2]]],pch=4,col="blue")
            if(extended){
                lines(c(df_region$xpos[i],df_region$xpos[i]),df_region[i,ors.ci[[3]]],lty="dotted",col="green",lw=2)
                points(c(df_region$xpos[i],df_region$xpos[i]),df_region[i,ors.ci[[3]]],pch=4,col="green")
            }
        }
        
    }

    # plot grid lines
    abline(v=(df_region$xpos)-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(df_region$xpos)+0.25,col="grey",lty="dotted",lw=2)
    abline(h=seq(-4,3,by=1),col="grey",lty="dotted",lw=2)
    abline(h=0,col="darkgrey",lty="dotted",lw=2)
    # axes
    axis(1, at=df_region$xpos, labels = FALSE)
    axis(2,at=seq(-4,3,by=1))
    # plot haplotypes
    if(length(start_hap1) > 0){
        plot_mean(df_region,start_hap1,stop_hap1,var=ors[1],col="red")
    }
    if(length(start_hap2) > 0){
        plot_mean(df_region,start_hap2,stop_hap2,var=ors[2],col="blue")
        if(extended){ plot_mean(df_region,start_hap2,stop_hap2,var=ors[3],col="green")}
    }
    text(x = df_region$xpos, par("usr")[3] - ycorr, labels = df_region$BPS, srt = 45, pos = 1, xpd = TRUE,cex=0.75)
    mtext(coords[1],side=1,line=3, xpd = TRUE,font=2)

    points((df_region$xpos)-0.15,df_region[,ors[2]],col="blue",cex=1,pch=16)
    ## points((df_region$xpos)-0.15,df_region[,"O_svd"],col="blue",cex=1,pch=17)
    ##  points((df_region$xpos)-0.15,df_region[,"O_sl"],col="blue",cex=0.75,pch=23)
    ## points((df_region$xpos)-0.15,df_region[,"O_sd"],col="blue",cex=0.75,pch=23,bg="cyan")
    
    points((df_region$xpos)+0.15,df_region[,ors[1]],col="red",cex=1,pch=16)
    if(extended){points((df_region$xpos),df_region[,ors[3]],col="green",cex=1,pch=16)}
    ## points((df_region$xpos)+0.15,df_region[,"O_el"],col="red",cex=1,pch=6)
    ## points((df_region$xpos)+0.15,df_region [,"O_ed"],col="red",cex=1,pch=17)
    #legend("bottomright",legend=c("very light","light","dark","very dark"),title="log(OR) to base",pch=c(6,23,23,17),col=c("black"),pt.bg=c("white","white","grey16","black"),pt.cex =c(1,0.75,0.75,1))
    print(ymax)
    return(df_region[,c("BPS","xpos")])

}

plot_mean = function(df_region,start_hap,stop_hap,var="lO_eu",col="red"){
    means= sapply(1:length(start_hap),function (x) mean(df_region[start_hap[x]:stop_hap[x],var],na.rm=T))
    segments(df_region$xpos[start_hap]-0.25,means,df_region$xpos[stop_hap]+0.25,means,col=col,lty=3,lw=2.5)
    print(means)
}

tan=c("X",8867539,9170829)
tan_p=c(16,12)
quartz(width=18,height=8)
x11(width=18,height=8)
#plot_OR_pV_haps("tan",tan,tan_p,df=eu_sa_odds_base,haplos=c("hap_eu","hap_sa"),pVs=c("lP_eu","lP_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"))
plot_OR_pV_haps("tan",tan,tan_p,df=eu.sa.odds,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"),ors.ci=list(c("lOa_eu","lOb_eu"),c("lOa_sa","lOb_sa")))

dev.copy2pdf(file="eu_sa_OR_haps_tan_16_12.pdf")

plot_OR_pV_haps("tan",tan,tan_p,df=eu.sa.xs.odds,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa","P_xs"),ors=c("lO_eu","lO_sa","lO_xs"),pVname=c("eu","sa"),ors.ci=list(c("lOa_eu","lOb_eu"),c("lOa_sa","lOb_sa"),c("lOa_xs","lOb_xs")))
dev.copy2pdf(file="eu_sa_xs_OR_haps_tan_16_12.pdf")

bab1=c("3L",1071314,1120551 )
bab1_p=c(10,16)
quartz(width=18,height=8)
#plot_OR_pV_haps("bab1",bab1,bab1_p,df=eu_sa_odds_base,haplos=c("hap_eu","hap_sa"),pVs=c("lP_eu","lP_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"))
plot_OR_pV_haps("bab1",bab1,bab1_p,df=eu.sa.odds,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa"),ors=c("O_eu","O_sa"),pVname=c("eu","sa"),ors.ci=list(c("Oa_eu","Ob_eu"),c("Oa_sa","Ob_sa")))
dev.copy2pdf(file="eu_sa_OR_haps_bab1_lpeu10_psa16.pdf")
bab1_p=c(10,18)
quartz(width=18,height=8)
#plot_OR_pV_haps("bab1",bab1,bab1_p,df=eu_sa_odds_base,haplos=c("hap_eu","hap_sa"),pVs=c("lP_eu","lP_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"))
plot_OR_pV_haps("bab1",bab1,bab1_p,df=eu.sa.odds,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa"),ors=c("O_eu","O_sa"),pVname=c("eu","sa"))
dev.copy2pdf(file="eu_sa_OR_haps_bab1_lpeu10_psa18.pdf")
plot_OR_pV_haps("bab1",bab1,bab1_p,df=eu.sa.xs.odds,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa","P_xs"),ors=c("lO_eu","lO_sa","lO_xs"),pVname=c("eu","sa"),ors.ci=list(c("lOa_eu","lOb_eu"),c("lOa_sa","lOb_sa"),c("lOa_xs","lOb_xs")))
dev.copy2pdf(file="eu_sa_xs_OR_haps_bab1_lpeu10_psa18.pdf")
bab1_p=c(12,22)
quartz(width=18,height=8)
#plot_OR_pV_haps("bab1",bab1,bab1_p,df=eu_sa_odds_base,haplos=c("hap_eu","hap_sa"),pVs=c("lP_eu","lP_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"))
plot_OR_pV_haps("bab1",bab1,bab1_p,df=eu.sa.odds,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa"),ors=c("O_eu","O_sa"),pVname=c("eu","sa"))
dev.copy2pdf(file="eu_sa_OR_haps_bab1_lpeu12_psa22.pdf")
bab1_p=c(14,24)
quartz(width=18,height=8)
#plot_OR_pV_haps("bab1",bab1,bab1_p,df=eu_sa_odds_base,haplos=c("hap_eu","hap_sa"),pVs=c("lP_eu","lP_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"))
plot_OR_pV_haps("bab1",bab1,bab1_p,df=eu.sa.odds,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa"),ors=c("O_eu","O_sa"),pVname=c("eu","sa"))
dev.copy2pdf(file="eu_sa_OR_haps_bab1_lpeu14_psa24.pdf")

ebony=c("3R",17055975,17069171)
ebony_p=c(6,6)
quartz(width=18,height=8)
plot_OR_pV_haps("ebony",ebony,ebony_p,df=eu_sa_odds_base,haplos=c("hap_eu","hap_sa"),pVs=c("lP_eu","lP_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"))
dev.copy2pdf(file="eu_sa_OR_haps_ebony_lpeu6_lpsa6.pdf")

bab2=c("3L",1128600,1181056)
bab2_p=c(8,8)
quartz(width=18,height=8)
plot_OR_pV_haps("bab2",bab2,bab2_p,df=eu_sa_odds_base,haplos=c("hap_eu","hap_sa"),pVs=c("lP_eu","lP_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"))
dev.copy2pdf(file="eu_sa_OR_haps_bab2_lpeu8_lpsa8.pdf")
bab2_p=c(6,6)
quartz(width=18,height=8)
plot_OR_pV_haps("bab2",bab2,bab2_p,df=eu_sa_odds_base,haplos=c("hap_eu","hap_sa"),pVs=c("lP_eu","lP_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"))
dev.copy2pdf(file="eu_sa_OR_haps_bab2_lpeu6_lpsa6.pdf")


pdm3=c("2R",4179149,4319221)
pdm3_p=c(6,6)
quartz(width=18,height=8)
plot_OR_pV_haps("pdm3",pdm3,pdm3_p,df=eu_sa_odds_base,haplos=c("hap_eu","hap_sa"),pVs=c("lP_eu","lP_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"))
dev.copy2pdf(file="eu_sa_OR_haps_pdm3_lpeu6_lpsa6.pdf")

pdm3=c("2R",4179149,4319221)
pdm3_p=c(6,8)
quartz(width=18,height=8)
plot_OR_pV_haps("pdm3",pdm3,pdm3_p,df=eu_sa_odds_base,haplos=c("hap_eu","hap_sa"),pVs=c("lP_eu","lP_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"))
dev.copy2pdf(file="eu_sa_OR_haps_pdm3_lpeu6_lpsa8.pdf")

bab1_p=c(12,22)
quartz(width=18,height=8)
plot_AF_pV_haps("bab1",bab1,bab1_p,df=eu.sa.freqs.p6,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa"),pVname=c("eu","sa"))
dev.copy2pdf(file="eu_sa_AFS_haps_bab1_lpeu12_psa22.pdf")


ebony=c("3R",17055975,17069171)
ebony=c("3R",17055561,17075171)
ebony_p=c(6,6)
#quartz(width=18,height=8)
x11(width=18,height=8)
a <- plot_OR_pV_haps("ebony",ebony,ebony_p,df=eu.sa.odds,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"),ors.ci=list(c("lOa_eu","lOb_eu"),c("lOa_sa","lOb_sa")),ycorr=0.125)
dev.copy2pdf(file="eu_sa_OR_haps_ebony_lpeu6_lpsa6.pdf")
quartz(width=18,height=8)
a <- plot_OR_pV_haps("ebony",ebony,ebony_p,df=eu.sa.odds,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"),ors.ci=list(c("lOa_eu","lOb_eu"),c("lOa_sa","lOb_sa")),ycorr=0.125)

dev.copy2pdf(file="eu_sa_OR_haps_ebony_lpeu6_lpsa6.pdf")

x11(width=18,height=8)
quartz(width=18,height=8)
a <- plot_AF_pV_haps("ebony",ebony,ebony_p,df=eu.sa.freqs.p6,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa"),afs=list(c("EUL","EUb","EUD"),c("SVL","Sb","SVD")),pVname=c("eu","sa"),ycorr=0.04,llpos="topright")
dev.copy2pdf(file="eu_sa_AF_haps_ebony_lpeu6_lpsa6.pdf")


bab1=c("3L",1071314,1120551 )
bab1_p=c(12,22)
quartz(width=18,height=8)
#plot_OR_pV_haps("bab1",bab1,bab1_p,df=eu_sa_odds_base,haplos=c("hap_eu","hap_sa"),pVs=c("lP_eu","lP_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"))
a <- plot_OR_pV_haps("bab1",bab1,bab1_p,df=eu.sa.xs.odds,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa","P_xs"),ors=c("lO_eu","lO_sa","lO_xs"),pVname=c("eu","sa"),ors.ci=list(c("lOa_eu","lOb_eu"),c("lOa_sa","lOb_sa"),c("lOa_xs","lOb_xs")))
a <- plot_OR_pV_haps("bab1",bab1,bab1_p,df=eu.sa.odds,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"),ors.ci=list(c("lOa_eu","lOb_eu"),c("lOa_sa","lOb_sa")))
ymax=3.55
bab1.x=get.xpos.from.coords(a,c(1099101,1101089))
dme.x=get.xpos.from.coords(a,c(1084814,1085477))
bab1.x=c(14.75,15.25)
lines(bab1i.x,c(ymax+0.1,ymax+0.1))
bab1i.x=get.xpos.from.coords(a,c(1040355,1099100))
rect(bab1.x[1],ymax,bab1.x[2],ymax+0.2,col="grey")
rect(dme.x[1],ymax,dme.x[2],ymax+0.2,col="white")
text(mean(bab1.x),ymax+0.085,labels="bab1")
text(mean(dme.x),ymax+0.085,labels="DME")
dev.copy2pdf(file="eu_sa_OR_haps_bab1_lpeu12_psa22.pdf")
dev.copy2pdf(file="eu_sa_OR_haps_bab1_lpeu12_psa22_noxs.pdf")
quartz(width=18,height=8)
a <- plot_AF_pV_haps("bab1",bab1,bab1_p,df=eu.sa.freqs.p6,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa","P_xs"),afs=list(c("EUL","EUb","EUD"),c("SVL","Sb","SVD"),c("XSL","Sb","XSD")),pVname=c("eu","sa"))
a <- plot_AF_pV_haps("bab1",bab1,bab1_p,df=eu.sa.freqs.p6,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa"),afs=list(c("EUL","EUb","EUD"),c("SVL","Sb","SVD")),pVname=c("eu","sa"))
ymax=1.0
bab1.x=get.xpos.from.coords(a,c(1099101,1101089))
dme.x=get.xpos.from.coords(a,c(1084814,1085477))
bab1.x=c(14.75,15.25)
bab1i.x=get.xpos.from.coords(a,c(1040355,1099100))
lines(bab1i.x,c(ymax+0.025,ymax+0.025))
rect(bab1.x[1],ymax,bab1.x[2],ymax+0.05,col="grey")
rect(dme.x[1],ymax,dme.x[2],ymax+0.05,col="white")
text(mean(bab1.x),ymax+0.025,labels="bab1")
text(mean(dme.x),ymax+0.025,labels="DME")

dev.copy2pdf(file="eu_sa_AFS_haps_bab1_lpeu12_psa22_noxs.pdf")




#tan=c("X",8867539,9170829)
tan=c("X",9111469,9127367)
tan_p=c(16,12)
quartz(width=18,height=8)
#plot_OR_pV_haps("tan",tan,tan_p,df=eu_sa_odds_base,haplos=c("hap_eu","hap_sa"),pVs=c("lP_eu","lP_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"))
#a=plot_OR_pV_haps("tan",tan,tan_p,df=eu.sa.xs.odds,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa","P_xs"),ors=c("lO_eu","lO_sa","lO_xs"),pVname=c("eu","sa"),ors.ci=list(c("lOa_eu","lOb_eu"),c("lOa_sa","lOb_sa"),c("lOa_xs","lOb_xs")))
a=plot_OR_pV_haps("tan",tan,tan_p,df=eu.sa.odds,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa"),ors=c("lO_eu","lO_sa"),pVname=c("eu","sa"),ors.ci=list(c("lOa_eu","lOb_eu"),c("lOa_sa","lOb_sa")))
ymax=3.65
tan.x=get.xpos.from.coords(a,c(9111688,9117290))
Gr8a.x=get.xpos.from.coords(a,c(9121170,9122860))
cg15370.x=get.xpos.from.coords(a,c(9119470,9120721))
CG12121.x=get.xpos.from.coords(a,c(9122820,9126406))
Ir8a.x=get.xpos.from.coords(a,c(9126759,9130618))
group1=get.xpos.from.coords(a,c(9119116,9119160))
group2=get.xpos.from.coords(a,c(9120683,9120730))
group3=get.xpos.from.coords(a,c(9123892,9123903))
rect(tan.x[1],ymax,tan.x[2],ymax+0.2,col="grey")
rect(Gr8a.x[1],ymax,Gr8a.x[2],ymax+0.2,col="grey")
rect(cg15370.x[1],ymax,cg15370.x[2],ymax+0.2,col="grey")
rect(CG12121.x[1],ymax,CG12121.x[2],ymax+0.2,col="grey")
rect(Ir8a.x[1],ymax,Ir8a.x[2],ymax+0.2,col="grey")
text(mean(Gr8a.x),ymax+0.085,labels="Gr8a")
text(mean(cg15370.x),ymax+0.085,labels="CG15370")
text(mean(Ir8a.x),ymax+0.085,labels="Ir8a")
text(mean(CG12121.x),ymax+0.085,labels="CG12121")
rect(Gr8a.x[1],ymax,cg15370.x[2],ymax+0.2,col="white")
text(mean(c(Gr8a.x[1],cg15370.x[2])),ymax+0.085,labels="MSE")
# add rectangle to special haplotypes
rect(group1[1],-0.5,group1[2],2.0,col=NA, border="darkgreen",lty="dashed",lwd=3)
rect(group2[1],-0.75,group2[2],1.5,col=NA, border="darkgreen",lty="dashed",lwd=3)
rect(group3[1],-0.75,group3[2],2.0,col=NA, border="darkgreen",lty="dashed",lwd=3)
dev.copy2pdf(file="eu_sa_OR_haps_tan_16_12_noxs_boxes.pdf")


quartz(width=18,height=8)
#a=plot_AF_pV_haps("tan",tan,tan_p,df=eu.sa.freqs.p6,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa","P_xs"),afs=list(c("EUL","EUb","EUD"),c("SVL","Sb","SVD"),c("XSL","Sb","XSD")),pVname=c("eu","sa"))
a=plot_AF_pV_haps("tan",tan,tan_p,df=eu.sa.freqs.p6,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa"),afs=list(c("EUL","EUb","EUD"),c("SVL","Sb","SVD")),pVname=c("eu","sa"))
tan.x=get.xpos.from.coords(a,c(9111688,9117290))
Gr8a.x=get.xpos.from.coords(a,c(9121170,9122860))
cg15370.x=get.xpos.from.coords(a,c(9119470,9120721))
CG12121.x=get.xpos.from.coords(a,c(9122820,9126406))
Ir8a.x=get.xpos.from.coords(a,c(9126759,9130618)) 
rect(tan.x[1],1.0,tan.x[2],1.05,col="grey")
rect(Gr8a.x[1],1.0,Gr8a.x[2],1.05,col="grey")
rect(cg15370.x[1],1.0,cg15370.x[2],1.05,col="grey")
rect(CG12121.x[1],1.0,CG12121.x[2],1.05,col="grey")
rect(Ir8a.x[1],1.0,Ir8a.x[2],1.05,col="grey")
text(mean(Gr8a.x),1.025,labels="Gr8a")
text(mean(cg15370.x),1.025,labels="CG15370")
text(mean(Ir8a.x),1.025,labels="Ir8a")
text(mean(CG12121.x),1.025,labels="CG12121")
rect(Gr8a.x[1],1.0,cg15370.x[2],1.05,col="white")
text(mean(c(Gr8a.x[1],cg15370.x[2])),1.025,labels="MSE")
dev.copy2pdf(file="eu_sa_AFs_haps_tan_16_12_noxs.pdf")

get.xpos.from.coords <- function(df,coords) {
    stretch <-  df$xpos[ coords[1] <= df$BPS & df$BPS <= coords[2]]
    return(c(min(stretch)-0.25,max(stretch)+0.25))
}

plot_AF_pV_haps <- function(region,coords,lpVals,ylab_pv=expression(-log[10](italic(P))),df=eu_sa_odds_base,haplos=c("hap_eu","hap_sa"),pVs=c("P_eu","P_sa"),afs=list(c("EUL","EUb","EUD"),c("SVL","SL","Sb","SD","SVD")),pVname=c("P_eu","P_sa"),ycorr=0.065,llpos="bottomright"){
    #length(region_idx)
    extended=F
    if(length(pVs)==3){
        extended=T
    }
    if (extended){
      region_idx = which(df$CHR == coords[1] & (coords[2] <= df$BPS & df$BPS <= coords[3] ) & (df[,pVs[1]] >= lpVals[1] | df[,pVs[2]] >= lpVals[2] | df[,pVs[3]] >= lpVals[2] ) )
    }
  else {
    region_idx = which(df$CHR == coords[1] & (coords[2] <= df$BPS & df$BPS <= coords[3] ) & 
                         (df[,pVs[1]] >= lpVals[1] | df[,pVs[2]] >= lpVals[2] ) )
    
  }
    
    df_region=df[region_idx,]
                                        #a=table(df_region[,c("hap_eu")])
    df_region$xpos=seq(1,length(df_region[,1]))
                                        #df_region$hap1=0
                                        #df_region$hap2=0
    start_hap1=c()
    stop_hap1=c()
    start_hap2=c()
    stop_hap2=c()
    hap1=F
    hap2=F
    for(i in seq(1,length(df_region[,1])-1)){
        if(df_region[i,haplos[1]] == df_region[i+1,haplos[1]]){
            if(! hap1){
                start_hap1=c(start_hap1,i)
                hap1=T
            }
        }
        else if (hap1 ){
            stop_hap1=c(stop_hap1,i)
            hap1=F
        }
        if(df_region[i,haplos[2]] == df_region[i+1,haplos[2]]){
            if(! hap2){
                start_hap2=c(start_hap2,i)
                hap2=T
            }
        }
        else if (hap2){
            stop_hap2=c(stop_hap2,i)
            hap2=F    
        }
        if(df_region[i,haplos[2]] == df_region[i+1,haplos[2]] | df_region[i,haplos[1]] == df_region[i+1,haplos[1]]){
            df_region$xpos[i+1]=df_region$xpos[i]+0.5
        }
        else{
            df_region$xpos[i+1]=df_region$xpos[i]+1
        }
    }
    i=length(length(df_region[,1]))
    if(hap1){
        stop_hap1=c(stop_hap1,i)
    }
    if(hap2){
        stop_hap2=c(stop_hap2,i)
    }    
    main_tit=bquote("SNPs around"~italic(.(region)) ~"with" ~ P[.(pVname[1])] <= 10^-.(lpVals[1]) ~ " or " ~  P[.(pVname[2])] <= 10 ^-.(lpVals[2]))
    xlimit=c(df_region$xpos[1]-0.25,df_region$xpos[length(df_region$xpos)]+0.15)
    layout(matrix(1:2, ncol = 1), widths = 1, heights = c(0.75,1.25), respect = FALSE)
    par(mar = c(0, 4.1, 4.1, 2.1))
    plot(df_region$xpos,df_region[,pVs[2]],col="white",cex=1,pch=NULL,ylim=c(0,max(df_region[,pVs])),ylab=ylab_pv, xaxt="n",main=main_tit,xlim=xlimit )
    print(start_hap1)
    print(start_hap2)
    print(stop_hap1)
    print(stop_hap2)
    abline(v=(df_region$xpos)-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(df_region$xpos)+0.25,col="grey",lty="dotted",lw=2)
    # plot haplotypes
    if(length(start_hap1) > 0){
        plot_mean(df_region,start_hap1,stop_hap1,var=pVs[1],col="red")
    }
    if(length(start_hap2) > 0){
        plot_mean(df_region,start_hap2,stop_hap2,var=pVs[2],col="blue")
        if (extended) {plot_mean(df_region,start_hap2,stop_hap2,var=pVs[3],col="green")}
    }
    # plot pV
    if (extended) {
        points(df_region$xpos-0.15,df_region[,pVs[2]],col="blue",cex=1,pch=16)
        points(df_region$xpos,df_region[,pVs[3]],col="green",cex=1,pch=16)
        points(df_region$xpos+0.15,df_region[,pVs[1]],col="red",cex=1,pch=16)
    }
    else{
        points(df_region$xpos,df_region[,pVs[2]],col="blue",cex=1,pch=16)
        points(df_region$xpos,df_region[,pVs[1]],col="red",cex=1,pch=16)
    }
    if (! extended) {
        legend("topright",legend=c("SA","Europe"),pch=c(16,16),col=c("blue","red"))
    }
    else{
        legend("topright",legend=c("SA","SA ext.","Europe"),pch=c(16,16),col=c("blue","green","red"))
    }
    par(mar = c(5.1, 4.1, 0, 2.1))
    # plot AFs
    ymin=min(df_region[,unlist(afs)],na.rm=T)
    plot((df_region$xpos)-0.15,df_region[,afs[[2]][3]],col="white",cex=1,pch=NULL,ylim=c(ymin,1.05),ylab="major AF",yaxt="n",xaxt="n",xlab="",xlim=xlimit )
    # plot grid lines
    abline(v=(df_region$xpos)-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(df_region$xpos)+0.25,col="grey",lty="dotted",lw=2)
    abline(h=seq(-4,3,by=1),col="grey",lty="dotted",lw=2)
    abline(h=0,col="darkgrey",lty="dotted",lw=2)
    # axes
    axis(1, at=df_region$xpos, labels = FALSE)
    axis(2,at=seq(0,1,by=0.2))
    # plot haplotypes
    if(length(start_hap1) > 0){
        plot_mean(df_region,start_hap1,stop_hap1,var=afs[[1]][2],col="red")
    }
    if(length(start_hap2) > 0){
        plot_mean(df_region,start_hap2,stop_hap2,var=afs[[2]][2],col="blue")
    }
    if (extended){
        points((df_region$xpos)-0.15,df_region[,afs[[2]][2]],col="blue",cex=1,pch=16)
        points((df_region$xpos)-0.15,df_region[,afs[[2]][1]],col="blue",cex=1,pch=6)
        points((df_region$xpos)-0.15,df_region[,afs[[2]][3]],col="blue",cex=1,pch=17)
        points((df_region$xpos),df_region[,afs[[3]][2]],col="green",cex=1,pch=16)
        points((df_region$xpos),df_region[,afs[[3]][1]],col="green",cex=1,pch=6)
        points((df_region$xpos),df_region[,afs[[3]][3]],col="green",cex=1,pch=17)
        
        points((df_region$xpos)+0.15,df_region[,afs[[1]][2]],col="red",cex=1,pch=16)
        points((df_region$xpos)+0.15,df_region[,afs[[1]][1]],col="red",cex=1,pch=6)
        points((df_region$xpos)+0.15,df_region [,afs[[1]][3]],col="red",cex=1,pch=17)
        legend(llpos,legend=c("light","base","dark"),title="average AF",pch=c(6,16,17),col=c("black"),pt.bg=c("white","black","black"),pt.cex =c(1,1,1))

    }
    else{
        points((df_region$xpos)-0.15,df_region[,afs[[2]][2]],col="blue",cex=1,pch=16)
        points((df_region$xpos)-0.15,df_region[,afs[[2]][1]],col="blue",cex=1,pch=6)
        points((df_region$xpos)-0.15,df_region[,afs[[2]][3]],col="blue",cex=1,pch=17)
        #points((df_region$xpos)-0.15,df_region[,afs[[2]][3]],col="blue",cex=1,pch=16)
#         points((df_region$xpos)-0.15,df_region[,afs[[2]][5]],col="blue",cex=1,pch=17)
#         points((df_region$xpos)-0.15,df_region[,afs[[2]][2]],col="blue",cex=0.75,pch=23)
#         points((df_region$xpos)-0.15,df_region[,afs[[2]][4]],col="blue",cex=0.75,pch=23,bg="cyan")
#         points((df_region$xpos)-0.15,df_region[,afs[[2]][1]],col="blue",cex=0.75,pch=6)
        
        points((df_region$xpos)+0.15,df_region[,afs[[1]][2]],col="red",cex=1,pch=16)
        points((df_region$xpos)+0.15,df_region[,afs[[1]][1]],col="red",cex=1,pch=6)
        points((df_region$xpos)+0.15,df_region [,afs[[1]][3]],col="red",cex=1,pch=17)
        #legend(llpos,legend=c("very light","light","dark","very dark"),title="average AF",pch=c(6,23,23,17),col=c("black"),pt.bg=c("white","white","grey16","black"),pt.cex =c(1,0.75,0.75,1))
        legend(llpos,legend=c("light","base","dark"),title="average AF",pch=c(6,16,17),col=c("black"),pt.bg=c("white","black","black"),pt.cex =c(1,1,1))
    }
    print( df_region$BPS)
    print( df_region$xpos)
    text(x = df_region$xpos, par("usr")[3]-ycorr, labels = df_region$BPS, srt = 45, pos = 1, xpd = TRUE,cex=0.85)
    mtext(coords[1],side=1,line=3, xpd = TRUE,font=2)
    return(df_region[,c("BPS","xpos")])

}


tan=c("X",8867539,9170829)
tan_p=c(20,12)
quartz(width=18,height=8)
plot_OR_pV("tan",tan,tan_p)
dev.copy2pdf(file="eu_sa_OR_comp_tan_20_12.pdf")

tan=c("X",8867539,9170829)
tan_p=c(12,12)
quartz(width=18,height=8)
plot_OR_pV("tan",tan,tan_p)
dev.copy2pdf(file="eu_sa_OR_comp_tan.pdf")

bab1=c("3L",1071314,1120551 )
bab1_p=c(10,16)
quartz(width=18,height=8)
plot_OR_pV("bab1",bab1,bab1_p)
dev.copy2pdf(file="eu_sa_OR_comp_bab1_lpeu10_psa16.pdf")
bab1_p=c(10,18)
quartz(width=18,height=8)
plot_OR_pV("bab1",bab1,bab1_p)
dev.copy2pdf(file="eu_sa_OR_comp_bab1_lpeu10_psa18.pdf")
bab1_p=c(12,22)
quartz(width=18,height=8)
plot_OR_pV("bab1",bab1,bab1_p)
dev.copy2pdf(file="eu_sa_OR_comp_bab1_lpeu12_psa22.pdf")

ebony=c("3R",17055975,17069171)
ebony_p=c(6,6)
quartz(width=18,height=8)
plot_OR_pV("ebony",ebony,ebony_p)
dev.copy2pdf(file="eu_sa_OR_comp_ebony_lpeu6_lpsa6.pdf")

bab2=c("3L",1128600,1181056)
bab2_p=c(8,8)
quartz(width=18,height=8)
plot_OR_pV("bab2",bab2,bab2_p)
dev.copy2pdf(file="eu_sa_OR_comp_bab2_lpeu8_lpsa8.pdf")
bab2_p=c(6,6)
quartz(width=18,height=8)
plot_OR_pV("bab2",bab2,bab2_p)
dev.copy2pdf(file="eu_sa_OR_comp_bab2_lpeu6_lpsa6.pdf")


pdm3=c("2R",4179149,4319221)
pdm3_p=c(6,6)
quartz(width=18,height=8)
plot_OR_pV("pdm3",pdm3,pdm3_p)
dev.copy2pdf(file="eu_sa_OR_comp_pdm3_lpeu6_lpsa6.pdf")

pdm3=c("2R",4179149,4319221)
pdm3_p=c(6,8)
quartz(width=18,height=8)
plot_OR_pV("pdm3",pdm3,pdm3_p)
dev.copy2pdf(file="eu_sa_OR_comp_pdm3_lpeu6_lpsa8.pdf")

# correlation plots lOR
library(car)
quartz()
pairs(~lO_vi+lO_ita+lO_eu+lO_sa,data=eu_sa_odds[eu_sa_odds$lP_eu >= 8,],
   main="log(P_eu)>8")
dev.copy2pdf(file="pair_scatter_lpeu_gt_8.pdf")
quartz()
pairs(~lO_vi+lO_ita+lO_eu+lO_sa,data=eu_sa_odds[eu_sa_odds$lP_ita >= 8,],
   main="log(P_ita)>8")
dev.copy2pdf(file="pair_scatter_lpsa_gt_8.pdf")
quartz()
pairs(~lO_vi+lO_ita+lO_eu+lO_sa,data=eu_sa_odds[eu_sa_odds$lP_sa >= 8,],
   main="log(P_sa)>8")
dev.copy2pdf(file="pair_scatter_lpsa_gt_8.pdf")
quartz()
pairs(~lO_vi+lO_ita+lO_eu+lO_sa,data=eu_sa_odds[eu_sa_odds$lP_vi >= 8,],
   main="log(P_vi)>8")
dev.copy2pdf(file="pair_scatter_lpvi_gt_8.pdf")


tan_snps=with(vie_ita_sa_odds, vie_ita_sa_odds[which(CHR == tan[1] & (tan[2] <= BPS & BPS <= tan[3]) & ( P_eu >= 14 | P_sa >= 10)),] )
bab_snps=with(vie_ita_sa_odds, vie_ita_sa_odds[which(CHR == bab1[1] & (bab1[2] <= BPS & BPS <= bab1[3]) & ( P_eu >= 10 | P_sa >= 14)),] )
inf_val=c("CHR","BPS","P_eu","lO_eu","P_sa","lO_sa","EUb","Sb","rO")


tan_snps_freq=merge(tan_snps[,colnames(tan_snps)],eu_sa_freqs[,colnames(eu_sa_freqs)[c(-3)]],by = c("CHR","BPS"))
bab_snps_freq=merge(bab_snps[,colnames(bab_snps)],eu_sa_freqs[,colnames(eu_sa_freqs)[c(-3)]],by = c("CHR","BPS"))
tan_snps_freq[,c("lO_eu","lO_sa")]=log(tan_snps_freq[,c("O_eu","O_sa")])
bab_snps_freq[,c("lO_eu","lO_sa")]=log(bab_snps_freq[,c("O_eu","O_sa")])
bab_snps_freq$rO=bab_snps_freq$lO_eu/bab_snps_freq$lO_sa
tan_snps_freq$rO=tan_snps_freq$lO_eu/tan_snps_freq$lO_sa

# write files
filename="best_5_tan_bab_eu_sa.txt"
write.table(tan_snps_freq[order(tan_snps$P_eu)[(length(tan_snps$P_sa)-5):length(tan_snps$P_sa)],inf_val],file=filename,append=T,quote=F)
write.table(bab_snps_freq[order(bab_snps$P_sa)[(length(bab_snps$P_sa)-5):length(bab_snps$P_sa)],inf_val],file=filename,append=T,quote=F)
write.table(colMeans(tan_snps_freq[order(tan_snps$P_eu)[(length(tan_snps$P_sa)-5):length(tan_snps$P_sa)],inf_val[c(-1,-2)]]),file=filename,append=T,quote=F)
write.table(colMeans(bab_snps_freq[order(bab_snps$P_sa)[(length(bab_snps$P_sa)-5):length(bab_snps$P_sa)],inf_val[c(-1,-2)]]),file=filename,append=T,quote=F)

sd(tan_snps_freq[order(tan_snps$P_eu)[(length(tan_snps$P_sa)-5):length(tan_snps$P_sa)],])
[1] 0.5442912
sd(bab_snps_freq[order(bab_snps$P_sa)[(length(bab_snps$P_sa)-5):length(bab_snps$P_sa)],"rO"])
[1] 0.2602508
## error propagation:
## s(x/y) = sqrt((1/y)^2*sx^2+(-x/y^2)^2*sy^2)
m.tan=colMeans(tan_snps_freq[order(tan_snps$P_eu)[(length(tan_snps$P_sa)-5):length(tan_snps$P_sa)],c("lO_sa","lO_eu")])
m.bab=colMeans(bab_snps_freq[order(bab_snps$P_sa)[(length(bab_snps$P_sa)-5):length(bab_snps$P_sa)],c("lO_sa","lO_eu")])
s.bab=unlist(lapply(bab_snps_freq[order(bab_snps$P_sa)[(length(bab_snps$P_sa)-5):length(bab_snps$P_sa)],c("lO_sa","lO_eu")],sd))
s.tan=unlist(lapply(tan_snps_freq[order(tan_snps$P_eu)[(length(tan_snps$P_sa)-5):length(tan_snps$P_sa)],c("lO_sa","lO_eu")],sd))
s.t.b=sqrt( (1/m.bab)^2*s.tan^2 + (m.tan/m.bab^2)^2*s.bab^2)
##   lO_sa     lO_eu 
##0.1336343 2.0817782
tan=c("X",9120721,9121170)
tan=c("X",9119470,9122860) # just CG.., MSE and GR8
ebony=c("3R",17055975,17069171)
bab1.3i=c("3L",1036369,1101089)
bab1.3i=c("3L",1071314,1120551 )
bab1.3i=c("3L",1080000,1090000 )
col_freqs=c("CHR","BPS","Allele","VL","VD","IL","ID","Vb","Ib","SVL","SL","SD","SVD","XSL","XSD","EUL","EUD","EUb","Sb")
col_odds=c("CHR","BPS","Allele","P_vi","P_ita","P_eu","P_sa","P_xs","lO_vi","lO_ita","lO_sa","lO_eu","lO_xs")
rank.cols=c("P_vi","P_ita","P_eu","P_sa","P_xs","lO_vi","lO_ita","lO_sa","lO_eu","lO_xs")
rank.names=c("rP_vi","rP_ita","rP_eu","rP_sa","rP_xs","rO_vi","rO_ita","rO_sa","rO_eu","rO_xs")
eu.sa.freqs=merge(eu_sa_freqs[,col_freqs],vie_ita_sa_odds[,col_odds])
eu.sa.freqs[,rank.names]=sapply(1:length(rank.cols), function (x) length(eu.sa.freqs[,rank.cols[x]])-rank(abs(eu.sa.freqs[,rank.cols[x]]))+1)

tan_snps=with(eu.sa.freqs, eu.sa.freqs[which(CHR == tan[1] & (tan[2] <= BPS & BPS <= tan[3]) & ( P_eu >= 3 | P_sa >= 3 |   P_ita >= 3 | P_vi >= 3 | P_xs >= 3  )),] )
bab_snps=with(eu.sa.freqs, eu.sa.freqs[which(CHR == bab1.3i[1] & (bab1.3i[2] <= BPS & BPS <= bab1.3i[3]) & ( P_eu >= 3 | P_sa >= 3 |   P_ita >= 3 | P_vi >= 3 | P_xs >= 3  )),] )
ebony_snps=with(eu.sa.freqs, eu.sa.freqs[which(CHR == ebony[1] & (ebony[2] <= BPS & BPS <= ebony[3]) & ( P_eu >= 3 | P_sa >= 3 |   P_ita >= 3 | P_vi >= 3 | P_xs >= 3  )),] )
ebony_snps[,rank.names]=sapply(1:length(rank.cols), function (x) length(ebony_snps[,rank.cols[x]])-rank(abs(ebony_snps[,rank.cols[x]]))+1)
ebony_snps[order(-ebony_snps$rP_eu),c("CHR","BPS","Vb","Ib","EUb","Sb",rank.names)]
bab_snps[,rank.names]=sapply(1:length(rank.cols), function (x) length(bab_snps[,rank.cols[x]])-rank(abs(bab_snps[,rank.cols[x]]))+1)
best.bab=with(bab_snps,bab_snps[rP_eu <= 50 & rP_sa <= 15 & rP_xs <= 15,])
tan_snps[,rank.names]=sapply(1:length(rank.cols), function (x) length(tan_snps[,rank.cols[x]])-rank(abs(tan_snps[,rank.cols[x]]))+1)
best.tan=with(tan_snps,tan_snps[rP_eu <= 5 & rP_sa <= 10 & rP_xs <= 10,])

best.tan[,c("CHR","BPS","EUb","Sb",rank.names)]
inf_val=c("CHR","BPS","P_eu","lO_eu","P_sa","lO_sa","P_xs","lO_xs","EUb","Sb","rOs","rOx")
best.bab$rOs=best.bab$lO_eu/best.bab$lO_sa
best.tan$rOs=best.tan$lO_eu/best.tan$lO_sa
best.bab$rOx=best.bab$lO_eu/best.bab$lO_xs
best.tan$rOx=best.tan$lO_eu/best.tan$lO_xs

# write files
filename="best_5_tan_bab_eu_sa.txt"
write.table(best.tan[,inf_val],file=filename,append=F,quote=F)
write.table(best.bab[,inf_val],file=filename,append=T,quote=F)
write.table(colMeans(best.tan[,inf_val[c(-1,-2)]]),file=filename,append=T,quote=F)
write.table(colMeans(best.bab[,inf_val[c(-1,-2)]]),file=filename,append=T,quote=F)

eu.sa.xs.odds=read.table("pV_lt_e_6/vie_pIel_pIIel_pIIIel_ita_pIhl_pIIhl_pIIIhl_sa_base_eu_LD_sa_VLVD_sa_xLD_pV_OR.light_alleles_eu_sa_xsa.ancestral_allele_filt.pVlte1e_6.gz",header=F)
colnames(eu.sa.xs.odds)=c("CHR","BPS","Alleles","av1","av2","av3","ab1","ab2","ab3","absa","P_eu","O_eu","Oa_eu","Ob_eu","P_sa","O_sa","Oa_sa","Ob_sa","P_xs","O_xs","Oa_xs","Ob_xs","leu","lsa","lxs","anc")
pvals=c("P_eu","P_sa","P_xs")
ors=c("O_eu","Oa_eu","Ob_eu","O_sa","Oa_sa","Ob_sa","O_xs","Oa_xs","Ob_xs")
eu.sa.xs.odds[,pvals]=-1*log10(eu.sa.xs.odds[,pvals])
eu.sa.xs.odds[,ors]=log(eu.sa.xs.odds[,ors])
haps=c("hap_eu","hap_sa")
eu.sa.xs.odds[,haps]=as.factor(1:length(eu.sa.xs.odds[,1]))
levels(eu.sa.xs.odds$hap_eu) <- c(levels(eu.sa.xs.odds$hap_eu),levels(hap_eu$HAP))
for(i in 1:length(hap_eu$CHR)){
eu.sa.xs.odds$hap_eu[eu.sa.xs.odds$CHR == as.character(hap_eu$CHR[i]) & eu.sa.xs.odds$BPS == hap_eu$BPS[i] ] = as.character(hap_eu$HAP[i])
}
levels(eu.sa.xs.odds$hap_sa) <- c(levels(eu.sa.xs.odds$hap_sa),levels(hap_sa$HAP))
for(i in 1:length(hap_sa$CHR)){
eu.sa.xs.odds$hap_sa[eu.sa.xs.odds$CHR == as.character(hap_sa$CHR[i]) & eu.sa.xs.odds$BPS == hap_sa$BPS[i] ] = as.character(hap_sa$HAP[i])
}
head(eu.sa.xs.odds)
eu.sa.xs.odds$LF=0
maxPs=apply(eu.sa.xs.odds[,c("P_eu","P_sa","P_xs")],1,function(x) max(x,na.rm = T))
eu.sa.xs.odds$LF[! is.na(eu.sa.xs.odds$P_eu) & eu.sa.xs.odds$P_eu >= maxPs ]=sign(eu.sa.xs.odds$O_eu[! is.na(eu.sa.xs.odds$P_eu) & eu.sa.xs.odds$P_eu >= maxPs])
eu.sa.xs.odds$LF[! is.na(eu.sa.xs.odds$O_xs) & ! is.na(eu.sa.xs.odds$P_xs) & eu.sa.xs.odds$P_xs  >= maxPs]=sign(eu.sa.xs.odds$O_xs[! is.na(eu.sa.xs.odds$O_xs) & ! is.na(eu.sa.xs.odds$P_xs) & eu.sa.xs.odds$P_xs  >= maxPs])
eu.sa.xs.odds$LF[! is.na(eu.sa.xs.odds$O_sa) & ! is.na(eu.sa.xs.odds$P_sa) & eu.sa.xs.odds$P_sa  >= maxPs]=sign(eu.sa.xs.odds$O_sa[! is.na(eu.sa.xs.odds$O_sa) & ! is.na(eu.sa.xs.odds$P_sa) & eu.sa.xs.odds$P_sa  >= maxPs])
eu.sa.xs.odds[,c("lO_sa","lOa_sa","lOb_sa","lO_xs","lOa_xs","lOb_xs","lO_eu","lOa_eu","lOb_eu")] = eu.sa.xs.odds[,c("O_sa","Oa_sa","Ob_sa","O_xs","Oa_xs","Ob_xs","O_eu","Oa_eu","Ob_eu")]*eu.sa.xs.odds$LF





all.joined=read.table("pV_lt_e_6/vie_pIel_pIIel_pIIIel_ita_pIhl_pIIhl_pIIIhl_sa_base.af.eu_L_sa_L.light_alleles_eu_LD_sa_VLVD_pVOR_filt.pV_ancest_altAll_eff_aa_gene_effrank_1e_6pV",na.strings = c("NA","nan"),header=F)
colnames(all.joined)=c("CHR","BPS","Alleles","av1","av2","av3","ab1","ab2","ab3","absa","leu","lsa","P_eu","O_eu","Oa_eu","Ob_eu","P_sa","O_sa","Oa_sa","Ob_sa","anc","altA","eff","aa","gene","effrank")
inf_val=c("CHR","BPS","maj","anc","leu","lsa","aEu","absa","P_eu","O_eu","P_sa","O_sa","rPeu","rPsa","rOeu","rOsa","eff","gene","aa")
inf_val_eu=c("CHR","BPS","leu","deu","cons","anc","alEu.Eu","alEu.Sa","nlP_eu","nlP_sa","rPeu","rPsa","l2euO_eu","l2euO_sa","rOeu","rOsa","eff","gene","aa")
inf_val_sa=c("CHR","BPS","lsa","dsa","cons","anc","alSa.Eu","alSa.Sa","nlP_eu","nlP_sa","rPeu","rPsa","l2saO_eu","l2saO_sa","rOeu","rOsa","eff","gene","aa")
tan=c("X",9120721,9121170)
tan=c("X",9000000,9170829)
bab1.3i=c("3L",1036369,1101089)
#bab1.3i=c("3L",1071314,1120551 )
#bab1.3i=c("3L",1080000,1090000 )
all.joined$aEu= rowMeans(all.joined[,c("av1","av2","av3","ab1","ab2","ab3")])
all.joined$rPeu= rank(all.joined$P_eu)
all.joined$rPsa= rank(all.joined$P_sa)
# get allele vectors correct and calcualte light allele frequencies
all.joined$anc=sapply(1:length(all.joined$leu), function (n) if (all.joined$anc[n] == "0" || all.joined$anc[n] == "-1") "-" else as.character(all.joined$anc[n]))
all.joined$maj=substr(all.joined$Alleles,1,1)
all.joined$cons=  sapply(1:length(all.joined$leu), function (n) if (all.joined$leu[n] == all.joined$lsa[n]) "y" else "n" )
all.joined$deu=  sapply(1:length(all.joined$leu), function (n) if (all.joined$leu[n] == all.joined$maj[n]) substr(all.joined$Alleles[n],2,2) else all.joined$maj[n] )
all.joined$dsa=  sapply(1:length(all.joined$lsa), function (n) if (all.joined$lsa[n] == all.joined$maj[n]) substr(all.joined$Alleles[n],2,2) else all.joined$maj[n] )
all.joined$alEu.Eu= abs( (all.joined$deu == all.joined$maj) - all.joined$aEu )
all.joined$alEu.Sa= abs( (all.joined$deu == all.joined$maj) - all.joined$absa )
all.joined$alSa.Eu= abs( (all.joined$dsa == all.joined$maj) - all.joined$aEu )
all.joined$alSa.Sa= abs( (all.joined$dsa == all.joined$maj) - all.joined$absa )
all.joined$O_eu=log2(all.joined$O_eu)
all.joined$O_sa=log2(all.joined$O_sa)
all.joined$nlP_eu=-log10(all.joined$P_eu)
all.joined$nlP_sa=-log10(all.joined$P_sa)

all.joined$l2euO_eu= ( 1- 2*(all.joined$deu == all.joined$maj)) * all.joined$O_eu
all.joined$l2euO_sa= ( 1- 2*(all.joined$deu == all.joined$maj)) * all.joined$O_sa
all.joined$l2saO_sa= ( 1- 2*(all.joined$dsa == all.joined$maj)) * all.joined$O_sa
all.joined$l2saO_eu= ( 1- 2*(all.joined$dsa == all.joined$maj)) * all.joined$O_eu

all.joined$rOsa=abs(rank(abs((all.joined$O_sa)),na.last=F,ties.method = "max")-length(all.joined$O_sa)-1)
all.joined$rOeu=abs(rank(abs((all.joined$O_eu)),na.last=F,ties.method = "max")-length(all.joined$O_eu)-1)

all.joined[( all.joined$absa >= 0.99 | all.joined$absa <= 0.01)  &  (all.joined$aEu >= 0.99 | all.joined$aEu <= 0.01),]
all.joined[ all.joined$BPS == 8888865,]
best100_eu=all.joined[order(all.joined$P_eu),][1:100,]
best100_eu$rOeu=abs(rank(abs((best100_eu$O_eu)),na.last=F,ties.method = "max")-length(best100_eu$O_eu))
best100_eu$rOsa=abs(rank(abs((best100_eu$O_sa)),na.last=F,ties.method = "max")-length(best100_eu$O_eu))
best100_sa=all.joined[order(all.joined$P_sa),][1:100,]
best100_sa$rOeu=abs(rank(abs((best100_sa$O_eu)),na.last=F,ties.method = "max")-length(best100_sa$O_sa))
best100_sa$rOsa=abs(rank(abs((best100_sa$O_sa)),na.last=F,ties.method = "max")-length(best100_sa$O_sa))

tan_snps=with(all.joined, all.joined[which(CHR == tan[1] & (as.numeric(tan[2]) <= BPS & BPS <= as.numeric(tan[3]))),] )
bab_snps=with(all.joined, all.joined[which(CHR == bab1.3i[1] & (as.numeric(bab1.3i[2]) <= BPS & BPS <= as.numeric(bab1.3i[3]))),] )

best50.tan.eu=tan_snps[order(tan_snps$P_eu),][1:50,]
best50.tan.sa=tan_snps[order(tan_snps$P_sa),][1:50,]
best50.tan.sa$rOeu=abs(rank(abs((best50.tan.sa$O_eu)),na.last=F,ties.method = "max")-length(best50.tan.sa$O_sa)-1)
best50.tan.sa$rOsa=abs(rank(abs((best50.tan.sa$O_sa)),na.last=F,ties.method = "max")-length(best50.tan.sa$O_sa)-1)
best50.tan.eu$rOeu=abs(rank(abs((best50.tan.eu$O_eu)),na.last=F,ties.method = "max")-length(best50.tan.eu$O_sa)-1)
best50.tan.eu$rOsa=abs(rank(abs((best50.tan.eu$O_sa)),na.last=F,ties.method = "max")-length(best50.tan.eu$O_sa)-1)

best50.bab.eu=bab_snps[order(bab_snps$P_eu),][1:50,]
best50.bab.sa=bab_snps[order(bab_snps$P_sa),][1:50,]
best50.bab.sa$rOeu=abs(rank(abs((best50.bab.sa$O_eu)),na.last=F,ties.method = "max")-length(best50.bab.sa$O_sa)-1)
best50.bab.sa$rOsa=abs(rank(abs((best50.bab.sa$O_sa)),na.last=F,ties.method = "max")-length(best50.bab.sa$O_sa)-1)
best50.bab.eu$rOeu=abs(rank(abs((best50.bab.eu$O_eu)),na.last=F,ties.method = "max")-length(best50.bab.eu$O_sa)-1)
best50.bab.eu$rOsa=abs(rank(abs((best50.bab.eu$O_sa)),na.last=F,ties.method = "max")-length(best50.bab.eu$O_sa)-1)

write.table(best100_eu[,inf_val_eu],file="best100_eu_maj_anc_leu_lsa_aEu_aSa_Peu_Oeu_Psa_Osa_rPeu_rPsa_rOeu_rOsa",append=F,quote=F,row.names=FALSE,sep="\t")
write.table(best100_sa[,inf_val_sa],file="best100_sa_maj_anc_leu_lsa_aEu_aSa_Peu_Oeu_Psa_Osa_rPeu_rPsa_rOeu_rOsa",append=F,quote=F,row.names=FALSE,sep="\t")
write.table(best50.tan.eu[,inf_val_eu],file="best50.tan.eu_maj_anc_leu_lsa_aEu_aSa_Peu_Oeu_Psa_Osa_rPeu_rPsa_rOeu_rOsa",append=F,quote=F,row.names=FALSE,sep="\t")
write.table(best50.tan.sa[,inf_val_sa],file="best50.tan.sa_maj_anc_leu_lsa_aEu_aSa_Peu_Oeu_Psa_Osa_rPeu_rPsa_rOeu_rOsa",append=F,quote=F,row.names=FALSE,sep="\t")
write.table(best50.bab.sa[,inf_val_sa],file="best50.bab.sa_maj_anc_leu_lsa_aEu_aSa_Peu_Oeu_Psa_Osa_rPeu_rPsa_rOeu_rOsa",append=F,quote=F,row.names=FALSE,sep="\t")
write.table(best50.bab.eu[,inf_val_eu],file="best50.bab.eu_maj_anc_leu_lsa_aEu_aSa_Peu_Oeu_Psa_Osa_rPeu_rPsa_rOeu_rOsa",append=F,quote=F,row.names=FALSE,sep="\t")

# test differences for significance:
median(best50.bab.sa$aEu)
#  0.740073
median(best50.bab.sa$absa)
# 0.5758068
wilcox.test(best50.bab.sa$absa,best50.bab.sa$aEu,alternative="less",paired=T) #  p-value = 0.0004107

median(best50.tan.sa$aEu)
#0.7484895
median(best50.tan.sa$absa)
# 0.7024557
wilcox.test(best50.tan.sa$aEu,best50.tan.sa$absa,alternative="less",paired=T) #  p-value = 0.9946

# get ratio of light/dark ancestral snps in different categories
bab.eu.light=sapply(1:length(best50.bab.eu$CHR),function (x) length(best50.bab.eu$CHR[c(with(best50.bab.eu, (deu == dsa) & anc == leu)[1:x],rep(F,length(best50.bab.eu$CHR)-x))]) )
bab.eu.dark=sapply(1:length(best50.bab.eu$CHR),function (x) length(best50.bab.eu$CHR[c(with(best50.bab.eu, (deu == dsa) & anc == deu)[1:x],rep(F,length(best50.bab.eu$CHR)-x))]) )
 bab.eu.lfr= bab.eu.light/(bab.eu.light+bab.eu.dark)
 bab.eu.dfr= bab.eu.dark/(bab.eu.light+bab.eu.dark)

bab.sa.light=sapply(1:length(best50.bab.sa$CHR),function (x) length(best50.bab.sa$CHR[c(with(best50.bab.sa, (deu == dsa) & anc == leu)[1:x],rep(F,length(best50.bab.sa$CHR)-x))]) )
bab.sa.dark=sapply(1:length(best50.bab.sa$CHR),function (x) length(best50.bab.sa$CHR[c(with(best50.bab.sa, (deu == dsa) & anc == deu)[1:x],rep(F,length(best50.bab.sa$CHR)-x))]) )
 bab.sa.lfr= bab.sa.light/(bab.sa.light+bab.sa.dark)
 bab.sa.dfr= bab.sa.dark/(bab.sa.light+bab.sa.dark)

tan.eu.light=sapply(1:length(best50.tan.eu$CHR),function (x) length(best50.tan.eu$CHR[c(with(best50.tan.eu, (deu == dsa) & anc == leu)[1:x],rep(F,length(best50.tan.eu$CHR)-x))]) )
tan.eu.dark=sapply(1:length(best50.tan.eu$CHR),function (x) length(best50.tan.eu$CHR[c(with(best50.tan.eu, (deu == dsa) & anc == deu)[1:x],rep(F,length(best50.tan.eu$CHR)-x))]) )
 tan.eu.lfr= tan.eu.light/(tan.eu.light+tan.eu.dark)
 tan.eu.dfr= tan.eu.dark/(tan.eu.light+tan.eu.dark)

tan.sa.light=sapply(1:length(best50.tan.sa$CHR),function (x) length(best50.tan.sa$CHR[c(with(best50.tan.sa, (deu == dsa) & anc == leu)[1:x],rep(F,length(best50.tan.sa$CHR)-x))]) )
tan.sa.dark=sapply(1:length(best50.tan.sa$CHR),function (x) length(best50.tan.sa$CHR[c(with(best50.tan.sa, (deu == dsa) & anc == deu)[1:x],rep(F,length(best50.tan.sa$CHR)-x))]) )
 tan.sa.lfr= tan.sa.light/(tan.sa.light+tan.sa.dark)
 tan.sa.dfr= tan.sa.dark/(tan.sa.light+tan.sa.dark)

bab.eu.dfr[bab.eu.dark+bab.eu.light == 10 ]
bab.sa.dfr[bab.sa.dark+bab.sa.light == 10 ]
tan.eu.dfr[tan.eu.dark+tan.eu.light == 10 ]
tan.sa.dfr[tan.sa.dark+tan.sa.light == 10 ]


x11()
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
  	widths=c(3,1), heights=c(1,2))
hist(mtcars$wt)
hist(mtcars$mpg)
hist(mtcars$disp)
dev.off()

setwd("~/vetgrid10/Lukas/Data/A7_pigmentation/Joined_SA_Ita")
sa_afs=read.table("sa_VL1_VL4_VL5_L1_L4_L5_DVD1_DVD4_DVD5_VVD1_VVD4_VVD5_sa_base_merged.af_best_3L_X",na.strings = c("NA","nan"),header=F)
#xs
colnames(sa_afs)=c("CHR","BPS","Allele","VL1","VL2","VL3","L1","L2","L3","D1","D2","D3","VD1","VD2","VD3","B")
vl=colnames(sa_afs)[4:6]
li=colnames(sa_afs)[7:9]
da=colnames(sa_afs)[10:12]
vd=colnames(sa_afs)[13:15]

sa_afs$mVL=rowMeans(sa_afs[,vl])
sa_afs$sdVL=apply(sa_afs[,vl],1,sd)
sa_afs$mL=rowMeans(sa_afs[,li])
sa_afs$sdL=apply(sa_afs[,li],1,sd)
sa_afs$mD=rowMeans(sa_afs[,da])
sa_afs$sdD=apply(sa_afs[,da],1,sd)
sa_afs$mVD=rowMeans(sa_afs[,vd])
sa_afs$sdVD=apply(sa_afs[,vd],1,sd)
sa_raising= (sa_afs$mVL > sa_afs$mVD)
all_freq=c(vl,li,da,vd,"mVL","mL","mD","mVD")

sa_afs[,all_freq]=sa_afs[,all_freq] + (1-sa_raising)*(1-2*sa_afs[,all_freq])

tans=with(sa_afs,sa_afs[CHR == "X" & BPS >= 9120922 & BPS <= 9121177,])
babs=with(sa_afs,sa_afs[CHR == "3L" & ( BPS %in% c(1084990,1086356,1085137,1085454)),])

avgs=c("mVD","mD","mL","mVL")
sds=c("sdVD","sdD","sdL","sdVL")
xs=c(1,2,3,4)
ymax = min(c(1, max(tans[,avgs])+ max(tans[,sds])))
ymin= max(c(0,min(tans[,avgs]) - max(tans[,sds])))
plot(c(1,2,3,4), tans[1,avgs],ylim=c(ymin,ymax),
     pch=NA, xlab=NA, ylab="AF light allele",
     main="tan locus", xaxt='n'
)
axis(1,at=xs,labels=c("VD","D","L","VL"), col="black", line=0.0)
colors=c("red","blue","green","cyan")
for(i in 1:4){
  arrows(xs,unlist(tans[i,avgs]-tans[i,sds]), xs, unlist(tans[i,avgs]+tans[i,sds]), length=0.025, angle=90, code=3, col=colors[i],lty=2)
  lines(xs,tans[i,avgs],pch=20,col=colors[i])
  points(xs,tans[i,avgs],pch=20,col=colors[i],lty=1)
} 
locs=sapply(1:4,function(x) paste(tans$CHR[x],tans$BPS[x]))
legend("topleft",locs,lwd=2,col=colors,title=c("SNPs"),cex=0.75)
dev.copy2pdf(file="tan_best4_afs.pdf")

## babs
ymax = min(c(1, max(babs[,avgs])+ max(babs[,sds])))
ymin= max(c(0,min(babs[,avgs]) - max(babs[,sds])))
plot(c(1,2,3,4), babs[1,avgs],ylim=c(ymin,ymax),
     pch=NA, xlab=NA, ylab="AF light allele",
     main="bab locus", xaxt='n'
)
axis(1,at=xs,labels=c("VD","D","L","VL"), col="black", line=0.0)
colors=c("red","blue","green","cyan")
for(i in 1:4){
  arrows(xs,unlist(babs[i,avgs]-babs[i,sds]), xs, unlist(babs[i,avgs]+babs[i,sds]), length=0.025, angle=90, code=3, col=colors[i],lty=2)
  lines(xs,babs[i,avgs],pch=20,col=colors[i])
  points(xs,babs[i,avgs],pch=20,col=colors[i],lty=1)
} 
locs=sapply(1:4,function(x) paste(babs$CHR[x],babs$BPS[x]))
legend("bottomright",locs,lwd=2,col=colors,title=c("SNPs"),cex=0.75)
dev.copy2pdf(file="bab_best4_afs.pdf")

best4=read.table("best_4_tan_bab.txt",header=T)
best4merged=merge(best4,eu_sa_freqs,by=c("CHR","BPS"))
rm(best4merg)
colnames(best4merged)
best4_af=best4merged[,c("CHR","BPS","EUb.x","Sb.x","EUL","EUD","SVL","SVD","lO_eu","lO_sa")]
afs=c("EUb.x","Sb.x","EUL","EUD","SVL","SVD")
best4_af[,afs]=abs(best4_af[,afs]-(best4_af$lO_eu < 0))
write.table(best4_af,file="best4_afs.txt",quote=F,row.names=F)

# look at afs for unique snps
av.afs=read.table("/Volumes/vetgrid10//Data/A7_pigmentation/Joined_SA_Ita/Exclusive_SNPs/eu_sabase2012_excleu_exclsa.joined_pos.af.gz",header = F)
colnames(av.afs)=c("CHR","BPS","ALLELE","EU","SA","xEU","xSA")
av.afs$EU[is.na(av.afs$EU)] = 1
av.afs$SA[is.na(av.afs$SA)] = 1
av.afs$EU[av.afs$EU > 0.5] = 1 - av.afs$EU[av.afs$EU > 0.5]
av.afs$SA[av.afs$SA > 0.5] = 1 - av.afs$SA[av.afs$SA > 0.5]
av.afs=av.afs[av.afs$EU > 0.001 | av.afs$SA > 0.001 , ]
av.afs$BOTH= ! (av.afs$xEU | av.afs$xSA ) 
head(av.afs)
hist(av.afs$EU[av.afs$BOTH],freq = FALSE,col=rgb(0,0,1,1/4), xlim=c(0,0.5))
hist(av.afs$EU[av.afs$xEU == 1],freq = FALSE,col=rgb(1,0,0,1/4), xlim=c(0,0.5), add=T)
hist(av.afs$SA[av.afs$BOTH],freq = FALSE,col=rgb(0,0,1,1/4), xlim=c(0,0.5))
hist(av.afs$SA[av.afs$xSA == 1],freq = FALSE,col=rgb(1,0,0,1/4), xlim=c(0,0.5), add=T)

wilcox.test(av.afs$EU[av.afs$BOTH]~av.afs$EU[av.afs$xEU == 1]) 
wilcox.test(av.afs$EU[av.afs$BOTH],av.afs$EU[av.afs$xEU == 1],alternative = c("greater"))
wilcox.test(av.afs$SA[av.afs$BOTH],av.afs$SA[av.afs$xSA == 1],alternative = c("greater"))

Wilcoxon rank sum test with continuity correction

data:  av.afs$EU[av.afs$BOTH] and av.afs$EU[av.afs$xEU == 1]
W = 1.0466e+12, p-value < 2.2e-16
alternative hypothesis: true location shift is greater than 0

> wilcox.test(av.afs$SA[av.afs$BOTH],av.afs$SA[av.afs$xSA == 1],alternative = c("greater"))

Wilcoxon rank sum test with continuity correction

data:  av.afs$SA[av.afs$BOTH] and av.afs$SA[av.afs$xSA == 1]
W = 3.47e+12, p-value < 2.2e-16
alternative hypothesis: true location shift is greater than 0

median(av.afs$EU[av.afs$BOTH])
median(av.afs$EU[av.afs$xEU == 1])
median(av.afs$SA[av.afs$BOTH])
median(av.afs$SA[av.afs$xSA == 1])
> median(av.afs$EU[av.afs$BOTH])
[1] 0.1179645
> median(av.afs$EU[av.afs$xEU == 1])
[1] 0.00641
> median(av.afs$SA[av.afs$BOTH])
[1] 0.1166667
> median(av.afs$SA[av.afs$xSA == 1])
[1] 0.01470588


boxplot(av.afs$EU[av.afs$BOTH],av.afs$EU[av.afs$xEU])
head(av.afs$EU[av.afs$xEU == 1])


