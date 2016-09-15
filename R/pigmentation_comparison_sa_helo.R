setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Joined_SA_Ita")
helo_sa_odds=read.table("vie_ita_SA_pV_or_orl_orh_joined.odds",na.strings = c("NA","nan"),header=F,colClasses = c(character(),numeric(),character(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric()))
colnames(helo_sa_odds)=c("CHR","BPS","Allele","P_vi","O_vi","Oa_vi","Ob_vi","P_sa","O_sa","Oa_sa","Ob_sa")
helo_sa_odds$lO_vi=log(helo_sa_odds$O_vi)
helo_sa_odds$lO_sa=log(helo_sa_odds$O_sa)
has_Pv = ! (is.na(helo_sa_odds$P_sa) | is.na(helo_sa_odds$P_vi))
summary(helo_sa_odds$lO_vi[is.finite(helo_sa_odds$lO_vi)]) # Min: -4.058, Max: 4.125
summary(helo_sa_odds$lO_sa[is.finite(helo_sa_odds$lO_sa)]) # Min: -4.2, Max: 3.827
# set inf: 5.0, -inf: -5.0 (for pearson and such)
helo_sa_odds$lO_vi_ni=helo_sa_odds$lO_vi
helo_sa_odds$lO_sa_ni=helo_sa_odds$lO_sa
helo_sa_odds$lO_vi_ni[helo_sa_odds$lO_vi_ni == Inf] = 5
helo_sa_odds$lO_vi_ni[helo_sa_odds$lO_vi_ni == -Inf] = -5
helo_sa_odds$lO_sa_ni[helo_sa_odds$lO_sa_ni == Inf] = 5
helo_sa_odds$lO_sa_ni[helo_sa_odds$lO_sa_ni == -Inf] = -5
# number of same alleles: 1843798
length(helo_sa_odds$lO_sa[has_Pv])
# number of equal sign on odds ratio: 922328 (or 909751 without inf)
length(helo_sa_odds$lO_sa[  sign(helo_sa_odds$lO_sa) == sign(helo_sa_odds$lO_vi) & has_Pv ]) # 921815
length(helo_sa_odds$lO_sa[  sign(helo_sa_odds$lO_sa) == sign(helo_sa_odds$lO_vi) & has_Pv &  (is.finite(helo_sa_odds$lO_vi) & is.finite(helo_sa_odds$lO_vi))]) # 909263
# how many are zero and equal?
length(helo_sa_odds$lO_sa[  sign(helo_sa_odds$lO_sa) == sign(helo_sa_odds$lO_vi) & sign(helo_sa_odds$lO_vi) == 0 & has_Pv ]) # 0

length(helo_sa_odds$lO_sa[  sign(helo_sa_odds$lO_sa) == sign(helo_sa_odds$lO_vi) & helo_sa_odds$P_vi <= 1e-8 & has_Pv] ) # 126
length(helo_sa_odds$lO_sa[  sign(helo_sa_odds$lO_sa) != sign(helo_sa_odds$lO_vi) & helo_sa_odds$P_vi <= 1e-8 & has_Pv] ) # 51

length(helo_sa_odds$lO_sa[  sign(helo_sa_odds$lO_sa) == sign(helo_sa_odds$lO_vi) & helo_sa_odds$P_sa <= 1e-8   & has_Pv] ) # 824
length(helo_sa_odds$lO_sa[  sign(helo_sa_odds$lO_sa) != sign(helo_sa_odds$lO_vi) & helo_sa_odds$P_sa <= 1e-8  & has_Pv ] ) # 719


length(helo_sa_odds$lO_sa[  sign(helo_sa_odds$lO_sa) == sign(helo_sa_odds$lO_vi) & helo_sa_odds$P_vi <= 1e-12  & has_Pv] ) # 27
length(helo_sa_odds$lO_sa[  sign(helo_sa_odds$lO_sa) != sign(helo_sa_odds$lO_vi) & helo_sa_odds$P_vi <= 1e-12  & has_Pv] ) # 7
length(helo_sa_odds$lO_sa[  sign(helo_sa_odds$lO_sa) == sign(helo_sa_odds$lO_vi) & helo_sa_odds$P_sa <= 1e-12  & has_Pv] ) # 61
length(helo_sa_odds$lO_sa[  sign(helo_sa_odds$lO_sa) != sign(helo_sa_odds$lO_vi) & helo_sa_odds$P_sa <= 1e-12  & has_Pv] ) # 16
length(helo_sa_odds$lO_sa[  sign(helo_sa_odds$lO_sa) == sign(helo_sa_odds$lO_vi) & helo_sa_odds$P_sa <= 1e-20  & has_Pv] ) # 12
length(helo_sa_odds$lO_sa[  sign(helo_sa_odds$lO_sa) != sign(helo_sa_odds$lO_vi) & helo_sa_odds$P_sa <= 1e-20  & has_Pv] ) # 5
# spearman between all lods: p=0.047, 0.001460245
cor.test(helo_sa_odds$lO_sa[has_Pv],helo_sa_odds$lO_vi[has_Pv],method="spearman")
# spearman between P_vi < 1e-3:  p-value = 4e-5; rho 0.02467096
cor.test(helo_sa_odds$lO_sa[helo_sa_odds$P_vi <= 1e-3 & has_Pv],helo_sa_odds$lO_vi[helo_sa_odds$P_vi <= 1e-3 & has_Pv],method="spearman")
# spearman between P_sa < 1e-3:  p-value  3e-4 rho 0.016
cor.test(helo_sa_odds$lO_sa[helo_sa_odds$P_sa <= 1e-3 & has_Pv],helo_sa_odds$lO_vi[helo_sa_odds$P_sa <= 1e-3 & has_Pv],method="spearman")
# spearman between P_vi < 1e-8:  p-value = 2.918e-08 rho=0.402
cor.test(helo_sa_odds$lO_sa[helo_sa_odds$P_vi <= 1e-8 & has_Pv],helo_sa_odds$lO_vi[helo_sa_odds$P_vi <= 1e-8 & has_Pv],method="spearman")
# spearman between P_sa < 1e-8:  p-value = 6.905e-08 rho 0.1963
cor.test(helo_sa_odds$lO_sa[helo_sa_odds$P_sa <= 1e-8 & has_Pv],helo_sa_odds$lO_vi[helo_sa_odds$P_sa <= 1e-8 & has_Pv],method="spearman")
# spearman between P_vi < 1e-12:  p-value = 1e-04 rho=0.615
cor.test(helo_sa_odds$lO_sa[helo_sa_odds$P_vi <= 1e-12 & has_Pv],helo_sa_odds$lO_vi[helo_sa_odds$P_vi <= 1e-12 & has_Pv],method="spearman")
# spearman between P_sa < 1e-12:  p-value = 8.905e-06 rho 0.56
cor.test(helo_sa_odds$lO_sa[helo_sa_odds$P_sa <= 1e-12 & has_Pv],helo_sa_odds$lO_vi[helo_sa_odds$P_sa <= 1e-12 & has_Pv],method="spearman")


#pearson between P_vi < 1e-8:  p-value = 3e-05 rho=0.3066184
cor.test(helo_sa_odds$lO_sa_ni[helo_sa_odds$P_vi <= 1e-8 & has_Pv],helo_sa_odds$lO_vi_ni[helo_sa_odds$P_vi <= 1e-8 & has_Pv],method="pearson")
#pearson between P_sa < 1e-8:  p-value = 2e-07 rho 0.22
cor.test(helo_sa_odds$lO_sa_ni[helo_sa_odds$P_sa <= 1e-8 & has_Pv],helo_sa_odds$lO_vi_ni[helo_sa_odds$P_sa <= 1e-8 & has_Pv])

#pearson between P_vi < 1e-12:  p-value = 4e-06 rho=0.615
cor.test(helo_sa_odds$lO_sa_ni[helo_sa_odds$P_vi <= 1e-12 & has_Pv],helo_sa_odds$lO_vi_ni[helo_sa_odds$P_vi <= 1e-12 & has_Pv],method="pearson")
#pearson between P_sa < 1e-12:  p-value = 7.7e-07 rho 0.528
cor.test(helo_sa_odds$lO_sa_ni[helo_sa_odds$P_sa <= 1e-12 & has_Pv],helo_sa_odds$lO_vi_ni[helo_sa_odds$P_sa <= 1e-12 & has_Pv])


setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Joined_SA_Ita")
eu_sa_freqs=read.table('vie_ita_a7_VL1_VL2_VD1_VD2_IL1_IL2_IL3_ID1_ID2_ID3_sa_a7_VL1_VL4_VL5_VVD1_VVD4_VVD5_vie10_pIel_pIIel_pIIIel_ita11_pIhl_pIIhl_pIIIhl_sa_L1_L4_L5_DVD1_DVD4_DVD5_base_2012_eu_p_or_sa_p_or.joined_pos.pV_lt_0.01.good.af.gz',header=F,na.strings=c("NA","NaN","nan","NAN","N"))
colnames(eu_sa_freqs)=c("CHR","BPS","ALLELES","VL1","VL2","VD1","VD2","IL1","IL2","IL3","ID1","ID2","ID3","SL1","SL4","SL5","SD1","SD4","SD5","Vb1","Vb2","Vb3","Ib1","Ib2","Ib3","SDL1","SDL4","SDL5","SLD1","SLD4","SLD5","Sb","Peu","Oeu","Psa","Osa")
#,"lPsa","lPeu")
eu_sa_freqs$lPeu=-1*log10(eu_sa_freqs$Peu)
eu_sa_freqs$lPsa=-1*log10(eu_sa_freqs$Psa)

plot_af_pV <- function(region,coords,lpVals,ylab_pv=expression(-log[10](italic(P))) ) {
    region_idx = which(eu_sa_freqs$CHR == coords[1] & (coords[2] <= eu_sa_freqs$BPS & eu_sa_freqs$BPS <= coords[3] ) & (eu_sa_freqs$lPeu >= lpVals[1] | eu_sa_freqs$lPsa >= lpVals[2]))
    length(region_idx)
    layout(matrix(1:2, ncol = 1), widths = 1, heights = c(0.75,1.25), respect = FALSE)
    par(mar = c(0, 4.1, 4.1, 2.1))
    main_tit=bquote("SNPs around" ~.(region) ~"with" ~ Peu <= 10^-.(lpVals[1]) ~ " or " ~  Psa <= 10 ^-.(lpVals[2]))
    xlimit=c(1-0.25,length(region_idx)+0.15)
    ymin=min(eu_sa_freqs[region_idx,c("SL1","SL4","SL5","Sb","SD1","SD4","SD5","Vb1","Vb2","Vb3","VL1","VL2","IL1","IL2","IL3","VD1","VD2","ID1","ID2","ID3","Ib1","Ib2","Ib3")])
    plot(1:length(region_idx),eu_sa_freqs$lPsa[region_idx],col="blue",cex=1,pch=16,ylim=c(0,max(eu_sa_freqs[region_idx,c("lPsa","lPeu")])),ylab=ylab_pv, xaxt="n",main=main_tit,xlim=xlimit )
    abline(v=(1:length(region_idx))-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+0.25,col="grey",lty="dotted",lw=2)
    points(1:length(region_idx),eu_sa_freqs$lPeu[region_idx],col="red",cex=1,pch=16)
    legend("topright",legend=c("SA","Europe"),pch=c(16,20),col=c("blue","red"))
    par(mar = c(4.1, 4.1, 0, 2.1))
    plot((1:length(region_idx))-0.15,eu_sa_freqs[region_idx,c("Sb")],col="blue",cex=1,pch=".",ylim=c(ymin,1),ylab="AF",xlab="SNPs",xlim=xlimit )
    abline(v=(1:length(region_idx))-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+0.25,col="grey",lty="dotted",lw=2)
    abline(h=seq(0,1,by=0.1),col="grey",lty="dotted",lw=2)
    points((1:length(region_idx))-0.15,eu_sa_freqs[region_idx,c("Sb")],col="blue",cex=1,pch=16)
    points((1:length(region_idx))-0.15,rowMeans(eu_sa_freqs[region_idx,c("SL1","SL4","SL5")]),col="blue",cex=1,pch=6)
    points((1:length(region_idx))-0.15,rowMeans(eu_sa_freqs[region_idx,c("SD1","SD4","SD5")]),col="blue",cex=1,pch=17)
    points((1:length(region_idx))+0.15,rowMeans(eu_sa_freqs[region_idx,c("Vb1","Vb2","Vb3","Ib1","Ib2","Ib3")]),col="red",cex=1,pch=16)
    points((1:length(region_idx))+0.15,rowMeans(eu_sa_freqs[region_idx,c("VL1","VL2","IL1","IL2","IL3")]),col="red",cex=1,pch=6)
    points((1:length(region_idx))+0.15,rowMeans(eu_sa_freqs[region_idx,c("VD1","VD2","ID1","ID2","ID3")]),col="red",cex=1,pch=17)
    legend("bottomright",legend=c("light","base","dark"),title="mean AF",pch=c(6,16,17),col=c("black"))
}

tan=c("X",8867539,9170829)
tan_p=c(12,12)
quartz(width=18,height=8)
plot_af_pV("tan",tan,tan_p)
dev.copy2pdf(file="vie_ita_base_freq_comp_tan.pdf")

tan_p=c(18,10)
quartz(width=18,height=8)
plot_af_pV("tan",tan,tan_p)
dev.copy2pdf(file="vie_ita_base_freq_comp_tan_p18_10.pdf")


bab1=c("3L",1071314,1120551 )
bab1_p=c(10,16)
quartz(width=18,height=8)
plot_af_pV("bab1",bab1,bab1_p)
dev.copy2pdf(file="vie_ita_base_freq_comp_bab1_lpeu10_psa16.pdf")
bab1_p=c(10,18)
quartz(width=18,height=8)
plot_af_pV("bab1",bab1,bab1_p)
dev.copy2pdf(file="vie_ita_base_freq_comp_bab1_lpeu10_psa18.pdf")

ebony=c("3R",17055975,17069171)
ebony_p=c(6,6)
plot_af_pV("ebony",ebony,ebony_p)
dev.copy2pdf(file="vie_ita_base_freq_comp_ebony_lpeu6_lpsa6.pdf")


plot_af_all_pV <- function(region,coords,lpVals,ylab_pv=expression(-log[10](italic(P))) ) {
    region_idx = which(eu_sa_freqs$CHR == coords[1] & (coords[2] <= eu_sa_freqs$BPS & eu_sa_freqs$BPS <= coords[3] ) & (eu_sa_freqs$lPeu >= lpVals[1] | eu_sa_freqs$lPsa >= lpVals[2]))
    length(region_idx)
    layout(matrix(1:2, ncol = 1), widths = 1, heights = c(0.75,1.25), respect = FALSE)
    par(mar = c(0, 4.1, 4.1, 2.1))
    main_tit=bquote("SNPs around" ~.(region) ~"with" ~ Peu <= 10^-.(lpVals[1]) ~ " or " ~  Psa <= 10 ^-.(lpVals[2]))
    xlimit=c(1-0.25,length(region_idx)+0.15)
    plot(1:length(region_idx),eu_sa_freqs$lPsa[region_idx],col="blue",cex=1,pch=16,ylim=c(0,max(eu_sa_freqs[region_idx,c("lPsa","lPeu")],na.rm=T)),ylab=ylab_pv, xaxt="n",main=main_tit,xlim=xlimit )
    abline(v=(1:length(region_idx))-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+0.25,col="grey",lty="dotted",lw=2)
    points(1:length(region_idx),eu_sa_freqs$lPeu[region_idx],col="red",cex=1,pch=16)
    legend("topright",legend=c("SA","Europe"),pch=c(16,20),col=c("blue","red"))
    par(mar = c(4.1, 4.1, 0, 2.1))
    ymin=min(eu_sa_freqs[region_idx,c("SL1","SL4","SL5","Sb","SD1","SD4","SD5","SDL1","SDL4","SDL5","SLD1","SLD4","SLD5","Vb1","Vb2","Vb3","VL1","VL2","IL1","IL2","Ib1","Ib2","Ib3","IL3","VD1","VD2","ID1","ID2","ID3")],na.rm=T)
    plot((1:length(region_idx))-0.15,eu_sa_freqs[region_idx,c("Sb")],col="blue",cex=1,pch=".",ylim=c(ymin,1),ylab="AF",xlab="SNPs",xlim=xlimit )
    abline(v=(1:length(region_idx))-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+0.25,col="grey",lty="dotted",lw=2)
    abline(h=seq(0,1,by=0.1),col="grey",lty="dotted",lw=2)
    points((1:length(region_idx))-0.15,eu_sa_freqs[region_idx,c("Sb")],col="blue",cex=1,pch=16)
    points((1:length(region_idx))-0.15,rowMeans(eu_sa_freqs[region_idx,c("SL1","SL4","SL5")]),col="blue",cex=1,pch=6)
    points((1:length(region_idx))-0.15,rowMeans(eu_sa_freqs[region_idx,c("SD1","SD4","SD5")]),col="blue",cex=1,pch=17)
     points((1:length(region_idx))-0.15,rowMeans(eu_sa_freqs[region_idx,c("SDL1","SDL4","SDL5")]),col="blue",cex=0.75,pch=23)
    points((1:length(region_idx))-0.15,rowMeans(eu_sa_freqs[region_idx,c("SLD1","SLD4","SLD5")]),col="blue",cex=0.75,pch=23,bg="cyan")

    points((1:length(region_idx))+0.15,rowMeans(eu_sa_freqs[region_idx,c("Ib1","Ib2","Ib3","Vb1","Vb2","Vb3")]),col="red",cex=1,pch=16)
    points((1:length(region_idx))+0.15,rowMeans(eu_sa_freqs[region_idx,c("VL1","VL2","IL1","IL2","IL3")]),col="red",cex=1,pch=6)
    points((1:length(region_idx))+0.15,rowMeans(eu_sa_freqs[region_idx,c("VD1","VD2","ID1","ID2","ID3")]),col="red",cex=1,pch=17)
    legend("bottomright",legend=c("very light","light","base","dark","very dark"),title="mean AF",pch=c(6,23,16,23,17),col=c("black"),pt.bg=c("white","white","white","grey16","black"),pt.cex =c(1,0.75,1,0.75,1))
}

tan=c("X",8867539,9170829)
tan_p=c(12,12)
quartz(width=18,height=8)
plot_af_all_pV("tan",tan,tan_p)
dev.copy2pdf(file="vie_ita_all_freqs_comp_tan.pdf")
tan_p=c(20,12)
quartz(width=18,height=12)
plot_af_all_pV("tan",tan,tan_p)
dev.copy2pdf(file="vie_ita_all_freqs_comp_tan_p20_12.pdf")


bab1=c("3L",1071314,1120551 )
bab1_p=c(10,16)
quartz(width=18,height=8)
plot_af_all_pV("bab1",bab1,bab1_p)
dev.copy2pdf(file="vie_ita_all_freq_comp_bab1_lpeu10_psa16.pdf")
bab1_p=c(12,22)
quartz(width=18,height=12)
plot_af_all_pV("bab1",bab1,bab1_p)
dev.copy2pdf(file="vie_ita_all_freq_comp_bab1_lpeu12_psa22.pdf")

ebony=c("3R",17055975,17069171)
ebony_p=c(6,6)
quartz(width=18,height=8)
plot_af_all_pV("ebony",ebony,ebony_p)
dev.copy2pdf(file="vie_ita_all_freq_comp_ebony_lpeu6_lpsa6.pdf")

bab2=c("3L",1128600,1181056)
bab2_p=c(8,10)
quartz(width=18,height=8)
plot_af_all_pV("bab2",bab2,bab2_p)
dev.copy2pdf(file="vie_ita_all_freq_comp_bab2_lpeu8_lpsa10.pdf")
bab2_p=c(6,6)
quartz(width=18,height=8)
plot_af_all_pV("bab2",bab2,bab2_p)
dev.copy2pdf(file="vie_ita_all_freq_comp_bab2_lpeu6_lpsa6.pdf")


pdm3=c("2R",4179149,4319221)
pdm3_p=c(6,6)
quartz(width=18,height=8)
plot_af_all_pV("pdm3",pdm3,pdm3_p)
dev.copy2pdf(file="vie_ita_all_freq_comp_pdm3_lpeu6_lpsa6.pdf")

pdm3_p=c(6,8)
quartz(width=18,height=8)
plot_af_all_pV("pdm3",pdm3,pdm3_p)
dev.copy2pdf(file="vie_ita_all_freq_comp_pdm3_lpeu6_lpsa8.pdf")


setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Joined_SA_Ita")
eu_sa_base_odds=read.table("eu_sa_a7_eu_VL_VD_sa_VVL_VVD_L_D_pV_odds._sa_or_eu_PVs_lt_1e_2.odds.gz",na.strings = c("NA","nan"),header=F,colClasses = c(character(),numeric(),character(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric()))
colnames(eu_sa_base_odds)=c("CHR","BPS","Allele","Pl_eu","Ol_eu","Pd_eu","Od_eu","Pvl_eu","Ovl_eu","Pvd_eu","Ovd_eu",)

setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Joined_SA_Ita")

fsts_joined=read.table('europe_pIel_pIIel_pIIIel_pIhl_pIIhl_pIIIhl_sa_2013_rep1_rep2_q20_joined_pos.fst1:2.1:3.2:3_4:5.4:6.5:6_1:4.1:5.1:6.2:4.2:5.2:6.3:4.3:5.3:6_1:7.2:7.3:7.1:8.2:8.3:8_4:7.5:7.6:7.4:8.5:8.6:8_1:7.2:7.3:7.4:7.5:7.6:7.1:8.2:8.3:8.4:8.5:8.6:8.avg.fsts',header=F,skip=1,na.strings = c("NA","nan"))
colnames(fsts_joined)=c("CHR","BPS","fst_vie","sd_vie","fst_ita","sd_ita","fst_eu","sd_eu","fst_vs","sd_vs","fst_is","sd_is","fst_es","sd_es")
summary(fsts_joined[fsts_joined$fst_es >= 0.000001,])

# calculate binned fsts for all chromosomes
max_bps=max(fsts_joined$BPS)
#28994946
#intervall vector
int_breaks=seq(1,max_bps,by=5000)
fsts_joined$Win5K=cut(fsts_joined$BPS,breaks=int_breaks)
library(foreach)
fst_win5K=fsts_joined[0,-c(15)]
for(i in levels(fsts_joined$CHR)){
 aa= aggregate(x=subset(fsts_joined, CHR==i)[,-c(1,15)], by=list(subset(fsts_joined, CHR ==i)$Win5K), FUN=mean)
 chrom = factor(i,levels(fsts_joined$CHR))
 aa$CHR=rep(chrom,nrow(aa))
 fst_win5K=rbind(fst_win5K,aa[,-1])
}
aggregate(fst_win5K$BPS,by=list(fst_win5K$CHR),FUN=max )
quartz(width=10,height=10)
chroms=c("X","2L","2R","3L","3R")
ylimit=c(0,0.1)
layout(matrix(c(1,6,2,3,4,5),byrow=T, ncol = 2), widths = 1, heights = 1, respect = FALSE)
#par(mar = c(4.1, 4.1, 4.1, 2.1))
for (chrom in chroms){
par(mar = c(4.1, 4.1, 4.1, 2.1))
aa=subset(fst_win5K,CHR==chrom)
fst_col=c("fst_eu","fst_es")
plot(aa$BPS, aa$fst_eu, pch="o", col="olivedrab1", cex=0.5, ylim=ylimit, xlab="base position", main=chrom, ylab=bquote(F[ST]))
points(aa$BPS, aa$fst_es, pch="x", col="rosybrown1", cex=0.5)
lines(aa$BPS,predict(loess(aa$fst_eu~aa$BPS,span=0.01),newdata=aa$BPS), col='green', lwd=1)
lines(aa$BPS,predict(loess(aa$fst_es~aa$BPS,span=0.01),newdata=aa$BPS), col='red', lwd=1)

lines(aa$BPS,predict(loess(aa$fst_vs~aa$BPS,span=0.01),newdata=aa$BPS)-0.001, col='cyan', lwd=1,lty=2)
lines(aa$BPS,predict(loess(aa$fst_is~aa$BPS,span=0.01),newdata=aa$BPS)-0.002, col='blue', lwd=1, lty=2)
lines(aa$BPS,predict(loess(aa$fst_vie~aa$BPS,span=0.01),newdata=aa$BPS), col="darkgoldenrod1", lwd=1,lty=2)
lines(aa$BPS,predict(loess(aa$fst_ita~aa$BPS,span=0.01),newdata=aa$BPS)-0.001, col="darkorange2", lwd=1, lty=2)

## points(aa$BPS, aa$fst_vie, pch=".", col=c("darkgoldenrod1"), cex=0.25)
## points(aa$BPS, aa$fst_ita, pch=".", col=c("darkorange2"), cex=0.25)
## points(aa$BPS, aa$fst_vs, pch=".", col=c("cyan"), cex=0.25)
## points(aa$BPS, aa$fst_is, pch=".", col=c("darkcyan"), cex=0.25)

}
plot(0:2, 0:2, pch = 1, lty = 1, ylim=c(0,2), type = "n", axes = FALSE, ann = FALSE)
#legend(0.75,1.15, legend=c(expression(F[ST] ~ Vienna ),expression(F[ST] ~  Italy),expression(F[ST] ~ Vienna ~ and ~ South ~ Africa),expression(F[ST] ~ Italy ~ and ~ South ~ Africa), expression(F[ST] ~ Vienna ~ and ~ Italy),expression(F[ST] ~ Europe ~ and ~ South ~ Africa)), pch=c(".",".",".",".","o","x"),col=c("darkgoldenrod1","darkorange2","cyan","darkcyan","green","red"),xpd=TRUE)
legend(0.75,1.15, legend=c(expression(F[ST] ~ Vienna ),expression(F[ST] ~  Italy),expression(F[ST] ~ Vienna ~ and ~ South ~ Africa),expression(F[ST] ~ Italy ~ and ~ South ~ Africa), expression(F[ST] ~ Vienna ~ and ~ Italy),expression(F[ST] ~ Europe ~ and ~ South ~ Africa)), lty=c(2,2,2,2,1,1),col=c("darkgoldenrod1","darkorange2","cyan","blue","green","red"),xpd=TRUE)
#dev.copy2pdf(file="fst_bases.pdf")
dev.copy2pdf(file="fst_all_comp_loess_bases.pdf")
a=subset(fst_win5K, CHR=='2R' & fst_es >= 0.04)
a=subset(fst_win5K, CHR=='3R' & fst_es >= 0.05)
a=subset(fst_win5K, CHR=='X' & fst_es >= 0.05)


fst_win5K=read.table("europe_pIel_pIIel_pIIIel_pIhl_pIIhl_pIIIhl_sa_2013_rep1_rep2_q20_joined_pos_win5000K.fst1:2.1:3.2:3_4:5.4:6.5:6_1:4.1:5.1:6.2:4.2:5.2:6.3:4.3:5.3:6_1:7.2:7.3:7.1:8.2:8.3:8_4:7.5:7.6:7.4:8.5:8.6:8_1:7.2:7.3:7.4:7.5:7.6:7.1:8.2:8.3:8.4:8.5:8.6:8.avg.fsts",header=F,skip=1,na.strings = c("NA","nan"))
colnames(fst_win5K)=c("CHR","BPS","fst_vie","sd_vie","fst_ita","sd_ita","fst_eu","sd_eu","fst_vs","sd_vs","fst_is","sd_is","fst_es","sd_es")
summary(fst_win5K[fst_win5K$fst_es >= 0.000001,])

setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Joined_SA_Ita")
eu_sa_odds=read.table("eu_LD_sa_VL_VVD_eu_L_B_eu_D_B_sa_VL_B_VVD_B_D_B_L_B.ident_pos.good.odds.gz",na.strings = c("NA","nan"),header=F,colClasses = c(character(),numeric(),character(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric()))
colnames(eu_sa_odds)=c("CHR","BPS","Allele","P_eu","O_eu","Oa_eu","Ob_eu","P_sa","O_sa","Oa_sa","Ob_sa","P_el","O_el","Oa_el","Ob_el","P_ed","O_ed","Oa_ed","Ob_ed","P_svl","O_svl","Oa_svl","Ob_svl","P_svd","O_svd","Oa_svd","Ob_svd","P_sd","O_sd","Oa_sd","Ob_sd","P_sl","O_sl","Oa_sl","Ob_sl")
drops=c("Oa_ed","Ob_ed","Oa_el","Ob_el","Oa_svl","Ob_svl","Oa_svd","Ob_svd","Oa_sd","Ob_sd","Oa_sl","Ob_sl")
eu_sa_odds <- eu_sa_odds[,!(names(eu_sa_odds) %in% drops)]
head (eu_sa_odds)
eu_sa_odds$O_eu=log(eu_sa_odds$O_eu)
eu_sa_odds$O_sa=log(eu_sa_odds$O_sa)
eu_sa_odds$Oa_eu=log(eu_sa_odds$Oa_eu)
eu_sa_odds$Oa_sa=log(eu_sa_odds$Oa_sa)
eu_sa_odds$Ob_eu=log(eu_sa_odds$Ob_eu)
eu_sa_odds$Ob_sa=log(eu_sa_odds$Ob_sa)
eu_sa_odds$O_el=log(eu_sa_odds$O_el)
eu_sa_odds$O_ed=log(eu_sa_odds$O_ed)
eu_sa_odds$O_sl=log(eu_sa_odds$O_sl)
eu_sa_odds$O_sd=log(eu_sa_odds$O_sd)
eu_sa_odds$O_svl=log(eu_sa_odds$O_svl)
eu_sa_odds$O_svd=log(eu_sa_odds$O_svd)

# get correlation of odds
has_Pv = ! (is.na(eu_sa_odds$P_sa) | is.na(eu_sa_odds$P_eu))
has_OR = ! (is.na(eu_sa_odds$O_sa) | is.na(eu_sa_odds$O_eu))
summary(eu_sa_odds$O_eu[is.finite(eu_sa_odds$O_eu)]) # Min: -4.101000, Max: 4.125000
summary(eu_sa_odds$O_sa[is.finite(eu_sa_odds$O_sa)]) # Min: -4.20200, Max: 3.89400

eu_sa_odds$O_eu_ni=eu_sa_odds$O_eu
eu_sa_odds$O_sa_ni=eu_sa_odds$O_sa
eu_sa_odds$O_eu_ni[eu_sa_odds$O_eu_ni == Inf] = 5
eu_sa_odds$O_eu_ni[eu_sa_odds$O_eu_ni == -Inf] = -5
eu_sa_odds$O_sa_ni[eu_sa_odds$O_sa_ni == Inf] = 5
eu_sa_odds$O_sa_ni[eu_sa_odds$O_sa_ni == -Inf] = -5

ident_snps=length(eu_sa_odds$O_sa[has_Pv]) # 1927905

length(eu_sa_odds$O_sa) # 1928842

# number of equal sign on odds ratio: 922328 (or 909751 without inf)
length(eu_sa_odds$O_sa[  sign(eu_sa_odds$O_sa) == sign(eu_sa_odds$O_eu) & has_Pv ]) # 964120
length(eu_sa_odds$O_sa[  sign(eu_sa_odds$O_sa) == sign(eu_sa_odds$O_eu) & has_Pv &  (is.finite(eu_sa_odds$O_eu) & is.finite(eu_sa_odds$O_eu))]) # 949954
length(eu_sa_odds$O_sa[  sign(eu_sa_odds$O_sa) == sign(eu_sa_odds$O_eu) & has_Pv ])/ident_snps #0.5000869
length(eu_sa_odds$O_sa[  sign(eu_sa_odds$O_sa) == sign(eu_sa_odds$O_eu) & has_Pv &  (is.finite(eu_sa_odds$O_eu) & is.finite(eu_sa_odds$O_eu))])/ident_snps #  0.492739
# how many are zero and equal?
length(eu_sa_odds$O_sa[  sign(eu_sa_odds$O_sa) == sign(eu_sa_odds$O_eu) & sign(eu_sa_odds$O_eu) == 0 & has_Pv ]) # 0

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
dev.copy2pdf(file="eu_sa_lOR_sign_quant.pdf")


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
dev.copy2pdf(file="eu_sa_lOR_cor_OR_quant.pdf")

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

plot_region = function(chrom,a,b, title,ylab=expression(-log[10](italic(P)))){
    region_sa = which(sa_plot$CHR == chrom & sa_plot$BP >= a  & sa_plot$BP <= b )
    region_eu = which(eu_plot$CHR == chrom & eu_plot$BP >= a  & eu_plot$BP <= b )
    ymax = max(-log10(sa_plot$P[region_sa]),-log10(eu_plot$P[region_eu]))
    plot(sa_plot$BP[region_sa],-log10(sa_plot$P[region_sa]), col="blue", pch=20, main=title, ylim=c(0,ymax),xlab=chrom,ylab=ylab,cex=0.75 )
    points(eu_plot$BP[region_eu],-log10(eu_plot$P[region_eu]), col="red", pch=20,cex=0.75)
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
plot_region("X",9116995,9121764,"tan")
rect(9111688,45,9117290,47,col="grey")
text(9117100,48,labels="tan")
text(9117050,46,labels="< < <")
rect(9121170,45,9122860,47,col="grey")
text(9121300,48,labels="Gr8a")
text(9121350,46,labels="> > >")
rect(9119470,45,9120721,47,col="grey")
text((9119470+9120721)/2,48,labels="CG15370")
text((9119470+9120721)/2,46,labels="> > >")
rect(9120721,45,9121170,47,col="white")
text((9120721+9121170)/2,48,labels="MSE")

dev.copy2pdf(file="eu_sa_tan_region.pdf")

quartz()
plot_region("3L",1039508,1111138,"bab1")
segments(1036369,36,1101089,36)
rect(c(1099101,1039796,1036369,1036369),35,c(1101089,1040223,1039215,1039215),37,col="grey")
rect(1079062,35.75,1079672,36.25,col="grey")
rect(1084814,35.5,1085477,36.5,col="white")
text((1084814+1085477)/2,34.5,labels="DME")

text(1060000,38,labels="bab1")
text(1060000,36,labels="< < <")
dev.copy2pdf(file="eu_sa_bab1_big_region.pdf")

quartz()
plot_region("3L",1070000,1111138,"bab1")
segments(1036369,36,1101089,36)
rect(c(1099101,1039796,1036369,1036369),35,c(1101089,1040223,1039215,1039215),37,col="grey")
rect(1079062,35.75,1079672,36.25,col="grey")
rect(1084814,35.5,1085477,36.5,col="white")
text((1084814+1085477)/2,34.5,labels="DME")
text(1075000,38,labels="bab1")
text(1075000,36,labels="< < <")
dev.copy2pdf(file="eu_sa_bab1_small_region.pdf")
#pdm3
quartz()
plot_region("2R",4210633,4235314,"pdm3")
segments(4214995,0.15,42837859,0.15)
rect(c(4214995,4235006),0,c(4215512,4235799),0.3,col="grey")
text(4225000,0.4,labels="pdm3")
text(4225000,0.14,labels="> > >")
dev.copy2pdf(file="eu_sa_pdm3_region.pdf")

# ebony
quartz()
plot_region("3R",17062142,17065683,"ebony")
rect(17055561,0,17062899,0.5,col="grey")
text(17062700,0.7,labels="e")
text(17062650,0.25,labels="< < <")
dev.copy2pdf(file="eu_sa_ebony_region.pdf")


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



indels_sa7=read.table("/Volumes/Temp/Lukas/Data/SA_A7/Polymorphisms/Indels/sa_a7_indels_all_filter_95_15_maxAF_gt_0.05_minDP_gt_15_allele_depths_VL_VVD.gwas",header=T)
off_sa_indels=get_offset_manh(indels_sa7,limitchromosomes=chroms)


quartz(height=6,width=12)
par(bg="white")
chroms=c("X","2L","2R","3L","3R","4")
manhattan(indels_sa7,colors=c("black","slategrey"),limitchromosomes=chroms,main="South Africa (2012), InDels")
rect(8867539+off_eu["X"],0,9170829+off_eu["X"],60,lty=3, border="red")
rect(1071314+off_eu["3L"],0,1120551+off_eu["3L"],60,lty=3, border="green")
rect(17055975+off_eu["3R"],0,17069171+off_eu["3R"],60,lty=3, border="blue")
rect(4179149+off_eu["2R"],0,4319221+off_eu["2R"],60,lty=3, border="yellow")

legend("topright",legend=c("tan","bab","ebony","pdm3"),lty=1,col=c("red","green","blue","yellow"))
dev.print(png,width=800,file="sa_indel_manhattan.png")

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

setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Joined_SA_Ita")
a7_inv_data =  read.table("vie_ita_a7_VL1_VL2_VD1_VD2_IL1_IL2_IL3_ID1_ID2_ID3_sa_a7_VL1_VL4_VL5_VVD1_VVD4_VVD5_vie10_pIel_pIIel_pIIIel_ita11_pIhl_pIIhl_pIIIhl_sa_L1_L4_L5_DVD1_DVD4_DVD5_base_2012.joined_pos.inversions",header=TRUE)
library(doBy)
head(a7_inv_data)
group_names = c( "IN(2L)t"  ,"IN(2R)NS" ,"IN(3L)P" , "IN(3R)C" , "IN(3R)Mo" ,"IN(3R)P")
colnames(a7_inv_data)=c("CHR","INV","BPS","REF","ALT","ViL1","ViL1c","ViL2","ViL2c","ViD1","ViD1c","ViD2","ViD2c","IL1","IL1C","IL2","IL2C","IL3","IL3C","ID1","ID1C","ID2","ID2C","ID3","ID3c","SaVL1","SaVL1c","SaVL4","SaVL4c","SaVL5","SaVL5C","SaVVD1","SaVVD1c","SaVVD4","SaVVD4c","SaVVD5","SaVVD5c","VipIel","VipIelc","VipIIel","VipIIelc","VipIIIel","VipIIIelc","IpIhl","IpIhlc","IpIIhl","IpIIhlc","IpIIIhl","IpIIIhlc","SaL1","SaL1c","SaL4","SaL4c","SaL5","SaL5c","SaDVD1","SaDVD1c","SaDVD4","SaDVD4c","SaDVD5","SaDVD5c","SAb","SAbc")


a7_base_vie=c("VipIel","VipIIel","VipIIIel")
a7_base_ita=c("IpIhl","IpIIhl","IpIIIhl")
a7l_vie=c("ViL1","ViL2")
a7d_vie=c("ViD1","ViD2")
a7l_ita=c("IL1","IL2","IL3")
a7d_ita=c("ID1","ID2","ID3")
a7l_sa=c("SaVL1","SaVL4","SaVL5")
a7d_sa=c("SaVVD1","SaVVD4","SaVVD5")
a7_base_sa=c("SAb")

a7_inv_list <- splitBy(formula=~ CHR + INV, data = a7_inv_data )
quartz()
boxplot(unlist(a7_inv_list[[6]][,a7l_vie]),unlist(a7_inv_list[[6]][,a7_base_vie]),unlist(a7_inv_list[[6]][a7d_vie]),unlist(a7_inv_list[[6]][,a7l_ita]),unlist(a7_inv_list[[6]][,a7_base_ita]),unlist(a7_inv_list[[6]][a7d_ita]),unlist(a7_inv_list[[6]][,a7l_sa]),unlist(a7_inv_list[[6]][,a7_base_sa]),unlist(a7_inv_list[[6]][a7d_sa]),data=a7_inv_list[[6]],notch=T,col=c("cornsilk","yellow","brown4"),at=c(0.5,1,1.5,2.5,3,3.5,4.5,5,5.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark", "light", "base","dark"), main="Inversion 3R(p)")
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


