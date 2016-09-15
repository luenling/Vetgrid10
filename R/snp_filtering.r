setwd("/Volumes/Temp/Lukas/Data/SA_A7/BGI_100_101/Polymorphims/CMH")
snp_sa7=read.table("sa_a7_metrics_ident_comb.tab.gz",header=F,na.strings = c("NA","nan"))
head(snp_sa7)
colnames(snp_sa7)=c("CHR","BPS","ALL","FS","SB","BSB","RTD","ATD","TDB","RPB","DD_i","DD_c")
snp_sa7$aRPB=abs(snp_sa7$RPB)

library(ggplot2)

xlimit=c(0,600)
quartz()
good=ggplot(data=snp_sa7) + xlim(xlimit) +
     stat_binhex(data=snp_sa7[snp_sa7$DD_i == 1,],aes(x=FS, y=SB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=snp_sa7) + xlim(xlimit) +
     stat_binhex(data=snp_sa7[snp_sa7$DD_c == 0,],aes(x=FS, y=SB),bins=40) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)  
dev.copy2pdf(file="sa_a7_FS_SB_double.pdf")

quartz()
good=ggplot(data=snp_sa7) +
     stat_binhex(data=snp_sa7[snp_sa7$DD_i == 1,],aes(x=aRPB, y=ATD),binwidth=c(0.5,1)) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=snp_sa7) +
     stat_binhex(data=snp_sa7[snp_sa7$DD_c == 0,],aes(x=aRPB, y=ATD),binwidth=c(0.5,1)) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)     
dev.copy2pdf(file="sa_a7_aRPB_ATD_double.pdf")


# for snp_sa7$DD_i == 1
good_snps= snp_sa7$DD_i == 1
# for snp_sa7$DD_c == 0
bad_snps= snp_sa7$DD_c == 0
num_good=length(snp_sa7$FS[good_snps])
num_bad= length(snp_sa7$FS[bad_snps])
num_tot= length(snp_sa7$FS)
num_tot #5646977
num_good #3429225
num_bad #1108450
length(snp_sa7$DD_i[( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 )]) # 26036
length(snp_sa7$DD_i[ ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 ) & bad_snps]) # 16822
length(snp_sa7$DD_i[( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 ) & good_snps]) # 4393
length(snp_sa7$DD_i[(snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) & bad_snps]) #105647
length(snp_sa7$DD_i[(snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) & good_snps]) #26457
length(snp_sa7$DD_i[(snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5)]) #159604
length(snp_sa7$DD_i[(snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) | ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 )  ]) # 169241
length(snp_sa7$DD_i[((snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) | ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 )) & good_snps]) # 28787
length(snp_sa7$DD_i[((snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) | ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 ) ) & bad_snps]) # 111327
length(snp_sa7$DD_i[! ((snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) | ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 ))]) # 5497692

# create ROC like curves for retrieval SNPs
quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
fs_vals=c(5,10,20,30,60,100,150,200,250)
good_vals = sapply(fs_vals,function (x) length(snp_sa7$FS[good_snps & snp_sa7$FS >= x] ))
bad_vals =  sapply(fs_vals,function (x) length(snp_sa7$FS[bad_snps & snp_sa7$FS >= x] ))
plot(fs_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.2),xlab="FS",ylab="fraction filtered")
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=1)
SB_vals=c(0.2,0.1,0.05,0.01)
ltypes=c(2,3,4,5)
for(i in 1:4){
good_vals = sapply(fs_vals,function (x) length(snp_sa7$FS[good_snps & snp_sa7$FS >= x & snp_sa7$SB <= SB_vals[i]] ))
bad_vals =  sapply(fs_vals,function (x) length(snp_sa7$FS[bad_snps & snp_sa7$FS >= x & snp_sa7$SB <=  SB_vals[i]] ))
lines(fs_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","SB <= 1.0","SB <= 0.2","SB <= 0.1","SB <= 0.05","SB <= 0.01"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,1,2,3,4,5))
dev.copy2pdf(file="sa_a7_FS_SB_SNPs_filtered.pdf")
quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
fs_vals=c(5,10,20,30,60,100,150,200,250)
good_vals = sapply(fs_vals,function (x) length(snp_sa7$FS[good_snps & snp_sa7$FS >= x] ))
bad_vals =  sapply(fs_vals,function (x) length(snp_sa7$FS[bad_snps & snp_sa7$FS >= x] ))
plot(fs_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.2),xlab="FS",ylab="fraction filtered")
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=1)
bsb_vals=c(25,50,100,150)
ltypes=c(2,3,4,5)
for(i in 1:4){
good_vals = sapply(fs_vals,function (x) length(snp_sa7$FS[good_snps & snp_sa7$FS >= x & snp_sa7$BSB >= bsb_vals[i]] ))
bad_vals =  sapply(fs_vals,function (x) length(snp_sa7$FS[bad_snps & snp_sa7$FS >= x & snp_sa7$BSB >=  bsb_vals[i]] ))
lines(fs_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","BSB >= 0.0","BSB >= 25","BSB >= 50","BSB >= 100","BSB >= 200"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,1,2,3,4,5))
dev.copy2pdf(file="sa_a7_FS_BSB_SNPs_filtered.pdf")

quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
aRPB_vals=c(0,1,2,3,3.5,4,4.5,5,7.5,10,15)
good_vals = sapply(aRPB_vals,function (x) length(snp_sa7$aRPB[good_snps & snp_sa7$aRPB >= x] ))
bad_vals =  sapply(aRPB_vals,function (x) length(snp_sa7$aRPB[bad_snps & snp_sa7$aRPB >= x] ))
plot(aRPB_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.25),xlab="aRPB",ylab="fraction filtered")
lines(aRPB_vals, bad_vals/num_bad,type="l", col="red", lty=1)
ATD_vals=c(15,10,7.5,5,2.5,1)
ltypes=c(2,3,4,5,6,2)
for(i in 1:4){
good_vals = sapply(aRPB_vals,function (x) length(snp_sa7$aRPB[good_snps & snp_sa7$aRPB >= x & snp_sa7$ATD <= ATD_vals[i]] ))
bad_vals =  sapply(aRPB_vals,function (x) length(snp_sa7$aRPB[bad_snps & snp_sa7$aRPB >= x & snp_sa7$ATD <=  ATD_vals[i]] ))
lines(aRPB_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(aRPB_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","ATD <= 15","ATD <= 10","ATD <= 7.5","ATD <= 5","ATD <= 2.5","ATD <= 1.0"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,2,3,4,5,6,2))
dev.copy2pdf(file="sa_a7_aRPB_ATD_SNPs_filtered.pdf")


quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
TDB_vals=c(0,10,25,50,75,100)
good_vals = sapply(TDB_vals,function (x) length(snp_sa7$TDB[good_snps & snp_sa7$TDB >= x] ))
bad_vals =  sapply(TDB_vals,function (x) length(snp_sa7$TDB[bad_snps & snp_sa7$TDB >= x] ))
plot(TDB_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.25),xlab="TDB",ylab="fraction filtered")
lines(TDB_vals, bad_vals/num_bad,type="l", col="red", lty=1)
ATD_vals=c(15,10,7.5,5,2.5,1)
ltypes=c(2,3,4,5,6,2)
for(i in 1:4){
good_vals = sapply(TDB_vals,function (x) length(snp_sa7$TDB[good_snps & snp_sa7$TDB >= x & snp_sa7$ATD <= ATD_vals[i]] ))
bad_vals =  sapply(TDB_vals,function (x) length(snp_sa7$TDB[bad_snps & snp_sa7$TDB >= x & snp_sa7$ATD <=  ATD_vals[i]] ))
lines(TDB_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(TDB_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","ATD <= 15","ATD <= 10","ATD <= 7.5","ATD <= 5","ATD <= 2.5","ATD <= 1.0"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,2,3,4,5,6,2))
rm("good_vals","bad_vals")
dev.copy2pdf(file="sa_a7_TDB_ATD_SNPs_filtered.pdf")
gc()

setwd("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/Polymorphisms/CMH")
snp_sa7=read.table("females_vie_ita_metrics_ident_comb.tab.gz",header=F,na.strings = c("NA","nan"))
head(snp_sa7)
colnames(snp_sa7)=c("CHR","BPS","ALL","FS","SB","BSB","RTD","ATD","TDB","RPB","DD_i","DD_c")
snp_sa7$aRPB=abs(snp_sa7$RPB)

library(ggplot2)

xlimit=c(0,600)
quartz()
good=ggplot(data=snp_sa7) + xlim(xlimit) +
     stat_binhex(data=snp_sa7[snp_sa7$DD_i == 1,],aes(x=FS, y=SB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=snp_sa7) + xlim(xlimit) +
     stat_binhex(data=snp_sa7[snp_sa7$DD_c == 0,],aes(x=FS, y=SB),bins=40) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)  
dev.copy2pdf(file="vie_ita_FS_SB_double.pdf")

quartz()
good=ggplot(data=snp_sa7) +
     stat_binhex(data=snp_sa7[snp_sa7$DD_i == 1,],aes(x=aRPB, y=ATD),binwidth=c(0.5,1)) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=snp_sa7) +
     stat_binhex(data=snp_sa7[snp_sa7$DD_c == 0,],aes(x=aRPB, y=ATD),binwidth=c(0.5,1)) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)     
dev.copy2pdf(file="vie_ita_aRPB_ATD_double.pdf")
rm("good","bad")
gc()

# for snp_sa7$DD_i == 1
good_snps= snp_sa7$DD_i == 1
# for snp_sa7$DD_c == 0
bad_snps= snp_sa7$DD_c == 0
num_good=length(snp_sa7$FS[good_snps])
num_bad= length(snp_sa7$FS[bad_snps])
num_tot= length(snp_sa7$FS)
num_tot #1914804
num_good #1482916
num_bad #202792
length(snp_sa7$DD_i[( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 )]) # 6652
length(snp_sa7$DD_i[ ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 ) & bad_snps]) # 4800
length(snp_sa7$DD_i[( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 ) & good_snps]) # 973
length(snp_sa7$DD_i[(snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) & bad_snps]) #31883
length(snp_sa7$DD_i[(snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) & good_snps]) #5955
length(snp_sa7$DD_i[(snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5)]) #44217
length(snp_sa7$DD_i[(snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) | ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 )  ]) # 46373
length(snp_sa7$DD_i[((snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) | ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 )) & good_snps]) # 6529
length(snp_sa7$DD_i[((snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) | ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 ) ) & bad_snps]) # 33019
length(snp_sa7$DD_i[! ((snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) | ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 ))]) # 1868431

# create ROC like curves for retrieval SNPs
quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
fs_vals=c(5,10,20,30,60,100,150,200,250)
good_vals = sapply(fs_vals,function (x) length(snp_sa7$FS[good_snps & snp_sa7$FS >= x] ))
bad_vals =  sapply(fs_vals,function (x) length(snp_sa7$FS[bad_snps & snp_sa7$FS >= x] ))
plot(fs_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.2),xlab="FS",ylab="fraction filtered")
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=1)
SB_vals=c(0.2,0.1,0.05,0.01)
ltypes=c(2,3,4,5)
for(i in 1:4){
good_vals = sapply(fs_vals,function (x) length(snp_sa7$FS[good_snps & snp_sa7$FS >= x & snp_sa7$SB <= SB_vals[i]] ))
bad_vals =  sapply(fs_vals,function (x) length(snp_sa7$FS[bad_snps & snp_sa7$FS >= x & snp_sa7$SB <=  SB_vals[i]] ))
lines(fs_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","SB <= 1.0","SB <= 0.2","SB <= 0.1","SB <= 0.05","SB <= 0.01"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,1,2,3,4,5))
dev.copy2pdf(file="vie_ita_FS_SB_SNPs_filtered.pdf")
quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
fs_vals=c(5,10,20,30,60,100,150,200,250)
good_vals = sapply(fs_vals,function (x) length(snp_sa7$FS[good_snps & snp_sa7$FS >= x] ))
bad_vals =  sapply(fs_vals,function (x) length(snp_sa7$FS[bad_snps & snp_sa7$FS >= x] ))
plot(fs_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.2),xlab="FS",ylab="fraction filtered")
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=1)
bsb_vals=c(25,50,100,150)
ltypes=c(2,3,4,5)
for(i in 1:4){
good_vals = sapply(fs_vals,function (x) length(snp_sa7$FS[good_snps & snp_sa7$FS >= x & snp_sa7$BSB >= bsb_vals[i]] ))
bad_vals =  sapply(fs_vals,function (x) length(snp_sa7$FS[bad_snps & snp_sa7$FS >= x & snp_sa7$BSB >=  bsb_vals[i]] ))
lines(fs_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","BSB >= 0.0","BSB >= 25","BSB >= 50","BSB >= 100","BSB >= 200"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,1,2,3,4,5))
dev.copy2pdf(file="vie_ita_FS_BSB_SNPs_filtered.pdf")

quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
aRPB_vals=c(0,1,2,3,3.5,4,4.5,5,7.5,10,15)
good_vals = sapply(aRPB_vals,function (x) length(snp_sa7$aRPB[good_snps & snp_sa7$aRPB >= x] ))
bad_vals =  sapply(aRPB_vals,function (x) length(snp_sa7$aRPB[bad_snps & snp_sa7$aRPB >= x] ))
plot(aRPB_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.25),xlab="aRPB",ylab="fraction filtered")
lines(aRPB_vals, bad_vals/num_bad,type="l", col="red", lty=1)
ATD_vals=c(15,10,7.5,5,2.5,1)
ltypes=c(2,3,4,5,6,2)
for(i in 1:4){
good_vals = sapply(aRPB_vals,function (x) length(snp_sa7$aRPB[good_snps & snp_sa7$aRPB >= x & snp_sa7$ATD <= ATD_vals[i]] ))
bad_vals =  sapply(aRPB_vals,function (x) length(snp_sa7$aRPB[bad_snps & snp_sa7$aRPB >= x & snp_sa7$ATD <=  ATD_vals[i]] ))
lines(aRPB_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(aRPB_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","ATD <= 15","ATD <= 10","ATD <= 7.5","ATD <= 5","ATD <= 2.5","ATD <= 1.0"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,2,3,4,5,6,2))
dev.copy2pdf(file="vie_ita_aRPB_ATD_SNPs_filtered.pdf")


quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
TDB_vals=c(0,10,25,50,75,100)
good_vals = sapply(TDB_vals,function (x) length(snp_sa7$TDB[good_snps & snp_sa7$TDB >= x] ))
bad_vals =  sapply(TDB_vals,function (x) length(snp_sa7$TDB[bad_snps & snp_sa7$TDB >= x] ))
plot(TDB_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.25),xlab="TDB",ylab="fraction filtered")
lines(TDB_vals, bad_vals/num_bad,type="l", col="red", lty=1)
ATD_vals=c(15,10,7.5,5,2.5,1)
ltypes=c(2,3,4,5,6,2)
for(i in 1:4){
good_vals = sapply(TDB_vals,function (x) length(snp_sa7$TDB[good_snps & snp_sa7$TDB >= x & snp_sa7$ATD <= ATD_vals[i]] ))
bad_vals =  sapply(TDB_vals,function (x) length(snp_sa7$TDB[bad_snps & snp_sa7$TDB >= x & snp_sa7$ATD <=  ATD_vals[i]] ))
lines(TDB_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(TDB_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","ATD <= 15","ATD <= 10","ATD <= 7.5","ATD <= 5","ATD <= 2.5","ATD <= 1.0"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,2,3,4,5,6,2))
rm("good_vals","bad_vals")
dev.copy2pdf(file="vie_ita_TDB_ATD_SNPs_filtered.pdf")
gc()


setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Polymorphisms/CMH")
snp_sa7=read.table("vie_ita_a7_metrics_ident_comb.tab.gz",header=F,na.strings = c("NA","nan"))
head(snp_sa7)
colnames(snp_sa7)=c("CHR","BPS","ALL","FS","SB","BSB","RTD","ATD","TDB","RPB","DD_i","DD_c")
snp_sa7$aRPB=abs(snp_sa7$RPB)

xlimit=c(0,600)
quartz()
good=ggplot(data=snp_sa7) + xlim(xlimit) +
     stat_binhex(data=snp_sa7[snp_sa7$DD_i == 1,],aes(x=FS, y=SB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=snp_sa7) + xlim(xlimit) +
     stat_binhex(data=snp_sa7[snp_sa7$DD_c == 0,],aes(x=FS, y=SB),bins=40) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)  
dev.copy2pdf(file="vie_ita_sa7_FS_SB_double.pdf")

quartz()
good=ggplot(data=snp_sa7) +
     stat_binhex(data=snp_sa7[snp_sa7$DD_i == 1,],aes(x=aRPB, y=ATD),binwidth=c(0.5,1)) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs in DGRP and DPGP")
bad=ggplot(data=snp_sa7) +
     stat_binhex(data=snp_sa7[snp_sa7$DD_c == 0,],aes(x=aRPB, y=ATD),binwidth=c(0.5,1)) + scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs not in DGRP or DPGP")
grid.arrange(good,bad,ncol=1)     
dev.copy2pdf(file="vie_ita_sa7_aRPB_ATD_double.pdf")
rm("good","bad")
gc()

# for snp_sa7$DD_i == 1
good_snps= snp_sa7$DD_i == 1
# for snp_sa7$DD_c == 0
bad_snps= snp_sa7$DD_c == 0
num_good=length(snp_sa7$FS[good_snps])
num_bad= length(snp_sa7$FS[bad_snps])
num_tot= length(snp_sa7$FS)
num_tot #4703410
num_good #3234890
num_bad #861256
length(snp_sa7$DD_i[( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 )]) # 61915
length(snp_sa7$DD_i[ ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 ) & bad_snps]) # 44778
length(snp_sa7$DD_i[( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 ) & good_snps]) # 9160
length(snp_sa7$DD_i[(snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) & bad_snps]) #167585
length(snp_sa7$DD_i[(snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) & good_snps]) #36256
length(snp_sa7$DD_i[(snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5)]) #238921
length(snp_sa7$DD_i[(snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) | ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 )  ]) # 255872
length(snp_sa7$DD_i[((snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) | ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 )) & good_snps]) # 40333
length(snp_sa7$DD_i[((snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) | ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 ) ) & bad_snps]) # 178041
length(snp_sa7$DD_i[! ((snp_sa7$ATD <= 8 & snp_sa7$aRPB >= 3.5) | ( snp_sa7$SB <= 0.01 & snp_sa7$FS >= 30 ))]) # 4447539

# create ROC like curves for retrieval SNPs
quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
fs_vals=c(5,10,20,30,60,100,150,200,250)
good_vals = sapply(fs_vals,function (x) length(snp_sa7$FS[good_snps & snp_sa7$FS >= x] ))
bad_vals =  sapply(fs_vals,function (x) length(snp_sa7$FS[bad_snps & snp_sa7$FS >= x] ))
plot(fs_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.2),xlab="FS",ylab="fraction filtered")
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=1)
SB_vals=c(0.2,0.1,0.05,0.01)
ltypes=c(2,3,4,5)
for(i in 1:4){
good_vals = sapply(fs_vals,function (x) length(snp_sa7$FS[good_snps & snp_sa7$FS >= x & snp_sa7$SB <= SB_vals[i]] ))
bad_vals =  sapply(fs_vals,function (x) length(snp_sa7$FS[bad_snps & snp_sa7$FS >= x & snp_sa7$SB <=  SB_vals[i]] ))
lines(fs_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","SB <= 1.0","SB <= 0.2","SB <= 0.1","SB <= 0.05","SB <= 0.01"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,1,2,3,4,5))
dev.copy2pdf(file="vie_ita_sa7_FS_SB_SNPs_filtered.pdf")
quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
fs_vals=c(5,10,20,30,60,100,150,200,250)
good_vals = sapply(fs_vals,function (x) length(snp_sa7$FS[good_snps & snp_sa7$FS >= x] ))
bad_vals =  sapply(fs_vals,function (x) length(snp_sa7$FS[bad_snps & snp_sa7$FS >= x] ))
plot(fs_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.2),xlab="FS",ylab="fraction filtered")
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=1)
bsb_vals=c(25,50,100,150)
ltypes=c(2,3,4,5)
for(i in 1:4){
good_vals = sapply(fs_vals,function (x) length(snp_sa7$FS[good_snps & snp_sa7$FS >= x & snp_sa7$BSB >= bsb_vals[i]] ))
bad_vals =  sapply(fs_vals,function (x) length(snp_sa7$FS[bad_snps & snp_sa7$FS >= x & snp_sa7$BSB >=  bsb_vals[i]] ))
lines(fs_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(fs_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.2,by=0.05),labels=round(seq(0,0.2,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","BSB >= 0.0","BSB >= 25","BSB >= 50","BSB >= 100","BSB >= 200"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,1,2,3,4,5))
dev.copy2pdf(file="vie_ita_sa7_FS_BSB_SNPs_filtered.pdf")

quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
aRPB_vals=c(0,1,2,3,3.5,4,4.5,5,7.5,10,15)
good_vals = sapply(aRPB_vals,function (x) length(snp_sa7$aRPB[good_snps & snp_sa7$aRPB >= x] ))
bad_vals =  sapply(aRPB_vals,function (x) length(snp_sa7$aRPB[bad_snps & snp_sa7$aRPB >= x] ))
plot(aRPB_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.25),xlab="aRPB",ylab="fraction filtered")
lines(aRPB_vals, bad_vals/num_bad,type="l", col="red", lty=1)
ATD_vals=c(15,10,7.5,5,2.5,1)
ltypes=c(2,3,4,5,6,2)
for(i in 1:4){
good_vals = sapply(aRPB_vals,function (x) length(snp_sa7$aRPB[good_snps & snp_sa7$aRPB >= x & snp_sa7$ATD <= ATD_vals[i]] ))
bad_vals =  sapply(aRPB_vals,function (x) length(snp_sa7$aRPB[bad_snps & snp_sa7$aRPB >= x & snp_sa7$ATD <=  ATD_vals[i]] ))
lines(aRPB_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(aRPB_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","ATD <= 15","ATD <= 10","ATD <= 7.5","ATD <= 5","ATD <= 2.5","ATD <= 1.0"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,2,3,4,5,6,2))
dev.copy2pdf(file="vie_ita_sa7_aRPB_ATD_SNPs_filtered.pdf")


quartz()
par(mar=c(5, 4, 4, 8) + 0.1)
TDB_vals=c(0,10,25,50,75,100)
good_vals = sapply(TDB_vals,function (x) length(snp_sa7$TDB[good_snps & snp_sa7$TDB >= x] ))
bad_vals =  sapply(TDB_vals,function (x) length(snp_sa7$TDB[bad_snps & snp_sa7$TDB >= x] ))
plot(TDB_vals, good_vals/num_good,type="l", col="blue", lty=1,ylim=c(0,0.25),xlab="TDB",ylab="fraction filtered")
lines(TDB_vals, bad_vals/num_bad,type="l", col="red", lty=1)
ATD_vals=c(15,10,7.5,5,2.5,1)
ltypes=c(2,3,4,5,6,2)
for(i in 1:4){
good_vals = sapply(TDB_vals,function (x) length(snp_sa7$TDB[good_snps & snp_sa7$TDB >= x & snp_sa7$ATD <= ATD_vals[i]] ))
bad_vals =  sapply(TDB_vals,function (x) length(snp_sa7$TDB[bad_snps & snp_sa7$TDB >= x & snp_sa7$ATD <=  ATD_vals[i]] ))
lines(TDB_vals, good_vals/num_good,type="l", col="blue", lty=ltypes[i])
lines(TDB_vals, bad_vals/num_bad,type="l", col="red", lty=ltypes[i])
}
axis(4,col.axis="red",at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_bad,digits=0))
axis(4,col.axis="blue",outer=T,at=seq(0,0.25,by=0.05),labels=round(seq(0,0.25,by=0.05)*num_good,digits=0),line=-3.5)
title(main="SNPs filtered out")
legend("topright",legend=c("in DGRP and DPGP","not in DGRP or DPGP","ATD <= 15","ATD <= 10","ATD <= 7.5","ATD <= 5","ATD <= 2.5","ATD <= 1.0"),col=c("blue","red","grey10","grey10","grey10","grey10","grey10","grey10"),lty=c(1,1,2,3,4,5,6,2))
rm("good_vals","bad_vals")
dev.copy2pdf(file="vie_ita_sa7_TDB_ATD_SNPs_filtered.pdf")
gc()


ram_tab=read.table("/Volumes/Temp/Lukas/12pop_NH_mc_5_mcov_5_all_snps_biases.tab",header=F,na.strings = c("NA","nan"))
head(ram_tab)
colnames(ram_tab)=c("CHR","BPS","ALL","FS","SB","BSB","RTD","ATD","TDB","aRPB")


library(ggplot2)
mybreaks=c(1,10,100,500,1000,5000)
xlimit=c(0,600)
quartz()
#par(mfrow=c(1,2))
ggplot(data=ram_tab) + xlim(xlimit) +
     stat_binhex(data=ram_tab,aes(x=FS, y=SB),bins=40) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs")

quartz()
ggplot(data=ram_tab) +
     stat_binhex(data=ram_tab,aes(x=aRPB, y=ATD),binwidth=c(0.5,1)) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs")

xlimit=c(0,60)

ylimit=c(-0.01,0.15)
mybreaks=c(1,50,250,500,1000)

quartz()
ggplot(data=ram_tab) + xlim(xlimit) + ylim(ylimit) +
     stat_binhex(data=ram_tab,aes(x=FS, y=SB),bins=20) +
    scale_fill_gradientn(colours=c("red","green","blue"),trans="log", name = "count",breaks= mybreaks, labels=mybreaks, na.value=NA) + ggtitle("SNPs")


gc()
