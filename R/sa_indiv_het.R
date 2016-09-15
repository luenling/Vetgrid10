inversions=read.table("/Volumes/Temp/Lukas/Inversions/inversion_breakpoints.bed",header=F)
inversions=read.table("/Volumes/vetgrid10/Inversions/inversion_breakpoints.bed",header=F)
colnames(inversions)=c("CHR","START","STOP","I_NAME")
inversions$START=inversions$START+1
inversions=inversions[order(inversions$CHR,(inversions$START+inversions$STOP)/2),]

setwd("/Volumes/Temp/Lukas/Data/SA_A7/BGI_135")
setwd("/Volumes/vetgrid10/Data/SA_A7/BGI_135")
sa_het=read.table("BGI_135.real.gatk.unif_genotyp_recal_99_het_b.sync", header=T)
#colnames(sa_het)=c("CHR","BPS","AL","tot_cal","homR","het","homA","Fi","fR","fA","HomExp","HetExp","fRt","fAt","Fit" )
#sa_het=sa_het[! (sa_het$fA == 0 | sa_het$fR == 0 |sa_het$fAt == 0 | sa_het$fRt == 0) ,]
sa_het=sa_het[! (sa_het$fAt <= 0.15 | sa_het$fRt <= 0.15) ,]
sa_het$het_fr = with(sa_het, het/tot_cal)
sa_het$He_S=with(sa_het,2*fR*fA)
sa_het$He_T=with(sa_het,2*fRt*fAt)
for (i in colnames(sa_het)[16:27]) {
     sa_het[,paste("H_",i,sep="")]=(sa_het[,i] == 0.5)*sa_het[,i]*2.0
 }


# calculate the dark and light indv. heterozygosity
sa_het$hetf_D=rowSums( (sa_het[,16:21] == 0.5)*sa_het[,16:21]*2,na.rm=T)/rowSums(! is.na(sa_het[,16:21]))
sa_het$hetf_L=rowSums( (sa_het[,22:27] == 0.5)*sa_het[,22:27]*2,na.rm=T)/rowSums(! is.na(sa_het[,22:27]))
sa_het$fR_D=( sa_het$hetf_D*0.5 + rowSums((sa_het[,16:21] == 0.0),na.rm=T)/rowSums(! is.na(sa_het[,16:21])))
sa_het$fR_D[sa_het$fR_D < 0.0001]=0.0001
sa_het$fR_D[sa_het$fR_D > 0.9999]=0.9999
sa_het$fR_L=( sa_het$hetf_L*0.5 + rowSums((sa_het[,22:27] == 0.0),na.rm=T)/rowSums(! is.na(sa_het[,22:27])))
sa_het$fR_L[sa_het$fR_L < 0.0001]=0.0001
sa_het$fR_L[sa_het$fR_L > 0.9999]=0.9999

sa_het$Fs_D=with(sa_het,1-hetf_D/(2*fR_D*(1-fR_D)))
sa_het$Fs_L=with(sa_het,1-hetf_L/(2*fR_L*(1-fR_L)))
sa_het$Ft_D=with(sa_het,1-hetf_D/(2*fRt*fAt))
sa_het$Ft_L=with(sa_het,1-hetf_L/(2*fRt*fAt))
# calculate binned fsts for all chromosomes
max_bps=max(sa_het$BPS)
#28994946
#intervall vector
int_breaks=seq(1,max_bps,by=100000)
sa_het$Win100K=cut(sa_het$BPS,breaks=int_breaks)
library(foreach)
cols=c("BPS","Fit","Fi","Fs_D","Fs_L","Ft_D","Ft_L","H_vd1","H_vd2","H_vd3","H_vd4","H_vd5","H_vd6","H_vvl1","H_vvl2","H_vvl3","H_vvl4","H_vvl5","H_vvl6","He_S","He_T")
fit_win100K=sa_het[0,c("CHR",cols)]
for(i in levels(sa_het$CHR)){
 aa= aggregate(x=subset(sa_het, CHR==i)[,cols], by=list(subset(sa_het, CHR ==i)$Win100K), FUN=mean, na.rm=TRUE, na.action="na.pass" )
 chrom = factor(i,levels(sa_het$CHR))
 aa$CHR=rep(chrom,nrow(aa))
 fit_win100K=rbind(fit_win100K,aa[,-1])
}

for(i in grep("H_",colnames(fit_win100K),value=T)){
    fit_win100K[,gsub("H_","Fs_",i)]=1-fit_win100K[,i]/fit_win100K$He_S
    fit_win100K[,gsub("H_","Ft_",i)]=1-fit_win100K[,i]/fit_win100K$He_T    
}

# only Fit
library(boot)
cols=c("BPS","Fit")
fit_win100K_m_ci<- data.frame(matrix(vector(), 0, 5, dimnames=list(c(), c("CHR",cols,"cia","cib"))), stringsAsFactors=F)
for(i in levels(sa_het$CHR)){
  data=subset(sa_het, CHR==i)[,cols]
  wins=subset(sa_het, CHR ==i)$Win100K  
  aa= aggregate(x=data, by=list(wins), FUN=mean, na.rm=TRUE, na.action="na.pass" )
  ab= aggregate(data$Fit, by=list(wins), FUN= get.boot.ci )
  chrom = factor(i,levels(sa_het$CHR))
  aa$CHR=rep(chrom,nrow(aa))
  alla=cbind(aa[,-1],ab[,-1])
  fit_win100K_m_ci=rbind(fit_win100K_m_ci,alla)
  
}
fit_win100K_m_sd <- data.frame(matrix(vector(), 0, 4, dimnames=list(c(), c("CHR",cols,"sd"))), stringsAsFactors=F)
for(i in levels(sa_het$CHR)){
  data=subset(sa_het, CHR==i)[,cols]
  wins=subset(sa_het, CHR ==i)$Win100K  
  aa= aggregate(x=data, by=list(wins), FUN=mean, na.rm=TRUE, na.action="na.pass" )
  ab= aggregate(data$Fit, by=list(wins), FUN= sd )
  chrom = factor(i,levels(sa_het$CHR))
  aa$CHR=rep(chrom,nrow(aa))
  alla=cbind(aa[,-1],ab[,-1])
  colnames(alla)=c(colnames(aa[,-1]),"sd")
  fit_win100K_m_sd=rbind(fit_win100K_m_sd,alla)
  
}




get.boot.ci <- function(x) {
  #bps=mean(x[,1],na.rm=TRUE)
  if (length(unique(x)) == 1){
    b <-c(mean(x),mean(x))
  }
  else{
    bfit=boot(x,statistic=boot.mean,R=1000)
    a <- boot.ci(bfit,conf=0.95,type="basic")
    b <- c(a$basic[4:5])
  }
  names(b)=c("cia","cib")
  return(b)
}

boot.mean <- function(data,indices){
  d <- data[indices]   
  return(mean(d, na.rm = TRUE))
}


# calculate binned fsts for all chromosomes
max_bps=max(sa_het$BPS)
#28994946
#intervall vector
int_breaks=seq(1,max_bps,by=25000)
sa_het$Win25K=cut(sa_het$BPS,breaks=int_breaks)
library(foreach)
cols=c("BPS","Fit","Fi","Fs_D","Fs_L","Ft_D","Ft_L","H_vd1","H_vd2","H_vd3","H_vd4","H_vd5","H_vd6","H_vvl1","H_vvl2","H_vvl3","H_vvl4","H_vvl5","H_vvl6","He_S","He_T")
fit_win25K=sa_het[0,c("CHR",cols)]
for(i in levels(sa_het$CHR)){
 aa= aggregate(x=subset(sa_het, CHR==i)[,cols], by=list(subset(sa_het, CHR ==i)$Win25K), FUN=mean, na.rm=TRUE, na.action="na.pass" )
 chrom = factor(i,levels(sa_het$CHR))
 aa$CHR=rep(chrom,nrow(aa))
 fit_win25K=rbind(fit_win25K,aa[,-1])
}

# test for signifcance of differences:
chroms=c("2L","2R","3L","3R","X")
fit.win25K = droplevels(fit_win25K[fit_win25K$CHR %in% chroms,])
fit.win100K = droplevels(fit_win100K[fit_win100K$CHR %in% chroms,])
library("boot")
boot.median <- function(data,indices){
    d <- data[indices]   
    return(median(d, na.rm = TRUE))
}
bfit.sa.100K=boot(fit.win100K$Fit,statistic=boot.median ,R=10000)
a=list()
for(chrom in chroms){
  bfit.sa=boot(fit.win100K$Fit[fit.win100K$CHR == chrom],statistic=boot.median ,R=1000)
  a[[chrom]]=boot.ci(bfit.sa) 
}


x11()
plot(bfit.sa.100K)
boot.ci(bfit.sa.100K)  
## Intervals : 
## Level      Normal              Basic         
## 95%   ( 0.3381,  0.3650 )   ( 0.3401,  0.3643 )  
## Level     Percentile            BCa          
## 95%   ( 0.3368,  0.3610 )   ( 0.3368,  0.3610 )  

kruskal.test(Fit ~ CHR, data = fit.win25K)
a=kruskal.test(Fit ~ CHR, data = fit.win100K)
library("dunn.test")
dunn25K=dunn.test(fit.win25K$Fit,fit.win25K$CHR, kw=TRUE, method="bonferroni")
dunn100K=dunn.test(fit.win100K$Fit,fit.win100K$CHR, kw=TRUE, method="bh")
aggregate(fit.win100K$Fit,by=list(chrom=fit.win100K$CHR),FUN=median)
f_avg_chrom=aggregate(x=fit_win100K[,c("Fit","Fi","Fs_D","Fs_L","Ft_D","Ft_L",grep("^Ft_v",colnames(fit_win100K),value=T) )], by=list(fit_win100K$CHR), FUN=mean, na.rm=TRUE, na.action="na.pass")
h_avg_chrom=aggregate(x=fit_win100K[,c(grep("^H_v",colnames(fit_win100K),value=T) )], by=list(fit_win100K$CHR), FUN=mean, na.rm=TRUE, na.action="na.pass")
colnames(f_avg_chrom)[1]="CHR"
colnames(h_avg_chrom)[1]="CHR"
write.table(f_avg_chrom,file="inbreeding_F_averages_chrom_sa.txt",sep="\t",quote=F,row.names=F)
write.table(h_avg_chrom,file="heterozygosity_H_averages_chrom_sa.txt",sep="\t",quote=F,row.names=F)

quartz(width=10,height=10)
chroms=c("X","2L","2R","3L","3R")
ylimit=c(-0.15,1.0)
layout(matrix(c(1,6,2,3,4,5),byrow=T, ncol = 2), widths = 1, heights = 1, respect = FALSE)
#par(mar = c(4.1, 4.1, 4.1, 2.1))
for (chrom in chroms){
  par(mar = c(4.1, 4.1, 4.1, 2.1))
  aa=subset(fit_win100K,CHR==chrom)
  #fit_col=c("fit_eu","fit_es")
  #ylimit=c(min(aa$Fit),max(aa$Fit))
  plot(aa$BPS, aa$Fit, pch=".", col="red", cex=1, ylim=ylimit, xlab="base position", main=chrom, ylab=bquote(F[IT]))
  lines(aa$BPS,aa$Fit, col='red', lwd=1)
  lines(aa$BPS,aa$Ft_L, col='cyan', lwd=1,lty=1)
  lines(aa$BPS,aa$Ft_D, col='darkblue', lwd=1,lty=1)
  #lines(aa$BPS,predict(loess(aa$Fit~aa$BPS,span=0.01),newdata=aa$BPS), col='red', lwd=1)
  #points(aa$BPS, aa$Fi, pch="o", col="red", cex=0.5)
  ## lines(aa$BPS,predict(loess(aa$Fi~aa$BPS,span=0.01),newdata=aa$BPS), col='red', lwd=1)
  ## lines(aa$BPS,predict(loess(aa$Fs_L~aa$BPS,span=0.01,na.action = na.exclude),newdata=aa$BPS), col='cyan', lwd=1)
  ## lines(aa$BPS,predict(loess(aa$Fs_D~aa$BPS,span=0.01,na.action = na.exclude),newdata=aa$BPS), col='darkblue', lwd=1)
  #points(aa$BPS, aa$Fs_L, col='cyan',pch="x")
  #points(aa$BPS, aa$Fs_D, col='darkblue',pch="x")
  inv_sub=subset(inversions,CHR==chrom)
  i=1
  while(i<=dim(inv_sub)[1]){
    inkr=i*0.08
    segments(inv_sub[i,]$START,-0.15+inkr,inv_sub[i,]$STOP,-0.15+inkr)
    text((inv_sub[i,]$START+inv_sub[i,]$STOP)/2.0,-0.125+inkr,labels=inv_sub[i,]$I_NAME)
    i = i+1
  }
  abline(h=mean(aa$Fit),col="red",lty=2)
}
plot(0:2, 0:2, pch = 1, lty = 1, ylim=c(0,2), type = "n", axes = FALSE, ann = FALSE,main="South Africa 2012, 12 Isofemale lines (6 light, 6 dark)")#
legend(0.75,1.15, legend=c(expression(F[IT]),expression(F[IT] ~ Dark),expression(F[IT] ~ Light)), lty=c(1,1,1,1),col=c("red","green","darkblue","cyan"),xpd=TRUE)
text(1,2.2,labels="South Africa 2012, 12 Isofemale lines (6 light, 6 dark)",cex=1.5,font=2,xpd=T)
#dev.copy2pdf(file="fst_bases.pdf")
dev.copy2pdf(file="inbreeding_FIT_sa_2010_a7_with_LD.pdf")

quartz(width=10,height=10)
chroms=c("X","2L","2R","3L","3R")
ylimit=c(-0.05,1.0)
layout(matrix(c(1,6,2,3,4,5),byrow=T, ncol = 2), widths = 1, heights = 1, respect = FALSE)
#par(mar = c(4.1, 4.1, 4.1, 2.1))
for (chrom in chroms){
  par(mar = c(4.1, 4.1, 4.1, 2.1))
  aa=subset(fit_win100K_m_sd,CHR==chrom)
  #fit_col=c("fit_eu","fit_es")
  #ylimit=c(min(aa$Fit),max(aa$Fit))
  plot(aa$BPS, aa$Fit, col="black", type="l",lwd=2,lty=1, ylim=ylimit, xlab="base position", main=chrom, ylab=bquote(F[IT]))
  lines(aa$BPS,aa$Fit+aa$sd, col='black', lwd=1,lty=3)
  lines(aa$BPS,aa$Fit-aa$sd, col='black', lwd=1,lty=3)
  inv_sub=subset(inversions,CHR==chrom)
  i=1
  while(i<=dim(inv_sub)[1]){
    inkr=i*0.08
    segments(inv_sub[i,]$START,-0.15+inkr,inv_sub[i,]$STOP,-0.15+inkr)
    text((inv_sub[i,]$START+inv_sub[i,]$STOP)/2.0,-0.125+inkr,labels=inv_sub[i,]$I_NAME)
    i = i+1
  }
  abline(h=mean(aa$Fit),col="red",lty=2)
}
plot(0:2, 0:2, pch = 1, lty = 1, ylim=c(0,2), type = "n", axes = FALSE, ann = FALSE,main="South Africa 2012 (12 individuals)") #
#legend(0.75,1.15, legend=c(expression(F[IT]),expression(F[IT] ~ Dark),expression(F[IT] ~ Light)), lty=c(1,1,1,1),col=c("red","green","darkblue","cyan"),xpd=TRUE)
text(1,2.2,labels="South Africa 2012, 12 Isofemale lines",cex=1.5,font=2,xpd=T)
#dev.copy2pdf(file="fst_bases.pdf")
dev.copy2pdf(file="inbreeding_FIT_sa_2010_a7.pdf")

aggregate(fit_win100K,by=,FUN=mean,na.action=)
# plot single individuals
quartz(width=10,height=10)
chroms=c("X","2L","2R","3L","3R")
ylimit=c(-0.5,1.0)
ltypes=c(1,1,1,1,1,1,4,4,4,4,4,4)
layout(matrix(c(1,6,2,3,4,5),byrow=T, ncol = 2), widths = 1, heights = 1, respect = FALSE)
indivs=c("Ft_vd1","Ft_vd2", "Ft_vd3", "Ft_vd4", "Ft_vd5", "Ft_vd6", "Ft_vvl1","Ft_vvl2", "Ft_vvl3","Ft_vvl4","Ft_vvl5","Ft_vvl6")
cols=c("green","red","blue","cyan","pink","brown","purple","orange","cornflowerblue","lightgreen","firebrick1")
for (chrom in chroms){
    par(mar = c(4.1, 4.1, 4.1, 2.1))
    aa=subset(fit_win100K,CHR==chrom)
    plot(aa$BPS,predict(loess(aa[,indivs[1]]~aa$BPS,span=0.01),newdata=aa$BPS) , col="green",type="l", ylim=ylimit, xlab="base position", main=chrom, lty=ltypes[1], ylab=bquote(F[I]))
    aa=subset(fit_win100K,CHR==chrom)
    for(i in 2:length(indivs)){
        lines(aa$BPS,predict(loess(aa[,indivs[i]]~aa$BPS,span=0.01),newdata=aa$BPS), col=cols[i], lty=ltypes[i], lwd=1)
#        lines(aa$BPS,aa[,i], col='red', lwd=1)
    }
    inv_sub=subset(inversions,CHR==chrom)
    i=1
    while(i<=dim(inv_sub)[1]){
        inkr=i*0.08
        segments(inv_sub[i,]$START,-0.45+inkr,inv_sub[i,]$STOP,-0.45+inkr)
        text((inv_sub[i,]$START+inv_sub[i,]$STOP)/2.0,-0.425+inkr,labels=inv_sub[i,]$I_NAME)
        i = i+1
    }
    
}
plot(0:2, 0:2, pch = 1, lty = 1, ylim=c(0,2), type = "n", axes = FALSE, ann = FALSE,main="South Africa 2012, 12 Isofemale lines (6 light, 6 dark)")
#legend(0.75,1.15, legend=c(expression(F[ST] ~ Vienna ),expression(F[ST] ~  Italy),expression(F[ST] ~ Vienna ~ and ~ South ~ Africa),expression(F[ST] ~ Italy ~ and ~ South ~ Africa), expression(F[ST] ~ Vienna ~ and ~ Italy),expression(F[ST] ~ Europe ~ and ~ South ~ Africa)), pch=c(".",".",".",".","o","x"),col=c("darkgoldenrod1","darkorange2","cyan","darkcyan","green","red"),xpd=TRUE)
legend(0.75,1.15, legend=c("ind vd 1","ind vd 2","ind vd 3","ind vd 4","ind vd 5","ind vd 6","ind vvl 1","ind vvl 2","ind vvl 3","ind vvl 4","ind vvl 5","ind vvl 6"), lty=ltypes,col=cols,xpd=TRUE)
text(1,2.2,labels="South Africa 2012, 12 Isofemale lines (6 light, 6 dark)",cex=1.5,font=2,xpd=T)
text(1,2.0,labels="100 kbp windows",cex=1.5,font=1,xpd=T)
dev.copy2pdf(file="inbreeding_sa_2010_a7_individuals.pdf")


# plot single individuals
quartz(width=10,height=10)
chroms=c("X","2L","2R","3L","3R")
ylimit=c(0.0,1.0)
ltypes=c(1,1,1,1,1,1,4,4,4,4,4,4)
layout(matrix(c(1,6,2,3,4,5),byrow=T, ncol = 2), widths = 1, heights = 1, respect = FALSE)
indivs=c("H_vd1","H_vd2", "H_vd3", "H_vd4", "H_vd5", "H_vd6", "H_vvl1","H_vvl2", "H_vvl3","H_vvl4","H_vvl5","H_vvl6")
cols=c("green","red","blue","cyan","pink","brown","purple","orange","cornflowerblue","lightgreen","firebrick1")
for (chrom in chroms){
    par(mar = c(4.1, 4.1, 4.1, 2.1))
    aa=subset(fit_win25K,CHR==chrom)
    plot(aa$BPS,predict(loess(aa[,indivs[1]]~aa$BPS,span=0.01),newdata=aa$BPS) , col="green",type="l", ylim=ylimit, xlab="base position", main=chrom,lty=ltypes[1], ylab=bquote(H[avg]))
    aa=subset(fit_win100K,CHR==chrom)
    for(i in 2:length(indivs)){
        lines(aa$BPS,predict(loess(aa[,indivs[i]]~aa$BPS,span=0.01),newdata=aa$BPS), col=cols[i],lty=ltypes[i], lwd=1)
#        lines(aa$BPS,aa[,i], col='red', lwd=1)
    }
    inv_sub=subset(inversions,CHR==chrom)
    i=1
    while(i<=dim(inv_sub)[1]){
        inkr=i*0.05
        segments(inv_sub[i,]$START,0.75+inkr,inv_sub[i,]$STOP,0.75+inkr)
        text((inv_sub[i,]$START+inv_sub[i,]$STOP)/2.0,0.775+inkr,labels=inv_sub[i,]$I_NAME)
        i = i+1
    }
    
}
plot(0:2, 0:2, pch = 1, lty = 1, ylim=c(0,2), type = "n", axes = FALSE, ann = FALSE,main="South Africa 2012, 12 Isofemale lines (6 light, 6 dark)")
#legend(0.75,1.15, legend=c(expression(F[ST] ~ Vienna ),expression(F[ST] ~  Italy),expression(F[ST] ~ Vienna ~ and ~ South ~ Africa),expression(F[ST] ~ Italy ~ and ~ South ~ Africa), expression(F[ST] ~ Vienna ~ and ~ Italy),expression(F[ST] ~ Europe ~ and ~ South ~ Africa)), pch=c(".",".",".",".","o","x"),col=c("darkgoldenrod1","darkorange2","cyan","darkcyan","green","red"),xpd=TRUE)
text(1,2.2,labels="South Africa 2012, 12 Isofemale lines (6 light, 6 dark)",cex=1.5,font=2,xpd=T)
legend(0.75,1.15, legend=c("ind vd 1","ind vd 2","ind vd 3","ind vd 4","ind vd 5","ind vd 6","ind vvl 1","ind vvl 2","ind vvl 3","ind vvl 4","ind vvl 5","ind vvl 6"), lty=ltypes,col=cols,xpd=TRUE)
text(1,2,labels="25 Kbp windows",cex=1.5,font=1,xpd=T)
dev.copy2pdf(file="heterozygosity_sa_2010_a7_individuals.pdf")

# plot average hetreozygosity
quartz(width=10,height=10)
chroms=c("X","2L","2R","3L","3R")
ylimit=c(0.0,1.0)
layout(matrix(c(1,6,2,3,4,5),byrow=T, ncol = 2), widths = 1, heights = 1, respect = FALSE)
indivs=c("H_vd1","H_vd2", "H_vd3", "H_vd4", "H_vd5", "H_vd6", "H_vvl1","H_vvl2", "H_vvl3","H_vvl4","H_vvl5","H_vvl6")
cols=c("green","red","blue","cyan","pink","brown","purple","orange","cornflowerblue","lightgreen","firebrick1")
for (chrom in chroms){
    par(mar = c(4.1, 4.1, 4.1, 2.1))
    aa=subset(fit_win25K,CHR==chrom)
    plot(aa$BPS,predict(loess(rowMeans(aa[,indivs],na.rm=T)~aa$BPS,span=0.01),newdata=aa$BPS) , col="red",type="l", ylim=ylimit, xlab="base position", main=chrom, ylab=bquote(H[avg]))
    lines(aa$BPS,predict(loess(rowMeans(aa[,indivs[1:6]],na.rm=T)~aa$BPS,span=0.01),newdata=aa$BPS), col="darkblue", lwd=1)
    lines(aa$BPS,predict(loess(rowMeans(aa[,indivs[7:12]],na.rm=T)~aa$BPS,span=0.01),newdata=aa$BPS), col="cyan", lwd=1)
    inv_sub=subset(inversions,CHR==chrom)
    i=1
    while(i<=dim(inv_sub)[1]){
        inkr=i*0.05
        segments(inv_sub[i,]$START,0.75+inkr,inv_sub[i,]$STOP,0.75+inkr)
        text((inv_sub[i,]$START+inv_sub[i,]$STOP)/2.0,0.775+inkr,labels=inv_sub[i,]$I_NAME)
        i = i+1
    }
    
}
plot(0:2, 0:2, pch = 1, lty = 1, ylim=c(0,2), type = "n", axes = FALSE, ann = FALSE,main="South Africa 2012, 12 Isofemale lines (6 light, 6 dark)")
#legend(0.75,1.15, legend=c(expression(F[ST] ~ Vienna ),expression(F[ST] ~  Italy),expression(F[ST] ~ Vienna ~ and ~ South ~ Africa),expression(F[ST] ~ Italy ~ and ~ South ~ Africa), expression(F[ST] ~ Vienna ~ and ~ Italy),expression(F[ST] ~ Europe ~ and ~ South ~ Africa)), pch=c(".",".",".",".","o","x"),col=c("darkgoldenrod1","darkorange2","cyan","darkcyan","green","red"),xpd=TRUE)
text(1,2.2,labels="South Africa 2012, 12 Isofemale lines (6 light, 6 dark)",cex=1.5,font=2,xpd=T)
text(1,2,labels="25 Kbp windows",cex=1.5,font=1,xpd=T)
legend(0.75,1.15, legend=c(expression(H~all),expression(H~Dark),expression(H ~ Light)), lty=c(1,1,1),col=c("red","darkblue","cyan"),xpd=TRUE)
dev.copy2pdf(file="heterozygosity_sa_2010_a7_averages.pdf")


setwd("/Volumes/Temp/Lukas/Data/Vienna_2011/Diploids")
ita_het=read.table("/Volumes/Temp/Lukas/Data/Italy_2011/Diploids/BGI_106_snps_unified_gt_recal99.het_b_tab.sync", header=T)
#colnames(ita_het)=c("CHR","BPS","AL","tot_cal","homR","het","homA","Fi","fR","fA","HomExp","HetExp","fRt","fAt","Fit" )
#ita_het=ita_het[! (ita_het$fA == 0 | ita_het$fR == 0 |ita_het$fAt == 0 | ita_het$fRt == 0) ,]
ita_het=ita_het[! (ita_het$fAt <= 0.15 | ita_het$fRt <= 0.15) ,]
ita_het$het_fr = with(ita_het, het/tot_cal)
ita_het$He_S=with(ita_het,2*fR*fA)
ita_het$He_T=with(ita_het,2*fRt*fAt)
ita_het=ita_het[! is.na(ita_het$BPS),]
for (i in grep("ind_",colnames(ita_het),value=T)) {
     ita_het[,paste("H_",i,sep="")]=(ita_het[,i] == 0.5)*ita_het[,i]*2.0
 }

# calculate binned fsts for all chromosomes
max_bps=max(ita_het$BPS)
#28994946
#intervall vector
int_breaks=seq(1,max_bps,by=100000)
ita_het$Win100K=cut(ita_het$BPS,breaks=int_breaks)
library(foreach)
cols=c("BPS","Fit","Fi",grep("^H_ind_",colnames(ita_het),value=T),"He_S","He_T")
fit_win100K=ita_het[0,c("CHR",cols)]
for(i in levels(ita_het$CHR)){
 aa= aggregate(x=subset(ita_het, CHR==i)[,cols], by=list(subset(ita_het, CHR ==i)$Win100K), FUN=mean, na.rm=TRUE, na.action="na.pass" )
 chrom = factor(i,levels(ita_het$CHR))
 aa$CHR=rep(chrom,nrow(aa))
 fit_win100K=rbind(fit_win100K,aa[,-1])
}

for(i in grep("H_",colnames(fit_win100K),value=T)){
    fit_win100K[,gsub("H_","Fs_",i)]=1-fit_win100K[,i]/fit_win100K$He_S
    fit_win100K[,gsub("H_","Ft_",i)]=1-fit_win100K[,i]/fit_win100K$He_T    
}

chroms=c("2L","2R","3L","3R","X")
fit.win100K_eu = droplevels(fit_win100K[fit_win100K$CHR %in% chroms,])
library("boot")
boot.median <- function(data,indices){
    d <- data[indices]   
    return(median(d, na.rm = TRUE))
}
bfit.eu.100K=boot(fit.win100K_eu$Fit,statistic=boot.median ,R=10000)
x11()
plot(bfit.eu.100K)
boot.ci(bfit.eu.100K)
## Intervals : 
## Level      Normal              Basic         
## 95%   (-0.0076,  0.0026 )   (-0.0069,  0.0029 )  
## Level     Percentile            BCa          
## 95%   (-0.0087,  0.0010 )   (-0.0088,  0.0008 )  


kruskal.test(Fit ~ CHR, data = fit.win100K_eu)
library("dunn.test")
dunn100K=dunn.test(fit.win100K_eu$Fit,fit.win100K_eu$CHR, kw=TRUE, method="bh")
aggregate(fit.win100K_eu$Fit,by=list(chrom=fit.win100K_eu$CHR),FUN=median)

quartz(width=10,height=10)
chroms=c("X","2L","2R","3L","3R")
ylimit=c(-0.25,0.75)
layout(matrix(c(1,6,2,3,4,5),byrow=T, ncol = 2), widths = 1, heights = 1, respect = FALSE)
#par(mar = c(4.1, 4.1, 4.1, 2.1))
for (chrom in chroms){
    par(mar = c(4.1, 4.1, 4.1, 2.1))
    aa=subset(fit_win100K,CHR==chrom)
                                        #fit_col=c("fit_eu","fit_es")
                                        #ylimit=c(min(aa$Fit),max(aa$Fit))
    plot(aa$BPS, aa$Fit, pch="o", col="green", cex=0.5, ylim=ylimit, xlab="base position", main=chrom, ylab=bquote(F[ST]~or~F[IS]))
    lines(aa$BPS,predict(loess(aa$Fit~aa$BPS,span=0.01),newdata=aa$BPS), col='green', lwd=1)
    points(aa$BPS, aa$Fi, pch="o", col="blue", cex=0.5)
    lines(aa$BPS,predict(loess(aa$Fi~aa$BPS,span=0.01),newdata=aa$BPS), col='blue', lwd=1)
    inv_sub=subset(inversions,CHR==chrom)
    i=1
    while(i<=dim(inv_sub)[1]){
        inkr=i*0.05
        segments(inv_sub[i,]$START,0.5+inkr,inv_sub[i,]$STOP,0.5+inkr)
        text((inv_sub[i,]$START+inv_sub[i,]$STOP)/2.0,0.52+inkr,labels=inv_sub[i,]$I_NAME)
        i = i+1
    }
    abline(h=mean(aa$Fit),col="red",lty=2)
}
plot(0:2, 0:2, pch = 1, lty = 1, ylim=c(0,2), type = "n", axes = FALSE, ann = FALSE,main="Italy, Bolzano (2011), 12 individuals chill coma")
legend(0.75,1.15, legend=c(expression(F[IS]),expression(F[IT])), lty=c(1,1),col=c("blue","green"),xpd=TRUE)
text(1,2.2,labels="Italy, Bolzano (2011) CCR, 12 individuals",cex=1.5,font=2,xpd=T)
#dev.copy2pdf(file="fst_bases.pdf")
dev.copy2pdf(file="inbreeding_ita_2011_chill_coma.pdf")


# plot single individuals
quartz(width=10,height=10)
chroms=c("X","2L","2R","3L","3R")
ylimit=c(-0.5,1.0)
layout(matrix(c(1,6,2,3,4,5),byrow=T, ncol = 2), widths = 1, heights = 1, respect = FALSE)
indivs=grep("^Ft_ind_",colnames(fit_win100K),value=T)
cols=c("green","red","blue","cyan","pink","brown","purple","orange","cornflowerblue","lightgreen","firebrick1")
for (chrom in chroms){
    par(mar = c(4.1, 4.1, 4.1, 2.1))
    aa=subset(fit_win100K,CHR==chrom)
    plot(aa$BPS,predict(loess(aa[,indivs[1]]~aa$BPS,span=0.01),newdata=aa$BPS) , col="green",type="l", ylim=ylimit, xlab="base position", main=chrom, ylab=bquote(F[I]))
    aa=subset(fit_win100K,CHR==chrom)
    for(i in 2:length(indivs)){
        lines(aa$BPS,predict(loess(aa[,indivs[i]]~aa$BPS,span=0.01),newdata=aa$BPS), col=cols[i], lwd=1)
#        lines(aa$BPS,aa[,i], col='red', lwd=1)
    }
    inv_sub=subset(inversions,CHR==chrom)
    i=1
    while(i<=dim(inv_sub)[1]){
        inkr=i*0.05
        segments(inv_sub[i,]$START,-0.25+inkr,inv_sub[i,]$STOP,-0.25+inkr)
        text((inv_sub[i,]$START+inv_sub[i,]$STOP)/2.0,-0.2+inkr,labels=inv_sub[i,]$I_NAME)
        i = i+1
    }
}
plot(0:2, 0:2, pch = 1, lty = 1, ylim=c(0,2), type = "n", axes = FALSE, ann = FALSE,main="Italy, Bolzano (2011) CCR, 12 individuals")
#legend(0.75,1.15, legend=c(expression(F[ST] ~ Vienna ),expression(F[ST] ~  Italy),expression(F[ST] ~ Vienna ~ and ~ South ~ Africa),expression(F[ST] ~ Italy ~ and ~ South ~ Africa), expression(F[ST] ~ Vienna ~ and ~ Italy),expression(F[ST] ~ Europe ~ and ~ South ~ Africa)), pch=c(".",".",".",".","o","x"),col=c("darkgoldenrod1","darkorange2","cyan","darkcyan","green","red"),xpd=TRUE)
legend(0.75,1.15, legend=c("ind 1","ind 2","ind 3","ind 4","ind 5","ind 6","ind 7","ind 8","ind 9","ind 10","ind 11","ind 12"), lty=c(1,1,1,1,1,1),col=cols,xpd=TRUE)
text(1,2.2,labels="Italy, Bolzano (2011) CCR, 12 individuals",cex=1.5,font=2,xpd=T)
dev.copy2pdf(file="inbreeding_ita_2011_ccr_individuals.pdf")

