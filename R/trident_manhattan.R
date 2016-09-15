library(hexbin)
qqplot_hex <- function(nl_oP,nl_dP1,title,alpha=FALSE,beta=" by coverage"){
	s_nl_oP <- sort(nl_oP,decreasing=F)
	s_nl_dP <- sort(nl_dP1,decreasing=F)
	bin <- hexbin(s_nl_dP,s_nl_oP,xbins=400,xlab=expression(paste("-log(P",scriptstyle(null),")")),ylab=expression(paste("-log(P",scriptstyle(obs),")")))
	if (alpha){
		main_tit=bquote(paste(.(title)," ",alpha,"=",.(alpha)," ", beta,.(beta)))
		}
	else{
		main_tit=bquote(.(title))
	}
	pp<-plot(bin,legend=FALSE,style="constant.col",main=main_tit)
	
	hvp=hexViewport(bin)
	hexVP.abline(hvp,0,1,col="red",lty=2)
}

source("/Volumes/Temp/Lukas/Tools/Scripts/R/manhattan_plot_all_chromosomes_function_qqplot.R")
setwd("/Volumes/Temp/Lukas/Data/Trident/Realigned/Polymorphisms/CMH")
data <- read.table("males_vie_pTL1_pTL2_pTD1_pTD2_q20_filt_mc1.5_chroms_only_mct8_mcv25.gwas",header=TRUE)
null_pV=scan("FDR/males_vie_pTL1_pTL2_pTD1_pTD2_q20_filt_mc1.5_chroms_only_mct8_mcv25.cmhout_a_24.0_beta_r1_17.78_r2_26.89_nullP.out")
null_pV=scan("FDR/males_vie_pTL1_pTL2_pTD1_pTD2_q20_filt_mc1.5_chroms_only_mct8_mcv25.cmhout_44.0_nullP.out")

data$ID <- paste(data$CHR,data$BP,sep="_")
quartz(height=6,width=12)
par(bg="white")
chroms=c("X","2L","2R","3L","3R","4")#,"XHet","2LHet","2RHet","3LHet","3RHet")
manhattan(data,colors=c("black","slategrey"),limitchromosomes=chroms,main="Vienna males trident (2010)", suggestiveline=c(-log10(3.126176223804206e-08)))
dev.print(png,width=800,file="vie2010_males_trident_manhattan.png")
chr_offset=get_offset_manh(data, ymax="max", limitchromosomes=chroms)
inv3r_p=c(12051982,20591144)
inv3r_p_coord=inv3r_p+chr_offset['3R']
lines(inv3r_p_coord,c(0,0),col=c("red"))
lines(inv3r_p_coord,c(-0.1,-0.1),lwd=3,col=c("red"))
lines(inv3r_p_coord,c(14.25,14.25),lwd=3,col=c("red"))
quartz()
nl_oP = -1 * log10(data$P)
nl_unif = -1 *log10(runif(length(nl_oP)))
nl_null=-1 *log10(null_pV)
qqplot_hex(nl_oP,nl_unif,title="Vienna 2010 males trident against uniform")
dev.copy2pdf(file="qq_unif_vie2010_males_trident.pdf")
quartz()
null_pV=scan("FDR/males_vie_pTL1_pTL2_pTD1_pTD2_q20_filt_mc1.5_chroms_only_mct8_mcv25.cmhout_60.0_nullP.out")
nl_null=-1 *log10(null_pV)
qqplot_hex(nl_oP,nl_null,title="Vienna 2010 males trident",alpha=60,beta="=60")
dev.copy2pdf(file="qq_a_60_b_60_vie2010_males_trident.pdf")
quartz(height=6,width=12)
manhattan(data,colors=c("black","slategrey"),limitchromosomes=c("3R:17050000-17071000"),main="Vienna males trident (2010)", suggestiveline=c(-log10(3.126176223804206e-08)))
dev.print(png,width=800,file="vie2010_males_trident_manhattan_3R_ebony.png")
data_nlp_g_1=subset(na.omit(data), (P>0 & P<=0.1))
manhattan(data_nlp_g_1,colors=c("black","slategrey"),limitchromosomes=chroms,main="Vienna males trident (2010)", suggestiveline=c(-log10(3.126176223804206e-08)))
rng<-par("usr")
#[1]  -4757442.56 124991315.56        -0.56        14.56
par(xpd=NA)
text(mean(inv3r_p_coord),rng[4],pos=3,labels=c("3R(p)"), col="red")
rect(inv3r_p_coord[1], rng[3] , inv3r_p_coord[2], rng[4], border=c("red"), lty = 3,lwd=2)
par(xpd=FALSE)
dev.print(png,width=800,file="vie2010_males_trident_manhattan_3Rp_lab.png")
# trident do inversion plotting
inv_data <- read.table("males_vie_pTL1_pTL2_pTD1_pTD2_q20_filt_mc1.5_chroms_only_mct8_mcv25.inversions",header=TRUE)
library(doBy)
summaryBy(pop1 + pop2 + pop3 + pop4 ~ CHR + INV, data = inv_data, FUN =c(median), keep.names = T)
aggregate(cbind(pop1,pop2,pop3,pop4)~CHR+INV,data=inv_data,FUN=median)
a <- splitBy(formula=~ CHR + INV, data = inv_data )
quartz()
boxplot(a[[6]]$pop1,a[[6]]$pop2,a[[6]]$pop3,a[[6]]$pop4,data=a[[6]],notch=F,col=c("cornsilk","cornsilk","brown4","brown4"),ylab="frequency of allele fixed in inversion",names=c("light I", "light II","dark I", "dark II"), main="Inversion 3R(p)")
dev.copy2pdf(file="inv_3Rp_vie2010_males_trident.pdf")
boxplot(a[[1]]$pop1,a[[1]]$pop2,a[[1]]$pop3,a[[1]]$pop4,data=a[[1]],notch=F,col=c("cornsilk","cornsilk","brown4","brown4"),ylab="frequency of allele fixed in inversion",names=c("light I", "light II","dark I", "dark II"), main="Inversion 2L(t)")
dev.copy2pdf(file="inv_2Lt_vie2010_males_trident.pdf")
boxplot(a[[2]]$pop1,a[[2]]$pop2,a[[2]]$pop3,a[[2]]$pop4,data=a[[2]],notch=F,col=c("cornsilk","cornsilk","brown4","brown4"),ylab="frequency of allele fixed in inversion",names=c("light I", "light II","dark I", "dark II"), main="Inversion 2R(ns)")
dev.copy2pdf(file="inv_2Rns_vie2010_males_trident.pdf")
boxplot(a[[3]]$pop1,a[[3]]$pop2,a[[3]]$pop3,a[[3]]$pop4,data=a[[3]],notch=F,col=c("cornsilk","cornsilk","brown4","brown4"),ylab="frequency of allele fixed in inversion",names=c("light I", "light II","dark I", "dark II"), main="Inversion 3L(p)")
dev.copy2pdf(file="inv_3Lp_vie2010_males_trident.pdf")
boxplot(a[[4]]$pop1,a[[4]]$pop2,a[[4]]$pop3,a[[4]]$pop4,data=a[[4]],notch=F,col=c("cornsilk","cornsilk","brown4","brown4"),ylab="frequency of allele fixed in inversion",names=c("light I", "light II","dark I", "dark II"), main="Inversion 3R(c)")
dev.copy2pdf(file="inv_3Rc_vie2010_males_trident.pdf")
boxplot(a[[5]]$pop1,a[[5]]$pop2,a[[5]]$pop3,a[[5]]$pop4,data=a[[5]],notch=F,col=c("cornsilk","cornsilk","brown4","brown4"),ylab="frequency of allele fixed in inversion",names=c("light I", "light II","dark I", "dark II"), main="Inversion 3R(mo)")
dev.copy2pdf(file="inv_3Rmo_vie2010_males_trident.pdf")

# qqplot fro masked 3R(p)
data_wo3rp =scan("males_vie_pTL1_pTL2_pTD1_pTD2_q20_filt_mc1.5_chroms_only_mct8_mcv25_wo_3R_p.pV")
null_wo3rp=scan("FDR/males_vie_pTL1_pTL2_pTD1_pTD2_q20_filt_mc1.5_chroms_only_mct8_mcv25_wo_3R_p.cmhout_78.0_nullP.out")
nl_oP_3rp = -1 * log10(data_wo3rp)
nl_unif_3rp = -1 *log10(runif(length(nl_oP_3rp)))
nl_null_3rp=-1 *log10(null_wo3rp)
quartz()
qqplot_hex(nl_oP_3rp,nl_unif_3rp,title="Vienna 2010 males trident against uniform, wo 3R(p)")
dev.copy2pdf(file="qq_uniform_wo_3Rp_vie2010_males_trident.pdf")
quartz()
qqplot_hex(nl_oP_3rp,nl_null_3rp,title="Vienna 2010 males trident, wo 3R(p)",alpha=78,beta="=78")
dev.copy2pdf(file="qq_a78_wo_3Rp_vie2010_males_trident.pdf")


# obtain FDR estimate
null_p=scan("FDR/males_vie_pTL1_pTL2_pTD1_pTD2_q20_filt_mc1.5_chroms_only_mct8_mcv25.cmhout_44.0_nullP.out10x")
null_p=scan("FDR/males_vie_pTL1_pTL2_pTD1_pTD2_q20_filt_mc1.5_chroms_only_mct8_mcv25.cmhout_60.0_nullP.out")
#null_p=null_pV
null_p=sort(null_p)
obs_pval=sort(data$P)
n_ratio=length(obs_pval)/length(null_p)

ranks = 1:50
for (r in ranks)
	{
	null_smaller =  n_ratio*sum(null_p < obs_pval[r])
	cat ("r=", r, obs_pval[r], null_smaller/r, "\n")
	}
#inversion
library(doBy)

boxplot_inv <- function(a,n,tit="",col1="cornsilk",col2="brown4",name1="L",name2="D")
{
	pops=(length(a[[n]])-5)/2
	names=1:(length(a[[n]])-5)
	for(i in 1:pops){names[i]=paste(name1,i,sep=""); names[i+pops]=paste(name2,i,sep="");}
	boxplot(a[[n]][6:length(a[[n]])],data=a[[n]],notch=F,col=c(rep(col1,pops),rep(col2,pops)),ylab="frequency of allele fixed in inversion",names=names, main=paste("Inversion",a[[n]]$CHR[1],a[[n]]$INV[1],tit))
	}


# trident 2011 do inversion plotting
setwd("/Volumes/Temp/Lukas/Data/Trident/Italy/Polymorphisms/CMH")
inv_data <- read.table("males_ita11_pTL1_pTL2_pTL3_pTD1_pTD2_pTD3_q20_filt_mc2_chroms_only_mct6_mcv15.cmhout.inv",header=TRUE)
summaryBy(pop1 + pop2 + pop3 + pop4 +pop5 + pop6 ~ CHR + INV, data = inv_data, FUN =c(median), keep.names = T)
#aggregate(cbind(pop1,pop2,pop3,pop4)~CHR+INV,data=inv_data,FUN=median)
a_11 <- splitBy(formula=~ CHR + INV, data = inv_data )
quartz()
boxplot_inv(a_11,1,tit=", Trident Italy 2011")
dev.copy2pdf(file="inv2L_tri_ita11.pdf")
boxplot_inv(a_11,2,tit=", Trident Italy 2011")
dev.copy2pdf(file="inv2Rns_tri_ita11.pdf")
boxplot_inv(a_11,3,tit=", Trident Italy 2011")
dev.copy2pdf(file="inv3Lp_tri_ita11.pdf")
boxplot_inv(a_11,4,tit=", Trident Italy 2011")
dev.copy2pdf(file="inv3Rc_tri_ita11.pdf")
boxplot_inv(a_11,5,tit=", Trident Italy 2011")
dev.copy2pdf(file="inv3Rmo_tri_ita11.pdf")
boxplot_inv(a_11,6,tit=", Trident Italy 2011")
dev.copy2pdf(file="inv3Rmo_tri_ita11.pdf")
boxplot_inv(a_11,5,tit=", Trident Italy 2011")
dev.copy2pdf(file="inv3Rmo_tri_ita11.pdf")
0:23:0:47:0:0   0:38:0:46:0:0   0:8:0:22:0:0    0:8:0:90:0:0    0:3:0:120:0:0   0:5:0:81:0:0    1.468466e-18

boxplot_inv_combine <- function(a,b,tit="",col1="lightblue",col2="orange",name1="l",name3="b",name2="d",pops1=1:2,pops2=3:4,pops3=1:3)
{	
	b1=list()
	b2=list()
	b3=list()
	group_names = c()
	pops1=pops1+5
	pops2=pops2+5
	for(i in names(a)){
		b1[[i]]=unlist(a[[i]][,pops1],use.names=F)
		b2[[i]]=unlist(a[[i]][,pops2],use.names=F)
		group_names = append(group_names,paste("IN","(",a[[i]]$CHR[1],")",toupper
(a[[i]]$INV[1]),sep ="")) 
		}
	group_names = c( "IN(2L)t"  ,"IN(2R)NS" ,"IN(3L)P" , "IN(3R)C" , "IN(3R)Mo" ,"IN(3R)P")
	num_all = length(names(a))
	boxplot(b1,at = 0:(num_all-1)*3 + 1, xlim = c(0,num_all*3), ylim = range(b1, b2), xaxt = "n", ylab="frequency of allele fixed in inversion",col=col1, main=tit,notch=T)
    boxplot(b2, at = 0:(num_all-1)*3 + 2, xaxt = "n",col=col2, add = TRUE,notch=T)
   legend("topright",c(name1,name2),fill=c(col1,col2))
  axis(1, at = 0:(num_all-1)*3 + 1.5, labels = group_names, tick = TRUE)
	}

setwd("/Volumes/Temp/Lukas/Data/Trident/BGI_111_112/Joined/CMH")
trid_inv_data =  read.table("males_join_pTLV1_pTLV2_pTDV1_pTDV2_pTLI1_pTLI2_pTLI3_pTDI1_pTDI2_pTDI3__q20_filt_mc2_chroms_only_mct8_mcv15.inversions",header=TRUE)
base_inv= read.table("base_vie_base_ita.inversions",header=TRUE)
library(doBy)
summaryBy(pop1 + pop2 + pop3 + pop4 +pop5 + pop6 +pop7 + pop8 + pop9 + pop10 ~ CHR + INV, data = trid_inv_data, FUN =c(median), keep.names = T)
head(trid_inv_data)
group_names = c( "IN(2L)t"  ,"IN(2R)NS" ,"IN(3L)P" , "IN(3R)C" , "IN(3R)Mo" ,"IN(3R)P")
colnames(trid_inv_data)=c("CHR","INV","BPS","REF","ALT","v1_L","v1_Lc","v2_L","v2_Lc","v1_D","v1_Dc","v2_D","v2_Dc","i1_L","i1_Lc","i2_L","i2_Lc","i3_L","i3_Lc","i1_D","i1_Dc","i2_D","i2_Dc","i3_D","i3_Dc")
colnames(base_inv)=c("CHR","INV","BPS","REF","ALT","v1_B","v1_Bc","v2_B","v2_Bc","v3_B","v3_Bc","i1_B","i1_Bc","i2_B","i2_Bc","i3_B","i3_Bc")

all_inv=merge(trid_inv_data,base_inv,by=c( "CHR","INV","BPS","REF","ALT"),all=T)

base_vie=c("v1_B","v2_B","v3_B")
base_ita=c("i1_B","i2_B","i3_B")
tl_vie=c("v1_L","v2_L")
td_vie=c("v1_D","v2_D")
tl_ita=c("i1_L","i2_L","i3_L")
td_ita=c("i1_D","i2_D","i3_D")

all_inv_list <- splitBy(formula=~ CHR + INV, data = all_inv )
quartz()
boxplot(unlist(all_inv_list[[6]][,tl_vie]),unlist(all_inv_list[[6]][,base_vie]),unlist(all_inv_list[[6]][td_vie]),unlist(all_inv_list[[6]][,tl_ita]),unlist(all_inv_list[[6]][,base_ita]),unlist(all_inv_list[[6]][td_ita]),data=all_inv_list[[6]],notch=T,col=c("cornsilk","yellow","brown4"),at=c(0.5,1,1.5,2.5,3,3.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark"), main="Inversion 3R(p)")
#text(c(1,3),-0.05,labels=c("Vienna 2010", "Italy 2011"))
mtext("Vienna 2010", at=1,line=2,side=1)
mtext("Italy 2011", at=3,line=2,side=1)
dev.copy2pdf(file="inv_3Rp_vie_ita_males_trident.pdf")
quartz()
boxplot(unlist(all_inv_list[[1]][,tl_vie]),unlist(all_inv_list[[1]][,base_vie]),unlist(all_inv_list[[1]][td_vie]),unlist(all_inv_list[[1]][,tl_ita]),unlist(all_inv_list[[1]][,base_ita]),unlist(all_inv_list[[1]][td_ita]),data=all_inv_list[[1]],notch=T,col=c("cornsilk","yellow","brown4"),at=c(0.5,1,1.5,2.5,3,3.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark"), main="Inversion 2L(t)")
#text(c(1,3),-0.05,labels=c("Vienna 2010", "Italy 2011"))
mtext("Vienna 2010", at=1,line=2,side=1)
mtext("Italy 2011", at=3,line=2,side=1)
dev.copy2pdf(file="inv_2Lt_vie_ita_males_trident.pdf")
quartz()
boxplot(unlist(all_inv_list[[2]][,tl_vie]),unlist(all_inv_list[[2]][,base_vie]),unlist(all_inv_list[[2]][td_vie]),unlist(all_inv_list[[2]][,tl_ita]),unlist(all_inv_list[[2]][,base_ita]),unlist(all_inv_list[[2]][td_ita]),data=all_inv_list[[2]],notch=T,col=c("cornsilk","yellow","brown4"),at=c(0.5,1,1.5,2.5,3,3.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark"), main="Inversion 2R(ns)")
#text(c(1,3),-0.05,labels=c("Vienna 2010", "Italy 2011"))
mtext("Vienna 2010", at=1,line=2,side=1)
mtext("Italy 2011", at=3,line=2,side=1)
dev.copy2pdf(file="inv_2Rns_vie_ita_males_trident.pdf")
quartz()
boxplot(unlist(all_inv_list[[3]][,tl_vie]),unlist(all_inv_list[[3]][,base_vie]),unlist(all_inv_list[[3]][td_vie]),unlist(all_inv_list[[3]][,tl_ita]),unlist(all_inv_list[[3]][,base_ita]),unlist(all_inv_list[[3]][td_ita]),data=all_inv_list[[3]],notch=T,col=c("cornsilk","yellow","brown4"),at=c(0.5,1,1.5,2.5,3,3.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark"), main="Inversion 3L(p)")
#text(c(1,3),-0.05,labels=c("Vienna 2010", "Italy 2011"))
mtext("Vienna 2010", at=1,line=2,side=1)
mtext("Italy 2011", at=3,line=2,side=1)
dev.copy2pdf(file="inv_3Lp_vie_ita_males_trident.pdf")
quartz()
boxplot(unlist(all_inv_list[[4]][,tl_vie]),unlist(all_inv_list[[4]][,base_vie]),unlist(all_inv_list[[4]][td_vie]),unlist(all_inv_list[[4]][,tl_ita]),unlist(all_inv_list[[4]][,base_ita]),unlist(all_inv_list[[4]][td_ita]),data=all_inv_list[[4]],notch=T,col=c("cornsilk","yellow","brown4"),at=c(0.5,1,1.5,2.5,3,3.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark"), main="Inversion 3R(c)")
#text(c(1,3),-0.05,labels=c("Vienna 2010", "Italy 2011"))
mtext("Vienna 2010", at=1,line=2,side=1)
mtext("Italy 2011", at=3,line=2,side=1)
dev.copy2pdf(file="inv_3Rc_vie_ita_males_trident.pdf")
quartz()
boxplot(unlist(all_inv_list[[5]][,tl_vie]),unlist(all_inv_list[[5]][,base_vie]),unlist(all_inv_list[[5]][td_vie]),unlist(all_inv_list[[5]][,tl_ita]),unlist(all_inv_list[[5]][,base_ita]),unlist(all_inv_list[[5]][td_ita]),data=all_inv_list[[5]],notch=T,col=c("cornsilk","yellow","brown4"),at=c(0.5,1,1.5,2.5,3,3.5),boxwex=0.4,ylab="frequency of allele fixed in inversion",names=c("light", "base","dark", "light", "base","dark"), main="Inversion 3R(mo)")
#text(c(1,3),-0.05,labels=c("Vienna 2010", "Italy 2011"))
mtext("Vienna 2010", at=1,line=2,side=1)
mtext("Italy 2011", at=3,line=2,side=1)
dev.copy2pdf(file="inv_3Rmo_vie_ita_males_trident.pdf")
setwd("/Volumes/Temp/Lukas/Data/Trident/BGI_111_112/Joined/CMH")
trid_pV=read.table("males_join_V_LD_I_LD_VI_LD_filt_mc2_chroms_only_mct8_mcv15.cmh_pV_OR.ple0.01.gz",header=F,na.strings = c("NA","nan"))
colnames(trid_pV)=c("CHR","BPS","Allele","P_vi","O_vi","Oa_vi","Ob_vi","P_ita","O_ita","Oa_ita","Ob_ita","P_all","O_all","Oa_all","Ob_all")
trid_pV$lO_vi=log(trid_pV$O_vi)
trid_pV$lO_ita=log(trid_pV$O_ita)
trid_pV$lO_all=log(trid_pV$O_all)
trid_pV$lP_vi=-1*log10(trid_pV$P_vi)
trid_pV$lP_ita=-1*log10(trid_pV$P_ita)
trid_pV$lP_all=-1*log10(trid_pV$P_all)

has_Pv = ! (is.na(trid_pV$P_ita) | is.na(trid_pV$P_vi))
summary(trid_pV$lO_vi[is.finite(trid_pV$lO_vi)]) # Min -4.1630 Max 3.6730 
summary(trid_pV$lO_ita[is.finite(trid_pV$lO_ita)]) # Min -4.2270 Max  4.6490
summary(trid_pV$lO_all[is.finite(trid_pV$lO_all)]) # Min -4.8030 Max  4.6970

trid_pV$lO_vi_ni=trid_pV$lO_vi
trid_pV$lO_ita_ni=trid_pV$lO_ita

trid_pV$lO_vi_ni[trid_pV$lO_vi_ni == Inf] = 5
trid_pV$lO_vi_ni[trid_pV$lO_vi_ni == -Inf] = -5
trid_pV$lO_ita_ni[trid_pV$lO_ita_ni == Inf] = 5
trid_pV$lO_ita_ni[trid_pV$lO_ita_ni == -Inf] = -5


plot_OR_all_pV <- function(region,coords,lpVals,ylab_pv=expression(-log[10](italic(P))) ){
    region_idx = which(trid_pV$CHR == coords[1] & (coords[2] <= trid_pV$BPS & trid_pV$BPS <= coords[3] ) & (trid_pV$lP_ita >= lpVals[1] | trid_pV$lP_vi >= lpVals[2]))
    length(region_idx)
    layout(matrix(1:2, ncol = 1), widths = 1, heights = c(0.75,1.25), respect = FALSE)
    par(mar = c(0, 4.1, 4.1, 2.1))
    main_tit=bquote("SNPs around" ~.(region) ~"with" ~ P_ita <= 10^-.(lpVals[1]) ~ " or " ~  P_vi <= 10 ^-.(lpVals[2]))
    xlimit=c(1-0.25,length(region_idx)+0.15)  
    plot((1:length(region_idx))-0.15,trid_pV$lP_vi[region_idx],col="blue",cex=1,pch=16,ylim=c(0,max(trid_pV[region_idx,c("lP_vi","lP_ita","lP_all")],na.rm=T)),ylab=ylab_pv, xaxt="n",main=main_tit,xlim=xlimit )
    abline(v=(1:length(region_idx))-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+0.25,col="grey",lty="dotted",lw=2)
    points(1:length(region_idx),trid_pV$lP_ita[region_idx],col="red",cex=1,pch=16)
    points((1:length(region_idx))+0.15,trid_pV$lP_all[region_idx],col="green",cex=1,pch=16)
    legend("topright",legend=c("Vie 2010","Ita 2011","both"),pch=c(16,16,16),col=c("blue","red","green"))
    par(mar = c(4.1, 4.1, 0, 2.1))
    ors=c("lO_vi_ni","lO_ita_ni","lO_all")
    vals=unlist(trid_pV[region_idx,ors])
    ymin=min(vals[is.finite(vals)],na.rm=T)
    ymax=max(vals[is.finite(vals)],na.rm=T)
    plot((1:length(region_idx))-0.15,trid_pV[region_idx,"lO_vi_ni"],col="blue",cex=1,pch=16,ylim=c(ymin,ymax),ylab="log(OR)",xlab="SNPs",xlim=xlimit)
    abline(v=(1:length(region_idx))-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(1:length(region_idx))+0.25,col="grey",lty="dotted",lw=2)
    abline(h=seq(-4,4,by=1),col="grey",lty="dotted",lw=2)

    points((1:length(region_idx)),trid_pV[region_idx,"lO_ita_ni"],col="red",cex=1,pch=16)
    points((1:length(region_idx))+0.15,trid_pV[region_idx,"lO_all"],col="green",cex=1,pch=16)
    #legend("bottomright",legend=c("very light","light","dark","very dark"),title="log(OR) to base",pch=c(6,23,23,17),col=c("black"),pt.bg=c("white","white","grey16","black"),pt.cex =c(1,0.75,0.75,1))
}

#tan=c("X",8867539,9170829)
tan=c("X",9111115,9124880)
tan_p=c(8,8)
quartz(width=18,height=8)
plot_OR_all_pV("tan",tan,tan_p)
dev.copy2pdf(file="trid_vie_ita_all_OR_comp_tan_8_8.pdf")

bab1=c("3L",1071314,1120551 )
bab1_p=c(6,6)
quartz(width=18,height=8)
plot_OR_all_pV("bab1",bab1,bab1_p)

ebony=c("3R",17050240,17068514)
ebony_p=c(8,8)
quartz(width=18,height=8)
plot_OR_all_pV("ebony",ebony,ebony_p)
dev.copy2pdf(file="tri_vie_ita_all_OR_comp_ebony_lpeu8_lpsa8.pdf")

vi_plot=trid.all[trid.all$P_tv <= 0.01 & ! is.na(trid.all$P_tv),c("CHR","BPS","P_tv")]
ita_plot=trid.all[trid.all$P_ti <= 0.01 & ! is.na(trid.all$P_ti),c("CHR","BPS","P_ti")]
all_plot=trid.all[trid.all$P_te <= 0.01 & ! is.na(trid.all$P_te),c("CHR","BPS","P_te")]
#sa_plot=eu_sa_freqs[eu_sa_freqs$Psa <= 0.01 & ! is.na(eu_sa_freqs$Psa),c("CHR","BPS","Psa")]
colnames(vi_plot)=c("CHR","BP","P")
colnames(ita_plot)=c("CHR","BP","P")
colnames(all_plot)=c("CHR","BP","P")
chroms=c("X","2L","2R","3L","3R","4")
off_vi=get_offset_manh(vi_plot,limitchromosomes=chroms)
off_ita=get_offset_manh(ita_plot,limitchromosomes=chroms)
off_all=get_offset_manh(all_plot,limitchromosomes=chroms)
in2lt=c(2204115,13212655)
in2rns=c(11319930,16163246)
in3lp=c(3160436,16362055)
in3rc=c(15922589,28234649)
in3rp=c(12401910,20580518)
in3rmo=c(17058538,24780914)

x11(height=6,width=12)
par(bg="white")
manhattan(vi_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="Vienna males trident (2010)", suggestiveline=c(-log10(3.126176223804206e-08)))
rect(8867539+off_vi["X"],0,9170829+off_vi["X"],60,lty=3, border="red")
rect(1071314+off_vi["3L"],0,1120551+off_vi["3L"],60,lty=3, border="green")
rect(17055975+off_vi["3R"],0,17069171+off_vi["3R"],60,lty=3, border="blue")
rect(17055975+off_vi["3R"],0,17069171+off_vi["3R"],60,lty=3, border="blue")
inv3r_p_coord=in3rp+off_vi['3R']
inv3rc_coord=in3rc+off_vi['3R']
inv2lt_coord=in2lt+off_vi['2L']
rng<-par("usr")
par(xpd=NA)
text(mean(inv3r_p_coord)-1500000,rng[4],pos=3,labels=c("In(3R)p"), col="black")
rect(inv3r_p_coord[1], rng[3] , inv3r_p_coord[2], rng[4], border=c("grey"), lty = 3,lwd=2)
text(mean(inv3rc_coord)+1750000,rng[4],pos=3,labels=c("In(3R)c"), col="black")
rect(inv3rc_coord[1], rng[3] , inv3rc_coord[2], rng[4], border=c("grey"), lty = 2,lwd=2)
text(mean(inv2lt_coord),rng[4],pos=3,labels=c("In(2L)t"), col="black")
rect(inv2lt_coord[1], rng[3] , inv2lt_coord[2], rng[4], border=c("grey"), lty = 3,lwd=2)
par(xpd=FALSE)
#rect(4179149+off_vi["2R"],0,4319221+off_vi["2R"],60,lty=3, border="yellow")
dev.print(png,width=800,file="vie_trid_manhattan.png")

x11(height=6,width=12)
par(bg="white")
manhattan(ita_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="Bolzano males trident (2010)", suggestiveline=c(-log10(3.126176223804206e-08)))
rect(8867539+off_ita["X"],0,9170829+off_ita["X"],60,lty=3, border="red")
rect(1071314+off_ita["3L"],0,1120551+off_ita["3L"],60,lty=3, border="green")
rect(17055975+off_ita["3R"],0,17069171+off_ita["3R"],60,lty=3, border="blue")
rect(17055975+off_ita["3R"],0,17069171+off_ita["3R"],60,lty=3, border="blue")
inv3r_p_coord=in3rp+off_ita['3R']
inv3rc_coord=in3rc+off_ita['3R']
inv2lt_coord=in2lt+off_ita['2L']
rng<-par("usr")
par(xpd=NA)
text(mean(inv3r_p_coord)-1500000,rng[4],pos=3,labels=c("In(3R)p"), col="black")
rect(inv3r_p_coord[1], rng[3] , inv3r_p_coord[2], rng[4], border=c("grey"), lty = 3,lwd=2)
text(mean(inv3rc_coord)+1750000,rng[4],pos=3,labels=c("In(3R)c"), col="black")
rect(inv3rc_coord[1], rng[3] , inv3rc_coord[2], rng[4], border=c("grey"), lty = 2,lwd=2)
text(mean(inv2lt_coord),rng[4],pos=3,labels=c("In(2L)t"), col="black")
rect(inv2lt_coord[1], rng[3] , inv2lt_coord[2], rng[4], border=c("grey"), lty = 3,lwd=2)
par(xpd=FALSE)
#rect(4179149+off_ita["2R"],0,4319221+off_ita["2R"],60,lty=3, border="yellow")
dev.print(png,width=800,file="ita_trid_manhattan.png")

x11(height=6,width=12)
par(bg="white")
manhattan(all_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="Vienna and Bolzano males trident (2010 & 2011)", suggestiveline=c(-log10(3.126176223804206e-08)))
rect(8867539+off_all["X"],0,9170829+off_all["X"],60,lty=3, border="red")
rect(1071314+off_all["3L"],0,1120551+off_all["3L"],60,lty=3, border="green")
rect(17055975+off_all["3R"],0,17069171+off_all["3R"],60,lty=3, border="blue")
rect(17055975+off_all["3R"],0,17069171+off_all["3R"],60,lty=3, border="blue")
inv3r_p_coord=in3rp+off_all['3R']
inv3rc_coord=in3rc+off_all['3R']
inv2lt_coord=in2lt+off_all['2L']
rng<-par("usr")
par(xpd=NA)
text(mean(inv3r_p_coord)-1500000,rng[4],pos=3,labels=c("In(3R)p"), col="black")
rect(inv3r_p_coord[1], rng[3] , inv3r_p_coord[2], rng[4], border=c("grey"), lty = 3,lwd=2)
text(mean(inv3rc_coord)+1750000,rng[4],pos=3,labels=c("In(3R)c"), col="black")
rect(inv3rc_coord[1], rng[3] , inv3rc_coord[2], rng[4], border=c("grey"), lty = 2,lwd=2)
text(mean(inv2lt_coord),rng[4],pos=3,labels=c("In(2L)t"), col="black")
rect(inv2lt_coord[1], rng[3] , inv2lt_coord[2], rng[4], border=c("grey"), lty = 3,lwd=2)
par(xpd=FALSE)
#rect(4179149+off_all["2R"],0,4319221+off_all["2R"],60,lty=3, border="yellow")
dev.print(png,width=800,file="all_trid_manhattan.png")

quartz(height=6,width=12)
par(bg="white")
manhattan(ita_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="Italy (2011) trident", suggestiveline=c(8.5))
rect(8867539+off_ita["X"],0,9170829+off_ita["X"],60,lty=3, border="red")
rect(1071314+off_ita["3L"],0,1120551+off_ita["3L"],60,lty=3, border="green")
rect(17055975+off_ita["3R"],0,17069171+off_ita["3R"],60,lty=3, border="blue")
#rect(4179149+off_ita["2R"],0,4319221+off_ita["2R"],60,lty=3, border="yellow")
dev.print(png,width=800,file="ita_trid_manhattan.png")


quartz(height=6,width=12)
par(bg="white")
manhattan(all_plot,colors=c("black","slategrey"),limitchromosomes=chroms,main="Vienna and Italy trident", suggestiveline=c(8.5))
rect(8867539+off_all["X"],0,9170829+off_all["X"],60,lty=3, border="red")
rect(1071314+off_all["3L"],0,1120551+off_all["3L"],60,lty=3, border="green")
rect(17055975+off_all["3R"],0,17069171+off_all["3R"],60,lty=3, border="blue")
#rect(4179149+off_all["2R"],0,4319221+off_all["2R"],60,lty=3, border="yellow")
dev.print(png,width=800,file="all_trid_manhattan.png")

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

plot_region = function(chrom,a,b, title,ylab=expression(-log[10](italic(P))),df=trid.all, minlp=2.5){
    #region_vi = which(vi_plot$CHR == chrom & vi_plot$BPS >= a  & vi_plot$BPS <= b )
    #region_ita = which(ita_plot$CHR == chrom & ita_plot$BPS >= a  & ita_plot$BPS <= b )
    region_all = (df$CHR == chrom & df$BPS >= a  & df$BPS <= b )
    ymax = max((df$lP_tv[region_all]),(df$lP_ti[region_all]),(df$lP_te[region_all]),na.rm=T)
    plot(df$BPS[region_all & df$lP_tv >= minlp],(df$lP_tv[region_all & df$lP_tv >= minlp]), col="blue", pch=20, main=title, ylim=c(0,ymax),xlab=chrom,ylab=ylab,cex=0.75 )
    points(df$BPS[region_all & df$lP_ti >= minlp],(df$lP_ti[region_all & df$lP_ti >= minlp]), col="red", pch=3,cex=0.95)
    points(df$BPS[region_all & df$lP_te >= minlp],(df$lP_te[region_all & df$lP_te >= minlp]), col="green", pch=5,cex=0.75)
}

require(grDevices)
quartz()
x11()
plot_region("X",9111115,9124880,"tan")
rect(9111688,0,9117290,-0.5,col="grey")
text(9117100,0.35,labels="tan")
text(9116500,-0.25,labels="< < <")
rect(9121170,0,9122860,-0.5,col="grey")
text(9122000,0.35,labels="Gr8a")
text(9121350,-0.25,labels=">")
rect(9119470,0,9120721,-0.5,col="grey")
text((9119470+9120721)/2-800,0.35,labels="CG15370")
text((9119470+9120721)/2,-0.25,labels=">")
rect(9120721,0,9121170,-0.5,col="white")
text((9120721+9121170)/2,0.35,labels="MSE")
dev.copy2pdf(file="trid_vie_ita_all_tan_region.pdf")

quartz()
plot_region("3R",17055354,17068822,"ebony")
segments(17055561,-0.25,17062899,-0.25)
rect(c(17055561,17056145,17056381,17056539,17057398,17057631,17062616),0,c(17056088,17056321,17056472,17057310,17057566,17058756,17062899),-0.5,col="grey")
text(17062700,0.35,labels="e")
text(17062650,-0.25,labels="<")
rect(17065936,0,17068318,-0.5,col="white")
text((17065936+17068318)/2,0.35,labels="aCRE")
dev.copy2pdf(file="trid_vie_ita_all_ebony_region.pdf")

setwd("/Volumes/Temp/Lukas/Data/Trident/BGI_111_112/Joined/CMH")
trid.all=read.table("males_trident_V_LD_I_LD_VI_LD.cmh_pV_OR.lightalleles.a7pig_females_vie_LD_ita_LD_eu_LD.pV_OR.lightalleles.afs.males_trid_pTLV1_pTLV2_pTDV1_pTDV2_pTLI1_pTLI2_pTLI3_pTDI1_pTDI2_pTDI3_fem_a7_pIel_pIIel_pIIIel_pIhl_pIIhl_pIIIhl.pV_lt_1e_2.af.gz",na.strings = c("NA","nan"),header=F)
colnames(trid.all)=c("CHR","BPS","Allele","P_tv","O_tv","Oa_tv","Ob_tv","P_ti","O_ti","Oa_ti","Ob_ti","P_te","O_te","Oa_te","Ob_te","Ltv","Lti","Lte","P_av","O_av","Oa_av","Ob_av","P_ai","O_ai","Oa_ai","Ob_ai","P_ae","O_ae","Oa_ae","Ob_ae","Lav","Lai","Lae","VL1","VL2","VD1","VD2","IL1","IL2","IL3","ID1","ID2","ID3","Vb1","Vb2","Vb3","Ib1","Ib2","Ib3")
pop.means=list(c("VL1","VL2"),c("VD1","VD2"),c("IL1","IL2","IL3"),c("ID1","ID2","ID3"),c("Vb1","Vb2","Vb3"),c("Ib1","Ib2","Ib3"),c("VL1","VL2","IL1","IL2","IL3"),c("VD1","VD2","ID1","ID2","ID3"),c("Vb1","Vb2","Vb3","Ib1","Ib2","Ib3"))
pop.means.names=c("VL","VD","IL","ID","Vb","Ib","EL","ED","Eb")
for(i in 1:length(pop.means.names)){
    trid.all[,pop.means.names[i]]=rowMeans(trid.all[,unlist(pop.means[i])])
}

P=c("P_tv","P_ti","P_te","P_av","P_ai","P_ae")
lP=c("lP_tv","lP_ti","lP_te","lP_av","lP_ai","lP_ae")
O=c("O_tv","Oa_tv","Ob_tv","O_ti","Oa_ti","Ob_ti","O_te","Oa_te","Ob_te","O_av","Oa_av","Ob_av","O_ai","Oa_ai","Ob_ai","O_ae","Oa_ae","Ob_ae")
lO=c("lO_tv","lOa_tv","lOb_tv","lO_ti","lOa_ti","lOb_ti","lO_te","lOa_te","lOb_te","lO_av","lOa_av","lOb_av","lO_ai","lOa_ai","lOb_ai","lO_ae","lOa_ae","lOb_ae")
trid.all[,O]=log(trid.all[,O])
trid.all[,lP]=-1*log10(trid.all[,P])
trid.all$LF=0
maxPs=apply(trid.all[,c("lP_tv","lP_ti","lP_te")],1,function(x) which.max(x))
a=c("O_tv","O_ti","O_te")
trid.all$LF=sapply(1:length(maxPs),function(x) sign(trid.all[x,a[maxPs[x]]]))
trid.all[,lO] = trid.all[,O]*trid.all$LF
for(i in  lO){
    #trid.all[is.na(trid.all[,i]),i] = 0
    trid.all[(trid.all[,i] == -Inf) & !(is.na(trid.all[,i])),i] = -5
    trid.all[(trid.all[,i] == Inf) & !(is.na(trid.all[,i])),i] = 5
}


tan=c("X",9114500,9122100)
tan_p=c(4,5,6)
x11(width=18,height=8)
a=plot_OR_pV("tan",tan,tan_p,df=trid.all,pVs=c("lP_tv","lP_ti","lP_te"),ors=c("lO_tv","lO_ti","lO_te"),pVname=c("vie","ita","eu"),ors.ci=list(c("lOa_tv","lOb_tv"),c("lOa_ti","lOb_ti"),c("lOa_te","lOb_te")))
ymax=4.0
tan.x=get.xpos.from.coords(a,c(9111688,9117290))
Gr8a.x=get.xpos.from.coords(a,c(9121170,9122860))
cg15370.x=get.xpos.from.coords(a,c(9119470,9120721))
CG12121.x=get.xpos.from.coords(a,c(9122820,9126406))
Ir8a.x=get.xpos.from.coords(a,c(9126759,9130618)) 
yh=0.225
rect(tan.x[1],ymax,tan.x[2],ymax+yh,col="grey")
rect(Gr8a.x[1],ymax,Gr8a.x[2],ymax+yh,col="grey")
rect(cg15370.x[1],ymax,cg15370.x[2],ymax+yh,col="grey")
rect(CG12121.x[1],ymax,CG12121.x[2],ymax+yh,col="grey")
rect(Ir8a.x[1],ymax,Ir8a.x[2],ymax+yh,col="grey")
text(mean(tan.x),ymax+0.1,labels="tan")
text(mean(Gr8a.x),ymax+0.1,labels="Gr8a")
text(mean(cg15370.x),ymax+0.1,labels="CG15370")
text(mean(Ir8a.x),ymax+0.1,labels="Ir8a")
text(mean(CG12121.x),ymax+0.1,labels="CG12121")
rect(Gr8a.x[1],ymax,cg15370.x[2],ymax+yh,col="white")
text(mean(c(Gr8a.x[1],cg15370.x[2])),ymax+0.1,labels="MSE")
#text(mean(c(Gr8a.x[1],cg15370.x[2])),1.025,labels="MSE")
dev.copy2pdf(file="males_trid_tan_pv_1e_4_pi_1e_5.pdf")

tan=c("X",9114500,9122100)
tan_p=c(4,5,6)
x11(width=18,height=8)
a=plot_AF_pV("tan",tan,tan_p,df=trid.all,pVs=c("lP_tv","lP_ti","lP_te"),afs=list(c("VL","Vb","VD"),c("IL","Ib","ID"),c("EL","Eb","ED")),pVname=c("vie","ita","eu"))
ymax=1.01
tan.x=get.xpos.from.coords(a,c(9111688,9117290))
Gr8a.x=get.xpos.from.coords(a,c(9121170,9122860))
cg15370.x=get.xpos.from.coords(a,c(9119470,9120721))
CG12121.x=get.xpos.from.coords(a,c(9122820,9126406))
Ir8a.x=get.xpos.from.coords(a,c(9126759,9130618)) 
yh=0.04
rect(tan.x[1],ymax,tan.x[2],ymax+yh,col="grey")
text(mean(tan.x),ymax+0.02,labels="tan")
rect(Gr8a.x[1],ymax,Gr8a.x[2],ymax+yh,col="grey")
rect(cg15370.x[1],ymax,cg15370.x[2],ymax+yh,col="grey")
rect(CG12121.x[1],ymax,CG12121.x[2],ymax+yh,col="grey")
rect(Ir8a.x[1],ymax,Ir8a.x[2],ymax+yh,col="grey")
text(mean(Gr8a.x),ymax+0.02,labels="Gr8a")
text(mean(cg15370.x),ymax+0.02,labels="CG15370")
text(mean(Ir8a.x),ymax+0.1,labels="Ir8a")
text(mean(CG12121.x),ymax+0.02,labels="CG12121")
rect(Gr8a.x[1],ymax,cg15370.x[2],ymax+yh,col="white")
text(mean(c(Gr8a.x[1],cg15370.x[2])),ymax+0.02,labels="MSE")
#text(mean(c(Gr8a.x[1],cg15370.x[2])),1.025,labels="MSE")
dev.copy2pdf(file="males_trid_afs_tan_pv_1e_4_pi_1e_5.pdf")

tan=c("X",9114500,9122100)
tan_p=c(8,14)
#x11(width=18,height=8)
a=plot_OR_pV("tan",tan,tan_p,df=trid.all,pVs=c("lP_te","lP_ae"),ors=c("lO_te","lO_ae"),pVname=c("trid","abd"),ors.ci=list(c("lOa_te","lOb_te"),c("lOa_ae","lOb_ae")))
ymax=4.0
tan.x=get.xpos.from.coords(a,c(9111688,9117290))
Gr8a.x=get.xpos.from.coords(a,c(9121170,9122860))
cg15370.x=get.xpos.from.coords(a,c(9119470,9120721))
CG12121.x=get.xpos.from.coords(a,c(9122820,9126406))
Ir8a.x=get.xpos.from.coords(a,c(9126759,9130618)) 
yh=0.225
if(max(tan.x) <= 0 ){
    tan.x=c(0,0.725)
}
rect(tan.x[1],ymax,tan.x[2],ymax+yh,col="grey")
rect(Gr8a.x[1],ymax,Gr8a.x[2],ymax+yh,col="grey")
rect(cg15370.x[1],ymax,cg15370.x[2],ymax+yh,col="grey")
rect(CG12121.x[1],ymax,CG12121.x[2],ymax+yh,col="grey")
rect(Ir8a.x[1],ymax,Ir8a.x[2],ymax+yh,col="grey")
text(mean(tan.x),ymax+0.1,labels="tan")
text(mean(Gr8a.x),ymax+0.1,labels="Gr8a")
text(mean(cg15370.x),ymax+0.1,labels="CG15370")
text(mean(Ir8a.x),ymax+0.1,labels="Ir8a")
text(mean(CG12121.x),ymax+0.1,labels="CG12121")
rect(Gr8a.x[1],ymax,cg15370.x[2],ymax+yh,col="white")
text(mean(c(Gr8a.x[1],cg15370.x[2])),ymax+0.1,labels="MSE")
#text(mean(c(Gr8a.x[1],cg15370.x[2])),1.025,labels="MSE")
dev.copy2pdf(file="comp_trid_abd_trid_tan_pt_1e_8_pi_1e_14.pdf")

ebony=c("3R",17055561,17075171)
ebony_p=c(10,12,14)
a=plot_AF_pV("ebony",ebony,ebony_p,df=trid.all,pVs=c("lP_tv","lP_ti","lP_te"),pVname=c("vie","ita","eu"),afs=list(c("VL","Vb","VD"),c("IL","Ib","ID"),c("EL","Eb","ED")))
ymax=1.01
yh=0.04
cre.x=get.xpos.from.coords(a,c(17065936,17068318))
ebony.x=get.xpos.from.coords(a,c(17055561,17062900))
rect(ebony.x[1],ymax,ebony.x[2],ymax +yh,col="grey")
text(mean(ebony.x),ymax+0.02,labels="e")
rect(cre.x[1],ymax,cre.x[2],ymax+yh,col="white")
text(mean(cre.x),ymax+0.02,labels="aCRE")
dev.copy2pdf(file="males_trid_afs_ebony_pv_1e_10_pi_1e_12.pdf")

ebony=c("3R",17055561,17075171)
ebony_p=c(10,12,14)
a=plot_OR_pV("ebony",ebony,ebony_p,df=trid.all,pVs=c("lP_tv","lP_ti","lP_te"),ors=c("lO_tv","lO_ti","lO_te"),pVname=c("vie","ita","eu"),ors.ci=list(c("lOa_tv","lOb_tv"),c("lOa_ti","lOb_ti"),c("lOa_te","lOb_te")))
ymax=4.0
yh=0.225
cre.x=get.xpos.from.coords(a,c(17065936,17068318))
ebony.x=get.xpos.from.coords(a,c(17055561,17062900))
rect(ebony.x[1],ymax,ebony.x[2],ymax +yh,col="grey")
text(mean(ebony.x),ymax+0.1,labels="e")
rect(cre.x[1],ymax,cre.x[2],ymax+yh,col="white")
text(mean(cre.x),ymax+0.1,labels="aCRE")
dev.copy2pdf(file="males_trid_ebony_pv_1e_10_pi_1e_12.pdf")

ebony=c("3R",17055561,17075171)
ebony_p=c(14,6)
a=plot_OR_pV("ebony",ebony,ebony_p,df=trid.all,pVs=c("lP_te","lP_ae"),ors=c("lO_te","lO_ae"),pVname=c("trid","abd"),ors.ci=list(c("lOa_te","lOb_te"),c("lOa_ae","lOb_ae")))
ymax=4.0
yh=0.225
cre.x=get.xpos.from.coords(a,c(17065936,17068318))
ebony.x=get.xpos.from.coords(a,c(17055561,17062900))
rect(ebony.x[1],ymax,ebony.x[2],ymax +yh,col="grey")
text(mean(ebony.x),ymax+0.1,labels="e")
rect(cre.x[1],ymax,cre.x[2],ymax+yh,col="white")
text(mean(cre.x),ymax+0.1,labels="aCRE")
dev.copy2pdf(file="comp_trid_abd_ebony_pt_1e_14_pa_1e_6.pdf")


get.xpos.from.coords <- function(df,coords) {
    stretch <-  df$xpos[ coords[1] <= df$BPS & df$BPS <= coords[2]]
    if (length(stretch) > 0){
        return(c(min(stretch)-0.25,max(stretch)+0.25))
    }
    else {
        return(c(-1000,-1000))
    }
}


plot_OR_pV <- function(region,coords,lpVals,ylab_pv=expression(-log[10](italic(P))),df=trid.all,pVs=c("lP_tv","lP_ti","lP_te"),ors=c("lO_tv","lO_ti","lO_te"),pVname=c("P_vie","P_ita","P_eu"),ors.ci=F,ycorr=0.35){
    extended=F
    if(length(pVs)==3){
        extended=T
        region_idx = which(df$CHR == coords[1] & (coords[2] <= df$BPS & df$BPS <= coords[3] ) & (df[,pVs[1]] >= lpVals[1] | df[,pVs[2]] >= lpVals[2] | ( extended &  df[,pVs[3]] >= lpVals[3]) ) )
    }
    else{
        region_idx = which(df$CHR == coords[1] & (coords[2] <= df$BPS & df$BPS <= coords[3] ) & (df[,pVs[1]] >= lpVals[1] | df[,pVs[2]] >= lpVals[2]  ) & ! ( is.na(df[,pVs[1]]) | is.na(df[,pVs[2]]) ) )
    }
    
    #length(region_idx)

    df_region=df[region_idx,]
    df_region$xpos=seq(1,length(df_region[,1]))
    main_tit=bquote("SNPs around" ~italic(.(region)) ~"with" ~ P[.(pVname[1])] <= 10^-.(lpVals[1]) ~ " or " ~  P[.(pVname[2])] <= 10 ^-.(lpVals[2]))
    xlimit=c(df_region$xpos[1]-0.25,df_region$xpos[length(df_region$xpos)]+0.15)
    layout(matrix(1:2, ncol = 1), widths = 1, heights = c(0.75,1.25), respect = FALSE)
    par(mar = c(0, 4.1, 4.1, 2.1))
    plot(df_region$xpos,df_region[,pVs[2]],col="white",cex=1,pch=NULL,ylim=c(0,max(df_region[,pVs],na.rm=T)),ylab=ylab_pv, xaxt="n",main=main_tit,xlim=xlimit )
    abline(v=(df_region$xpos)-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(df_region$xpos)+0.25,col="grey",lty="dotted",lw=2)
    # plot pV    
    if (extended) {
        points(df_region$xpos-0.15,df_region[,pVs[1]],col="red",cex=1,pch=16)
        points(df_region$xpos,df_region[,pVs[3]],col="green",cex=1,pch=16)
        points(df_region$xpos+0.15,df_region[,pVs[2]],col="blue",cex=1,pch=16)
    }
    else{
        points(df_region$xpos,df_region[,pVs[2]],col="blue",cex=1,pch=16)
        points(df_region$xpos,df_region[,pVs[1]],col="red",cex=1,pch=16)
    }
    if (! extended) {
        legend("topright",legend=c("Trident","A7"),pch=c(16,16),col=c("red","blue"))
    }
    else{
        legend("topright",legend=c("Vienna","Bolzano","Vie & Bolz"),pch=c(16,16),col=c("red","blue","green"))
    }
    par(mar = c(5.1, 4.1, 0, 2.1))
    if (length(ors.ci) > 1){
        vals=c(unlist(df_region[,c(ors,unlist(ors.ci))]))
    }
    else{
        vals=c(unlist(df_region[,ors])) 
    }
    ymin=min(c(vals[is.finite(vals)],-0.25),na.rm=T)
    ymax=max(vals[is.finite(vals)],na.rm=T)    
    plot((df_region$xpos)+0.15,df_region[,ors[2]],col="blue",cex=1,pch=16,ylim=c(ymin,ymax+0.15),ylab="log(OR)",yaxt="n",xaxt="n",xlab="",xlim=xlimit )
    # plot errorbars
    if (length(ors.ci) > 1){
        for (i in 1:length(df_region$xpos)) {
            lines(c(df_region$xpos[i]-0.15,df_region$xpos[i]-0.15),df_region[i,ors.ci[[1]]],lty="dotted",col="red",lw=2)
            lines(c(df_region$xpos[i]+0.15,df_region$xpos[i]+0.15),df_region[i,ors.ci[[2]]],lty="dotted",col="blue",lw=2)
            points(c(df_region$xpos[i]-0.15,df_region$xpos[i]-0.15),df_region[i,ors.ci[[1]]],pch=4,col="red")
            points(c(df_region$xpos[i]+0.15,df_region$xpos[i]+0.15),df_region[i,ors.ci[[2]]],pch=4,col="blue")
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
    text(x = df_region$xpos, par("usr")[3] - ycorr, labels = df_region$BPS, srt = 45, pos = 1, xpd = TRUE,cex=0.75)
    mtext(coords[1],side=1,line=3, xpd = TRUE,font=2)

    points((df_region$xpos)+0.15,df_region[,ors[2]],col="blue",cex=1,pch=16)
    
    points((df_region$xpos)-0.15,df_region[,ors[1]],col="red",cex=1,pch=16)
    if(extended){points((df_region$xpos),df_region[,ors[3]],col="green",cex=1,pch=16)}
    print(ymax)
    return(df_region[,c("BPS","xpos")])

}

plot_AF_pV <- function(region,coords,lpVals,ylab_pv=expression(-log[10](italic(P))),df=trid.all,pVs=c("lP_tv","lP_ti","lP_te"),afs=list(c("VL","Vb","VD"),c("IL","Ib","ID"),c("EL","Eb","ED")),pVname=c("P_vie","P_ita","P_eu"),ycorr=0.065,llpos="bottomright"){
    #length(region_idx)
    extended=F
    if(length(pVs)==3){
        extended=T
        region_idx = which(df$CHR == coords[1] & (coords[2] <= df$BPS & df$BPS <= coords[3] ) & (df[,pVs[1]] >= lpVals[1] | df[,pVs[2]] >= lpVals[2] | ( extended &  df[,pVs[3]] >= lpVals[3]) ) )
    }
    else{
        region_idx = which(df$CHR == coords[1] & (coords[2] <= df$BPS & df$BPS <= coords[3] ) & (df[,pVs[1]] >= lpVals[1] | df[,pVs[2]] >= lpVals[2]  ) & ! ( is.na(df[,pVs[1]]) | is.na(df[,pVs[2]]) ) )
    }
    

    df_region=df[region_idx,]
                                        #a=table(df_region[,c("hap_eu")])
    df_region$xpos=seq(1,length(df_region[,1]))
                                        #df_region$hap1=0
                                        #df_region$hap2=0
    main_tit=bquote("SNPs around"~italic(.(region)) ~"with" ~ P[.(pVname[1])] <= 10^-.(lpVals[1]) ~ " or " ~  P[.(pVname[2])] <= 10 ^-.(lpVals[2]))
    xlimit=c(df_region$xpos[1]-0.25,df_region$xpos[length(df_region$xpos)]+0.15)
#print(xlimit)
    layout(matrix(1:2, ncol = 1), widths = 1, heights = c(0.75,1.25), respect = FALSE)
    par(mar = c(0, 4.1, 4.1, 2.1))
    plot(df_region$xpos,df_region[,pVs[2]],col="white",cex=1,pch=NULL,ylim=c(0,max(df_region[,pVs],na.rm=T)),ylab=ylab_pv, xaxt="n",main=main_tit,xlim=xlimit )
    abline(v=(df_region$xpos)-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(df_region$xpos)+0.25,col="grey",lty="dotted",lw=2)
    # plot pV
    if (extended) {
        points(df_region$xpos-0.15,df_region[,pVs[1]],col="red",cex=1,pch=16)
        points(df_region$xpos,df_region[,pVs[3]],col="green",cex=1,pch=16)
        points(df_region$xpos+0.15,df_region[,pVs[2]],col="blue",cex=1,pch=16)
    }
    else{
        points(df_region$xpos,df_region[,pVs[2]],col="blue",cex=1,pch=16)
        points(df_region$xpos,df_region[,pVs[1]],col="red",cex=1,pch=16)
    }
    if (! extended) {
        legend("topright",legend=c("Trident","A7"),pch=c(16,16),col=c("blue","red"))
    }
    else{
        legend("topright",legend=c("Vienna","Bolzano","Vie & Bolz"),pch=c(16,16),col=c("red","blue","green"))
    }
    par(mar = c(5.1, 4.1, 0, 2.1))
    # plot AFs
    ymin=min(df_region[,unlist(afs)],na.rm=T)
    ymin=max(c(ymin,0))
    plot((df_region$xpos)-0.15,df_region[,afs[[2]][3]],col="white",cex=1,pch=NULL,ylim=c(ymin,1.05),ylab="major AF",yaxt="n",xaxt="n",xlab="",xlim=xlimit )
    # plot grid lines
    abline(v=(df_region$xpos)-0.25,col="grey",lty="dotted",lw=2)
    abline(v=(df_region$xpos)+0.25,col="grey",lty="dotted",lw=2)
    abline(h=seq(-4,3,by=1),col="grey",lty="dotted",lw=2)
    abline(h=0,col="darkgrey",lty="dotted",lw=2)
    # axes
    axis(1, at=df_region$xpos, labels = FALSE)
    axis(2,at=seq(0,1,by=0.2))
    if (extended){
        points((df_region$xpos)+0.15,df_region[,afs[[2]][2]],col="blue",cex=1,pch=16)
        points((df_region$xpos)+0.15,df_region[,afs[[2]][1]],col="blue",cex=1,pch=6)
        points((df_region$xpos)+0.15,df_region[,afs[[2]][3]],col="blue",cex=1,pch=17)
        points((df_region$xpos),df_region[,afs[[3]][2]],col="green",cex=1,pch=16)
        points((df_region$xpos),df_region[,afs[[3]][1]],col="green",cex=1,pch=6)
        points((df_region$xpos),df_region[,afs[[3]][3]],col="green",cex=1,pch=17)
        
        points((df_region$xpos)-0.15,df_region[,afs[[1]][2]],col="red",cex=1,pch=16)
        points((df_region$xpos)-0.15,df_region[,afs[[1]][1]],col="red",cex=1,pch=6)
        points((df_region$xpos)-0.15,df_region [,afs[[1]][3]],col="red",cex=1,pch=17)
        legend(llpos,legend=c("light","base","dark"),title="average AF",pch=c(6,16,17),col=c("black"),pt.bg=c("white","black","black"),pt.cex =c(1,1,1))

    }
    else{
        points((df_region$xpos)-0.15,df_region[,afs[[2]][3]],col="blue",cex=1,pch=16)
        points((df_region$xpos)-0.15,df_region[,afs[[2]][5]],col="blue",cex=1,pch=17)
        points((df_region$xpos)-0.15,df_region[,afs[[2]][2]],col="blue",cex=0.75,pch=23)
        points((df_region$xpos)-0.15,df_region[,afs[[2]][4]],col="blue",cex=0.75,pch=23,bg="cyan")
        points((df_region$xpos)-0.15,df_region[,afs[[2]][1]],col="blue",cex=0.75,pch=6)
        
        points((df_region$xpos)+0.15,df_region[,afs[[1]][2]],col="red",cex=1,pch=16)
        points((df_region$xpos)+0.15,df_region[,afs[[1]][1]],col="red",cex=1,pch=6)
        points((df_region$xpos)+0.15,df_region [,afs[[1]][3]],col="red",cex=1,pch=17)
        legend(llpos,legend=c("light","base","dark"),title="average AF",pch=c(6,16,17),col=c("black"),pt.bg=c("white","white","black"),pt.cex =c(1,1,1))
    }
    print( df_region$BPS)
    print( df_region$xpos)
    text(x = df_region$xpos, par("usr")[3]-ycorr, labels = df_region$BPS, srt = 45, pos = 1, xpd = TRUE,cex=0.85)
    mtext(coords[1],side=1,line=3, xpd = TRUE,font=2)
    return(df_region[,c("BPS","xpos")])

}
