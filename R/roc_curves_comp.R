# this file loads two GWAS format like files (created with cmh2gwas.pl from popoolation2), and gives ROC like curves comparing the significant SNPs found in both files.

# read the pValues
setwd=("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/Comparison_Vie2010_2011")
popA=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.gwas",header=T)
popB=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.gwas",header=T)
# create a unique ID giving for each SNP position
popA$SNP <- paste(popA$CHR,popA$BP,sep=":")
popB$SNP <- paste(popB$CHR,popB$BP,sep=":")
# merging all SNPs
popAB=merge(popA[,c(3,4)],popB[,c(3,4)],by="SNP",all=T)
names(popAB) = c("SNP","PA","PB")
# getting the numbers of all common and exclusive SNPs
popB_only=length(popAB$PA[is.na(popAB$PA)])
popA_only=length(popAB$PB[is.na(popAB$PB)])
popAB_common=length(popAB$PA)-popA_only-popB_only
popA_all=length(popAB$PA)-popB_only
popB_all=length(popAB$PA)-popA_only
# setting the significance cutoff using the FDR 0.05 correction
sigA=4.134235e-07
sigB=2.324404e-09
# setting the significance cutoff using the Bonferroni correction
sigB=0.05/popB_all
sigA=0.05/popA_all
# get logical vector for significant SNPs in A and B and common SNPs
popA_sig = popAB$PA < sigA & ! is.na(popAB$PA)
popB_sig = popAB$PB < sigB & ! is.na(popAB$PB)
common_AB = !(is.na(popAB$PA) | is.na(popAB$PB))
num_A_sig = length(popAB$PA[popA_sig])
num_B_sig = length(popAB$PB[popB_sig])
num_B_inA_sig = length(popAB$PB[popB_sig  & common_AB])
num_A_inB_sig = length(popAB$PA[popA_sig  & common_AB])
num_A_sig_B_sig_common = length(popAB$PA[ popA_sig  &  popB_sig & common_AB])
# Venn diagrams for overlaps
# needs package VennDiagram
# not very nice looking, as the empty fields are also displayed, but gives an overiew
library(VennDiagram)
quartz()
venn.plot <- draw.quad.venn(popA_all,popB_all,num_A_sig,num_B_sig,popAB_common,num_A_sig,num_B_inA_sig,num_A_inB_sig,num_B_sig,num_A_sig_B_sig_common,num_A_inB_sig,num_B_inA_sig,num_A_sig_B_sig_common,num_A_sig_B_sig_common, num_A_sig_B_sig_common, c("pop A", "pop B","sigA","sigB"),fill=c("blue","green","red","yellow"))
# ROC like curves
# with quantiles
# creating the quantile ranges for the comparison
popB_q = quantile(popAB$PB[common_AB],probs=c(seq(0,0.01,0.0001),seq(0.01,1,0.01)),na.rm=T )
popA_q = quantile(popAB$PA[common_AB],probs=c(seq(0,0.01,0.0001),seq(0.01,1,0.01)),na.rm=T )
# popA against popB
pop_PA_all = na.omit(popAB$PA[common_AB])
popA_B_x <- sapply(popA_q,function(x) {length(pop_PA_all[pop_PA_all < x ])})
popA_B_x <- popA_B_x / length(pop_PA_all)
popA_B_y <- sapply(popA_q,function(x) {length(popAB$PB[ popAB$PA < x & popB_sig & common_AB ])})
popA_B_y <- popA_B_y/num_B_inA_sig
quartz()
plot(popA_B_x,popA_B_y,type="l",col="red",lwd=2,main="population A against B",xlab="fraction of common SNPs in population A below pValue threshold", ylab="fract. of common SNPs sig. in B and below threshold in A")
# plot both in one:
quartz()
plot(popA_B_x,popA_B_y,type="l",col="red",lwd=2,xlab="fraction of SNPs in A with lowest P values", ylab="common SNPs found significant in B")
lines(popB_A_x,popB_A_y, col="blue",lwd=2)
lines(popB_A_x,popB_Ar_y,type="l",col="grey",lwd=1)
legend("bottomright",c("A:Vienna 2010, B:Bolzano 2011","A:Bolzano 2011, B:Vienna 2010"),col=c("red","blue"))
dev.copy2pdf(file="sign_overlap_fdr0.05.pdf")
plot(steps,pop_rA_By,type="l",col="red",lwd=2,xlab="SNPs in A with rank < x", ylab="common SNPs found significant in B")
lines(steps,pop_rB_Ay, col="blue",lwd=2)
legend("topleft",c("A:Vienna 2010, B:Bolzano 2011","A:Bolzano 2011, B:Vienna 2010"),col=c("red","blue"))
dev.copy2pdf(file="sign_overlap_fdr05_logscale.pdf")

# popB against popA
pop_PB_all = na.omit(popAB$PB[common_AB])
popB_A_x <- sapply(popB_q,function(x) {length(pop_PB_all[pop_PB_all < x ])})
popB_A_x <- popB_A_x / length(pop_PB_all)
popB_A_y <- sapply(popB_q,function(x) {length(popAB$PA[ popAB$PB < x & popA_sig & common_AB ])})
popB_A_y <- popB_A_y/num_A_inB_sig
quartz()
plot(popB_A_x,popB_A_y,type="l",col="red",lwd=2,main="population B against A",xlab="fraction of common SNPs in population B below pValue threshold", ylab="fract. of common SNPs sig. in A and below threshold in B",log="x")
# create random vector of SNPs from common ones 
rcom_A_idx = sample(which(common_AB),num_A_inB_sig )
rcom_A=rep(FALSE,length(common_AB))
rcom_A[rcom_A_idx]=TRUE
rcom_B_idx = sample(which(common_AB),num_B_inA_sig)
rcom_B=rep(FALSE,length(common_AB))
rcom_B[rcom_B_idx]=TRUE
# for B against A
popA_Br_y <- sapply(popA_q,function(x) length(popAB$PB[ popAB$PA < x & rcom_B ]))
popA_Br_y <- popA_Br_y/num_B_inA_sig
lines(popA_B_x,popA_Br_y,type="l",col="grey",lwd=1)
# for A against B
popB_Ar_y <- sapply(popB_q,function(x) length(popAB$PA[ popAB$PB < x & rcom_A ]))
popB_Ar_y <- popB_Ar_y/num_A_inB_sig
lines(popB_A_x,popB_Ar_y,type="l",col="grey",lwd=1)
# with ranks of the pValues against significance
# getting the ranks of each p Value
rank_pA = rank(popAB$PA,ties.method="first")
rank_pB = rank(popAB$PB,ties.method="first")
# set steps for going through rank cutoffs
steps= c(seq(1,100,10),seq(101,1000,100),seq(1001,10000,1000),seq(10001,100000,10000),seq(100001,2.2e6,1e6))
# B against A
pop_rB_Ay <- sapply(steps,function(x) {length(popAB$PA[ rank_pB < x & popA_sig & common_AB ]) })
pop_rB_Ay <- pop_rB_Ay/num_A_inB_sig
# random sample
pop_rB_Ary <- sapply(steps,function(x) {length(popAB$PA[ rank_pB < x & rcom_A  ]) })
pop_rB_Ary <- pop_rB_Ary/num_A_inB_sig
quartz()
plot(steps,pop_rB_Ay,type="l",col="red",lwd=2,main="population B against A",xlab="SNPs in B with rank < x", ylab="common SNPs found significant in A", log="x")
plot(steps,pop_rB_Ay,type="l",col="red",lwd=2,main="population B against A",xlab="SNPs in B with rank < x", ylab="common SNPs found significant in A")
# add random sample line
lines(steps,pop_rB_Ary,type="l",col="blue",lwd=1)
# A against B
pop_rA_By <- sapply(steps,function(x) {length(popAB$PB[ rank_pA < x & popB_sig & common_AB ]) })
pop_rA_By <- pop_rA_By/num_B_inA_sig
quartz()
plot(steps,pop_rA_By,type="l",col="red",lwd=2,main="population A against B",xlab="SNPs in A with rank < x", ylab="common SNPs found significant in B", log="x")
# plot both in one:
quartz()
plot(steps,pop_rA_By,type="l",col="red",lwd=2,xlab="SNPs in A with rank < x", ylab="common SNPs found significant in B", log="x")
lines(steps,pop_rB_Ay, col="blue",lwd=2)
legend("topleft",c("A:Vienna 2010, B:Bolzano 2011","A:Bolzano 2011, B:Vienna 2010"),col=c("red","blue"))
dev.copy2pdf(file="sign_overlap_fdr05_logscale.pdf")

# read the mckay data
mcKay=read.table("/Volumes/Temp/Lukas/Data/Huang2012/sd01.txt",header=T,skip=1)
mcKay$SNP <- paste(mcKay$CHR,mcKay$BP,sep=":")
#overlap with A, B and sigA, sigB
popABM=merge(popAB[,c(1,2,3)],mcKay[,c(10,9)],by="SNP",all=T)
popA_sigM = popABM$PA < sigA & ! is.na(popABM$PA)
popB_sigM = popABM$PB < sigB & ! is.na(popABM$PB)
common_AB_M = !(is.na(popABM$PA) | is.na(popABM$PB))
common_ABM =  !(is.na(popABM$PA) | is.na(popABM$PB) | is.na(popABM$P))
shared_ABM = popABM[common_ABM,]
dim(shared_ABM)
#[1] 154   4
min(shared_ABM$PA)
#[1] 0.002933016
min(shared_ABM$PB)
#[1] 1.254746e-06
sig_ABM = popABM$P[common_ABM & popA_sigM & popB_sigM]
# length 0
sig_AM = popABM[ popA_sigM & !(is.na(popABM$PA) | is.na(popABM$P)),]
# length 0
sig_BM = popABM[ popB_sigM & !(is.na(popABM$PB) | is.na(popABM$P)),]
# length 0
sig_ABM_nonoverlap = popABM[ popA_sigM | popB_sigM | ! is.na(popABM$P),]

