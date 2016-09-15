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
popB_q = quantile(popAB$PB,probs=c(seq(0,0.01,0.0001),seq(0.01,1,0.01)),na.rm=T )
popA_q = quantile(popAB$PA,probs=c(seq(0,0.01,0.0001),seq(0.01,1,0.01)),na.rm=T )
# popA against popB
pop_PA_all = na.omit(popAB$PA)
popA_B_x <- sapply(popA_q,function(x) {length(pop_PA_all[pop_PA_all < x ])})
popA_B_x <- popA_B_x / popA_all
popA_B_y <- sapply(popA_q,function(x) {length(popAB$PB[ popAB$PA < x & popB_sig & common_AB ])})
popA_B_y <- popA_B_y/num_B_inA_sig
quartz()
plot(popA_B_x,popA_B_y,type="l",col="red",lwd=2,main="population A against B",xlab="fraction of SNPs in population A below pValue threshold", ylab="fract. of common SNPs sig. in B and below threshold in A")

# popA against popB
pop_PB_all = na.omit(popAB$PB)
popB_A_x <- sapply(popB_q,function(x) {length(pop_PB_all[pop_PB_all < x ])})
popB_A_x <- popB_A_x / popB_all
popB_A_y <- sapply(popB_q,function(x) {length(popAB$PA[ popAB$PB < x & popA_sig & common_AB ])})
popB_A_y <- popB_A_y/num_A_inB_sig
quartz()
plot(popB_A_x,popB_A_y,type="l",col="red",lwd=2,main="population B against A",xlab="fraction of SNPs in population B below pValue threshold", ylab="fract. of common SNPs sig. in A and below threshold in B")

# with ranks of the pValues against significance
# getting the ranks of each p Value
rank_pA = rank(popAB$PA,ties.method="first")
rank_pB = rank(popAB$PB,ties.method="first")
# set steps for going through rank cutoffs
steps= c(seq(1,100,10),seq(101,1000,100),seq(1001,10000,1000),seq(10001,100000,10000),seq(100001,2.2e6,1e6))
# B against A
pop_rB_Ay <- sapply(steps,function(x) {length(popAB$PA[ rank_pB < x & popA_sig & common_AB ]) })
pop_rB_Ay <- pop_rB_Ay/num_A_inB_sig
quartz()
plot(steps,pop_rB_Ay,type="l",col="red",lwd=2,main="population B against A",xlab="SNPs in B with rank < x", ylab="common SNPs found significant in A", log="x")
# A against B
pop_rA_By <- sapply(steps,function(x) {length(popAB$PB[ rank_pA < x & popB_sig & common_AB ]) })
pop_rA_By <- pop_rA_By/num_B_inA_sig
quartz()
plot(steps,pop_rA_By,type="l",col="red",lwd=2,main="population A against B",xlab="SNPs in A with rank < x", ylab="common SNPs found significant in A", log="x")


 
