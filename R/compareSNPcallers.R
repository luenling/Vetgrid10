setwd("/Volumes/Temp/Lukas/Data/SNPcall_Test/Results_AF")
library("VennDiagram")

551815 filtered_haps_2L_only.afs
  170900 test_recal_pl50.recal.VQSLODmin3.af
  362565 real.recal.filter.af
 160312 hap_gatk_50_joined_polarised.af
  164078 popoo_recal_gatk_50_joined_polarised.af
  337276 hap_popoo_joined_polarised.af
  154946 hap_gatk_50_popoo_recal_joined_polarised.af


  362281 real.filter.af
  361577 popoo_recal_nonrecal_joined_polarised.af
quartz()
draw.triple.venn(551815,170900,362565, 160312,164078,337276,154946,category=c("Reference","GATK","PoPoolation"),col=c("black"),fill=c("red","green","yellow"),cex=1.5,cat.cex=1.7)
dev.copy2pdf(file="venn_comparison.pdf")

haps=read.table("filtered_haps_2L_only_repmasked.afs",header=F)
colnames(haps)=c("CHR","BPS","A","AF","N")
haps=subset(haps,select=-c(CHR))
haps_NR=read.table("haps_nonref_repmasked.af",header=F)
colnames(haps_NR)=c("CHR","BPS","An","AnF","Nn")
haps_NR=subset(haps_NR,select=-c(CHR))
haps_gatk=read.table("hap_gatk_50_joined_polarised_repmasked.af",header=F)
colnames(haps_gatk)=c("CHR","BPS","A","AF","N","G","GF")
haps_gatk=subset(haps_gatk,select=-c(CHR))
haps_gatk_popoo=read.table("hap_gatk_50_popoo_recal_joined.af",header=F)
colnames(haps_gatk_popoo)=c("CHR","BPS","A","AF","N","G","GF","P","PF","PC")
haps_gatk_popoo=subset(haps_gatk_popoo,select=-c(CHR))
haps_popoo=read.table("hap_popoo_joined_polarised.af",header=F)
colnames(haps_popoo)=c("CHR","BPS","A","AF","N","P","PF","PC")
haps_popoo=subset(haps_popoo,select=-c(CHR))
popoo_gatk=read.table("popoo_recal_gatk_50_joined_polarised.af",header=F)
colnames(popoo_gatk)=c("CHR","BPS","P","PF","PC","G","GF")
popoo_gatk=subset(popoo_gatk,select=-c(CHR))
haps_gatk$Ratio=abs(haps_gatk$GF - haps_gatk$AF)/haps_gatk$AF
haps_popoo$Ratio=abs(haps_popoo$PF - haps_popoo$AF)/haps_popoo$AF
haps_gatk$Dev=abs(haps_gatk$GF - haps_gatk$AF)
haps_popoo$Dev=abs(haps_popoo$PF - haps_popoo$AF)
haps_popoo$Bial = nchar(as.character(haps_popoo$A)) == 2 & nchar(as.character(haps_popoo$P)) == 2
haps_gatk$Bial = nchar(as.character(haps_gatk$A)) == 2 & nchar(as.character(haps_gatk$G)) == 2
popoo=read.table("real.recal.filter.af",header=F)
colnames(popoo)=c("CHR","BPS","P","PF","PC")
popoo=subset(popoo,select=-c(CHR))
gatk=read.table("test_recal_pl50.recal.VQSLODmin3_repmasked.af",header=F)
colnames(gatk)=c("CHR","BPS","G","GF")
gatk=subset(gatk,select=-c(CHR))

merge_all=merge(haps,haps_gatk[,c("BPS","G","GF")],by = c("BPS"), all = T)
merge_all=merge(merge_all,haps_popoo[,c("BPS","P","PF")],by = c("BPS"), all = T)
merge_all=merge(merge_all,haps_gatk_popoo[,c("BPS","P","G")],by = c("BPS"), all = T, suffixes = c("",".all"))
merge_all=merge(merge_all,popoo_gatk[,c("BPS","P","PF","G","GF")],by = c("BPS"), all = T, suffixes = c("",".fp"))
merge_all=merge(merge_all,popoo[,c("BPS","P","PF")],by = c("BPS"), all = T, suffixes = c("",".tot"))
merge_all=merge(merge_all,gatk[,c("BPS","G","GF")],by = c("BPS"), all = T, suffixes = c("",".tot"))
merge_all=merge(merge_all,haps_NR[,c("BPS","An","AnF","Nn")],by = c("BPS"), all = T)

hap_NA = is.na(merge_all$A)
hapNR_NA = is.na(merge_all$An)
hap_tot_NA = hap_NA & hapNR_NA
gatk_NA = is.na(merge_all$G.tot)
popoo_NA = is.na(merge_all$P.tot)
hgp_NA = is.na(merge_all$P.all)
gp_NA = is.na(merge_all$P.fp)
hg_NA= is.na(merge_all$G)
hp_NA= is.na(merge_all$P)
hap_trial = ! hap_NA & nchar(as.character(merge_all$A)) > 2
hapNR_trial = ! hapNR_NA & nchar(as.character(merge_all$An)) > 2
hap_tot_trial = hap_trial | hapNR_trial
g_trial =  ! gatk_NA & nchar(as.character(merge_all$G.tot)) > 2
p_trial =  ! popoo_NA & nchar(as.character(merge_all$P.tot)) > 2
all_trial = g_trial | p_trial | hap_tot_trial
merge_all$N[is.na(merge_all$N)] = -1
merge_all$Nn[is.na(merge_all$Nn)] = -1
excl_pop_NR = ( hap_NA | hp_NA ) & ! ( hapNR_NA | hap_tot_trial | p_trial | popoo_NA)  
excl_gatk_NR = ( hap_NA | hg_NA ) & ! ( hapNR_NA | hap_tot_trial |  g_trial | gatk_NA )   
excl_gatk_pop_NR =  ( hap_NA | hp_NA | hg_NA ) & ! ( hapNR_NA | all_trial | gatk_NA |  popoo_NA)   

num = 100
maf_thresh = 1.0

hap_low = ! hap_NA & merge_all$AF > maf_thresh # | ( hap_NA & ( merge_all$PF > maf_thresh |  merge_all$GF > maf_thresh ) ) 
num_ok = merge_all$N <= num & merge_all$Nn <= num
all_accept = ! hap_tot_trial &  num_ok & ! hap_low & ! all_trial


A=length(merge_all$BPS[ ! hap_tot_NA & all_accept ])
C_=length(merge_all$P.tot[ ! popoo_NA &  ! p_trial  & all_accept  ])
B=length(merge_all$G.tot[ ! gatk_NA &  ! g_trial & all_accept])
AB=length(merge_all$A[ ( ! hg_NA | excl_gatk_NR ) &  ! g_trial  & all_accept  ])
BC=length(merge_all$P.fp[ ! gp_NA & ! excl_gatk_pop_NR  & ! g_trial & ! p_trial & all_accept ])
AC=length(merge_all$A[ (! hp_NA  | excl_pop_NR ) &  ! p_trial & all_accept ])
ABC=length(merge_all$A[ (! hgp_NA | excl_gatk_pop_NR ) & ! g_trial & ! p_trial & all_accept ])
c(A,B,C_,AB,BC,AC,ABC)
FP_popoo=C_-(BC+AC-ABC)
FP_gatk=B-(BC+AB-ABC)
FP_popoo_perc=FP_popoo/C_
FP_gatk_perc=FP_gatk/B
Rest_hap=A-(AB+AC-ABC)
Res_hap_perc=(A-(AB+AC-ABC))/A
prettyNum(c(FP_popoo,FP_popoo_perc,FP_gatk,FP_gatk_perc,Rest_hap,Res_hap_perc,BC-ABC,(BC-ABC)/B))
quartz()
draw.triple.venn(A,B,C_,AB,BC,AC,ABC,category=c("Reference","GATK","PoPoolation"),col=c("black"),fill=c("red","green","yellow"),cex=1.5,cat.cex=1.7)
dev.copy2pdf(file="venn_comparison_N_100_MAF_1.pdf")


merge_all


# alternative:
bia <- levels(haps_popoo$A)[nchar(levels(haps_popoo$A)) == 2]
biag <- levels(haps_popoo$G)[nchar(levels(haps_popoo$G)) == 2]
hpp2 <- subset(haps_popoo,A%in%bia)
# or 
haps_gatk$Bial = haps_gatk$A %in% bia &  haps_gatk$G %in% biag
quartz()
boxplot(haps_gatk$Ratio[haps_gatk$N <= 5 & haps_gatk$Bial],haps_popoo$Ratio[haps_popoo$N <= 5 & haps_popoo$Bial],names=c("GATK","PoPoolation"),notch=F,ylab="rel. deviation of MAF", main="SNP Calling, relative deviation")
dev.copy2pdf(file="maf_deviation_relative.pdf")
quartz()
boxplot(haps_gatk$Dev[haps_gatk$N <= 5 & haps_gatk$Bial],haps_popoo$Dev[haps_popoo$N <= 5 & haps_popoo$Bial],names=c("GATK","PoPoolation"),notch=F,ylab="abs. deviation of MAF", main="SNP Calling, absolute deviation")
dev.copy2pdf(file="maf_deviation_absolute.pdf")
  