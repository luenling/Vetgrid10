setwd("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/CMH")
cmh_joined=read.table("fem_vie2010__pI25_pII25_pIII25_pIel_pIIel_pIIIel_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_masked_mct20_mcv25.gwas",header=T)
fst_joined=read.table("Fst/fem_vie2010__pI25_pII25_pIII25_pIel_pIIel_pIIIel_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_masked.fst1:4_2:5_3:6_7:10_8:11_9:12.avg.fsts",header=T)
cmh_joined$SNP= paste(cmh_joined$CHR,cmh_joined$BP,sep=":")
fst_joined$SNP= paste(fst_joined$CHR,fst_joined$BP,sep=":")
popAB=merge(cmh_joined[,c(3,4)],fst_joined[,c(3,4,5)],by="SNP",all=T)
rank_fst = rank(popAB$MEAN,ties.method="first")
rank_fst = max(rank_fst[ ! is.na(popAB$MEAN) ]) - rank_fst
rank_cmh= rank(popAB$P,ties.method="first")

rank_cmh[ which( rank_fst >=0 & rank_fst <= 100 )]
steps= c(seq(1,100,1),seq(101,200,10),seq(201,1000,100),seq(1001,10000,1000),seq(10001,100000,10000),seq(100001,2.2e6,1e6))
# B against A
pop_rB_Ay <- sapply(steps,function(x) {length(popAB$P[ rank_fst < x & rank_cmh < x & rank_fst >= 0 & rank_cmh >= 0])/x })
quartz()
plot(steps[1:120],pop_rB_Ay[1:120],type="l",col="red",lwd=2,main="rank CMH against fst",xlab="SNPs cmh with rank < x", ylab=" fraction of SNPs with rank  < x")

setwd("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/CMH")
cmh_joined=read.table("fem_vie2010__pI25_pII25_pIII25_pIel_pIIel_pIIIel_boz2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_realigned_q20_masked_mct16_mcv25.gwas",header=T)
pV1=read.delim("shuffles/fem_vie2010__pI25_pII25_pIII25_pIel_pIIel_pIIIel_boz2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_realigned_q20_masked_mct16_mcv25.cmhout_base_comp__1_2_3_7_8_9_4_5_6_10_11_12.gwas",header=T,colClasses=c('NULL', 'NULL', 'NULL', 'numeric'))
pV2=read.delim("shuffles/fem_vie2010__pI25_pII25_pIII25_pIel_pIIel_pIIIel_boz2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_realigned_q20_masked_mct16_mcv25.cmhout_base_comp__1_3_2_9_7_8_4_6_5_12_10_11.gwas",header=T,colClasses=c('NULL', 'NULL', 'NULL', 'numeric'))
pV3=read.delim("shuffles/fem_vie2010__pI25_pII25_pIII25_pIel_pIIel_pIIIel_boz2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_realigned_q20_masked_mct16_mcv25.cmhout_base_comp__1_7_2_8_3_9_4_10_5_11_6_12.gwas",header=T,colClasses=c('NULL', 'NULL', 'NULL', 'numeric'))
pV4=read.delim("shuffles/fem_vie2010__pI25_pII25_pIII25_pIel_pIIel_pIIIel_boz2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_realigned_q20_masked_mct16_mcv25.cmhout_base_comp__1_8_2_9_3_7_4_11_5_12_6_10.gwas",header=T,colClasses=c('NULL', 'NULL', 'NULL', 'numeric'))
pV5=read.delim("shuffles/fem_vie2010__pI25_pII25_pIII25_pIel_pIIel_pIIIel_boz2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_realigned_q20_masked_mct16_mcv25.cmhout_base_comp__1_9_2_7_3_8_4_12_5_10_6_11.gwas",header=T,colClasses=c('NULL', 'NULL', 'NULL', 'numeric'))
pV_shuf=scan(file="shuffles/all_shuffle_ps.dat")
quartz()
qqplot(-1*log10(pV_shuf),-1 * log10(cmh_joined$P),xlab="null P (10 shuffles)", ylab="observed P", main="Vienna 2010 & Italy 2011")
abline(a=0,b=1,col="red")






