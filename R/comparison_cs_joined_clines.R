setwd("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/Joined_Analysis")
lat_joined=read.table("lat_FLOR_PENN_MAINE_vie2010_bolz2011.af",header=F)
colnames(lat_joined) <- unlist(strsplit("CHR_BPS_Alleles_Flo_nFlo_Penn_nPenn_Maine_nMaine_Cvie_nCvie_Nvie_nNvie_Pvie_Cboz_nCboz_Nboz_nNboz","_"))
seas_joined=read.table("seasonal_July_Sept_Oct_Nov_vie2010_bolz2011.af",header=F)
lat_ranks = (length(lat_joined$Pvie)-rank(lat_joined$Pvie))
lat_cols=c("Flo","Penn","Maine","Cvie","Nvie","Cboz","Nboz")
lat_af = lat_joined[,lat_cols]
lat_af = abs((0.5<lat_af$Cvie)-lat_af)
colnames(seas_joined) <- unlist(strsplit("chrom_pos_snp_July_nJuly_Sept_nSept_Oct_nOct_Nov_nNov_Cvie_nCvie_Nvie_nNvie_Pvie_Cboz_nCboz_Nboz_nNboz","_"))
seas_ranks = (length(seas_joined$Pvie)-rank(seas_joined$Pvie))
seas_cols=c("July","Sept","Oct","Nov","Cvie","Nvie","Cboz","Nboz")
seas_af = seas_joined[,seas_cols]
seas_af = abs((0.5<seas_af$Cvie)-seas_af)

ms_joined=read.table("m_schaeffer_KEN_ETH_ATH_BAR_DIJ_ESX_BOW_LNZ_vie2010_bolz2011.af",header=F)
colnames(ms_joined) <- unlist(strsplit("CHR_BPS_Alleles_KEN_nKEN_ETH_nETH_ATH_nATH_BAR_nBAR_DIJ_nDIJ_ESX_nESX_BOW_nBOW_LNZ_nLNZ_Cvie_nCvie_Nvie_nNvie_Pvie_Cboz_nCboz_Nboz_nNboz","_"))
ms_ranks = (length(ms_joined$Pvie)-rank(ms_joined$Pvie))
ms_cols=c("KEN","ETH","ATH","BAR","DIJ","ESX","BOW","LNZ","Cvie","Nvie","Cboz","Nboz")
ms_af = ms_joined[,ms_cols]
ms_af = abs((0.5<ms_af$Cvie)-ms_af)



quartz()
boxplot(lat_af[lat_ranks <= 150,lat_cols],names=lat_cols,notch=F,ylab="minor allele freq", main="Latitudinal cline")
dev.copy2pdf(file="latitudinal_joined_comp_box.pdf")
quartz()
boxplot(seas_af[,seas_cols],names=seas_cols,notch=F,ylab="minor allele freq", main="Seasonal data")
dev.copy2pdf(file="seasonal_joined_comp_box.pdf")
quartz()
bp <-boxplot(ms_af[,ms_cols],names=ms_cols,notch=F,xaxt = "n", ylab="minor allele freq", main="Schaeffer geographical data")
## Set up x axis with tick marks alone
axis(1, at=1:length(ms_cols) ,labels = FALSE)
## Plot x axis labels at default tick marks
text(1:length(ms_cols), par("usr")[3] - 0.025, srt = 45, adj = 1, labels = ms_cols, xpd = TRUE)
dev.copy2pdf(file="m_schaeffer_joined_comp_box.pdf")
# ranks 50
quartz()
boxplot(lat_af[lat_ranks <= 50,lat_cols],names=lat_cols,notch=F,ylab="minor allele freq", main="Latitudinal cline, best 50 in vie")
dev.copy2pdf(file="latitudinal_joined_comp_best_50_box.pdf")
quartz()
boxplot(seas_af[seas_ranks <= 50,seas_cols],names=seas_cols,notch=F,ylab="minor allele freq", main="Seasonal data, best 50 in vie")
dev.copy2pdf(file="seasonal_joined_comp_best_50_box.pdf")
quartz()
bp <-boxplot(ms_af[ms_ranks <= 50,ms_cols],names=ms_cols,notch=F,xaxt = "n", ylab="minor allele freq", main="Schaeffer geographical data, best 50 in vie")
## Set up x axis with tick marks alone
axis(1, at=1:length(ms_cols) ,labels = FALSE)
## Plot x axis labels at default tick marks
text(1:length(ms_cols), par("usr")[3] - 0.025, srt = 45, adj = 1, labels = ms_cols, xpd = TRUE)
dev.copy2pdf(file="m_schaeffer_joined_best_50_comp_box.pdf")


ms_cols=c("KEN","ETH","ATH","BAR","DIJ","ESX","BOW","LNZ","Cm","Nm")
lat_cols=c("Flo","Penn","Maine","Nm","Cm")
lat_cline=c("Flo","Penn","Maine")

lat_vie =read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/vie_2010_fdr_lat_FLOR_PENN_MAINE_only_vie_fdr.sync.af",header=F)
colnames(lat_vie) <- unlist(strsplit("CHR_BPS_Alleles_C1_C2_C3_N1_N2_N3_P_Flo_nFlo_Penn_nPenn_Maine_nMaine","_"))
lat_ita =read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita_2011_fdr_lat_FLOR_PENN_MAINE_only_ita_fdr.sync.af",header=F)
colnames(lat_ita) <- unlist(strsplit("CHR_BPS_Alleles_C1_C2_C3_N1_N2_N3_P_Flo_nFlo_Penn_nPenn_Maine_nMaine","_"))
cs_cols=c("C1","C2","C3")
b_cols=c("N1","N2","N3")
lat_vie$Cm=rowMeans(lat_vie[,cs_cols])
lat_vie$Nm=rowMeans(lat_vie[,b_cols])
lat_ita$Cm=rowMeans(lat_ita[,cs_cols])
lat_ita$Nm=rowMeans(lat_ita[,b_cols])
dir_lat_vie=lat_vie[,cs_cols]-lat_vie[,b_cols]
ms_cols=c("ATH","BAR","DIJ","LNZ","ESX","Nm","Cm")
ms_cline=c("ATH","BAR","DIJ","ESX","LNZ")
ms_vie =read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/vie_2010_fdr_ms_KEN_ETH_ATH_BAR_DIJ_ESX_BOW_LNZ_only_vie_fdr.sync.af",header=F)
colnames(ms_vie) <- unlist(strsplit("CHR_BPS_Allele_C1_C2_C3_N1_N2_N3_P_KEN_nKEN_ETH_nETH_ATH_nATH_BAR_nBAR_DIJ_nDIJ_ESX_nESX_BOW_nBOW_LNZ_nLNZ","_"))
ms_ita =read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita_2011_fdr_ms_KEN_ETH_ATH_BAR_DIJ_ESX_BOW_LNZ_only_ita_fdr_no_N.sync.af",header=F)
colnames(ms_ita) <- unlist(strsplit("CHR_BPS_Alleles_C1_C2_C3_N1_N2_N3_P_KEN_nKEN_ETH_nETH_ATH_nATH_BAR_nBAR_DIJ_nDIJ_ESX_nESX_BOW_nBOW_LNZ_nLNZ","_"))
cs_cols=c("C1","C2","C3")
b_cols=c("N1","N2","N3")
ms_vie$Cm=rowMeans(ms_vie[,cs_cols])
ms_vie$Nm=rowMeans(ms_vie[,b_cols])
ms_ita$Cm=rowMeans(ms_ita[,cs_cols])
ms_ita$Nm=rowMeans(ms_ita[,b_cols])
dim(lat_vie[rowMeans(dir_lat_vie) > 0,])
# in 33 the minor alleles are decreasing instead of increasing


# nonsyn
ms_cols=c("KEN","ETH","ATH","BAR","DIJ","ESX","BOW","LNZ","Cm","Nm")
lat_cols=c("Flo","Penn","Maine","Nm","Cm")
lat_cline=c("Flo","Penn","Maine")

lat_vie_ns =read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/vie_2010_fdr_lat_FLOR_PENN_MAINE_only_vie_fdr.sync.non_syn.af",header=F)
colnames(lat_vie_ns) <- unlist(strsplit("CHR_BPS_Alleles_C1_C2_C3_N1_N2_N3_P_Flo_nFlo_Penn_nPenn_Maine_nMaine","_"))
lat_ita_ns =read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita_2011_fdr_lat_FLOR_PENN_MAINE_only_ita_fdr.sync.non_syn.af",header=F)
colnames(lat_ita_ns) <- unlist(strsplit("CHR_BPS_Alleles_C1_C2_C3_N1_N2_N3_P_Flo_nFlo_Penn_nPenn_Maine_nMaine","_"))
cs_cols=c("C1","C2","C3")
b_cols=c("N1","N2","N3")
lat_vie_ns$Cm=rowMeans(lat_vie_ns[,cs_cols])
lat_vie_ns$Nm=rowMeans(lat_vie_ns[,b_cols])
lat_ita_ns$Cm=rowMeans(lat_ita_ns[,cs_cols])
lat_ita_ns$Nm=rowMeans(lat_ita_ns[,b_cols])
dir_lat_vie_ns=lat_vie_ns[,cs_cols]-lat_vie_ns[,b_cols]
ms_cols=c("ATH","BAR","DIJ","LNZ","ESX","Nm","Cm")
ms_cline=c("ATH","BAR","DIJ","ESX","LNZ")
ms_vie_ns =read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/vie_2010_fdr_ms_KEN_ETH_ATH_BAR_DIJ_ESX_BOW_LNZ_only_vie_fdr.sync.non_syn.af",header=F)
colnames(ms_vie_ns) <- unlist(strsplit("CHR_BPS_Alleles_C1_C2_C3_N1_N2_N3_P_KEN_nKEN_ETH_nETH_ATH_nATH_BAR_nBAR_DIJ_nDIJ_ESX_nESX_BOW_nBOW_LNZ_nLNZ","_"))
ms_ita_ns =read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita_2011_fdr_ms_KEN_ETH_ATH_BAR_DIJ_ESX_BOW_LNZ_only_ita_fdr_no_N.sync.non_syn.af",header=F)
colnames(ms_ita_ns) <- unlist(strsplit("CHR_BPS_Alleles_C1_C2_C3_N1_N2_N3_P_KEN_nKEN_ETH_nETH_ATH_nATH_BAR_nBAR_DIJ_nDIJ_ESX_nESX_BOW_nBOW_LNZ_nLNZ","_"))
cs_cols=c("C1","C2","C3")
b_cols=c("N1","N2","N3")
ms_vie_ns$Cm=rowMeans(ms_vie_ns[,cs_cols])
ms_vie_ns$Nm=rowMeans(ms_vie_ns[,b_cols])
ms_ita_ns$Cm=rowMeans(ms_ita_ns[,cs_cols])
ms_ita_ns$Nm=rowMeans(ms_ita_ns[,b_cols])
dim(lat_vie_ns[rowMeans(dir_lat_vie_ns) > 0,])

quartz()
boxplot(lat_vie[,"Flo"] - lat_vie[,lat_cols],names=lat_cols,notch=F,ylab="AF difference to Flo", main="North American cline, Vienna")
dev.copy2pdf(file="lat_fdr_vienna_diff_box.pdf")
quartz()
boxplot(ms_vie[,"ATH"] - ms_vie[,ms_cols],names=ms_cols,notch=F,ylab="AF difference to ATH", main="European cline, Vienna")
dev.copy2pdf(file="ms_fdr_vienna_diff_box.pdf")
quartz()
boxplot(lat_ita[,"Flo"] - lat_ita[,lat_cols],names=lat_cols,notch=F,ylab="AF difference to Flo", main="North American cline, Italy")
dev.copy2pdf(file="lat_fdr_italy_diff_box.pdf")
quartz()
boxplot(ms_ita[,"ATH"] - ms_ita[,ms_cols],names=ms_cols,notch=F,ylab="AF difference to ATH", main="European cline, Italy")
dev.copy2pdf(file="ms_fdr_italy_diff_box.pdf")

quartz()
boxplot(1-lat_vie[,lat_cols],names=lat_cols,notch=F,ylab="minor allele freq", main="North American cline, Vienna")
dev.copy2pdf(file="lat_fdr_vienna_comp_box.pdf")
quartz()
boxplot(1-ms_vie[,ms_cols],names=ms_cols,notch=F,ylab="minor allele freq", main="European cline, Vienna")
dev.copy2pdf(file="ms_fdr_vienna_comp_box.pdf")
quartz()
boxplot(1-lat_ita[,lat_cols],names=lat_cols,notch=F,ylab="minor allele freq", main="North American cline, Italy")
dev.copy2pdf(file="lat_fdr_italy_comp_box.pdf")
quartz()
boxplot(1-ms_ita[,ms_cols],names=ms_cols,notch=F,ylab="minor allele freq", main="European cline, Italy")
dev.copy2pdf(file="ms_fdr_italy_comp_box.pdf")

latitudes=c(Flo=25.53,Penn=39.883,Maine=42.3,KEN=-0.4252,ETH=-9.4969,ATH=37.9778,BAR=41.3857,DIJ=47.3239,ESX=51.87,BOW=42.3,LNZ=48.3067)

fill_lat_frame = function(lat_tab,cline,latitudes,poly_f=0.95){
        lat_tab_poly=lat_tab[(rowSums(lat_tab[,cline] > poly_f) < length(cline)),]
	df_ret=data.frame(MAJ=as.integer(),MIN=as.integer(),LAT=as.numeric(),LOC=factor())
	for(i in cline){
		n_i=paste("n",i,sep="")
		df=matrix(c(round(lat_tab_poly[,i] * lat_tab_poly[,n_i],digits=0),lat_tab_poly[,n_i]-round(lat_tab_poly[,i] * lat_tab_poly[,n_i],digits=0), rep(latitudes[[i]],length(lat_tab_poly[,i]))), ncol=3)
		colnames(df)=c("MAJ","MIN","LAT")
		df=as.data.frame(df)
		df$LOC=as.character(rep(i,length(lat_tab_poly[,i])))
		df_ret=rbind(df_ret,df)
	}
	df_ret$Freq=with(df_ret, (MIN/(MAJ+MIN)))
	return(df_ret)
}
data_ms_ita=fill_lat_frame(ms_ita,ms_cline,latitudes)
data_ms_ita$LOC=factor(data_ms_ita$LOC)
data_ms_ita$days_above_10=as.numeric(d_a_10[data_ms_ita$LOC])
data_ms_ita$m_min_mean=as.numeric(m_min_mean[data_ms_ita$LOC])
data_ms_ita$y_mean=as.numeric(y_mean[data_ms_ita$LOC])

#only take polymorph data
ms_ita_poly = ms_ita[! rowSums(ms_ita[,ms_cline] > 0.95) == length(ms_cline) & rowMeans(ms_ita[,b_cols] - ms_ita[,cs_cols]) > 0,]
dim(ms_ita)
dim(ms_ita_poly)
data_ms_ita=fill_lat_frame(ms_ita_poly,ms_cline,latitudes)
data_ms_ita$LOC=factor(data_ms_ita$LOC)
data_ms_ita$days_above_10=as.numeric(d_a_10[data_ms_ita$LOC])
data_ms_ita$m_min_mean=as.numeric(m_min_mean[data_ms_ita$LOC])
data_ms_ita$y_mean=as.numeric(y_mean[data_ms_ita$LOC])
ms_vie_poly = ms_vie[! rowSums(ms_vie[,ms_cline] > 0.95) == length(ms_cline) & rowMeans(ms_vie[,b_cols] - ms_vie[,cs_cols]) > 0,]
dim(ms_vie)
dim(ms_vie_poly)
data_ms_vie=fill_lat_frame(ms_vie_poly,ms_cline,latitudes)
data_ms_vie$LOC=factor(data_ms_vie$LOC)
data_ms_vie$days_above_10=as.numeric(d_a_10[data_ms_vie$LOC])
lat_ita_poly = lat_ita[! rowSums(lat_ita[,lat_cline] > 0.95) == length(ms_cline) & rowMeans(lat_ita[,b_cols] - lat_ita[,cs_cols]) > 0,]
lat_ita_poly = lat_ita[! rowSums(lat_ita[,lat_cline] > 0.95) == length(ms_cline),]
dim(lat_ita)
dim(lat_ita_poly)
data_lat_ita=fill_lat_frame(lat_ita_poly,lat_cline,latitudes)
data_lat_ita$LOC=factor(data_lat_ita$LOC)
lat_vie_poly = lat_vie[! rowSums(lat_vie[,lat_cline] > 0.95) == length(lat_cline) & rowMeans(lat_vie[,b_cols] - lat_vie[,cs_cols]) > 0,]
lat_vie_poly = lat_vie[! rowSums(lat_vie[,lat_cline] > 0.95) == length(lat_cline),]
dim(lat_vie)
dim(lat_vie_poly)
data_lat_vie=fill_lat_frame(lat_vie_poly,lat_cline,latitudes)
data_lat_vie$LOC=factor(data_lat_vie$LOC)
data_lat_vie$Freq=with(data_lat_vie, (MIN/(MAJ+MIN)))
data_lat_ita$Freq=with(data_lat_ita, (MIN/(MAJ+MIN)))
lm_lat_ita=lm(sqrt(Freq) ~ 1 + LAT, data=data_lat_ita)
lm_lat_vie=lm(sqrt(Freq) ~ 1 + LAT, data=data_lat_vie)

#adding random factor similat to latitude:
rand_lat=runif(length(latitudes),min=0,max=100)
names(rand_lat)=names(latitudes)
ms_levels= 1:length(ms_cline)
names(ms_levels)=c("ATH","BAR","DIJ", "LNZ","ESX")
data_ms_ita$LocLev=as.numeric(ms_levels[data_ms_ita$LOC])
data_ms_ita$Rand=as.numeric(rand_lat[data_ms_ita$LOC])
data_ms_ita$Freq =with(data_ms_ita, (MIN/(MAJ+MIN))) 
data_ms_vie$LocLev=as.numeric(ms_levels[data_ms_vie$LOC])
data_ms_vie$Rand=as.numeric(rand_lat[data_ms_vie$LOC])
data_ms_vie$Freq =with(data_ms_vie, (MIN/(MAJ+MIN))) 
quartz(width=12,height=12)
scatterplotMatrix(~ LAT + days_above_10 + y_mean + m_min_mean + Freq | LOC, data=data_ms_ita, diagonal="density",smoother=F,spread=T)
dev.copy2pdf(file="ms_ita_scatter_matrix_temperatures.pdf")

quartz(width=12,height=12)
par(mfrow=c(2,2))
lm_ms_ita=lm(sqrt(Freq) ~ 1 + LAT, data=data_ms_ita)
lm_ms_vie=lm(sqrt(Freq) ~ 1 + LAT, data=data_ms_vie)
summary(lm_ita)$"coefficients"["LAT",]
summary(lm_ita)$"adj.r.squared"
a=summary(lm_ita)$"fstatistic"
pf(a[1],a[2],a[3],lower.tail=F)
anova(lm_ita)
plot(lm(sqrt(Freq) ~ 1 + LAT, data=data_ms_ita))
summary(lm(sqrt(Freq) ~ 1 + LAT, data=data_ms_ita))
plot(lm(sqrt(Freq) ~ 1 + LocLev, data=data_ms_ita))
summary(lm(sqrt(Freq) ~ 1 + LocLev, data=data_ms_ita))
summary(lm(sqrt(Freq) ~ 1 + y_mean, data=data_ms_ita))
summary(lm(sqrt(Freq) ~ 1 + days_above_10, data=data_ms_ita))
plot(lm(sqrt(Freq) ~ 1 + LocLev, data=data_ms_ita))
summary(lm(sqrt(Freq) ~ 1 + LocLev, data=data_ms_ita))
plot(lm(sqrt(Freq) ~ 1 + LOC, data=data_ms_ita))
summary(lm(sqrt(Freq) ~ 1 + LOC, data=data_ms_ita))
plot(lm(asin(Freq) ~ 1 + LAT, data=data_ms_ita))

library("lme4")

glm_ms_ita<-glm(cbind(MIN,MAJ) ~ 1 + LAT, data=data_ms_ita, family=binomial)
glm_ms_ita_nolat<-glm(cbind(MIN,MAJ) ~ 1 , data=data_ms_ita, family=binomial)
anova(glm_ms_ita_nolat,glm_ms_ita)

summary(glm_ms_ita)
glm_ms_ita_rand<-glm(cbind(MIN,MAJ) ~ 1 + Rand + LAT , data=data_ms_ita, family=binomial)
summary(glm_ms_ita)
summary(aov(glm_ms_ita_rand))

quartz(width=12,height=12)
par(mfrow=c(2,2))
plot(glm_ms_ita)
anova(glm_ms_ita)
glm_ms_vie<-glm(cbind(MIN,MAJ) ~ 1 + LAT, data=data_ms_vie, family=binomial)
summary(glm_ms_vie)
glm_ms_ita_da10<-glm(cbind(MIN,MAJ) ~ 1 +  days_above_10 , data=data_ms_ita, family=binomial)
summary(glm_ms_ita_da10)
glm_ms_ita_mmm<-glm(cbind(MIN,MAJ) ~ 1 +  m_min_mean , data=data_ms_ita, family=binomial)
summary(glm_ms_ita_mmm)
glm_ms_ita_ym<-glm(cbind(MIN,MAJ) ~ 1 +  y_mean , data=data_ms_ita, family=binomial)
summary(glm_ms_ita_ym)
glm_ms_vie_da10<-glm(cbind(MIN,MAJ) ~ 1 + days_above_10 , data=data_ms_vie, family=binomial)
summary(glm_ms_vie_da10)
glm_lat_ita<-glm(cbind(MIN,MAJ) ~ 1 + LAT, data=data_lat_ita, family=binomial)
summary(glm_lat_ita)
glm_lat_vie<-glm(cbind(MIN,MAJ) ~ 1 + LAT, data=data_lat_vie, family=binomial)
summary(glm_lat_vie)

library(car)
quartz(width=12,height=12)
data_ms_ita$Freq =with(data_ms_ita, (MIN/(MAJ+MIN))) 
scatterplotMatrix(~ LAT + days_above_10 + y_mean + m_min_mean + Freq | LOC, data=data_ms_ita, diagonal="density",smoother=F,spread=T)
dev.copy2pdf(file="ms_ita_scatter_matrix_temperatures.pdf")

quartz(width=12,height=12)
par(mfrow=c(2,2))
plot(lm(sqrt(Freq) ~ 1 + days_above_10, data=data_ms_ita))
summary(lm(sqrt(Freq) ~ 1 + days_above_10, data=data_ms_ita))
plot(lm(asin(Freq) ~ 1 + days_above_10, data=data_ms_ita))

plot(glm(Freq ~ 1 + LAT, data=data_ms_ita, family=poisson))



# reading temperature data 
get_days_above = function(temp_tab,temp){
	a=with(temp_tab, tapply(QCP,as.Date(Date,format="%m"), min) )
	return(length(a[a>=temp]))	
}

get_monthly_mean_min = function(temp_tab){
	a=tapply(tapply(temp_tab$QCP, as.Date(temp_tab$Date), min), months.Date(as.Date(names(tapply(temp_tab$QCP, as.Date(temp_tab$Date), min)))),mean)
	return(a)	
}

temp=10
temp_locs = c("athens_temp","barc_temp","dijon_temp","linz_temp","essex_temp")
names(temp_locs)=c("ATH","BAR","DIJ","LNZ","ESX")
d_a_10=list()
m_min_mean=list()
y_mean=list()
for(x in names(temp_locs) ){
	d_a_10[[x]]=get_days_above(get(temp_locs[[x]]),temp)
	m_min_mean[[x]]=mean(get_monthly_mean_min(get(temp_locs[[x]])))
	y_mean[[x]]=mean(with(get(temp_locs[[x]]), tapply(QCP, months.Date(Date), mean) ))
}
d_a_10=as.vector(d_a_10)
m_min_mean=as.vector(m_min_mean)
y_mean=as.vector(y_mean)
athens_temp=read.table("/Volumes/Temp/Lukas/Data/M_Schaeffer_PE_bams/Temperatures/temp_athens_helik_2005_2009.dat", header=T, skip=0,sep="\t",colClasses=c(NA,"character","character",NA))
athens_temp$Date=as.POSIXlt(paste(athens_temp$Date,athens_temp$HrMn),format="%Y%m%d %H%M")
with(athens_temp, tapply(QCP, months.Date(athens_temp$Date), min) )
with(athens_temp, tapply(QCP, months.Date(athens_temp$Date), mean) )
quartz()
boxplot(athens_temp$QCP~format(athens_temp$Date,"%m"),main="Athens",ylim=c(-15,40))
dev.copy2pdf(file="/Volumes/Temp/Lukas/Data/M_Schaeffer_PE_bams/Temperatures/temp_box_athens.pdf")

essex_temp=read.table("/Volumes/Temp/Lukas/Data/M_Schaeffer_PE_bams/Temperatures/tem_essex_ANDREWSFIELD_2006_2009.dat", header=T, skip=0,sep="\t",colClasses=c(NA,"character","character",NA))
essex_temp$Date=as.POSIXlt(paste(essex_temp$Date,essex_temp$HrMn),format="%Y%m%d %H%M")
with(essex_temp, tapply(QCP, months.Date(essex_temp$Date), min) )
with(essex_temp, tapply(QCP, months.Date(essex_temp$Date), mean) )
quartz()
boxplot(essex_temp$QCP~format(essex_temp$Date,"%m"),main="Essex",ylim=c(-15,40))
dev.copy2pdf(file="/Volumes/Temp/Lukas/Data/M_Schaeffer_PE_bams/Temperatures/temp_box_essex.pdf")

linz_temp=read.table("/Volumes/Temp/Lukas/Data/M_Schaeffer_PE_bams/Temperatures/temp_linz_stadt_2006_2009.dat", header=T, skip=0,sep="\t",colClasses=c(NA,"character","character",NA))
linz_temp$Date=as.POSIXlt(paste(linz_temp$Date,linz_temp$HrMn),format="%Y%m%d %H%M")
with(linz_temp, tapply(QCP, months.Date(linz_temp$Date), min) )
with(linz_temp, tapply(QCP, months.Date(linz_temp$Date), mean) )
quartz()
boxplot(linz_temp$QCP~format(linz_temp$Date,"%m"),main="Linz",ylim=c(-15,40))
dev.copy2pdf(file="/Volumes/Temp/Lukas/Data/M_Schaeffer_PE_bams/Temperatures/temp_box_linz.pdf")

dijon_temp=read.table("/Volumes/Temp/Lukas/Data/M_Schaeffer_PE_bams/Temperatures/temp_dijon_2006_2009.dat", header=T, skip=0,sep="\t",colClasses=c(NA,"character","character",NA))
dijon_temp$Date=as.POSIXlt(paste(dijon_temp$Date,dijon_temp$HrMn),format="%Y%m%d %H%M")
with(dijon_temp, tapply(QCP, months.Date(dijon_temp$Date), min) )
with(dijon_temp, tapply(QCP, months.Date(dijon_temp$Date), mean) )
quartz()
boxplot(dijon_temp$QCP~format(dijon_temp$Date,"%m"),main="Dijon",ylim=c(-15,40))
dev.copy2pdf(file="/Volumes/Temp/Lukas/Data/M_Schaeffer_PE_bams/Temperatures/temp_box_dijon.pdf")

barc_temp=read.table("/Volumes/Temp/Lukas/Data/M_Schaeffer_PE_bams/Temperatures/temp_barcelona_aerop_2007_2010.dat", header=T, skip=0,sep="\t",colClasses=c(NA,"character","character",NA))
barc_temp$Date=as.POSIXlt(paste(barc_temp$Date,barc_temp$HrMn),format="%Y%m%d %H%M")
with(barc_temp, tapply(QCP, months.Date(barc_temp$Date), min) )
with(barc_temp, tapply(QCP, months.Date(barc_temp$Date), mean) )
quartz()
boxplot(barc_temp$QCP~format(barc_temp$Date,"%m"),main="Barcelona",ylim=c(-15,40))
dev.copy2pdf(file="/Volumes/Temp/Lukas/Data/M_Schaeffer_PE_bams/Temperatures/temp_box_barc.pdf")

# sample model fit
#load fdr files
setwd("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/Joined_Analysis/")
maf_files=c("CHR","BPS","Alleles","MAF")
colnames(over_maf) = maf_files
vie2010_fdr_af = read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25_fdr_0.05.cmhout.af",header=F)
ita2011_fdr_af = read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20_filt_mc8_chroms_only_mct8_mcv25_fdr_0.05.cmhout.af",header=F)
over_fdr_af = read.table("vie2010_ita2011_avgMAF_fdr0.05.af",header=F)
colnames(over_fdr_af) = maf_files
af_files=c("CHR","BPS","Alleles","CS1","CS2","CSI3","B1","B2","B3","P")
colnames(vie2010_fdr_af) = af_files
colnames(ita2011_fdr_af) = af_files
vie2010_fdr_af$MAF=apply(vie2010_fdr_af[,7:9],1,mean)
#vie2010_fdr_af$MAF[vie2010_fdr_af$MAF < 0.5] = 1 - vie2010_fdr_af$MAF[vie2010_fdr_af$MAF < 0.5]
ita2011_fdr_af$MAF=apply(ita2011_fdr_af[,7:9],1,mean)
#ita2011_fdr_af$MAF[ita2011_fdr_af$MAF < 0.5] = 1 - ita2011_fdr_af$MAF[ita2011_fdr_af$MAF < 0.5] 
keeps=c("CHR","BPS","Alleles","MAF")
vie2010_fdr_af=vie2010_fdr_af[,keeps]
ita2011_fdr_af=ita2011_fdr_af[,keeps]

# creating histograms for fdr files, for each chrom indiviually
create_sampling_hist = function(snp_tab,fdr_tab,breaks=seq(0,1.0,0.02)){
	snp_afbins_chrom = list()
	fdr_hist_chrom = list()
	chroms=levels(fdr_tab$CHR)
	for(j in chroms){
		fdr_hist_chrom[[j]]=hist(fdr_tab$MAF[fdr_tab$CHR == j],breaks=breaks,plot=F)$counts
		# initialise lists for chromosomes
		snp_afbins_chrom[[j]]=list()
		for(i in 1:(length(breaks)-1) ) {
			snp_afbins_chrom[[j]][i]=list(which(snp_tab$MAF < breaks[i+1] & snp_tab$MAF >= breaks[i]& snp_tab$CHR == j))
		}	
	}
	ret_list=list("fdr"=fdr_hist_chrom,"all"=snp_afbins_chrom)
	return(ret_list)
}
create_fdr_sampling_hist = function(fdr_tab,breaks=seq(0,1.0,0.02)){
	fdr_hist_chrom = list()
	chroms=levels(fdr_tab$CHR)
	for(j in chroms){
            fdr_hist_chrom[[j]]=hist(fdr_tab$MAF[fdr_tab$CHR == j],breaks=breaks,plot=F)$counts
	}
	return(fdr_hist_chrom)
}

chroms=c("2L","2R","3L","3R","X")


# load snp files for MS
vie_ita_ms_names= unlist(strsplit("CHR_BPS_Alleles_Cf_Cn_MAF_Bn_KEN_nKEN_ETH_nETH_ATH_nATH_BAR_nBAR_DIJ_nDIJ_ESX_nESX_BOW_nBOW_LNZ_nLNZ","_"))
keeps=unlist(strsplit("CHR_BPS_Alleles_MAF_ATH_nATH_BAR_nBAR_DIJ_nDIJ_ESX_nESX_LNZ_nLNZ","_"))
vie2010_snps_ms = read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/vie_2010_avg_cs_b_KEN_ETH_ATH_BAR_DIJ_ESX_BOW_LNZ.sync.af",header=F,na.strings = c("NA","N","NaN"))
colnames(vie2010_snps_ms) = vie_ita_ms_names
vie2010_snps_ms=vie2010_snps_ms[,keeps]
ita2011_snps_ms = read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita_2011_avg_cs_b_KEN_ETH_ATH_BAR_DIJ_ESX_BOW_LNZ.sync.af",header=F,na.strings = c("NA","N","NaN"))
colnames(ita2011_snps_ms) = vie_ita_ms_names
ita2011_snps_ms=ita2011_snps_ms[,keeps]

vie_hist_list=create_sampling_hist(vie2010_snps_ms,vie2010_fdr_af)
#vie_hist_list[["fdr"]]=create_fdr_sampling_hist(vie2010_fdr_af)
vie_hist_list[["fdr"]]=vie_hist_list[["fdr"]][chroms]
ita_hist_list=create_sampling_hist(ita2011_snps_ms,ita2011_fdr_af)
#ita_hist_list[["fdr"]]=create_fdr_sampling_hist(ita2011_fdr_af)
ita_hist_list[["fdr"]]=ita_hist_list[["fdr"]][chroms]

# load snp files for lat
vie_ita_lat_names= unlist(strsplit("CHR_BPS_Alleles_Cf_Cn_MAF_Bn_Flo_nFlo_Penn_nPenn_Maine_nMaine","_"))
keeps=unlist(strsplit("CHR_BPS_Alleles_MAF_Flo_nFlo_Penn_nPenn_Maine_nMaine","_"))
vie2010_snps_lat = read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/vie_2010_cs_b_FLOR_PENN_MAINE.sync.af",header=F,na.strings = c("NA","N","NaN"))
colnames(vie2010_snps_lat) = vie_ita_lat_names
vie2010_snps_lat=vie2010_snps_lat[,keeps]
ita2011_snps_lat = read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/ita_2011_cs_b_FLOR_PENN_MAINE.sync.af",header=F,na.strings = c("NA","N","NaN"))
colnames(ita2011_snps_lat) = vie_ita_lat_names
ita2011_snps_lat=ita2011_snps_lat[,keeps]

vie_hist_list_lat=create_sampling_hist(vie2010_snps_lat,vie2010_fdr_af)
#vie_hist_list_lat[["fdr"]]=create_fdr_sampling_hist(vie2010_fdr_af)
vie_hist_list_lat[["fdr"]]=vie_hist_list_lat[["fdr"]][chroms]
ita_hist_list_lat=create_sampling_hist(ita2011_snps_lat,ita2011_fdr_af)
#ita_hist_list_lat[["fdr"]]=create_fdr_sampling_hist(ita2011_fdr_af)
ita_hist_list_lat[["fdr"]]=ita_hist_list_lat[["fdr"]][chroms]

# sample
get_samples_chrom= function(chrom_counts,afbins){
	# takes a vector with counts and a list with indeces in af bins and returns a random sample of indeces with the same size and af spectrum
	vec0=c()
	for(j in names(chrom_counts)){
		counts=chrom_counts[[j]]
	vec=foreach(i=1:length(counts), .combine='c') %:% when(counts[i] > 0) %do% sample(afbins[[j]][[i]],counts[i])
	vec0=c(vec0,vec)
	}
	return(vec0)	
}

iters=10
# load the foreach library and set it up for parallelisation
library(parallel)
library(doParallel)
registerDoParallel(cores=4)
library(foreach) 
save.image("rimage_260913_2pm.rData")
get_ret_vals = function(ms_poly,lats,cline,lm_model=T){
    ms_poly = ms_poly[! rowSums(ms_poly[,cline] > 0.95) == length(cline),]
    names(lats)=cline
    data_ms=fill_lat_frame(ms_poly,cline,lats)
    data_ms$Freq =with(data_ms, (MIN/(MAJ+MIN)))
    if (lm_model){
        lm_ms=lm(sqrt(Freq) ~ 1 + LAT, data=data_ms)
        sum_lm=summary(lm_ms)
        a=sum_lm$"fstatistic"
        ret_val=c(sum_lm$"coefficients"["LAT",],sum_lm$"adj.r.squared",pf(a[1],a[2],a[3],lower.tail=F),extractAIC(lm_ms)[2])
    }
    else {
        glm_ms <-glm(cbind(MIN,MAJ) ~ 1 + LAT, data=data_ms, family=binomial)
        sum_glm=summary(glm_ms)
        ret_val=c(sum_glm$"coefficients"["LAT",],sum_glm$"aic")
    }
    return(ret_val)
}

library(combinat)
lat_perms=permn(latitudes[ms_cline])
iters=10000
# parallel foreach
ms_vie_samples = foreach(i=1:iters, .combine=rbind) %dopar% {
    ms_poly=vie2010_snps_ms[get_samples_chrom(vie_hist_list[["fdr"]],vie_hist_list[["all"]]),]
    get_ret_vals(ms_poly,latitudes,ms_cline,lm_model=T)
}
ms_vie_glm_samples = foreach(i=1:iters, .combine=rbind) %dopar% {
    ms_poly=vie2010_snps_ms[get_samples_chrom(vie_hist_list[["fdr"]],vie_hist_list[["all"]]),]
    get_ret_vals(ms_poly,latitudes,ms_cline,lm_model=F)
}
quartz()
hist(ms_vie_samples[,"Estimate"],breaks=100,xlab=expression("slope [(fract)/degree latitude]"),xlim=c(min(ms_vie_samples[,"Estimate"]),lm_ms_vie$coef[2]))
abline(v=median(ms_vie_samples[,"Estimate"]),col="blue",lwd=1,lty=2)
abline(v=lm_ms_vie$coef[2],col="red",lwd=2,lty=2)

quartz()
hist(ms_vie_glm_samples[,"Estimate"],breaks=100,xlab=expression("slope [(fract)/degree latitude]"),xlim=c(min(ms_vie_glm_samples[,"Estimate"]),lm_ms_vie$coef[2]))
abline(v=median(ms_vie_glm_samples[,"Estimate"]),col="blue",lwd=1,lty=2)
abline(v=glm_ms_vie$coef[2],col="red",lwd=2,lty=2)

quartz()
hist(ms_vie_samples[,5])
# parallel foreach
ms_ita_samples = foreach(i=1:iters, .combine=rbind) %dopar% {
    ms_poly=ita2011_snps_ms[get_samples_chrom(ita_hist_list[["fdr"]],ita_hist_list[["all"]]),]
    get_ret_vals(ms_poly,latitudes,ms_cline,lm_model=T)
}
quartz()
hist(ms_ita_samples[,"Estimate"])
quartz()
hist(ms_vie_samples[,5])

ms_vie_samp_lat = foreach(i=1:length(lat_perms), .combine=rbind) %do% {
    ms_poly = ms_vie[rowMeans(ms_vie[,b_cols] - ms_vie[,cs_cols]) > 0,]
    lat_samp=lat_perms[[i]]
    get_ret_vals(ms_poly,lat_samp,ms_cline,lm_model=T)
    }
quartz()
hist(ms_vie_samp_lat[-1:-3,"Estimate"],breaks=10)
max(ms_vie_samp_lat[,"Estimate"])
max(ms_vie_samp_lat[,5])
sum

lat_vie_samples = foreach(i=1:iters, .combine=rbind) %dopar% {
    ms_poly=vie2010_snps_lat[get_samples_chrom(vie_hist_list_lat[["fdr"]],vie_hist_list_lat[["all"]]),]
    get_ret_vals(ms_poly,latitudes,ms_cline,lm_model=T)
}
quartz()
hist(lat_vie_samples[,"Estimate"])
quartz()
hist(ms_vie_samples[,5])

# parallel foreach
ms_ita_samples = foreach(i=1:iters, .combine=rbind) %dopar% {
    ms_poly=ita2011_snps_ms[get_samples_chrom(ita_hist_list_lat[["fdr"]],ita_hist_list_lat[["all"]]),]
    get_ret_vals(ms_poly,latitudes,ms_cline,lm_model=T)
}
quartz()
hist(ms_ita_samples[,"Estimate"])


quartz()
boxplot(sqrt(1-ms_vie_poly[,ms_cline]),xaxt="n",at=latitudes[ms_cline],notch=F,ylab=expression(sqrt( "AF" )), main="European cline, Vienna",xlab="latitude")
axis(1)
text(x =  latitudes[ms_cline], y = par("usr")[3] - 0.05, srt = 45, adj = 1,labels = ms_cline, xpd = TRUE)
abline(lm_vie,lty=2,col="red",lwd=2)
dev.copy2pdf(file="ms_vie_fdr_fit_lm.pdf")
quartz()
boxplot(sqrt(1-ms_ita_poly[,ms_cline]),xaxt="n",at=latitudes[ms_cline],notch=F,ylab=expression(sqrt( "AF" )), main="European cline, Italy",xlab="latitude")
axis(1)
text(x =  latitudes[ms_cline], y = par("usr")[3] - 0.05, srt = 45, adj = 1,labels = ms_cline, xpd = TRUE)
abline(lm_ita,lty=2,col="red",lwd=2)
dev.copy2pdf(file="ms_ita_fdr_fit_lm.pdf")

quartz()
boxplot(sqrt(1-lat_vie_poly[,lat_cline]),xaxt="n",at=latitudes[lat_cline],notch=F,ylab=expression(sqrt( "AF" )), main="North American cline, Vienna",xlab="latitude")
axis(1)
text(x =  latitudes[lat_cline], y = par("usr")[3] - 0.08, srt = 0, adj = 1,labels = lat_cline, xpd = TRUE)
abline(lm_lat_vie,lty=2,col="red",lwd=2)
dev.copy2pdf(file="lat_vie_fdr_fit_lm.pdf")
quartz()
boxplot(sqrt(1-lat_ita_poly[,lat_cline]),xaxt="n",at=latitudes[lat_cline],notch=F,ylab=expression(sqrt( "AF" )), main="North American cline, Italy",xlab="latitude")
axis(1)
text(x =  latitudes[lat_cline], y = par("usr")[3] - 0.08, srt = 0, adj = 1,labels = lat_cline, xpd = TRUE)
abline(lm_lat_ita,lty=2,col="red",lwd=2)
dev.copy2pdf(file="lat_ita_fdr_fit_lm.pdf")


#for ns data:
ms_ita_nns=merge(ms_ita_ns[,c("CHR","BPS","C1")],ms_ita,all.y=T,by=c("CHR","BPS"))
ms_ita_nns=ms_ita_nns[is.na(ms_ita_nns$C1.x),]
lat_ita_nns=merge(lat_ita_ns[,c("CHR","BPS","C1")],lat_ita,all.y=T,by=c("CHR","BPS"))
lat_ita_nns=lat_ita_nns[is.na(lat_ita_nns$C1.x),]
ms_vie_nns=merge(ms_vie_ns[,c("CHR","BPS","C1")],ms_vie,all.y=T,by=c("CHR","BPS"))
ms_vie_nns=ms_vie_nns[is.na(ms_vie_nns$C1.x),]
lat_vie_nns=merge(lat_vie_ns[,c("CHR","BPS","C1")],lat_vie,all.y=T,by=c("CHR","BPS"))
lat_vie_nns=lat_vie_nns[is.na(lat_vie_nns$C1.x),]


data_ms_ita_ns=fill_lat_frame(ms_ita_ns,ms_cline,latitudes)
data_lat_ita_ns=fill_lat_frame(lat_ita_ns,lat_cline,latitudes)
data_ms_vie_ns=fill_lat_frame(ms_vie_ns,ms_cline,latitudes)
data_lat_vie_ns=fill_lat_frame(lat_vie_ns,lat_cline,latitudes)

lm_lat_ita_ns=lm(sqrt(Freq) ~ 1 + LAT, data=data_lat_ita_ns)
lm_lat_vie_ns=lm(sqrt(Freq) ~ 1 + LAT, data=data_lat_vie_ns)
lm_ms_ita_ns=lm(sqrt(Freq) ~ 1 + LAT, data=data_ms_ita_ns)
lm_ms_vie_ns=lm(sqrt(Freq) ~ 1 + LAT, data=data_ms_vie_ns)
glm_lat_ita_ns=glm(cbind(MIN,MAJ) ~ 1 + LAT, data=data_lat_ita_ns, family=binomial)
glm_lat_vie_ns=glm(cbind(MIN,MAJ) ~ 1 + LAT, data=data_lat_vie_ns, family=binomial)
glm_ms_ita_ns=glm(cbind(MIN,MAJ) ~ 1 + LAT, data=data_ms_ita_ns, family=binomial)
glm_ms_vie_ns=glm(cbind(MIN,MAJ) ~ 1 + LAT, data=data_ms_vie_ns, family=binomial)

data_ms_ita_nns=fill_lat_frame(ms_ita_nns,ms_cline,latitudes)
data_lat_ita_nns=fill_lat_frame(lat_ita_nns,lat_cline,latitudes)
data_ms_vie_nns=fill_lat_frame(ms_vie_nns,ms_cline,latitudes)
data_lat_vie_nns=fill_lat_frame(lat_vie_nns,lat_cline,latitudes)
lm_lat_ita_nns=lm(sqrt(Freq) ~ 1 + LAT, data=data_lat_ita_nns)
lm_lat_vie_nns=lm(sqrt(Freq) ~ 1 + LAT, data=data_lat_vie_nns)
lm_ms_ita_nns=lm(sqrt(Freq) ~ 1 + LAT, data=data_ms_ita_nns)
lm_ms_vie_nns=lm(sqrt(Freq) ~ 1 + LAT, data=data_ms_vie_nns)
glm_lat_ita_nns=glm(cbind(MIN,MAJ) ~ 1 + LAT, data=data_lat_ita_nns, family=binomial)
glm_lat_vie_nns=glm(cbind(MIN,MAJ) ~ 1 + LAT, data=data_lat_vie_nns, family=binomial)
glm_ms_ita_nns=glm(cbind(MIN,MAJ) ~ 1 + LAT, data=data_ms_ita_nns, family=binomial)
glm_ms_vie_nns=glm(cbind(MIN,MAJ) ~ 1 + LAT, data=data_ms_vie_nns, family=binomial)

summary(lm_lat_ita_ns)
summary(lm_lat_ita_nns)
summary(lm_lat_vie_ns)
summary(lm_lat_vie_nns)
summary(lm_ms_ita_ns)
summary(lm_ms_ita_nns)
summary(lm_ms_vie_ns)
summary(lm_ms_vie_nns)
summary(glm_lat_ita_ns)
summary(glm_lat_ita_nns)
summary(glm_lat_vie_ns)
summary(glm_lat_vie_nns)
summary(glm_ms_ita_ns)
summary(glm_ms_ita_nns)
summary(glm_ms_vie_ns)
summary(glm_ms_vie_nns)

ms_ita_nns=merge(ms_ita_ns[,c("CHR","BPS","P")],ms_ita,all.y=T,by=c("CHR","BPS"))
ms_ita_nns=ms_ita_nns[is.na(ms_ita_not_ns$P.x),]
    
quartz()
boxplot(sqrt(1-ms_vie_ns[,ms_cline]),xaxt="n",at=latitudes[ms_cline],notch=F,ylab=expression(sqrt( "AF" )), main="North American cline, Italy, non. syn.",xlab="latitude")
axis(1)
text(x =  latitudes[ms_cline], y = par("usr")[3] - 0.08, srt = 0, adj = 1,labels = ms_cline, xpd = TRUE)
abline(lm_ms_vie_ns,lty=2,col="red",lwd=2)
abline(lm_ms_vie,lty=2,col="blue",lwd=2)
