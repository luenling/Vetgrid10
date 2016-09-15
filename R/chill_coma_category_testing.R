setwd("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/Joined_Analysis/")
chill_snps=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/Joined_Analysis/chill_coma_snps_comparisons.txt", header=T,sep = "\t")
fisher_ns = function(t_a,t_b,a,b){
return(fisher.test(matrix(c(t_a-a,a,t_b-b,b),nrow=2)))
}
chi2_ns = function(t_a,t_b,a,b){
return(chisq.test(matrix(c(t_a-a,a,t_b-b,b),nrow=2)))
}

test_cat=c( "Non.Syn_Coding","Syn_Coding","Intron","IN.2L.t","IN.2R.NS","IN.3L.P","IN.3R.C","IN.3R.P","IN.3R.Mo")
comp_cat=c("Total","Total","Total","X2L","X2R","X3L","X3R","X3R","X3R")
pvals=data.frame(rbind(1:length(test_cat)))
names(pvals)=test_cat
odds=data.frame(rbind(1:length(test_cat)))
names(odds)=test_cat
row_num=seq(1,5,2)
for(i in 1:length(row_num)){
	for (j in 1:length(test_cat)){ 
		x=row_num[i]
		y_test=test_cat[j]
		y_tot=comp_cat[j]
		pvals[i,j]=fisher_ns(chill_snps[x,y_tot],chill_snps[x+1,y_tot],chill_snps[x,y_test],chill_snps[x+1,y_test])$p.value
		odds[i,j]=fisher_ns(chill_snps[x,y_tot],chill_snps[x+1,y_tot],chill_snps[x,y_test],chill_snps[x+1,y_test])$estimate
	}
}
row.names(odds)=c("Vie10","Ita11","Overlap")
row.names(pvals)=c("Vie10","Ita11","Overlap")


> signif(pvals,3)
        Non.Syn_Coding Syn_Coding Intron  IN.2L.t IN.2R.NS  IN.3L.P  IN.3R.C  IN.3R.P IN.3R.Mo
Vie10        2.90e-105   2.71e-31  0.267 1.16e-06    0.117 1.00e-08 2.03e-09 6.43e-08 5.46e-05
Ita11         4.09e-20   2.02e-14  0.851 1.14e-02    0.085 1.06e-02 8.07e-03 8.16e-02 2.14e-01
Overlap       1.48e-05   1.72e-05  0.631 6.19e-03    0.251 5.67e-02 1.57e-03 7.94e-03 1.53e-01

> signif(odds,3)
        Non.Syn_Coding Syn_Coding Intron IN.2L.t IN.2R.NS IN.3L.P IN.3R.C IN.3R.P IN.3R.Mo
Vie10             4.07      0.385  0.956   0.631    0.840   0.587   0.617   0.620    0.697
Ita11             2.79      0.350  0.987   0.684    0.728   0.677   0.707   0.780    0.830
Overlap           2.87      0.270  0.933   0.433    0.620   0.516   0.367   0.368    0.578

for(i in seq(1,5,2)){
	for (j in test_cat){ 
		pvals[i,j]=chi2_ns(chill_snps$Total[i],chill_snps$Total[i+1],chill_snps[i,j],chill_snps[i+1,j])$p.value
		odds[i,j]=fisher_ns(chill_snps$Total[i],chill_snps$Total[i+1],chill_snps[i,j],chill_snps[i+1,j])$
estimate
	}
}

setwd("/Volumes/Temp/Lukas/Recombination_rates/")
bed_invs=read.table("Dmel_recomb_rates.invs.bed",sep="\t",header=F)
bed_noinvs=read.table("Dmel_recomb_rates.noinvs.bed",sep="\t",header=F)
quartz()
boxplot(V4~V5,data=bed_invs,main="Recombination rates")
dev.copy2pdf(file="rec_rates_inversions_boxplot.pdf")
boxplot(V4~V1,data=bed_noinvs,main="Recombination rates")
dev.copy2pdf(file="rec_rates_outside_inversions_boxplot.pdf")

inv_rates=aggregate(V4~V5,bed_invs,mean,na.rm=TRUE)
inv_rates=merge(inv_rates,aggregate(V4~V5,bed_invs,sd),all=T,by=c("V5"))
invs_rec=aggregate(V4~V5,bed_invs,summary)
chroms=c("2L","2R","3L","3R","3R","3R")
for(i in 1:6){
	a = bed_noinvs[bed_noinvs$V1 == chroms[i] & bed_noinvs$V4 <= invs_rec[i,2][5] &  bed_noinvs$V4 >= invs_rec[i,2][2], ]
	assign(gsub("\\((.*)\\)", "_\\1_", invs_rec[i,1]),a)

}
file_nam = ls(pattern="IN.*")
for(i in file_nam){
	write.table(eval(as.name(i)),file=paste("complement_",i,".bed",sep=""),quote=F,col.names=F,row.names=F,sep="\t")
	}


