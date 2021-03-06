setwd("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/CMH")
pop_fsts=read.table("fem_vie2010__pI25_pII25_pIII25_pIel_pIIel_pIIIel_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_masked.fst.tab.bases", header=T)
aggregate(X4.10 ~ CHR, pop_fsts,mean)
aggregate(X5.11 ~ CHR, pop_fsts,mean)
aggregate(X6.12 ~ CHR, pop_fsts,mean)
m4_10=mean(pop_fsts$X4.10)
m5_11=mean(pop_fsts$X5.11)
m6_12=mean(pop_fsts$X6.12)
fst_average=mean(m4_10,m5_11,m6_12)
v4_10=var(pop_fsts$X4.10)
v5_11=var(pop_fsts$X5.11)
v6_12=var(pop_fsts$X6.12)
v_tot=var(c(pop_fsts$X4.10,pop_fsts$X5.11,pop_fsts$X6.12))
#fsttot = 0.01448959+/-0.01653625
# creating histograms with colored quantiles
values=rnorm(15000,100,15)
quant = 100/1500
h=hist(values,breaks=25,plot=F)
clr = rep(rgb(0,0,1,1/4), length(h$counts))
clr[which(h$breaks <= quantile(values,quant))] = rgb(1,0,0,1/4)
clr[which(h$breaks >= quantile(values,1-quant))] = rgb(1,0,0,1/4)
plot(h, col=clr,main="",xlab="",ylab="",axes = F)
axis(1, label=F, tick=T)
dev.copy2pdf(filename="histogram_15.pdf")
