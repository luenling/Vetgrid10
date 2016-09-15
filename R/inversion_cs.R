library(doBy)

boxplot_inv <- function(a,n,tit="",col1="lightblue",col2="orange",name1="CS",name2="B")
{
	pops=(length(a[[n]])-5)/2
	names=1:(length(a[[n]])-5)
	for(i in 1:pops){names[i]=paste(name1,i,sep=""); names[i+pops]=paste(name2,i,sep="");}
	boxplot(a[[n]][6:length(a[[n]])],data=a[[n]],notch=F,col=c(rep(col1,pops),rep(col2,pops)),ylab="frequency of allele fixed in inversion",names=names, main=paste("Inversion",a[[n]]$CHR[1],a[[n]]$INV[1],tit))
	}

boxplot_inv_combine <- function(a,tit="",col1="lightblue",col2="orange",name1="CS",name2="B",pops1=1:3,pops2=4:6)
{	
	b1=list()
	b2=list()
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


# cs inversion plotting
setwd("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/Joined_Analysis")
#setwd("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH")
inv_data_2011 <- read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.inversions",header=TRUE)
inv_count_2011 <- read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.counts.inversions",header=TRUE)
summaryBy(pop1 + pop2 + pop3 + pop4 +pop5 + pop6 ~ CHR + INV, data = inv_data_2011, FUN =c(median), keep.names = T)
#aggregate(cbind(pop1,pop2,pop3,pop4)~CHR+INV,data=inv_data,FUN=median)
a_11 <- splitBy(formula=~ CHR + INV, data = inv_data_2011 )
quartz()
boxplot_inv(a_11,1,tit=", Bolzano 2011")
dev.copy2pdf(file="inv_cs_ita2011_fem.pdf")
#setwd("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_Poly/CMH")
inv_data_2010 <- read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_Poly/CMH/females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.inversions",header=TRUE)
inv_count_2010 <- read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_Poly/CMH/females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.counts.inversions",header=TRUE)
c_10 <- splitBy(formula=~ CHR + INV, data = inv_count_2010 )
c_11 <- splitBy(formula=~ CHR + INV, data = inv_count_2011 )


summaryBy(pop1 + pop2 + pop3 + pop4 +pop5 + pop6 ~ CHR + INV, data = inv_data_2010, FUN =c(median), keep.names = T)
#aggregate(cbind(pop1,pop2,pop3,pop4)~CHR+INV,data=inv_data,FUN=median)
a_10 <- splitBy(formula=~ CHR + INV, data = inv_data_2010 )
c_10 <- splitBy(formula=~ CHR + INV, data = inv_count_2010 )
for(i in names(c_10)){c_10[[i]] = sapply(c_10[[i]][,-(1:5)],median)}
for(i in c_10){ a = c(i[seq(from=1,to=length(i),by=2)]*i[seq(from=2,to=length(i),by=2)],(1-i[seq(from=1,to=length(i),by=2)])*i[seq(from=2,to=length(i),by=2)] ) }
for(i in c_10){ print(mantelhaen.test( array(matrix( c(i[seq(from=1,to=length(i),by=2)]*i[seq(from=2,to=length(i),by=2)],(1-i[seq(from=1,to=length(i),by=2)])*i[seq(from=2,to=length(i),by=2)]) ,nrow=4,byrow=T),dim=c(2,2,3)) )$p.value) }
c_10 <- splitBy(formula=~ CHR + INV, data = inv_count_2010 )
c_11 <- splitBy(formula=~ CHR + INV, data = inv_count_2011 )
for(i in names(c_11)){c_11[[i]] = sapply(c_11[[i]][,-(1:5)],median)}
for(i in c_11){ print(mantelhaen.test( array(matrix( c(i[seq(from=1,to=length(i),by=2)]*i[seq(from=2,to=length(i),by=2)],(1-i[seq(from=1,to=length(i),by=2)])*i[seq(from=2,to=length(i),by=2)]) ,nrow=4,byrow=T),dim=c(2,2,3)) )$p.value) }

quartz()
boxplot_inv(a_10,1,tit=", Vienna 2010")

quartz()
boxplot_inv(a_10,1,tit=", Vienna 2010")
dev.copy2pdf(file="inv_cs_Vie2010.pdf")
quartz()
boxplot_inv_combine(a_10,tit="Vienna 2010")
dev.copy2pdf(file="inv_cs_all_Vie2010.pdf")
quartz()
boxplot_inv_combine(a_11,tit="Bolzano 2011")
dev.copy2pdf(file="inv_cs_all_Ita2011.pdf")

inv_data_over = read.table("vie2010_2011.merged.inversions",header=T)
quartz()
boxplot_inv(a_10,6,tit=", Vienna 2010")
dev.copy2pdf(file="inv_3Rp_Vie2010.pdf")
for(i in 1:6){
	quartz()
boxplot_inv(a_11,i,tit=", Vienna 2010")
	}
# transforming list entries into matrices for marlies:
# count cov rep b_cs

fill_data = function(raw_list){
	data_list = list()
	for(inversion in names(raw_list)) {
		int_tab=raw_list[[inversion]]
		for(i in seq(from=6,to=length(int_tab[1,]),by=2)){int_tab[,i]=int_tab[,i]*int_tab[,i+1]}
		marl_mat=matrix(nrow=0,ncol=4)
	colnames(marl_mat)=c("count","cov","rep","cs")
	for(i in 1:3){
		j=(i-1)*2
		loc_mat=cbind(int_tab[,(6+j):(6+j+1)],i,1)
		colnames(loc_mat)=c("count","cov","rep","cs")
		marl_mat=rbind(marl_mat,loc_mat) 
		loc_mat=cbind(int_tab[,(12+j):(12+j+1)],i,0)
		colnames(loc_mat)=c("count","cov","rep","cs")
		marl_mat=rbind(marl_mat,loc_mat)
		}
	rownames(marl_mat)=NULL
	data=as.data.frame(marl_mat)	
	data$freq<-round((data$count/data$cov)*100,digits=0)
	data$inv_count<-round(data$count,digits=0)
	data$not_inv_count<- data$cov-data$inv_count
	data$experiment<-as.factor(data$rep*10+data$cs)
	data$rep_fac<-as.factor(data$rep)
	data$cs_fac<-as.factor(data$cs)
	data_list[[inversion]] = data		
	}		
	return(data_list)	
}	

data_11=fill_data(c_11)	
data_10=fill_data(c_10)	
str(data_10)
library("lme4")
lmer_result_11=list()
for (inversion in names(data_11)){
	lmer_result_11[[inversion]]<-lmer(freq ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data_11[[inversion]], family=poisson)
}
lmer_result_10=list()
for (inversion in names(data_10)){
	lmer_result_10[[inversion]]<-lmer(freq ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data_10[[inversion]], family=poisson)
}
summary(lmer_result_10[[1]])
lmer_binom_result_11=list()
for (inversion in names(data_11)){
	lmer_binom_result_11[[inversion]]<-lmer(cbind(inv_count, not_inv_count) ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data_11[[inversion]], family=binomial)
}
summary(lmer_binom_result_11[[1]])
lmer_binom_result_10=list()
for (inversion in names(data_10)){
	lmer_binom_result_10[[inversion]]<-lmer(cbind(inv_count, not_inv_count) ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data_10[[inversion]], family=binomial)
}
summary(lmer_binom_result_10[[1]])

# only significant result in Bolzano, no significant result for in(2L)t in Vienna
# try 10 2L_t without replicate 2
data_10_2Lt_wo_rep2=data_10[["2L|t"]][data_10[["2L|t"]]$rep != 2, ]
summary(lmer(freq ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data_10_2Lt_wo_rep2, family=poisson))
summary(lmer(cbind(inv_count, not_inv_count) ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data_10_2Lt_wo_rep2, family=binomial))
# install.packages("coefplot2",repos="http://www.math.mcmaster.ca/bolker/R", type="source")
library(coefplot2)
coefplot2(test_lmer)
coeftab(test_lmer)
library(influence.ME)
x11()
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
f <- fitted(fm1)
r <- residuals(fm1)
plot(f,r)
sm <- loess(r~f)
v <- seq(min(f),max(f),length=101)
lines(v,predict(sm,data.frame(f=v)),col=2)
test_lmer=lmer(freq ~ 1+ (1|rep_fac)+ (1|rep_fac/experiment) + cs_fac, data=data_10_2Lt_wo_rep2,family=poisson)
f <- fitted(test_lmer)
r <- residuals(test_lmer)
plot(f,r)
sm <- loess(r~f)
v <- seq(min(f),max(f),length=101)
lines(v,predict(sm,data.frame(f=v)),col=2)
library(ggplot2)
ggplot(data = sleepstudy, aes(x = Days, y = Reaction, color = factor(Subject))) + geom_line(aes(group = Subject)) + geom_point()
