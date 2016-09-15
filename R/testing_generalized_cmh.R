# get data from best SNP, two alleles, four levels of pigmentation (very light, light, dark, very dark)
tab = as.table(array(c(7,55,9,47,14,103,
11,53,23,46,29,50,
21,19,65,38,65,31,
59,23,64,13,76,30), dim=c(2,3,4),dimnames=list(Alleles=c("A1","A2"), Reps=c("rep1","rep2","rep3"),Pigmentation=c("VL","L","D","VD"))))
# get it into a nice format
tab_reps=as.table(array(c(tab[,1,], tab[,2,],tab[,3,] ),dim=c(2,4,3),dimnames=list(Alleles=c("A1","A2"), Pigmentation=c("VL","L","D","VD"),Reps=c("rep1","rep2","rep3")) ))
# same for a less significant snp
tab2 = as.table(array(c(67,5,54,11,132,7,55,6,51,18,72,15,45,12,87,30,70,33,60,21,53:41,87,21), dim=c(2,3,4),dimnames=list(Alleles=c("A1","A2"), Reps=c("rep1","rep2","rep3"),Pigmentation=c("VL","L","D","VD"))))
tab2_reps=as.table(array(c(tab2[,1,], tab2[,2,],tab2[,3,] ),dim=c(2,4,3),dimnames=list(Alleles=c("A1","A2"), Pigmentation=c("VL","L","D","VD"),Reps=c("rep1","rep2","rep3")) ))
# look at stuff
ftab=ftable(. ~ Alleles + Pigmentation + Reps,tab_reps)
mantelhaen.test(tab_reps[,c(1,2,3,4),])$p.value
# vl againts vd (should be p ~ 3.879513e-43)
mantelhaen.test(tab_reps[,c(1,4),])$p.value
# shuffling the different pigmentation levels
mantelhaen.test(tab_reps[,c(3,4,1,2),])$p.value
library("vcdExtra")
# gen CMH for extremes - not smae result as no continuum correction - gives same asn mantel-haens with cor=FALSE
CMHtest(tab_reps[,c(1,4),],types=c("ALL"),overall=TRUE)
# gen with correct ranked columns
CMHtest(tab_reps[c(1,2),c(1,2,3,4),],types=c("ALL"),overall=TRUE)
# shuffled columns
CMHtest(tab_reps[,c(4,2,3,1),],types=c("ALL"),overall=TRUE)
CMHtest(tab_reps[,c(4,1,3,2),],types=c("ALL"),overall=TRUE)

CMHtest(tab_reps[,c(1,4),],types=c("ALL"),overall=TRUE)

CMHtest(tab2_reps[c(2,1),c(1,2,3,4),],types=c("ALL"),overall=TRUE)
CMHtest(tab2_reps[,c(3,4,1,2),],types=c("ALL"),overall=TRUE)
CMHtest(tab_reps[c(2,1),c(2,4,1,3),c(1,2)],types=c("ALL"),cscores="midrank",overall=TRUE)


tab_reps_test=tab2_reps
# try glm
r1=as.data.frame(tab_reps_test[,,"rep1"])
r2=as.data.frame(tab_reps_test[,,"rep2"])
r3=as.data.frame(tab_reps_test[,,"rep3"])
r1$REP=1
#r1$ALL=c(1,2)
r2$REP=2
#r2$ALL=c(1,2)
r3$REP=3
#r3$ALL=c(1,2)
r_tot2=df_ret=data.frame(A1=as.integer(),A2=as.integer(),Rep=as.integer(),Lev=as.integer())
for(i in 1:4){
	r_tot2=rbind(r_tot2, c(r1[1,i],r1[2,i],r1$REP[1],i))
	r_tot2=rbind(r_tot2, c(r2[1,i],r2[2,i],r2$REP[1],i))
	r_tot2=rbind(r_tot2, c(r3[1,i],r3[2,i],r3$REP[1],i))		
}
colnames(r_tot2)=c("A1","A2","Rep","Lev")
r_tot2$Rep=factor(r_tot2$Rep)
r_tot2$Lev=factor(r_tot2$Lev)
r_tot2$Freq=r_tot2$A1/(r_tot2$A1+r_tot2$A2)
r_tot2$REP=factor(10*as.integer(r_tot2$Rep)+as.integer(r_tot2$Lev))
library("lme4")
glm_test<-glm(cbind(A1,A2) ~ 1 + Lev, data=r_tot, family=binomial)
summary(glm_test)
summary(lmer_test<-lmer(cbind(A1,A2) ~ 1 + Lev + (1|Rep), data=r_tot, family=binomial))
summary(lmer_lin<-lmer(sqrt(Freq) ~ 1 + Lev + (1|Rep), data=r_tot))
summary(lmer_test<-lmer(cbind(A1,A2) ~ 1 + (1|Lev) , data=r_tot, family=binomial))
summary(lmer_null<-lmer(cbind(A1,A2) ~ 1 + (1|Rep) , data=r_tot, family=binomial))
plot(r_tot$Lev,sqrt(r_tot$Freq))
summary(lm_test<-lmer(sqrt(Freq) ~  1 + Lev , data=r_tot))
abline(lm_test,col="red")
par(mfrow=c(2,2))
plot(lm_test)
par(mfrow=c(1,1))
plot(r_tot2$Lev,sqrt(r_tot2$Freq))
summary(lmer_lin<-lmer(sqrt(Freq) ~ 1 + Lev + (-1|Rep), data=r_tot2))
