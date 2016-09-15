setwd("/Volumes/vetgrid10/Data/A7_pigmentation/Simulation")
setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Simulation")
library(stats)
library(reshape2,ggplot2)
library("RColorBrewer")
library(abc)
### use uniformly distributed parameter values
load("simres.unif.eps.2mio.RData")
load("simres.unif.500K.special.eps.1.RData")
uniform.sim.eps <- simres.eps
uniform.sim.eps = as.data.frame(uniform.sim.eps)
load("simres.unif.500K.special.eps.2a.RData")
uniform.sim.eps <- rbind(uniform.sim.eps,simres.eps)
load("simres.unif.500K.special.eps.2b.RData")
uniform.sim.eps <- rbind(uniform.sim.eps,simres.eps)
load("simres.unif.500K.special.eps.2.RData")
uniform.sim.eps <- rbind(uniform.sim.eps,simres.eps)


## get 2 mio simulations
rm(simres.eps)
#load("simres.unif.noeps.RData.gz")
#uniform.sim.noeps <- simulation_results
uniform.sim.noeps <- readRDS("simres.unif.noeps.2mio.rds")
rm(simulation_results)

load("thur07_08_2015.RSession")

# uniform.sim.eps = as.data.frame(uniform.sim.eps[,-c(6:7,10:17,20:25)])
# uniform.sim.noeps = as.data.frame(uniform.sim.noeps[,-c(6:7,10:17,20:25)])
# uniform.sim.eps[,c("or.eu.tan","or.eu.bab","or.sa.tan","or.sa.bab")]=log(uniform.sim.eps[,c("or1Eu","or2Eu","or1Sa","or2Sa")])
# uniform.sim.noeps[,c("or.eu.tan","or.eu.bab","or.sa.tan","or.sa.bab")]=log(uniform.sim.noeps[,c("or1Eu","or2Eu","or1Sa","or2Sa")])
# eps int
uniform.sim.eps = as.data.frame(uniform.sim.eps)
uniform.sim.eps[,c("or.eu.tan","or.eu.bab","or.sa.tan","or.sa.bab")]=log(uniform.sim.eps[,c("or.eu.tan","or.eu.bab","or.sa.tan","or.sa.bab")])
# no eps
uniform.sim.noeps = as.data.frame(uniform.sim.noeps)
uniform.sim.noeps[,c("or.eu.tan","or.eu.bab","or.sa.tan","or.sa.bab")]=log(uniform.sim.noeps[,c("or.eu.tan","or.eu.bab","or.sa.tan","or.sa.bab")])



cln=colnames(uniform.sim.eps)
cln[3]="rel.eff"
colnames(uniform.sim.eps)=cln
colnames(uniform.sim.noeps)=cln




#saveRDS(uniform.sim.eps,file="uniform.sim.eps.rds",compress=TRUE)

criteria=c(2.51,0.57,1.89,0.44)
names(criteria)=c("eu","sa","tan_eusa","bab_eusa")
# calculate criteria

uniform.sim.eps = calculate_criteria(uniform.sim.eps)
uniform.sim.noeps = calculate_criteria(uniform.sim.noeps)

calculate_criteria <- function(unidf) {
  for (i in c("or.eu.tan","or.eu.bab","or.sa.tan","or.sa.bab")){
    minim= min(unidf[is.finite(unidf[,i]),i])-1
    maxim= max(unidf[is.finite(unidf[,i]),i])+1
    unidf[! is.na(unidf[,i]) & unidf[,i] == Inf, i ] = maxim
    unidf[! is.na(unidf[,i]) & unidf[,i] == -Inf, i ] = minim
    #unidf[unidf[,i] == 0, i ] = 1e-16
  }
  
  unidf$eu=unidf$or.eu.tan/unidf$or.eu.bab
  unidf$sa=unidf$or.sa.tan/unidf$or.sa.bab
  unidf$tan_eusa=unidf$or.eu.tan/unidf$or.sa.tan
  unidf$bab_eusa=unidf$or.eu.bab/unidf$or.sa.bab
  return(unidf)
}

# compare best eps to no eps
load("simres.unif.eps.2.25.500K.RData")
uniform.sim.eps.max <- simres.eps
rm(simres.eps)
uniform.sim.eps.max = as.data.frame(uniform.sim.eps.max)
uniform.sim.eps.max[,c("or.eu.tan","or.eu.bab","or.sa.tan","or.sa.bab")]=log(uniform.sim.eps.max[,c("or.eu.tan","or.eu.bab","or.sa.tan","or.sa.bab")])
cln=colnames(uniform.sim.eps.max)
cln[3]="rel.eff"
colnames(uniform.sim.eps.max)=cln
uniform.sim.eps.max = calculate_criteria(uniform.sim.eps.max)

best.sim=rbind(uniform.sim.noeps[1:nrow(uniform.sim.eps.max),c(names(criteria))],uniform.sim.eps.max[,c(names(criteria))])
best.sim$model=factor(c(rep("no.eps",nrow(uniform.sim.eps.max)),rep("eps",nrow(uniform.sim.eps.max))))
modelsel.max=postpr(target = criteria,index=as.character(best.sim$model),sumstat=best.sim[,c(names(criteria))],tol=0.1,method="rejection")
summary(modelsel.max)
# Proportion of accepted simulations (rejection):
#   eps no.eps 
# 0.9917 0.0083 
# 
# Bayes factors:
#   eps   no.eps
# eps      1.0000 119.1634
# no.eps   0.0084   1.0000

# do abc for eps
tol = 0.1
method = "rejection"
unif.eps.rej = abc(target = criteria, param=uniform.sim.eps[,c("h1","h2","rel.eff","eps.int")], sumstat = uniform.sim.eps[,c("eu","sa","tan_eusa","bab_eusa")], tol = tol, method = method)
#quartz()
hist(unif.eps.rej)
tol = 0.05
method = "rejection"
unif.eps.rej.05 = abc(target = criteria, param=uniform.sim.eps[,c("h1","h2","rel.eff","eps.int")], sumstat = uniform.sim.eps[,c("eu","sa","tan_eusa","bab_eusa")], tol = tol, method = method)
#quartz()
hist(unif.eps.rej.05)
tol = 0.01
method = "rejection"
unif.eps.rej.01 = abc(target = criteria, param=uniform.sim.eps[,c("h1","h2","rel.eff","eps.int")], sumstat = uniform.sim.eps[,c("eu","sa","tan_eusa","bab_eusa")], tol = tol, method = method)
#quartz()
hist(unif.eps.rej.01)

# do abc for no eps
tol = 0.1
method = "rejection"
unif.noeps.rej =abc(target = c(criteria), param=uniform.sim.noeps[,c("h1","h2","rel.eff")], sumstat = uniform.sim.noeps[,c(names(criteria))], tol = tol, method = method)
#quartz()
hist(unif.noeps.rej)
str(unif.noeps.rej)
tol = 0.05
unif.noeps.rej.05 =abc(target = c(criteria), param=uniform.sim.noeps[,c("h1","h2","rel.eff")], sumstat = uniform.sim.noeps[,c(names(criteria))], tol = tol, method = method)
#quartz()
hist(unif.noeps.rej)
tol = 0.01
unif.noeps.rej.01 =abc(target = c(criteria), param=uniform.sim.noeps[,c("h1","h2","rel.eff")], sumstat = uniform.sim.noeps[,c(names(criteria))], tol = tol, method = method)
#quartz()
hist(unif.noeps.rej.05)
hist(unif.noeps.rej.01)
str(unif.noeps.rej.01[["ss"]])
# percent criteria  
sum(unif.noeps.rej.01[["ss"]][,1] > 1 & unif.noeps.rej.01[["ss"]][,2] < 1)/length(unif.noeps.rej.01[["ss"]][,2])
# 0.917
sum(unif.eps.rej.01[["ss"]][,1] > 1 & unif.eps.rej.01[["ss"]][,2] < 1)/length(unif.eps.rej.01[["ss"]][,2])

# plot accepted and all values for abc
cat.names=c(expression("EU(lOR"[tan]*"/lOR"[bab1]*")"),expression("SA(lOR"[tan]*"/lOR"[bab1]*")"),
            expression("EU(lOR"[tan]*")/SA(lOR"[tan]*")"),expression("EU(lOR"[bab1]*")/SA(lOR"[bab1]*")"))
plot_acc_hist(unif.eps.rej,uniform.sim.eps,criteria,cat.names=cat.names,filename="uni_values_witheps_spec")
plot_acc_hist(unif.noeps.rej,uniform.sim.noeps,criteria,cat.names=cat.names,filename="uni_values_noeps")

# tol 0.05
plot_acc_hist(unif.eps.rej.05,uniform.sim.eps,criteria,cat.names=cat.names,filename="uni_values_witheps_0.05_spec")
plot_acc_hist(unif.noeps.rej.05,uniform.sim.noeps,criteria,cat.names=cat.names,filename="uni_values_noeps_0.05")
# tol 0.01
plot_acc_hist(unif.eps.rej.01,uniform.sim.eps,criteria,cat.names=cat.names,filename="uni_values_witheps_0.01_spec")
plot_acc_hist(unif.noeps.rej.01,uniform.sim.noeps,criteria,cat.names=cat.names,filename="uni_values_noeps_0.01")

get.mean.dists <- function(abc.obj,criteria) {
  sapply(names(criteria), function (x) c( mean(abs(abc.obj[["ss"]][,x] - criteria[x])),
                                          sd(abs(abc.obj[["ss"]][,x] - criteria[x]))/sqrt(length(abc.obj[["ss"]][,x])) ) )
}

dev.off()

mdist.eps=get.mean.dists(unif.eps.rej,criteria)
mdist.noeps=get.mean.dists(unif.noeps.rej,criteria)
mdist = rbind(mdist.noeps[1,],mdist.eps[1,])
barplot(mdist, beside=T,legend=c("no eps.","with eps. interaction"),ylab="mean distance",names.arg=cat.names,cex.names=0.8)
dev.copy2pdf(file="uni_mean_dist_comp_spec.pdf")


plot(unif.noeps.rej)
# cross validation
all.sim=rbind(uniform.sim.noeps[,c(names(criteria))],uniform.sim.eps[,c(names(criteria))])
all.sim$model=factor(c(rep("no.eps",nrow(uniform.sim.noeps)),rep("eps",nrow(uniform.sim.eps))))
modelsel=postpr(target = criteria,index=as.character(all.sim$model),sumstat=all.sim[,c(names(criteria))],tol=0.1,method="rejection")
summary(modelsel)
#  old eps
# Proportion of accepted simulations (rejection):
#   eps no.eps 
# 0.7406 0.2594 
# 
# Bayes factors:
#   eps no.eps
# eps    1.0000 2.8548

# Proportion of accepted simulations (rejection):
#   eps no.eps 
# 0.8603 0.1397 
# 
# Bayes factors:
#   eps no.eps
# eps    1.0000 6.1577
# no.eps 0.1624 1.0000
modelsel05=postpr(target = criteria,index=as.character(all.sim$model),sumstat=all.sim[,c(names(criteria))],tol=0.05,method="rejection")
summary(modelsel05)
# old eps:
# Proportion of accepted simulations (rejection):
#   eps no.eps 
# 0.8949 0.1051 
# 
# Bayes factors:
#   eps no.eps
# eps    1.0000 8.5147

# Proportion of accepted simulations (rejection):
#   eps no.eps 
# 0.9584 0.0416 
# 
# Bayes factors:
#   eps  no.eps
# eps     1.0000 23.0616
# no.eps  0.0434  1.0000

modelsel01=postpr(target = criteria,index=as.character(all.sim$model),sumstat=all.sim[,c(names(criteria))],tol=0.01,method="rejection")
summary(modelsel01)
## for old eps: 1.5mio eps against 2 mio noeps
# Proportion of accepted simulations (rejection):
#   eps no.eps 
# 0.9984 0.0016 
# 
# Bayes factors:
#   eps   no.eps
# eps      1.0000 644.1613


# Proportion of accepted simulations (rejection):
#   eps no.eps 
# 0.9986 0.0014 
# 
# Bayes factors:
#   eps   no.eps
# eps      1.0000 688.6552
# no.eps   0.0015   1.0000

crossval=cv4postpr(index=all.sim$model,all.sim[,c(names(criteria))],nval=10,tols=0.1,method="rejection")

# plot accepted and all values for abc
plot_acc_hist <- function(abc.obj,uniform.sim.eps,criteria,cat.names=NULL,filename=FALSE,xlims=list(c(0,5),c(0,3),c(0,3),c(0,3)),ylims=list(c(0,1.5),c(0,1.5),c(0,2.0),c(0,1.5))){
  quartz()
  par(mfrow=c(2,2))
  cols=c(rgb(0.2,0.2,0.2,0.5),rgb(0.9,0.9,0.9,0.5))
  if (length(cat.names) == 0) {
    cat.names = names(criteria)
  }
  for(i in 1:length(criteria)){
    hist(abc.obj[["ss"]][,names(criteria)[i]],xlab="value",main=cat.names[i],breaks=c(-Inf,seq(-10,10,by=0.25),Inf),
         xlim=xlims[[i]],ylim=ylims[[i]],col=cols[1]) 
    hist(uniform.sim.eps[,names(criteria)[i]],breaks=c(-Inf,seq(-10,10,by=0.25),Inf),col=cols[2],add=TRUE) 
    
    abline(v=criteria[i],col="red")
  }
  legend("topright", c("accepted", "all"), title="Simulations",fill=cols)
  if (filename != FALSE){
    dev.copy2pdf(file=paste0("abc_accepted_rej_values",filename,".pdf"))
  }
  quartz()
  par(mfrow=c(2,2))
  for(i in 1:length(criteria)){
    hist(abs(abc.obj[["ss"]][,names(criteria)[i]] -criteria[i]),xlab="distance",main=cat.names[i],breaks=c(seq(0,10,by=0.25),Inf),
         xlim=xlims[[i]],ylim=ylims[[i]],col=cols[1]) 
    hist(abs(uniform.sim.eps[,names(criteria)[i]] - criteria[i]),breaks=c(seq(0,10,by=0.25),Inf),col=cols[2],add=TRUE) 
  }
  legend("topright", c("accepted", "all"), title="Simulations",fill=cols)
  if (filename != FALSE){
    dev.copy2pdf(file=paste0("abc_accepted_rej_distances",filename,".pdf"))
  }
}

par.old=par
plot_hists <- function (abc.obj,cols=c("h1","h2","rel.eff","eps.int"),
                        titles=c(expression("dominance"~italic("tan")),expression("dominance"~italic("bab1")),expression("allelic eff. of"~italic("bab1")~"relative to"~italic("tan")),"epistatic interaction strength"), xlims=list(c(0,1),c(0,1),c(0.5,4),c(-1,3)), 
                        fileprefix="abc_hist") {
  # plots separate histograms for each parameter
  pdf(file=paste0(fileprefix,"accepted_sims.pdf"),width=7,height=7)
  #quartz()
  layout(matrix(1:length(cols),nrow=2,byrow=T))
  if(length(cols) == 3)
  layout(matrix(c(1,1,1,1,2,2,2,2,0,0,3,3,3,3,0,0), nrow=2, 8, byrow = TRUE))
  for (i in 1:length(cols)){
    par(mar=c(5,4,2,2)+0.1)
    hist(abc.obj[["unadj.values"]][,cols[i]],xlab=titles[i], xlim=xlims[[i]], main=NULL,freq=FALSE)
  }
  dev.off()
}

plot_hists(unif.eps.rej,cols=c("h1","h2","rel.eff","eps.int"),fileprefix="abc_hist_uni_eps_spec")
plot_hists(unif.noeps.rej,cols=c("h1","h2","rel.eff"),fileprefix="abc_hist_uni_noeps_")

plot_hists(unif.eps.rej.05,cols=c("h1","h2","rel.eff","eps.int"),fileprefix="abc_hist_uni_eps_05_spec")
plot_hists(unif.noeps.rej.05,cols=c("h1","h2","rel.eff"),fileprefix="abc_hist_uni_noeps_05")

plot_hists(unif.eps.rej.01,cols=c("h1","h2","rel.eff","eps.int"),fileprefix="abc_hist_uni_eps_01_spec")
plot_hists(unif.noeps.rej.01,cols=c("h1","h2","rel.eff"),fileprefix="abc_hist_uni_noeps_01")

save.image(file="mon_02_08_2015.RSession",compress = T)


load("mon_02_08_2015.RSession")
# bin values
cut.h=seq(-0.125,1.125,by=0.25)
cut.eps= seq(-1.5,3.5,by=1.0)
cut.eff=c(0.25,0.75,seq(1.5,4.5,by=1.0))
uniform.sim.eps$eps.cut=cut(uniform.sim.eps$eps.int,breaks=cut.eps)
levels(uniform.sim.eps$eps.cut)=as.character(-1:3)
uniform.sim.eps$eff.cut=cut(uniform.sim.eps$rel.eff,breaks=cut.eff)
levels(uniform.sim.eps$eff.cut)=as.character(c(0.5,1,2,3,4))
uniform.sim.eps$h1.cut=cut(uniform.sim.eps$h1,breaks=cut.h)
levels(uniform.sim.eps$h1.cut)=as.character(seq(0,1,by=0.25))
uniform.sim.eps$h2.cut=cut(uniform.sim.eps$h2,breaks=cut.h)
levels(uniform.sim.eps$h2.cut)=as.character(seq(0,1,by=0.25))
uniform.sim.eps$h1xh2xeffxeps=interaction(uniform.sim.eps[,c("h1.cut","h2.cut","eff.cut","eps.cut")],sep=":")
# 875 combinations of levels
# calculate max values for bins
unif.eps.counts=summary(uniform.sim.eps$h1xh2xeffxeps,maxsum=1500)
hist(unif.eps.counts,breaks=25)
# no eps
uniform.sim.noeps$eff.cut=cut(uniform.sim.noeps$rel.eff,breaks=cut.eff)
levels(uniform.sim.noeps$eff.cut)=as.character(c(0.5,1,2,3,4))
uniform.sim.noeps$h1.cut=cut(uniform.sim.noeps$h1,breaks=cut.h)
levels(uniform.sim.noeps$h1.cut)=as.character(seq(0,1,by=0.25))
uniform.sim.noeps$h2.cut=cut(uniform.sim.noeps$h2,breaks=cut.h)
levels(uniform.sim.noeps$h2.cut)=as.character(seq(0,1,by=0.25))
uniform.sim.noeps$h1xh2xeff=interaction(uniform.sim.noeps[,c("h1.cut","h2.cut","eff.cut")],sep=":")
# 125 categories
unif.noeps.counts=summary(uniform.sim.noeps$h1xh2xeff,maxsum=1500)
unif.noeps.countsb = sum.cuts(uniform.sim.noeps,eps=F)
hist(unif.noeps.counts,breaks=25)
# construct category vectors
h1h2=expand.grid(seq(0,1,by=0.25),seq(0,1,by=0.25))
colnames(h1h2)=c("h1","h2")
rownames(h1h2)=apply(h1h2,1,paste0,collapse=":")
effxeps.u=expand.grid(c(0.5,1,2,3,4),-1:3)
rownames(effxeps.u)=apply(effxeps.u,1,paste0,collapse=":")
colnames(effxeps.u)=c("rel.eff","eps.int")

fill.unif.matrix <- function(counts,rows,cols){
  count.mat=matrix(0,nrow=length(rows),ncol=length(cols))
  rownames(count.mat)=rows
  colnames(count.mat)=cols
  for (cmb in names(counts) ){
    cmbs=strsplit(cmb,':',fixed=FALSE,perl=TRUE)[[1]]
    row=paste0(cmbs[1:2],collapse=":")
    col=paste0(cmbs[3:length(cmbs)],collapse=":")
    count.mat[row,col]=counts[cmb]
  }
  return(count.mat)
}

# fill into matrix
all.unif.eps.fq=fill.unif.matrix(unif.eps.counts,rownames(h1h2),rownames(effxeps.u))
all.unif.noeps.fq = fill.unif.matrix(unif.noeps.counts,rownames(h1h2),levels(uniform.sim.noeps$eff.cut))

sum.cuts <- function(reja,eff.cut=c(0.5,1,2,3,4),h.cut=seq(0,1,by=0.25),eps.cut=-1:3
                     ,cut.eff.=cut.eff,cut.h.=cut.h,cut.eps.=cut.eps,eps=T) {
  rej = data.frame(reja)
  rej$eff.cut=cut(rej$rel.eff,breaks=cut.eff)
  levels(rej$eff.cut)=as.character(eff.cut)
  rej$h1.cut=cut(rej$h1,breaks=cut.h)
  levels(rej$h1.cut)=as.character(h.cut)
  rej$h2.cut=cut(rej$h2,breaks=cut.h)
  levels(rej$h2.cut)=as.character(h.cut)
  if (eps == TRUE){
    rej$eps.cut=cut(rej$eps.int,breaks=cut.eps)
    levels(rej$eps.cut)=as.character(eps.cut)
    rej$h1xh2xeffxeps=interaction(rej[,c("h1.cut","h2.cut","eff.cut","eps.cut")],sep=":")
    return(summary(rej$h1xh2xeffxeps,maxsum=15000))
  }
  else {
    rej$h1xh2xeff=interaction(rej[,c("h1.cut","h2.cut","eff.cut")],sep=":")
    # calculate max values for bins
    return(summary(rej$h1xh2xeff,maxsum=15000))
  }
}


eff.cut=c(0.5,1,2,3,4)
# tol 0.05
rej.unif.noeps.05.counts=sum.cuts(unif.noeps.rej.05[["unadj.values"]],eps=F)
rej.unif.eps.05.counts=sum.cuts(unif.eps.rej.05[["unadj.values"]],eps=T)
rej.unif.eps.05.fq=fill.unif.matrix(rej.unif.eps.05.counts,rownames(h1h2),rownames(effxeps.u))
rej.unif.noeps.05.fq = fill.unif.matrix(rej.unif.noeps.05.counts,rownames(h1h2),as.character(eff.cut))
rej.unif.eps.05.fq=rej.unif.eps.05.fq/all.unif.eps.fq
rej.unif.noeps.05.fq=rej.unif.noeps.05.fq/all.unif.noeps.fq
# 0.01
rej.unif.noeps.01.counts=sum.cuts(unif.noeps.rej.01[["unadj.values"]],eps=F)
rej.unif.eps.01.counts=sum.cuts(unif.eps.rej.01[["unadj.values"]],eps=T)
rej.unif.eps.01.fq=fill.unif.matrix(rej.unif.eps.01.counts,rownames(h1h2),rownames(effxeps.u))
rej.unif.noeps.01.fq = fill.unif.matrix(rej.unif.noeps.01.counts,rownames(h1h2),as.character(eff.cut))
rej.unif.eps.01.fq=rej.unif.eps.01.fq/all.unif.eps.fq
rej.unif.noeps.01.fq=rej.unif.noeps.01.fq/all.unif.noeps.fq

quartz(width=10,height=8)
max.scale=as.integer(max(rej.unif.eps.05.fq)*10+1)*10
plot.heatmap.fq(rej.unif.eps.05.fq*100,doms=doms,breaknum=50,effxeps=effxeps.u,hs2=h1h2,
                print.breaks=seq(0,max.scale,by=10),mybreaks=seq(0,max.scale,length.out=50+1),best.tresh=1,scale.name="Percent")
dev.copy2pdf(file="heatmap_rej_method_unif_eps_05.pdf")

quartz(width=9,height=4.0)
max.scale=as.integer(max(rej.unif.noeps.05.fq)*10+1)*10
plot.heatmap.fq(rej.unif.noeps.05.fq*100,breaknum=50,no.eps=T,heights=c(4,2),effxeps=effxeps.u,hs2=h1h2,
                print.breaks=seq(0,max.scale,by=10),mybreaks=seq(0,max.scale,length.out=50+1),best.tresh=1,scale.name="Percent")
dev.copy2pdf(file="heatmap_rej_method_unif_noeps_05.pdf")

quartz(width=10,height=8)
max.scale=as.integer(max(rej.unif.eps.01.fq)*10+1)*10
plot.heatmap.fq(rej.unif.eps.01.fq*100,doms=doms,breaknum=50,effxeps=effxeps.u,hs2=h1h2,
                print.breaks=seq(0,max.scale,by=10),
                mybreaks=seq(0,max.scale,length.out=50+1),best.tresh=1,scale.name="Percent")
dev.copy2pdf(file="heatmap_rej_method_unif_eps_01_spec.pdf")

quartz(width=9,height=4.0)
max.scale=as.integer(max(rej.unif.noeps.01.fq)*10+1)*10
plot.heatmap.fq(rej.unif.noeps.01.fq*100,breaknum=50,no.eps=T,heights=c(4,2),
                effxeps=effxeps.u,hs2=h1h2,print.breaks=seq(0,max.scale,by=10),
                mybreaks=seq(0,max.scale,length.out=50+1),best.tresh=1,scale.name="Percent")
dev.copy2pdf(file="heatmap_rej_method_unif_noeps_01.pdf")

# get the rejected counts
rej.unif.noeps = data.frame(unif.noeps.rej[["unadj.values"]])
rej.unif.noeps$eff.cut=cut(rej.unif.noeps$rel.eff,breaks=cut.eff)
levels(rej.unif.noeps$eff.cut)=as.character(c(0.5,1,2,3,4))
rej.unif.noeps$h1.cut=cut(rej.unif.noeps$h1,breaks=cut.h)
levels(rej.unif.noeps$h1.cut)=as.character(seq(0,1,by=0.25))
rej.unif.noeps$h2.cut=cut(rej.unif.noeps$h2,breaks=cut.h)
levels(rej.unif.noeps$h2.cut)=as.character(seq(0,1,by=0.25))
rej.unif.noeps$h1xh2xeff=interaction(rej.unif.noeps[,c("h1.cut","h2.cut","eff.cut")],sep=":")
# calculate max values for bins
rej.unif.noeps.counts=summary(rej.unif.noeps$h1xh2xeff,maxsum=1500)

# get the rejected counts
rej.unif.eps = data.frame(unif.eps.rej[["unadj.values"]])
rej.unif.eps$eps.cut=cut(rej.unif.eps$eps.int,breaks=cut.eps)
levels(rej.unif.eps$eps.cut)=as.character(-1:3)
rej.unif.eps$eff.cut=cut(rej.unif.eps$rel.eff,breaks=cut.eff)
levels(rej.unif.eps$eff.cut)=as.character(c(0.5,1,2,3,4))
rej.unif.eps$h1.cut=cut(rej.unif.eps$h1,breaks=cut.h)
levels(rej.unif.eps$h1.cut)=as.character(seq(0,1,by=0.25))
rej.unif.eps$h2.cut=cut(rej.unif.eps$h2,breaks=cut.h)
levels(rej.unif.eps$h2.cut)=as.character(seq(0,1,by=0.25))
rej.unif.eps$h1xh2xeffxeps=interaction(rej.unif.eps[,c("h1.cut","h2.cut","eff.cut","eps.cut")],sep=":")
# 875 combinations of levels
# calculate max values for bins
rej.unif.eps.counts=summary(rej.unif.eps$h1xh2xeffxeps,maxsum=1500)

rej.unif.eps.fq=fill.unif.matrix(rej.unif.eps.counts,rownames(h1h2),rownames(effxeps.u))
rej.unif.noeps.fq = fill.unif.matrix(rej.unif.noeps.counts,rownames(h1h2),levels(uniform.sim.noeps$eff.cut))

rej.unif.eps.fq=rej.unif.eps.fq/all.unif.eps.fq
rej.unif.noeps.fq=rej.unif.noeps.fq/all.unif.noeps.fq

doms=unique(sort(unlist(h1h2[,1])))
eff2=as.numeric(unique(effxeps.u[,1]))
eps.int=as.numeric(unique(sort(effxeps.u[,2])))

save.image(file="thur07_08_2015.RSession",compress=T)
#load("thur07_08_2015.RSession")
quartz(width=10,height=8)
plot.heatmap.fq(rej.unif.eps.fq,doms=doms,breaknum=50,effxeps=effxeps.u,hs2=h1h2)
dev.copy2pdf(file="heatmap_rej_method_unif_eps.pdf")

quartz(width=9,height=4.0)
plot.heatmap.fq(rej.unif.noeps.fq,breaknum=50,no.eps=T,heights=c(4,2),effxeps=effxeps.u,hs2=h1h2)
dev.copy2pdf(file="heatmap_rej_method_unif_noeps.pdf")

# eps.int > -2
eps.not=grep("-[23]",colnames(rej.unif.eps.fq))
quartz(width=10,height=8)
plot.heatmap.fq(rej.unif.eps.fq[,- eps.not],breaknum=50,effxeps=effxeps.u[effxeps.u$eps.int > -2,],hs2=h1h2)
dev.copy2pdf(file="heatmap_rej_method_unif_eps.pdf")

library(reshape2,ggplot2)
source("/Volumes/Temp/Lukas/Tools/Scripts/R/image.scale.R")
plot.heatmap.fq <- function(fq.array,doms=NULL,breaknum=20, no.eps=FALSE,widths = c(7,1)
                            , heights = c(2,1), hs2=hs2, eff2=NULL, effxeps=effxeps
                            ,eps.int=NULL, max=FALSE,print.breaks=seq(0,1,by=0.1)
                            ,mybreaks=NULL,best.tresh = 0.05,scale.name="Fraction" )
{
  if(length(doms) == 0){
    doms <- unique(hs2[,1])
  }
  if (length(eff2) == 0){ 
    eff2 <- unique(effxeps[,1])
  }
  if (length(eps.int) == 0){ 
    eps.int <- unique(effxeps[,2])
  }
  if (max == TRUE) {
    max = max(fq.array)
  }
  else {
    max = 1.0
  }
  plot.dat=fq.array/max
  ## cols=c(colorRampPalette(c("lightblue","darkblue","red"))(breaknum))    
  cols=rev(heat.colors(breaknum))
  ## cols=rev(terrain.colors(breaknum))
  ## cols=rev(topo.colors(breaknum))
  ## cols=colorRampPalette(brewer.pal(9,"Greens"))(50)
  if (length(mybreaks) == 0) {
    mybreaks=seq(0,1,length.out=breaknum+1)
  }
  #print.ends=12
  horizontal = FALSE 
  if (no.eps){
    horizontal = TRUE
    layout(matrix(1:2, ncol = 1,byrow=T), heights = heights, respect = FALSE)
    par(mar=c(3,4,4,1))        
    # horizontal image
    image(1:length(rownames(hs2)),1:length(eff2),plot.dat,col=cols, breaks=sort(mybreaks),xaxt='n',yaxt='n',xlab="",ylab="")
    abline(h=(1:length(eff2))+0.5)
    abline(v=(1:length(rownames(hs2))+0.5))
    # draw lines
    abline(v=seq(1,length(hs2[,1])/length(doms)-1)*length(doms)+0.5,lw=3)  
    axis(2,at=1:length(eff2),labels=abs(eff2),cex.axis=0.80)
    mtext(expression("allelic eff. of"~italic("bab1")~"rel. to"~italic("tan")), side=2.5, line=2,cex=1.1)
    mtext(expression("dominance ("*italic("h")*") of the light allele"),side=3, line=2.5, cex=1.1)
    mtext(expression(italic("bab1")~"locus"), side=3, line=1.5,cex=1.1)
    mtext(expression(italic("tan")~"locus"),side=1, line=2.0, cex=1.1)
    # tan dominances
    mtext(rep(doms,length(hs2[,1])/length(doms)),side=1,at=1:length(plot.dat[,1]),line=1.0,las=3,adj=c(0.5,0.5),cex=0.85)
    # bab1 dominances
    mtext(doms,side=3,at=(1:length(doms))*length(hs2[,1])/length(doms)-length(doms)/2,line=0.5,las=1,adj=c(0.5,0.5),cex=0.85)
    iscale.par = c(4,4,1.75,1)
    iscale.side = 1
  }
  else{
    layout(matrix(1:2, ncol = 2,byrow=T), widths = widths, respect = FALSE)
    par(mar=c(5,7,5,1))    
    image(1:length(rownames(effxeps)),1:length(rownames(hs2)),t(plot.dat),col=cols, breaks=sort(mybreaks),xaxt='n',yaxt='n',xlab="",ylab="")
    abline(v=(1:length(rownames(effxeps))+0.5))
    abline(v=seq(1,length(rownames(effxeps))/length(eff2)-1)*length(eff2)+0.5,lw=3)
    axis(1,at=1:length(rownames(effxeps)),labels=rep(abs(eff2),length(eps.int)),cex.axis=0.80)
    mtext(expression("effect of"~italic("bab1")~"relative to"~italic("tan")), side=1, line=2.5,cex=1.5)
    mtext("dominance",side=2, line=5.0,cex=1.25)
    mtext("locus", at=length(hs2[,1])+2,side=2, line=2.75,las=1,adj=0.5, cex=1.25)
    mtext(expression(italic(" tan")~":"~italic("bab1")), at=length(hs2[,1])+1,side=2, line=2.75,las=1,adj=0.5,cex=1.25)
    #mtext(c(0,0.25,0.5,0.75,1),side=2,at=seq(0,length(hs2[,1])/5-1)*5+3,las=1,line=4)
    mtext(rep(doms,length(hs2[,1])/length(doms)),side=2,at=1:length(plot.dat[,1]),line=4,las=1,adj=c(0.5,0.5),cex=0.85)
    mtext(":",side=2,at=1:length(plot.dat[,1]),las=1,line=2.75,cex=0.85)
    #axis(2,at=1:length(hs2[,1]),labels=rep(c(0,0.25,0.5,0.75,1),5),las=1)
    axis(2,at=1:length(hs2[,1]),labels=NA,las=1)
    mtext(sapply(doms,function(x) rep(x,length(doms))),side=2,at=1:length(hs2[,1]),las=1,line=1.5,adj=c(0.5,0.5),cex=0.85)
    #    mtext(sapply(c(1,0.75,0.5,0.25,0),function(x) rep(x,5)),side=2,at=1:length(hs2[,1]),las=1,line=1.5,adj=c(0.5,0.5))
    mtext(eps.int,side=3,at=seq(0,length(eps.int)-1)*length(eff2)+2.5,las=1,line=0.5,cex=1.25)
    mtext("epistatic interaction strength",side=3,at=length(rownames(effxeps))/2+0.5,line=2.0,cex=1.5)
    # draw lines
    abline(h=(1:length(hs2[,1])+0.5))
    abline(h=seq(1,length(hs2[,1])/length(doms)-1)*length(doms)+0.5,lw=3)   
    iscale.par = c(5,0,5,5.5)
    iscale.side=4 
  }
  #axis(2, at=seq(0,1,,dim(db)[2]), labels=colnames(db))
  plot.rank=array(length(plot.dat)-rank(plot.dat),dim=dim(plot.dat))
  for(i in 1:length(plot.dat[,1])){
    for(j in 1:length(plot.dat[1,])){
      x = j
      y = i
      if (no.eps) {
        x=i
        y=j
      }
      if (plot.rank[i,j] < 5 & plot.dat[i,j] >= best.tresh){
        text(x,y,prettyNum(plot.dat[i,j],digits=2,drop0trailing=T),col="black",cex=0.6,font=2)
      }
      else if (plot.rank[i,j] < 10 & plot.dat[i,j] >= best.tresh){  #(plot.dat[i,j] >= 1e-9){
        #if plot.dat[i,j] >= plot.dat
        text(x,y,prettyNum(plot.dat[i,j],digits=2,drop0trailing=T),col="black",cex=0.6,font=1)
      }
      else if (plot.rank[i,j] < 50 & plot.dat[i,j] >= best.tresh){
        text(x,y,prettyNum(plot.dat[i,j],digits=2,drop0trailing=T),col="black",cex=0.6,font=1)
      }
    }
  }
  par(mar=iscale.par)
  image.scale(c(max(mybreaks) * 1.01,0.01), col=cols, breaks=sort(mybreaks),horiz=horizontal, yaxt="n",xaxt="n",xlab="",ylab="")  #print.breaks=c(0,-1,seq(-3,-24,by=-3))
  #  axis(iscale.side,at=print.breaks,labels=print.breaks, las=2)
  axis(iscale.side,at=print.breaks,labels=print.breaks)
  mtext(scale.name, side=iscale.side, line=2.0,cex=1.1)
}

# do everything with preselected eu > 1 and sa < 1
uniform.sim.eps$pre = uniform.sim.eps$eu > 1.0 & uniform.sim.eps$sa < 1.0
uniform.sim.eps$pre[is.na(uniform.sim.eps$pre)] = FALSE
uniform.sim.noeps$pre = uniform.sim.noeps$eu > 1.0 & uniform.sim.noeps$sa < 1.0
uniform.sim.noeps$pre[is.na(uniform.sim.noeps$pre)] = FALSE


tol = 0.1
method = "rejection"
unif.pre.eps.rej = abc(target = criteria, param=uniform.sim.eps[,c("h1","h2","rel.eff","eps.int")], subset=uniform.sim.eps$pre, sumstat = uniform.sim.eps[,c("eu","sa","tan_eusa","bab_eusa")], tol = tol, method = method)
#quartz()
hist(unif.pre.eps.rej)

# do abc for no eps
tol = 0.1
method = "rejection"
unif.pre.noeps.rej =abc(target = c(criteria), param=uniform.sim.noeps[,c("h1","h2","rel.eff")], subset=uniform.sim.noeps$pre, sumstat = uniform.sim.noeps[,c(names(criteria))], tol = tol, method = method)
#quartz()
hist(unif.pre.noeps.rej)


plot_acc_hist(unif.pre.eps.rej,uniform.sim.eps,criteria,cat.names=cat.names,filename="uni_pre_values_witheps")
plot_acc_hist(unif.pre.noeps.rej,uniform.sim.noeps,criteria,cat.names=cat.names,filename="uni_pre_values_noeps")

plot_hists(unif.pre.eps.rej,cols=c("h1","h2","rel.eff","eps.int"),fileprefix="abc_pre_hist_uni_eps_")
plot_hists(unif.pre.noeps.rej,cols=c("h1","h2","rel.eff"),fileprefix="abc_pre_hist_uni_noeps_")

# bin values
unif.pre.eps.counts=summary(uniform.sim.eps$h1xh2xeffxeps[uniform.sim.eps$pre],maxsum=1500)
hist(unif.pre.eps.counts,breaks=25)
# no eps
unif.pre.noeps.counts=summary(uniform.sim.noeps$h1xh2xeff[uniform.sim.noeps$pre],maxsum=1500)
hist(unif.pre.noeps.counts,breaks=25)
# fill into matrix
all.unif.pre.eps.fq=fill.unif.matrix(unif.pre.eps.counts,rownames(h1h2),rownames(effxeps.u))
all.unif.pre.noeps.fq = fill.unif.matrix(unif.pre.noeps.counts,rownames(h1h2),levels(uniform.sim.noeps$eff.cut))

# get the rejected counts
rej.pre.unif.noeps = data.frame(unif.pre.noeps.rej[["unadj.values"]])
rej.pre.unif.noeps$eff.cut=cut(rej.pre.unif.noeps$rel.eff,breaks=cut.eff)
levels(rej.pre.unif.noeps$eff.cut)=as.character(c(0.5,1,2,3,4))
rej.pre.unif.noeps$h1.cut=cut(rej.pre.unif.noeps$h1,breaks=cut.h)
levels(rej.pre.unif.noeps$h1.cut)=as.character(seq(0,1,by=0.25))
rej.pre.unif.noeps$h2.cut=cut(rej.pre.unif.noeps$h2,breaks=cut.h)
levels(rej.pre.unif.noeps$h2.cut)=as.character(seq(0,1,by=0.25))
rej.pre.unif.noeps$h1xh2xeff=interaction(rej.pre.unif.noeps[,c("h1.cut","h2.cut","eff.cut")],sep=":")
# calculate max values for bins
rej.pre.unif.noeps.counts=summary(rej.pre.unif.noeps$h1xh2xeff,maxsum=1500)

# get the rejected counts
rej.pre.unif.eps = data.frame(unif.pre.eps.rej[["unadj.values"]])
rej.pre.unif.eps$eps.cut=cut(rej.pre.unif.eps$eps.int,breaks=cut.eps)
levels(rej.pre.unif.eps$eps.cut)=as.character(-1:3)
rej.pre.unif.eps$eff.cut=cut(rej.pre.unif.eps$rel.eff,breaks=cut.eff)
levels(rej.pre.unif.eps$eff.cut)=as.character(c(0.5,1,2,3,4))
rej.pre.unif.eps$h1.cut=cut(rej.pre.unif.eps$h1,breaks=cut.h)
levels(rej.pre.unif.eps$h1.cut)=as.character(seq(0,1,by=0.25))
rej.pre.unif.eps$h2.cut=cut(rej.pre.unif.eps$h2,breaks=cut.h)
levels(rej.pre.unif.eps$h2.cut)=as.character(seq(0,1,by=0.25))
rej.pre.unif.eps$h1xh2xeffxeps=interaction(rej.pre.unif.eps[,c("h1.cut","h2.cut","eff.cut","eps.cut")],sep=":")
# 875 combinations of levels
# calculate max values for bins
rej.pre.unif.eps.counts=summary(rej.pre.unif.eps$h1xh2xeffxeps,maxsum=1500)

rej.pre.unif.eps.fq=fill.unif.matrix(rej.pre.unif.eps.counts,rownames(h1h2),rownames(effxeps.u))
rej.pre.unif.noeps.fq = fill.unif.matrix(rej.pre.unif.noeps.counts,rownames(h1h2),levels(uniform.sim.noeps$eff.cut))

#rej.pre.unif.eps.fq=rej.pre.unif.eps.fq/all.unif.pre.eps.fq
#rej.pre.unif.eps.fq[is.nan(rej.pre.unif.eps.fq)]=0
#rej.pre.unif.noeps.fq=rej.pre.unif.noeps.fq/all.unif.pre.noeps.fq
#rej.pre.unif.noeps.fq[is.nan(rej.pre.unif.noeps.fq)]=0

quartz(width=10,height=8)
plot.heatmap.fq(rej.pre.unif.eps.fq,doms=doms,breaknum=50,effxeps=effxeps.u,hs2=h1h2,max=TRUE)
dev.copy2pdf(file="heatmap_rej_method_unif_pre_eps.pdf")

quartz(width=9,height=4.0)
plot.heatmap.fq(rej.pre.unif.noeps.fq,breaknum=50,no.eps=T,heights=c(4,2),effxeps=effxeps.u,hs2=h1h2,max=TRUE)
dev.copy2pdf(file="heatmap_rej_method_unif_pre_noeps.pdf")


all.pre=c(uniform.sim.noeps$pre,uniform.sim.eps$pre)
modelsel_pre=postpr(target = criteria,index=as.character(all.sim$model),sumstat=all.sim[,c(names(criteria))],tol=0.1,
                method="rejection",subset=all.pre)
summary(modelsel_pre)

# Proportion of accepted simulations (rejection):
#   eps no.eps 
# 0.7354 0.2646 
# 
# Bayes factors:
#   eps no.eps
# eps    1.0000 2.7796
# no.eps 0.3598 1.0000

##### looking at the accepted simulations and selecting the right ones


str(unif.eps.rej)
post.eps=with(unif.eps.rej[[ss]], 



