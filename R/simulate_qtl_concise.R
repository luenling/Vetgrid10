library(msm)
library(foreach)
library(doMC)
library(reshape2)

registerDoMC(18)
getDoParWorkers()
setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Simulation")
gt.vals = c(2,1,0) ### value for each gt: aa:2, aA: 1, AA: 0

### get gt for locus
get.gt.locus <- function(gt.freq,gt.vals,total.ind){
    ### gets genotype frequencies and gt values for a locus and returns random genotypes for all ind. for this locus
    return(colSums(gt.vals*rmultinom(total.ind,1,gt.freq)))
}

### locus effects
loc.eff <- function(h.vals,g.type) { # takes the dominance and genotype vectors and returns a vector of gt effects (aa:1,aA:h,AA:0)
    return(g.type/2 + (g.type==1)*(h.vals-0.5))
}

### environmental variance of phenotype
v.E = 0.1
afs_cosm=c(0.2,0.83) # cosmopolitan allele frequencies
h2=0.3
get.V.E <- function(afs_loc,eff_loc,h_loc,h2=0.5,F.in=0,n=1000,eps.int=0.0){
    gt.freqs = sapply(afs_loc,function (p) c(p^2 + (p-p^2)*F.in,2*(p-p^2)*(1-F.in),(1-p)^2 + (p-p^2)*F.in))
    pos <- c()
    pos <- apply(gt.freqs,2,function (x) get.gt.locus(x,c(2,1,0),n))
    pos <- matrix(pos,nrow=length(afs_loc),ncol=n,byrow=T)
    ## pigmentation phenotypes resulting from genotypes
    ## additive effects (eff centered around 0)
    ind.eff <- apply(pos,2,function(x) loc.eff(h_loc,x)) 
    ## interaction component
    intxn <- ( eps.int*ind.eff[1,]*ind.eff[2,] - eps.int/2.0  )*mean(abs(eff_loc))
    ## additive effects (eff centered around 0)
    pheno <- eff_loc %*% ( ind.eff - 0.5 ) 
    pheno <- pheno + intxn
    return((1-h2)/h2*var(as.vector(pheno)))
}

run_simulations.eps= function(n.runs=50,v.E=0.0,total.ind =2100,n.replicates=3,n.pools=100,eps.int=0.0) {
    
#######  Initialize
    gt.freqs = sapply(n.freqs,function (p) c(p^2 + (p-p^2)*F.inbreed,2*(p-p^2)*(1-F.inbreed),(1-p)^2 + (p-p^2)*F.inbreed)) ### genotype frequencies for each locus, aa, aA, AA
    mean_eff=sum(gt.vals * gt.freqs %*% eff)
    d.replicates=total.ind/n.replicates # ind. per replicate
    results = foreach (i=1:n.runs,.combine="rbind") %dopar% {
        ## generate genotypes randomly (independent loci)
                                        # alleles in sample by binomial sampling from population. alleles
                                        # pos <- matrix(rbinom(n.loci*total.ind,2,p),nrow=n.loci,ncol=total.ind)    ### alleles A a, a:light allele, A: dark, gt.value: 2:aa 1:aA 0:AA
        pos <- c()
                                        #pos <- replicate(total.ind,get.gt(gt.freqs,gt.vals))
                                        #pos <- matrix(pos,nrow=n.loci,ncol=total.ind,byrow=T)
        pos <- apply(gt.freqs,2,function (x) get.gt.locus(x,gt.vals,total.ind))
        pos <- matrix(pos,nrow=n.loci,ncol=total.ind,byrow=T)
        # ind.eff only holds the dominance corrected effect of each locus, goes from 0 to 1 (or higher/lower in case of under/overdominance)
        ind.eff <- apply(pos,2,function(x) loc.eff(h.loci,x)) 
        ## pigmentation phenotypes resulting from genotypes
        ## interaction component, mean(abs(eff)) should always be total effect of both loci
        ## centered around zero roughly too
        int <- (eps.int*ind.eff[1,]*ind.eff[2,] - eps.int/2 )*mean(abs(eff))
        ## additive effects (ind.eff centered around 0)
        pheno <- eff %*% (ind.eff - 0.5)
        ## environmental effect
        pheno <- pheno + int + rnorm(total.ind,0,sqrt(v.E))
        ## get the randomized array of indices 
        indices=sample(1:total.ind)
        p.val <-  rep(0,n.loci)
        odds.ratio_ld <- rep(0,n.loci)
        odds.ratio_l <- rep(0,n.loci)
        odds.ratio_d <- rep(0,n.loci)
                                        # log.odds.ratio_ajb <- matrix(0,nrow=n.runs,ncol=n.loci)
                                        #     cors <- matrix(0,nrow=n.runs,ncol=4)
        af_avg <- rep(0,n.loci)
        af_light <- rep(0,n.loci)
        af_dark <- rep(0,n.loci)    
        
        dat = array(dim=c(n.loci,4,n.replicates)) # table for light dark cmh
        dat_l =  array(dim=c(n.loci,4,n.replicates)) # light base cmh
        dat_d =  array(dim=c(n.loci,4,n.replicates)) # dark base cmh
        alldark = c()
        alllight = c()
        ## get replicates
        for (j in 1:n.replicates) {
            ## indviduals in replicates
            sample.indices = indices[(1+(j-1)*d.replicates):(j*d.replicates)]
            ##  phenotypes in repl.		
            sampled.pheno = pheno[sample.indices]
            ## order indices by pheno and take extreme ones
            ## light allele high gt value: reverse order as light have hihger phenotypic value!
            o.pheno <- sample.indices[order(sampled.pheno,decreasing=T)]
            light   <- o.pheno[1:n.pools]   # take "n.pools" lightest colored individuals for light pool
            dark    <- o.pheno[((d.replicates-n.pools)+1):d.replicates] # take "n.pools" darkest individuals for dark pool
            alllight=c(alllight, light)
            alldark=c(alldark, dark)
                                        # par(mfrow=c(1,1)); b.points=seq(min(pheno[sample.indices]),max(pheno[sample.indices]),by=0.5); hist(pheno[sample.indices],breaks=b.points); hist(pheno[light], col="yellow",add=T,breaks=b.points); hist(pheno[dark], col="blue",add=T,breaks=b.points)
           
            ## genotypes for the pools
            geno.light  <-  pos[,light]
            geno.dark   <-  pos[,dark]
            ## allele counts       
            ac.light=rowSums(geno.light)
            ac.dark =rowSums(geno.dark)
            ac.base =rowSums(pos[,sample.indices])
            ac.rep=rowSums(pos[,sample.indices])
            dat[,,j]=array(c(ac.light,ac.dark, n.pools*2-ac.light, n.pools*2-ac.dark), dim=c( n.loci,4))
            dat_l[,,j]=array(c(ac.light,ac.base, n.pools*2-ac.light,d.replicates*2-ac.base), dim=c( n.loci,4))
            dat_d[,,j]=array(c(ac.dark,ac.base, n.pools*2-ac.dark,d.replicates*2-ac.base), dim=c( n.loci,4))
        }
        af_avg=rowSums(pos)/(2*total.ind)
        ## calculate OR and AFs
        for (x in 1:n.loci){
            mhtable = array(dat[x,,],dim=c(2,2,n.replicates))
            cmh.res= mantelhaen.test(mhtable)
            odds.ratio_ld[x] = cmh.res$estimate
            p.val[x]=cmh.res$p.value
            mhtable = array(dat_l[x,,],dim=c(2,2,n.replicates))
            odds.ratio_l[x] = mantelhaen.test(mhtable)$estimate
            mhtable = array(dat_d[x,,],dim=c(2,2,n.replicates))
            odds.ratio_d[x] = mantelhaen.test(mhtable)$estimate
            af_light[x]=sum(dat[x,1,])/(2*n.replicates*n.pools)
            af_dark[x]=sum(dat[x,2,])/(2*n.replicates*n.pools)
        }
        return(c(p.val,odds.ratio_ld,odds.ratio_l,odds.ratio_d,af_avg,af_light,af_dark))
    }
    
    ## return all tables in list structure
    l.names=c("pv","or","ol","od","af.m","af.l","af.d")
    ret.list=list()
    for (i in 1:length(l.names)){
        ret.list[[l.names[i]]]=results[,((i-1)*n.loci+1):(i*n.loci)]
    }
    return(ret.list)    
}


### Runs for assessing parameter combinations
### Fixed afs
afs.eu=c(0.2,0.83)
afs.sa=c(0.17,0.47)
### Fis
###F.sa=0.25
F.sa=0.35
F.eu=0.0
### base effect
### light allele effect positive
eff1=1.0
eps.eff=c(1)
### loci
n.loci=2
### heritability
h2=0.30
afs.cosm=c(0.2,0.83)
### Experimental setup
n.pools.sa = 65 # size of light and dark pools
n.pools.xs = 150 # size of light and dark pools
n.rep.sa = 3 #number of replicates
total.ind.sa =700*n.rep.sa # total individuals
n.pools.eu = 100 # size of light and dark pools
n.rep.eu = 6 #number of replicates
total.ind.eu =1500*n.rep.eu # total individuals
n.runs=5000
#n.runs=10
# parameters to go over
#eff2=c(3,2,1,0.5)
eff2=c(4,2,1,0.5) # relative effect 
#hs2=c(0.0,0.25,0.5,0.75,1.0)
hs2=c(0.0,0.25,0.5,0.75,1.0) #,1.25,1.5) # for dominance
hs2=expand.grid(hs2,hs2)
rownames(hs2)=apply(hs2,1,function (x) paste(x,collapse=":"))
eps.int=c(-1,-0.5,0,0.5,1.0,2.0,3.0) 
res.fields2=c("pv","or","ol","od","af.m","af.l","af.d")
two.loc.eps.sa=array(0,dim=c(length(hs2[,1]),length(eff2),length(eps.int),length(res.fields2),2,n.runs),dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"eps.int"=eps.int,"resfields"=res.fields2,"locus"=c(1,2),"run"=1:n.runs))
two.loc.eps.eu=array(0,dim=c(length(hs2[,1]),length(eff2),length(eps.int),length(res.fields2),2,n.runs),dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"eps.int"=eps.int,"resfields"=res.fields2,"locus"=c(1,2),"run"=1:n.runs))
#two.loc.eps.xs=array(0,dim=c(length(hs2[,1]),length(eff2),length(eps.int),length(res.fields2),2,n.runs),dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"eps.int"=eps.int,"resfields"=res.fields2,"locus"=c(1,2),"run"=1:n.runs))

#eff=c(0.5,-0.5)
for(j in 1:length(hs2[,1])){
    for(k in 1:length(eff2)){
        for(l in 1:length(eps.int)){
            eff=c(eff1/(1+abs(eff2[k])),eff1*eff2[k]/(1+abs(eff2[k])))
            h.loci=unlist(hs2[j,])
            v.E=get.V.E(afs.cosm,eff,h.loci,h2,eps.int=eps.int[l])
            F.inbreed=F.sa
            n.freqs=afs.sa
            run.result=run_simulations.eps(n.runs,v.E,total.ind=total.ind.sa,n.replicates=n.rep.sa,n.pools=n.pools.sa,eps.int=eps.int[l])
            two.loc.eps.sa[j,k,l,,1,]=t(sapply(res.fields2, function (x) run.result[[x]][,1]))
            two.loc.eps.sa[j,k,l,,2,]=t(sapply(res.fields2, function (x) run.result[[x]][,2]))
            ## run.result=run_simulations.eps(n.runs,v.E,total.ind=total.ind.sa,n.replicates=n.rep.sa,n.pools=n.pools.xs,eps.int=eps.int[l])
            ## two.loc.eps.xs[j,k,l,,1,]=t(sapply(res.fields2, function (x) run.result[[x]][,1]))
            ## two.loc.eps.xs[j,k,l,,2,]=t(sapply(res.fields2, function (x) run.result[[x]][,2]))

            F.inbreed=F.eu
            n.freqs=afs.eu
            run.result=run_simulations.eps(n.runs,v.E,total.ind=total.ind.eu,n.replicates=n.rep.eu,n.pools=n.pools.eu,eps.int=eps.int[l])
            two.loc.eps.eu[j,k,l,,1,]=t(sapply(res.fields2, function (x) run.result[[x]][,1]))
            two.loc.eps.eu[j,k,l,,2,]=t(sapply(res.fields2, function (x) run.result[[x]][,2]))
        }
    }
}

save(two.loc.eps.eu,two.loc.eps.sa,file="two.loc.eu.sa_d.RData")
load("two.loc.eu.sa_d.RData")

### calculate variances
criteria=c(2.51,0.57,1.89,0.44)
names(criteria)=c("eu","sa","tan_eusa","bab_eusa")


### calculate variances
get.var <- function (two.loc.eps.1,two.loc.eps.2,citeria,drop.eff=1)
{
  var.crit.eps =  array(0,dim=c(length(hs2[,1]),length(eff2)-drop.eff,length(eps.int),length(criteria)), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2[1:(length(eff2)-drop.eff)],"eps.int"=eps.int,"criteria"=names(criteria)))
  for(j in 1:length(hs2[,1])){
    for(k in 1:(length(eff2)-drop.eff)){
      for(l in 1:length(eps.int)){
        a=log(two.loc.eps.1[j,k,l,"or",1,])/log(two.loc.eps.1[j,k,l,"or",2,])
        var.crit.eps[j,k,l,1]=var((a[is.finite(a) & ! is.na(a)]-criteria[1]))
        a=log(two.loc.eps.2[j,k,l,"or",1,])/log(two.loc.eps.2[j,k,l,"or",2,])
        var.crit.eps[j,k,l,2]=var((a[is.finite(a) & ! is.na(a)]-criteria[2]))
        a=log(two.loc.eps.1[j,k,l,"or",1,])/log(two.loc.eps.2[j,k,l,"or",1,])
        var.crit.eps[j,k,l,3]=var((a[is.finite(a) & ! is.na(a)]-criteria[3]))
        a=log(two.loc.eps.1[j,k,l,"or",2,])/log(two.loc.eps.2[j,k,l,"or",2,])
        var.crit.eps[j,k,l,4]=var((a[is.finite(a) & ! is.na(a)]-criteria[4]))
      }
    }
  }
  return(var.crit.eps)
}


### use upper limit for OR
min.or=1.0
max.or=1.0
a=two.loc.eps.eu[,,,"or",,]
min.or=min(1.0,min(a[which(a != 0)]))
max.or=max(1.0,max(a[which(a != 0)]))
min.or=min.or/2
max.or=max.or*2
a=two.loc.eps.eu[,,,"or",,]
a[a <= min.or] = min.or; a[a >= max.or] = max.or
two.loc.eps.eu[,,,"or",,]=a
a=two.loc.eps.sa[,,,"or",,]
min.or=min(1.0,min(a[which(a != 0)]))
max.or=max(1.0,max(a[which(a != 0)]))
### set the infinite or 0 ORs to min/2 or 2*max 
min.or=min.or/2
max.or=max.or*2
a=two.loc.eps.sa[,,,"or",,]
a[a <= min.or] = min.or; a[a >= max.or] = max.or
two.loc.eps.sa[,,,"or",,]=a

var.crit.eps =get.var(two.loc.eps.eu,two.loc.eps.sa,citeria,drop.eff=0)

## test.array=melt(var.crit.eps[,,,1],id=c("rel.eff"))
## boxplot(value~hs,data=melt(test.array))
med.var.crit=sapply(1:4,function(x) median(var.crit.eps[,,,x],na.rm=T))
mean.var.crit=sapply(1:4,function(x) mean(var.crit.eps[,,,x],na.rm=T))


### calculate LHs
kernel.dens <- function(a,criteria,crit,med.var.crit.int){
  n=length(a)
  bw=tryCatch(bw.nrd0((a-criteria[crit])/med.var.crit.int[crit]),
              error=function (e) return(NA) )
  #return(mean((a-criteria[crit])^2))
  return(1/(2*pi*n*bw)*sum(exp(-(a-criteria[crit])^2/(2*bw*med.var.crit.int[crit]))))
  #return(density((a-criteria[crit])/med.var.crit.int[crit],bw=bw.nrd0((a-criteria[crit])/med.var.crit.int[crit]),kernel="gaussian"))
}



get.lh <- function (two.loc.eps.1,two.loc.eps.2,med.var.crit.int,criteria,drop.eff=0)
{
  lh.two.loc.eps=array(0,dim=c(length(hs2[,1]),length(eff2)-drop.eff,length(eps.int),length(criteria)), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2[1:(length(eff2)-drop.eff)],"eps.int"=eps.int,"criteria"=names(criteria)))
  for(j in 1:length(hs2[,1])){
    for(k in 1:(length(eff2)-drop.eff)){
      for(l in 1:length(eps.int)){
        a=log(two.loc.eps.1[j,k,l,"or",1,])/log(two.loc.eps.1[j,k,l,"or",2,])
        attr(a,"values")=paste(j,k,l)
        a=a[! is.na(a) & is.finite(a)]; #a[a>=1e4]=1e4; a[a<=-1e4]=-1e4
        crit=1
        lh.two.loc.eps[j,k,l,crit]=kernel.dens(a,criteria,crit,med.var.crit.int)
        a=log(two.loc.eps.2[j,k,l,"or",1,])/log(two.loc.eps.2[j,k,l,"or",2,])
        a=a[! is.na(a) & is.finite(a)]; #a[a>=1e4]=1e4; a[a<=-1e4]=-1e4
        crit=2
        lh.two.loc.eps[j,k,l,crit]=kernel.dens(a,criteria,crit,med.var.crit.int)
        a=log(two.loc.eps.1[j,k,l,"or",1,])/log(two.loc.eps.2[j,k,l,"or",1,])
        a=a[! is.na(a) & is.finite(a)]; #a[a>=1e4]=1e4; a[a<=-1e4]=-1e4
        crit=3
        lh.two.loc.eps[j,k,l,crit]=kernel.dens(a,criteria,crit,med.var.crit.int)
        a=log(two.loc.eps.1[j,k,l,"or",2,])/log(two.loc.eps.2[j,k,l,"or",2,])
        a=a[! is.na(a) & is.finite(a)]; #a[a>=1e4]=1e4; a[a<=-1e4]=-1e4  
        crit=4
        lh.two.loc.eps[j,k,l,crit]=kernel.dens(a,criteria,crit,med.var.crit.int)
      }
    }
  }
  return(lh.two.loc.eps)
}

###lh.two.loc.eps=get.lh(two.loc.eps.eu,two.loc.eps.sa,med.var.crit,criteria)
lh.two.loc.eps=get.lh(two.loc.eps.eu,two.loc.eps.sa,med.var.crit,criteria,drop.eff=0)

# for overdominance
# var.crit.eps =get.var(two.loc.eps.eu,two.loc.eps.sa,citeria,drop.eff=0)
# med.var.crit=sapply(1:4,function(x) median(var.crit.eps[,,,x]))
# mean.var.crit=sapply(1:4,function(x) mean(var.crit.eps[,,,x]))
# lh.two.loc.eps=get.lh(two.loc.eps.eu,two.loc.eps.sa,med.var.crit,criteria,drop.eff=0)


### multiply LHs
lh.two.loc.eps.c12=lh.two.loc.eps[,,,1]*lh.two.loc.eps[,,,2]
lh.two.loc.eps.c34=lh.two.loc.eps[,,,3]*lh.two.loc.eps[,,,4]
lh.two.loc.eps.alls=lh.two.loc.eps[,,,1]*lh.two.loc.eps[,,,2]*lh.two.loc.eps[,,,3]*lh.two.loc.eps[,,,4]


### create plots
library(stats)



#x11(width=10,height=9)
#quartz(width=10,height=9)
library(reshape2,ggplot2)
source("/Volumes/Temp/Lukas/Tools/Scripts/R/image.scale.R")
plot.heatmap.lh <- function(lh.array,doms=c(0,0.25,0.5,0.75,1))
{
  plot.dat=array(lh.array[,3:1,],dim=c(length(hs2[,1]),length(effxeps)), dimnames=list("hs"=rownames(hs2),"effxeps"=effxeps))
  cols=c(colorRampPalette(c("beige","bisque","khaki","yellow","greenyellow","green1","green3","lightblue","blue"))(19))    
  ## cols=rev(heat.colors(18))
  ## cols=rev(terrain.colors(50))
  ## cols=rev(topo.colors(50))
  ## cols=colorRampPalette(brewer.pal(9,"Greens"))(50)
  mybreaks=c(log(10^seq(0,-15,length.out=15)),log(1e-20),log(1e-30),log(1e-50),log(1e-99),0)
  print.breaks=c(0,-1,seq(-3,-24,by=-3))
  print.ends=12
  layout(matrix(1:2, ncol = 2,byrow=T), widths = rep(7,1), respect = FALSE)
  par(mar=c(5,7,5,1))
  image(1:length(effxeps),1:length(rownames(hs2)),log(t(plot.dat)),col=cols, breaks=sort(mybreaks),xaxt='n',yaxt='n',xlab="",ylab="")
  abline(v=(1:length(effxeps)+0.5))
  abline(v=seq(1,length(effxeps)/length(eff2)-1)*length(eff2)+0.5,lw=3)
  abline(h=(1:length(hs2[,1])+0.5))
  abline(h=seq(1,length(hs2[,1])/length(doms)-1)*length(doms)+0.5,lw=3)
  axis(1,at=1:length(effxeps),labels=rep(abs(eff2),length(eps.int)))
  mtext("effect of bab1 relative to tan", side=1, line=2.5,cex=1.5)
  mtext("dominance",side=2, line=5.0,cex=1.25)
  mtext("locus", at=-0.75,side=3, line=1.0,cex=1.25)
  mtext("tan : bab1", at=-0.75,side=3, line=0.0,cex=1.25)
  #mtext(c(0,0.25,0.5,0.75,1),side=2,at=seq(0,length(hs2[,1])/5-1)*5+3,las=1,line=4)
  mtext(rep(doms,length(hs2[,1])/length(doms)),side=2,at=1:length(plot.dat[,1]),line=4,las=1,adj=c(0.5,0.5))
  mtext(":",side=2,at=1:length(plot.dat[,1]),las=1,line=2.75)
  #axis(2,at=1:length(hs2[,1]),labels=rep(c(0,0.25,0.5,0.75,1),5),las=1)
  axis(2,at=1:length(hs2[,1]),labels=NA,las=1)
  mtext(sapply(doms,function(x) rep(x,length(doms))),side=2,at=1:length(hs2[,1]),las=1,line=1.5,adj=c(0.5,0.5))
  #    mtext(sapply(c(1,0.75,0.5,0.25,0),function(x) rep(x,5)),side=2,at=1:length(hs2[,1]),las=1,line=1.5,adj=c(0.5,0.5))
  mtext(eps.int,side=3,at=seq(0,length(eps.int)-1)*length(eff2)+2.5,las=1,line=0.5,cex=1.25)
  mtext("epistatic interaction strength",side=3,at=length(effxeps)/2+0.5,line=2.0,cex=1.5)
  #axis(2, at=seq(0,1,,dim(db)[2]), labels=colnames(db))
  plot.rank=array(length(plot.dat)-rank(plot.dat,na.last=F),dim=dim(plot.dat))
  for(i in 1:length(plot.dat[,1])){
    for(j in 1:length(plot.dat[1,])){
      if (plot.rank[i,j] < 5){
        text(j,i,prettyNum(plot.dat[i,j],digits=2,drop0trailing=T),col="red",cex=0.5,font=2)
      }
      else if (plot.rank[i,j] < 10){  #(plot.dat[i,j] >= 1e-9){
        #if plot.dat[i,j] >= plot.dat
        text(j,i,prettyNum(plot.dat[i,j],digits=2,drop0trailing=T),col="red",cex=0.5,font=1)
      }
      else if (plot.rank[i,j] < 50 & plot.dat[i,j] >= 1e-9){
        text(j,i,prettyNum(plot.dat[i,j],digits=2,drop0trailing=T),col="black",cex=0.5,font=1)
      }
    }
  }
  par(mar=c(5,0,5,5.5))
  image.scale(c(1.01,0.01), col=cols[(length(cols)-print.ends):length(cols)], breaks=sort(mybreaks[1:(print.ends+2)]),horiz=FALSE, yaxt="n",xaxt="n",xlab="")
  #print.breaks=c(0,-1,seq(-3,-24,by=-3))
  axis(4,at=log(10^print.breaks),labels=(format(10^print.breaks,digits=2)), las=2)
  mtext("Likelihood", side=4, line=3.5,cex=1.5)
}


#quartz(width=10,height=9)
eff2=rev(eff2)
eps.int=eps.int[-1]
effxeps=expand.grid(eff2,eps.int)
effxeps=apply(effxeps,1,function (x) paste("Eff:EI",x,collapse=":"))
seff=as.character(eff2)
seps=as.character(eps.int)

quartz(width=10,height=9)
pdf(file="lh_kerneldensity_parameters_d.pdf",width=10,height=9)
plot.heatmap.lh(lh.two.loc.eps[,seff,seps,1]*lh.two.loc.eps[,seff,seps,2]*lh.two.loc.eps[,seff,seps,3]*lh.two.loc.eps[,seff,seps,4],doms=sort(unique(unlist(hs2[1]))))
dev.off()
pdf(file="lh_kerneldensity_crit_eu_parameters_d.pdf",width=10,height=9)
plot.heatmap.lh(lh.two.loc.eps[,seff,seps,1],doms=sort(unique(unlist(hs2[1]))))
dev.off()
pdf(file="lh_kerneldensity_crit_sa_parameters_d.pdf",width=10,height=9)
plot.heatmap.lh(lh.two.loc.eps[,seff,seps,2],doms=sort(unique(unlist(hs2[1]))))
dev.off()
pdf(file="lh_kerneldensity_crit_eusa1_parameters_d.pdf",width=10,height=9)
plot.heatmap.lh(lh.two.loc.eps[,seff,seps,3],doms=sort(unique(unlist(hs2[1]))))
dev.off()
pdf(file="lh_kerneldensity_crit_eusa2_parameters_d.pdf",width=10,height=9)
plot.heatmap.lh(lh.two.loc.eps[,seff,seps,4],doms=sort(unique(unlist(hs2[1]))))
dev.off()
pdf(file="lh_kerneldensity_crit_eusa_parameters_d.pdf",width=10,height=9)
plot.heatmap.lh(lh.two.loc.eps[,seff,seps,3]*lh.two.loc.eps[,seff,seps,4],doms=sort(unique(unlist(hs2[1]))))
dev.off()
pdf(file="lh_kerneldensity_crit12_eu_sa_parameters_d.pdf",width=10,height=9)
plot.heatmap.lh(lh.two.loc.eps[,seff,seps,1]*lh.two.loc.eps[,seff,seps,2],doms=sort(unique(unlist(hs2[1]))))
dev.off()


