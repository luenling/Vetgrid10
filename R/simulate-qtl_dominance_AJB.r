setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Simulation")
### environmental variance of phenotype
v.E = 0.1
afs_cosm=c(0.75,0.8) # cosmopolitan allele frequencies
h2=0.3
#get.V.E <- function(afs_loc,eff_loc,h_loc,h2=0.5,F.in=0,n=1000,eps.int=0.0){
#    gt.freqs = sapply(afs_loc,function (p) c(p^2 + (p-p^2)*F.in,2*(p-p^2)*(1-F.in),(1-p)^2 + (p-p^2)*F.in))
#    pos <- c()
#    pos <- apply(gt.freqs,2,function (x) get.gt.locus(x,c(2,1,0),n))
#    pos <- matrix(pos,nrow=length(afs_loc),ncol=n,byrow=T)
#    ## pigmentation phenotypes resulting from genotypes
#    ## additive effects (eff centered around 0)
#    ind.eff <- (apply(pos,2,function(x) loc.eff(h_loc,x)) - 0.5)
#    ## interaction component
#    intxn <- eps.int*ind.eff[1]*ind.eff[2]*mean(abs(eff_loc))
#    ## additive effects (eff centered around 0)
#    pheno <- eff_loc %*% ind.eff
#    pheno <- pheno + intxn
#    return((1-h2)/h2*var(as.vector(pheno)))
#}
get.V.E <- function(afs_loc,eff_loc,h_loc,h2=0.5,F.in=0,n=1000,eps.int=0.0){
  gt.freqs = sapply(afs_loc,function (p) c(p^2 + (p-p^2)*F.in,2*(p-p^2)*(1-F.in),(1-p)^2 + (p-p^2)*F.in))
  pos <- c()
  pos <- apply(gt.freqs,2,function (x) get.gt.locus(x,c(2,1,0),n))
  pos <- matrix(pos,nrow=length(afs_loc),ncol=n,byrow=T)
  ## pigmentation phenotypes resulting from genotypes
  ## additive effects (eff centered around 0)
  ind.eff <- (apply(pos,2,function(x) loc.eff(h_loc,x)) - 0.5)
  ## interaction component
  #int <- eps.int*ind.eff[1,]*ind.eff[2,]*mean(abs(eff_loc))
  ## additive effects (eff centered around 0)
  pheno <- eff_loc * ind.eff
  #pheno <- pheno + int
  return((1-h2)/h2*var(as.vector(pheno)))
}

library(msm)

### get individual genotype
get.gt <- function(gt.freqs,gt.vals){
    ### gets genotype frequencies and gt values for each locus and returns random genotype for one individual
    return(apply(gt.freqs,2,function (x) colSums(gt.vals*rmultinom(1,1,x))))
}

### get gt for locus
get.gt.locus <- function(gt.freq,gt.vals,total.ind){
    ### gets genotype frequencies and gt values for a locus and returns random genotypes for all ind. for this locus
    return(colSums(gt.vals*rmultinom(total.ind,1,gt.freq)))
}

### locus effects
loc.eff <- function(h.vals,g.type) { # takes the dominance and genotype vectors and returns a vector of gt effects (aa:1,aA:h,AA:0)
    return(g.type/2 + (g.type==1)*(h.vals-0.5))
}


library(foreach)
library(doMC)
registerDoMC(18)
getDoParWorkers()

run_simulations = function(n.runs=50,v.E=0.0,total.ind =2100,n.replicates=3,n.pools=100) {
    
#######  Initialize
    gt.freqs = sapply(n.freqs,function (p) c(p^2 + (p-p^2)*F.inbreed,2*(p-p^2)*(1-F.inbreed),(1-p)^2 + (p-p^2)*F.inbreed)) 
    mean_eff=sum(gt.vals * gt.freqs %*% eff)
    d.replicates=total.ind/n.replicates # ind. per replicate
    results = foreach (i=1:n.runs,.combine="rbind") %dopar% {
        pos <- c()
        pos <- apply(gt.freqs,2,function (x) get.gt.locus(x,gt.vals,total.ind))
        pos <- matrix(pos,nrow=n.loci,ncol=total.ind,byrow=T)  
        pheno <- eff %*% (apply(pos,2,function(x) loc.eff(h.loci,x)) - 0.5)
        ## environmental effect
        pheno <- pheno + rnorm(total.ind,0,sqrt(v.E))
        ## get the randomized array of indices 
        indices=sample(1:total.ind)
         ## get replicates
        for (j in 1:n.replicates) {
            ## indviduals in replicates
            sample.indices = indices[(1+(j-1)*d.replicates):(j*d.replicates)]
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


run_simulations.eps= function(n.runs=50,v.E=0.0,total.ind =2100,n.replicates=3,n.pools=100,eps.int=0.0) {
    
#######  Initialize
    gt.freqs = sapply(n.freqs,function (p) c(p^2 + (p-p^2)*F.inbreed,2*(p-p^2)*(1-F.inbreed),(1-p)^2 + (p-p^2)*F.inbreed)) ### genotype frequencies for each locus, aa, aA, AA
    mean_eff=sum(gt.vals * gt.freqs %*% eff)
    d.replicates=total.ind/n.replicates # ind. per replicate
    ## odds.ratio_ld <- matrix(0,nrow=n.runs,ncol=n.loci)
    ## odds.ratio_l <- matrix(0,nrow=n.runs,ncol=n.loci)
    ## odds.ratio_d <- matrix(0,nrow=n.runs,ncol=n.loci)
                                    #log.odds.ratio_ajb <- matrix(0,nrow=n.runs,ncol=n.loci)
                                        #cors <- matrix(0,nrow=n.runs,ncol=4)
    ## af_avg <- matrix(0,nrow=n.runs,ncol=n.loci)
    ## af_light <- matrix(0,nrow=n.runs,ncol=n.loci)
    ## af_dark <- matrix(0,nrow=n.runs,ncol=n.loci)    

    results = foreach (i=1:n.runs,.combine="rbind") %dopar% {
        ## generate genotypes randomly (independent loci)
                                        # alleles in sample by binomial sampling from population. alleles
                                        # pos <- matrix(rbinom(n.loci*total.ind,2,p),nrow=n.loci,ncol=total.ind)    ### alleles A a
        pos <- c()
                                        #pos <- replicate(total.ind,get.gt(gt.freqs,gt.vals))
                                        #pos <- matrix(pos,nrow=n.loci,ncol=total.ind,byrow=T)
        pos <- apply(gt.freqs,2,function (x) get.gt.locus(x,gt.vals,total.ind))
        pos <- matrix(pos,nrow=n.loci,ncol=total.ind,byrow=T)
        ind.eff <- (apply(pos,2,function(x) loc.eff(h.loci,x)) - 0.5)
        ## pigmentation phenotypes resulting from genotypes
        ## interaction component
        int <- eps.int*ind.eff[1,]*ind.eff[2,]*mean(abs(eff))
        ## additive effects (eff centered around 0)
        pheno <- eff %*% ind.eff
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
            o.pheno <- sample.indices[order(sampled.pheno)]
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

##### simulate diploid individuals
### set total number of ind., number of replicates, and pool sizes
gt.vals <- c(2,1,0)
total.ind <- 700*3 # total individuals
n.replicates <- 3 #number of replicates
d.replicates <- total.ind/n.replicates # ind. per replicate
n.pools <- 100 # size of light and dark pools
n.runs <- 1000
# set loci and heritability
n.loci <- 1
h2 <- 0.25
# set dominance, frequencies, inbreeding and effects
h.loci <- rep(0.5,n.loci)
n.freqs <- rep(0.25,n.loci)
eff <- rep(0.1,n.loci)
F.inbreed <- 0.0
res.fields=c("pv","or")
afs=c(0.1,0.25,0.5)
effs=seq(-1,1,by=0.5)
hs=c(0.0,0.25,0.5,0.75,1.0)
Fs=c(0.0,0.25,0.5)

eff.afs=array(0,dim=c(length(afs),length(effs),length(res.fields),n.runs),dimnames=list("afs"=afs,"effs"=effs,"resfields"=res.fields,"run"=1:n.runs))
# calculate env. variance
v.E=get.V.E(n.freqs,eff,h.loci,h2)
#result=run_simulations(n.runs=50,v.E=v.E)
for(i in 1:length(afs)){
    for(l in 1:length(effs)){
        n.freqs[1]=afs[i]
        eff[1]=effs[l]
        v.E=get.V.E(n.freqs,eff,h.loci,h2)
        run.result=run_simulations(n.runs,v.E=v.E)
        eff.afs[i,l,,]=t(sapply(res.fields, function (x) as.vector(run.result[[x]][,1])))
    }
}

x11(width=10)
x11(width=10,height=6)
layout(matrix(1:2, ncol = 2,byrow=T), widths = rep(1,2), respect = FALSE)
par(oma=c(5.1,5.1,4.1,4.1),mar=c(0,0,0,0))
cols=c("blue","cyan","red","lightgreen","green")
cols=c("cyan","green","red")
shift = (length(afs)+1)/2
bxwd=0.05
for (i in 1:length(afs)){
    if (i == 1){
        boxplot(t(log(eff.afs[i,order(effs),2,])),boxwex=bxwd,at=effs[order(effs)]+bxwd*(i-shift),xlim=c(min(effs)-bxwd*shift,max(effs)+bxwd*shift), xaxt='n',col=cols[i],outline=F,xlab="effect of locus",ylab=expression(log(OR)),add=F)
        axis(1)
        mtext("Odds Ratio",side=3,font=2,line=1.0)
        mtext(bquote("effect of locus"),side=1, line=2.5, outer=T)
        mtext(bquote(log(OR)), side=2, line=2.0)
    }
    else{
        boxplot(t(log(eff.afs[i,order(effs),2,])),boxwex=bxwd,at=effs[order(effs)]+bxwd*(i-shift),xlim=c(min(effs),max(effs)),col=cols[i],outline=F,add=T,xaxt='n')
    }
}
legend("bottomleft",c("10%","25%","50%"),title="MAF",fill=c("cyan","green","red"))
for(i in 1:length(afs)){
    if (i == 1){
        boxplot(t(-log10(eff.afs[i,order(effs),1,])),boxwex=bxwd,at=effs[order(effs)]+bxwd*(i-shift),xlim=c(min(effs)-bxwd*shift,max(effs)+bxwd*shift), ylim=c(0,max(-log10(eff.afs[1:length(afs),order(effs),1,]))), xaxt='n',yaxt='n',col=cols[i],outline=F,xlab="effect of locus",add=F)
        mtext("P-value",side=3,font=2,line=1.0)
        axis(4)
        axis(1)
        mtext(bquote(-log[10](P)), side=4, line=2.0)

    }
    else{
        boxplot(t(-log10(eff.afs[i,order(effs),1,])),boxwex=bxwd,at=effs[order(effs)]+bxwd*(i-shift),xlim=c(min(effs),max(effs)),col=cols[i],outline=F,add=T,xaxt='n',yaxt='n')
    }
}

dev.copy2pdf(file="/Volumes/Temp/Lukas/Data/A7_pigmentation/Simulation/eff_afs.pdf")





eff.afs=array(0,dim=c(length(afs),length(effs),length(res.fields),n.runs),dimnames=list("afs"=afs,"effs"=effs,"resfields"=res.fields,"run"=1:n.runs))
# calculate env. variance
v.E=get.V.E(n.freqs,eff,h.loci,h2)
#result=run_simulations(n.runs=50,v.E=v.E)
for(i in 1:length(afs)){
    for(j in 1:length(hs)){
        for(k in 1:length(Fs)){
            for(l in 1:length(effs)){
                n.freqs=c(afs[i],0.5)
                h.loci=c(hs[j],0.5)
                F.inbreed=Fs[k]
                eff=c(effs[l],0)
                v.E=get.V.E(n.freqs,eff,h.loci,h2)
                run.result=run_simulations(n.runs,v.E=0.75)
            }
        }
    }
}

### number of loci
n.loci <- 20   # number of loci

### population allele frequencies randomly chosen for each locus: either 0.4 or 0.5 or 0.6
fr.vec <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8)
freqs<-runif(n.loci,(length(fr.vec)-1),0.5)
n.freqs<- sample(fr.vec,n.loci,replace=T)
### Alternative: set all pop. allele freq = 0.5
# n.freqs <- 0.5
### Alternative have 2 different loci with different frequencies
# locus 1
fr.vec1<- c(0.5,0.6,0.7)
fr.vec2<- c(0.5,0.6,0.7)
n.freqs<- c(fr.vec1[rbinom(n.loci/2,(length(fr.vec1)-1),0.5)+1],fr.vec2[rbinom(n.loci/2,(length(fr.vec2)-1),0.5)+1])

### set inbreeding
F.inbreed = 0 # inbreeding coefficient
gt.freqs = sapply(n.freqs,function (p) c(p^2 + (p-p^2)*F.inbreed,2*(p-p^2)*(1-F.inbreed),(1-p)^2 + (p-p^2)*F.inbreed)) ### genotype frequencies for each locus, aa, aA, AA
gt.vals = c(2,1,0) ### value for each gt: aa:2, aA: 1, AA: 0

### set dominance: normally distributed around mean
mean.h <- 0.5 # dominance 0: second allele dominant, 0.5: codominant, 1.0: second allele recessive
sd.h <- 0.0 # for sampling dominance from normal dist. (min,max:0,1 )
h.loci <- rtnorm(n.loci,mean.h,sd.h,lower=0.0,upper=1.0)
                                        #h.loci[h.loci <= 0] <- 0.0
                                        #h.loci[h.loci >= 1.0] <- 1.0


#eff     <- 2*runif(n.loci, min=0.4)-1   ### uniformly distributed effect sizes
#eff     <- c(rep(-0.2,10),rep(0.5,10))   ### uniformly distributed effect sizes
n.runs  <- 100  #### 500 simulation runs
eff = runif(n.loci,min=0.5, max =1)*sample(c(1,-1), size=n.loci,replace=T)
### for two loci with different effects
eff1=c(0.0,1.0)
eff2=c(0.0,1.0)
eff = c( runif(n.loci/2,min=eff1[1], max =eff1[2])*sample(c(1,-1), size=n.loci/2,replace=T), runif(n.loci-n.loci/2,min=eff2[1], max =eff2[2])*sample(c(1,-1), size=n.loci-n.loci/2,replace=T))

n.freqs=rep(c(0.5,0.75,0.9),3)

n.loci=20
h.loci=rep(0.5,n.loci)

log_odds_ratio=run_simulations.eps(100,v.E=0.0)
log_odds_ratio_v.E1.0=run_simulations(100,v.E=1.0)
log_odds_ratio_v.E2.0=run_simulations(100,v.E=2.0)
n.runs=100
quartz()
x11()
boxplot(log(log_odds_ratio$or[,order(eff)]),boxwex=0.02,at=eff[order(eff)],xlim=c(-1,1),xaxt='n',col=rgb(0,0,1),outline=F,xlab="effect of locus",ylab=expression(log(OR)),add=F)
axis(1)
boxplot(log(log_odds_ratio_v.E1.0$or[,order(eff)]),boxwex=0.02,at=eff[order(eff)],xlim=c(-1,1),add=T,col=rgb(1,0,0),xaxt='n',yaxt='n',outline=F)
boxplot(log(log_odds_ratio_v.E2.0$or[,order(eff)]),boxwex=0.02,at=eff[order(eff)],xlim=c(-1,1),add=T,col=rgb(0,1,0),xaxt='n',yaxt='n',outline=F)
abline(lm(as.vector(t(log(log_odds_ratio$or)))~rep(eff,n.runs)),lty=2,col="blue")
abline(lm(as.vector(t(log(log_odds_ratio_v.E1.0$or)))~rep(eff,n.runs)),lty=2,col="red")
abline(lm(as.vector(t(log(log_odds_ratio_v.E2.0$or)))~rep(eff,n.runs)),lty=2,col="green")
legend("topright",c(expression(paste(V[E],"= 0.0")),expression(paste(V[E],"= 1.0")),expression(paste(V[E],"= 2.0"))),fill=c("blue","red","green"))
dev.copy2pdf(file="logOR_100x_700_ind_3_r_VE.pdf")
quartz()
boxplot(-log10(log_odds_ratio$pv[,order(eff)]),boxwex=0.02,at=eff[order(eff)],xlim=c(-1,1),xaxt='n',col=rgb(0,0,1),outline=F,xlab="effect of locus",ylab=expression(-log[10](P)),add=F)
axis(1)
boxplot(-log10(log_odds_ratio_v.E1.0$pv[,order(eff)]),boxwex=0.02,at=eff[order(eff)],xlim=c(-1,1),add=T,col=rgb(1,0,0),xaxt='n',yaxt='n',outline=F)
boxplot(-log10(log_odds_ratio_v.E2.0$pv[,order(eff)]),boxwex=0.02,at=eff[order(eff)],xlim=c(-1,1),add=T,col=rgb(0,1,0),xaxt='n',yaxt='n',outline=F)
## abline(lm(as.vector(t(-log10(log_odds_ratio$pv)))~rep(eff,n.runs)),lty=2,col="blue")
## abline(lm(as.vector(t(-log10(log_odds_ratio_v.E1.0$pv)))~rep(eff,n.runs)),lty=2,col="red")
## abline(lm(as.vector(t(-log10(log_odds_ratio_v.E2.0$pv)))~rep(eff,n.runs)),lty=2,col="green")
legend("topright",c(expression(paste(V[E],"= 0.0")),expression(paste(V[E],"= 1.0")),expression(paste(V[E],"= 2.0"))),fill=c("blue","red","green"))
dev.copy2pdf(file="logPV_100x_700_ind_3_r_VE.pdf")


n.freqs=rep(c(0.5,0.75,0.9),3)
n.loci=9
eff=c(0.25,0.25,0.25,0.5,0.5,0.5,1.0,1.0,1.0)
F.inbreed=0
h.loci=rep(0.5,9)

nine.loc=run_simulations(100,v.E=1.0)


#plot(colMeans(log(a$or)),pch=0,col="green",ylim=c(-2,2))


boxplot(log(a$ol),boxwex=0.2,xaxt='n',col=rgb(0,0,1),outline=F,xlab="effect of locus",ylab=expression(log(OR)),ylim=c(-2,2))
axis(1)
boxplot(log(a$od),boxwex=0.2,add=T,col=rgb(1,0,0),xaxt='n',yaxt='n',outline=F)
boxplot(log(a$ol)+log(a$od),boxwex=0.2,add=T,col=rgb(0,1,0),xaxt='n',yaxt='n',outline=F)

plot(colMeans(log(a$ol)),pch=1,col="red",ylim=c(-3,3))
points(colMeans(log(a$od)),pch=2,col="blue")
points(colMeans(log(a$ol))+colMeans(log(a$od)),pch=3,col="black")
abline(0,0,lty=2)
plot(colMeans(a$af.m),pch=0,col="green",ylim=c(0,1))
points(colMeans(a$af.l),pch=1,col="red")
points(colMeans(a$af.d),pch=2,col="blue")

#######  Initialize
odds.ratio_ajb <- matrix(0,nrow=n.runs,ncol=n.loci)
log.odds.ratio_ajb <- matrix(0,nrow=n.runs,ncol=n.loci)
cors       <- matrix(0,nrow=n.runs,ncol=4)

total.ind =700*3 # total individuals
n.replicates = 3 #number of replicates
d.replicates = total.ind/n.replicates # ind. per replicate
n.pools <- 100 # size of light and dark pools

#### experiments to run

### one locus, big variation, different afs, F and h
# should see AF changes depending on afs, F and h
afs=c(0.5,0.75,0.9)
eff=c(0.5,0)
hs=c(0.0,0.25,0.5,0.75,1.0)
Fs=c(0.0,0.25,0.5)
n.loci=2
n.runs=100
res.fields=c("pv","or","ol","od","af.m","af.l","af.d")
one.loc=array(0,dim=c(length(afs),length(hs),length(Fs),length(res.fields),n.runs),dimnames=list("afs"=afs,"hs"=hs,"Fs"=Fs,"resfields"=res.fields,"run"=1:n.runs))
##one.loc=array(0,dim=c(length(afs),length(hs),length(Fs)),dimnames=list(afs,hs,Fs))
for(i in 1:length(afs)){
    for(j in 1:length(hs)){
        for(k in 1:length(Fs)){
            v.E=0.5
            n.loci=2
            n.freqs=c(afs[i],0.5)
            h.loci=c(hs[j],0.5)
            F.inbreed=Fs[k]
            run.result=run_simulations(n.runs,v.E=0.75)
            #one.loc[i,j,k,]=sapply(res.fields,function (x) mean(run.result[[x]][,1]))
            one.loc[i,j,k,,]=t(sapply(res.fields, function (x) as.vector(run.result[[x]][,1])))
        }
    }
}

### plot one locus
quartz(width=12,height=7)
layout(matrix(1:3, ncol = 3), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
axis.t=c("t","n","n")
hs=c(0.0,0.25,0.5,0.75,1.0)
for (i in 1:length(Fs)){
    d.s=0.2
    ylims=c(0,-3)
    xlims=c(1-d.s,length(hs)+d.s)

    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    for(j in 1:3){
        boxplot(log(t(one.loc[j,,i,"or",])), at=1:length(hs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",cex.axis=1.5,xaxt=axis.t[j],add=bool.add[j])
    }
    mtext(bquote(F[is]~"="~.(Fs[i])), side=3, line=1.0)
    ## boxplot(log(t(one.loc["0.75",,i,"or",]))), at=hs,boxwex=0.03, col="green",ylim=c(0,6),xlim=c(-0.2,1.2), yaxt="n",xaxt=axis.t[j],add=T)
    ## boxplot(log(t(one.loc["0.9",,i,"or",]))), at=hs+d.s,boxwex=0.03, col="deepskyblue",ylim=c(0,6),xlim=c(-0.2,1.2), yaxt="n",xaxt="n",add=T)
    if (i==1){
        axis(2,cex.axis=1.5)
    }
}
    mtext("h",side=1,outer=T,line=2.5,cex=1.25)
legend("bottomright",c("0.5","0.75","0.9"),title=c("avg. AFs"),fill=c("red","green","deepskyblue"),cex=1.5)
mtext(expression(log(OR)),side=2,outer=T,line=2.0,cex=1.25)
dev.copy2pdf(file="logOR_100x700_3r_1locus_h_F.pdf")

quartz(width=12,height=7)
layout(matrix(1:3, ncol = 3), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
axis.t=c("t","n","n")
hs=c(0.0,0.25,0.5,0.75,1.0)
for (i in 1:length(Fs)){
    d.s=0.2
    ylims=c(0,50)
    xlims=c(1-d.s,length(hs)+d.s)

    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    for(j in 1:3){
        boxplot(-log10(t(one.loc[j,,i,"pv",])), at=1:length(hs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",cex.axis=1.5,xaxt=axis.t[j],add=bool.add[j])
    }
    mtext(bquote(F[is]~"="~.(Fs[i])), side=3, line=1.0)
    ## boxplot(log(t(one.loc["0.75",,i,"or",]))), at=hs,boxwex=0.03, col="green",ylim=c(0,6),xlim=c(-0.2,1.2), yaxt="n",xaxt=axis.t[j],add=T)
    ## boxplot(log(t(one.loc["0.9",,i,"or",]))), at=hs+d.s,boxwex=0.03, col="deepskyblue",ylim=c(0,6),xlim=c(-0.2,1.2), yaxt="n",xaxt="n",add=T)
    if (i==1){
        axis(2,cex.axis=1.5)
    }
}
    mtext("h",side=1,outer=T,line=2.5,cex=1.25)
## legend("bottomright",c("0.5","0.75","0.9"),title=c("avg. AFs"),fill=c("red","green","deepskyblue"),cex=1.5)
mtext(expression(-log[10](P)),side=2,outer=T,line=2.0,cex=1.25)
dev.copy2pdf(file="logPV_100x700_3r_1locus_h_F.pdf")

quartz(width=12,height=7)
layout(matrix(1:3, ncol = 3), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
axis.t=c("t","n","n")
hs=c(0.0,0.25,0.5,0.75,1.0)
for (i in 1:length(Fs)){
    d.s=0.2
    ylims=c(0.3,1.0)
    xlims=c(1-d.s,length(hs)+d.s)

    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    for(j in 1:3){
        boxplot(t(one.loc[j,,i,"af.m",]), at=1:length(hs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",cex.axis=1.5,xaxt=axis.t[j],add=bool.add[j])
         boxplot(t(one.loc[j,,i,"af.l",]), at=1:length(hs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",cex.axis=1.5,xaxt="n",add=T)
         boxplot(t(one.loc[j,,i,"af.d",]), at=1:length(hs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",cex.axis=1.5,xaxt="n",add=T)
   }
    mtext(bquote(F[is]~"="~.(Fs[i])), side=3, line=1.0)
    ## boxplot(log(t(one.loc["0.75",,i,"or",]))), at=hs,boxwex=0.03, col="green",ylim=c(0,6),xlim=c(-0.2,1.2), yaxt="n",xaxt=axis.t[j],add=T)
    ## boxplot(log(t(one.loc["0.9",,i,"or",]))), at=hs+d.s,boxwex=0.03, col="deepskyblue",ylim=c(0,6),xlim=c(-0.2,1.2), yaxt="n",xaxt="n",add=T)
    if (i==1){
        axis(2,cex.axis=1.5)
    }
}
    mtext("h",side=1,outer=T,line=2.5,cex=1.25)
## legend("bottomright",c("0.5","0.75","0.9"),title=c("avg. AFs"),fill=c("red","green","deepskyblue"),cex=1.5)
mtext(expression(AFS),side=2,outer=T,line=2.0,cex=1.25)
dev.copy2pdf(file="afs_100x700_3r_1locus_h_F.pdf")


### two loci, different effect
## one locus less of an effect

### set total number of ind., number of replicates, and pools sizes
### SA

total.ind =700*3 # total individuals
n.replicates = 3 #number of replicates
d.replicates = total.ind/n.replicates # ind. per replicate
n.pools <- 65 # size of light and dark pools

afs2=c(0.5,0.75,0.8,0.9)
afs2=expand.grid(afs2,afs2)
rownames(afs2)=apply(afs2,1,function (x) paste(x,collapse=":"))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
afs2b=afs2[afs_val,]
eff1=c(1.0)
eff1=1.0
eff2=c(-5,-3,-1,1,3,5)
eff2b=c(-3,-2,-1)
hs2=c(0.0,0.5,1.0)
hs2=expand.grid(hs2,hs2)
rownames(hs2)=apply(hs2,1,function (x) paste(x,collapse=":"))
H_levs=c("0.5:0","0.5:0.5","0.5:1","0:1","1:1")
hs2b=hs2[H_levs,]
Fs2=c(0.0,0.25,0.5,0.75)
Fs2b=c(0.0,0.25,0.5)
n.loci2=2
n.runs2=100
eff2.eu=c(-5,-3,-1)
eps.int=c(-1,0,1)
eps.eff=c(-1)
Fs2.eu=c(0)
res.fields2=c("pv","or","ol","od","af.m","af.l","af.d")
two.loc.b=array(0,dim=c(length(afs2b[,1]),length(hs2b[,1]),length(Fs2b),length(eff2b),length(res.fields2),2,n.runs2),dimnames=list("afs"=rownames(afs2b),"hs"=rownames(hs2b),"Fs"=Fs2b,"rel.eff"=eff2b,"resfields"=res.fields2,"locus"=c(1,2),"run"=1:n.runs2))
h2=0.15

for(j in 1:length(hs2b[,1])){
    for(k in 1:length(Fs2b)){
        for(l in 1:length(eff2b)){
            eff=c(eff1/(1+abs(eff2b[l])),eff1*eff2b[l]/(1+abs(eff2b[l])))
            n.loci=n.loci2
            h.loci=unlist(hs2b[j,])
            F.inbreed=Fs2b[k]
            v.E=get.V.E(afs_cosm,eff,h.loci,h2)
            for(i in 1:length(afs2b[,1])){
                n.freqs=unlist(afs2b[i,])
                run.result=run_simulations(n.runs2,v.E,total.ind =2100,n.replicates=3,n.pools=65)
                two.loc.b[i,j,k,l,,1,]=t(sapply(res.fields2, function (x) run.result[[x]][,1]))
                two.loc.b[i,j,k,l,,2,]=t(sapply(res.fields2, function (x) run.result[[x]][,2]))
            }
        }
    }
}
#str(two.loc[i,j,k,l,,1,])
eff2.eub=c(-3,-2,-1)
two.loc.eu.b=array(0,dim=c(length(afs2b[,1]),length(hs2b[,1]),length(Fs2.eu),length(eff2.eub),length(res.fields2),2,n.runs2),dimnames=list("afs"=rownames(afs2b),"hs"=rownames(hs2b),"Fs"=Fs2.eu,"rel.eff"=eff2.eub,"resfields"=res.fields2,"locus"=c(1,2),"run"=1:n.runs2))


for(j in 1:length(hs2b[,1])){
    for(k in 1:length(Fs2.eu)){
        for(l in 1:length(eff2.eub)){
            eff=c(eff1/(1+abs(eff2b[l])),eff1*eff2b[l]/(1+abs(eff2b[l])))
            n.loci=n.loci2
            h.loci=unlist(hs2b[j,])
            F.inbreed=Fs2[k]
            v.E=get.V.E(afs_cosm,eff,h.loci,h2)
            #v.E=1.0
            for(i in 1:length(afs2b[,1])){
                n.freqs=unlist(afs2b[i,])
                run.result=run_simulations(n.runs2,v.E,total.ind =9000,n.replicates=6,n.pools=100)
                two.loc.eu.b[i,j,k,l,,1,]=t(sapply(res.fields2, function (x) run.result[[x]][,1]))
                two.loc.eu.b[i,j,k,l,,2,]=t(sapply(res.fields2, function (x) run.result[[x]][,2]))
            }
        }
    }
}

### epistatic
two.loc.eps.b=array(0,dim=c(length(afs2b[,1]),length(hs2b[,1]),length(Fs2b),length(eps.int),length(res.fields2),2,n.runs2),dimnames=list("afs"=rownames(afs2b),"hs"=rownames(hs2b),"Fs"=Fs2b,"eps.int"=eps.int,"resfields"=res.fields2,"locus"=c(1,2),"run"=1:n.runs2))
h2=0.15

for(j in 1:length(hs2b[,1])){
    for(k in 1:length(Fs2b)){
        for(l in 1:length(eps.int)){
            eff=c(0.5,-0.5)
            n.loci=n.loci2
            h.loci=unlist(hs2b[j,])
            F.inbreed=Fs2b[k]
            v.E=get.V.E(afs_cosm,eff,h.loci,h2)
            for(i in 1:length(afs2b[,1])){
                n.freqs=unlist(afs2b[i,])
                run.result=run_simulations.eps(n.runs2,v.E,total.ind =2100,n.replicates=3,n.pools=65,eps.int=eps.int[l])
                two.loc.eps.b[i,j,k,l,,1,]=t(sapply(res.fields2, function (x) run.result[[x]][,1]))
                two.loc.eps.b[i,j,k,l,,2,]=t(sapply(res.fields2, function (x) run.result[[x]][,2]))
            }
        }
    }
}
#str(two.loc[i,j,k,l,,1,])
eff2.eub=c(-3,-2,-1)
two.loc.eu.eps.b=array(0,dim=c(length(afs2b[,1]),length(hs2b[,1]),length(Fs2.eu),length(eps.int),length(res.fields2),2,n.runs2),dimnames=list("afs"=rownames(afs2b),"hs"=rownames(hs2b),"Fs"=Fs2.eu,"eps.int"=eps.int,"resfields"=res.fields2,"locus"=c(1,2),"run"=1:n.runs2))


for(j in 1:length(hs2b[,1])){
    for(k in 1:length(Fs2.eu)){
        for(l in 1:length(eps.int)){
            eff=c(0.5,-0.5)
            n.loci=n.loci2
            h.loci=unlist(hs2b[j,])
            F.inbreed=Fs2[k]
            v.E=get.V.E(afs_cosm,eff,h.loci,h2)
            #v.E=1.0
            for(i in 1:length(afs2b[,1])){
                n.freqs=unlist(afs2b[i,])
                run.result=run_simulations.eps(n.runs2,v.E,total.ind =9000,n.replicates=6,n.pools=100,eps.int=eps.int[l])
                two.loc.eu.eps.b[i,j,k,l,,1,]=t(sapply(res.fields2, function (x) run.result[[x]][,1]))
                two.loc.eu.eps.b[i,j,k,l,,2,]=t(sapply(res.fields2, function (x) run.result[[x]][,2]))
            }
        }
    }
}
### plot two locus effects
quartz(width=12,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1","0:1","1:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val=1
for (i in as.character(eff2[1:3])){
    d.s=0.2
    ylims=c(0,-4.5)
    xlims=c(1-d.s,length(H_levs)+d.s)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(log(t(two.loc[afs_val[1],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-d.s,boxwex=d.s, col="red",ylim=ylims,xlim=xlims, yaxt="n",xaxt="n")
    boxplot(log(t(two.loc[afs_val[2],H_levs,F.val,i,"or",1,])), at=1:length(H_levs),boxwex=d.s, col="green",ylim=c(0,-3),xlim=xlims, yaxt="n",cex.axis=1.5,xaxt="n",add=T)
    #mtext("h",side=1,outer=F,line=2.5)
    mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)
    boxplot(log(t(two.loc[afs_val[3],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)+d.s,boxwex=d.s, col="deepskyblue",ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    if (i=="-5"){
        axis(2,cex.axis=1.5)
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[1]~ "="~.(afs2[x,1])~","~loc[2]~"="~.(afs2[x,2]))))
legend("bottomleft",legend=leg_expr,title=c("avg. AFs"),fill=c("red","green","deepskyblue"),cex=1.5)
for (i in as.character(eff2[1:3])){
    d.s=0.2
    xlims=c(1-d.s,length(H_levs)+d.s)
    ylims=c(0,8)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(abs(log(t(two.loc[afs_val[1],H_levs,F.val,i,"or",2,]))), at=1:length(H_levs)-d.s,boxwex=d.s, col="red",ylim=ylims,xlim=xlims, yaxt="n",xaxt="n")
    axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    boxplot(abs(log(t(two.loc[afs_val[2],H_levs,F.val,i,"or",2,]))), at=1:length(H_levs),boxwex=d.s, col="green",ylim=ylims,xlim=xlims, yaxt="n",cex.axis=1.5,xaxt="n",add=T)
    mtext("h",side=1,outer=F,line=2.5)
    #mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)
    boxplot(abs(log(t(two.loc[afs_val[3],H_levs,F.val,i,"or",2,]))), at=1:length(H_levs)+d.s,boxwex=d.s, col="deepskyblue",ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    if (i=="-5"){
        axis(2,cex.axis=1.5)
    }
}
mtext(expression(log(OR)),side=2,outer=T,line=2.0)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0)
dev.copy2pdf(file="logOR_100x700_3r_2loc_h_afs_loc1_0.75_neg_eff_h2_0.3.pdf")
### plot two locus effects restricted (b values)
quartz(width=12,height=7)

layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1","0:1","1:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val=2
for (i in as.character(eff2b[1:3])){
    d.s=0.2
    ylims=c(0,-3.75)
    xlims=c(1-d.s,length(H_levs)+d.s)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(log(t(two.loc.b[afs_val[1],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-d.s,boxwex=d.s, col="red",ylim=ylims,xlim=xlims, yaxt="n",xaxt="n")
    boxplot(log(t(two.loc.b[afs_val[2],H_levs,F.val,i,"or",1,])), at=1:length(H_levs),boxwex=d.s, col="green",ylim=c(0,-3),xlim=xlims, yaxt="n",cex.axis=1.5,xaxt="n",add=T)
    #mtext("h",side=1,outer=F,line=2.5)
    mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)
    boxplot(log(t(two.loc.b[afs_val[3],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)+d.s,boxwex=d.s, col="deepskyblue",ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    if (i=="-3"){
        axis(2,cex.axis=1.5)
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[1]~ "="~.(afs2[x,1])~","~loc[2]~"="~.(afs2[x,2]))))
legend("bottomleft",legend=leg_expr,title=c("avg. AFs"),fill=c("red","green","deepskyblue"),cex=1.5)
for (i in as.character(eff2b[1:3])){
    d.s=0.2
    xlims=c(1-d.s,length(H_levs)+d.s)
    ylims=c(0,6.5)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(abs(log(t(two.loc.b[afs_val[1],H_levs,F.val,i,"or",2,]))), at=1:length(H_levs)-d.s,boxwex=d.s, col="red",ylim=ylims,xlim=xlims, yaxt="n",xaxt="n")
    axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    boxplot(abs(log(t(two.loc.b[afs_val[2],H_levs,F.val,i,"or",2,]))), at=1:length(H_levs),boxwex=d.s, col="green",ylim=ylims,xlim=xlims, yaxt="n",cex.axis=1.5,xaxt="n",add=T)
    mtext("h",side=1,outer=F,line=2.5)
    #mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)
    boxplot(abs(log(t(two.loc.b[afs_val[3],H_levs,F.val,i,"or",2,]))), at=1:length(H_levs)+d.s,boxwex=d.s, col="deepskyblue",ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    if (i=="-3"){
        axis(2,cex.axis=1.5)
    }
}
mtext(expression(log(OR)),side=2,outer=T,line=2.0)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0)
dev.copy2pdf(file="logOR_65x700_3r_2loc_h_afs_loc1_0.8_neg_eff_0.5_h2_0.2_F_0.25.pdf")


mean(log(two.loc[afs_val[3],"0.5:1","0.25","-5","or",1,]))
#[1] -1.276026
mean(log(two.loc[afs_val[1],"0.5:1","0.25","-5","or",1,]))
#[1] -0.8224399
#> >  -1.049521/-0.8224399
#[1] 1.276107
mean(log(two.loc[afs_val[1],"0.5:1","0","-1","or",1,]))
#[1] -0.6581223
mean(log(two.loc[afs_val[3],"0.5:1","0","-1","or",1,]))
#[1] -0.7878388
mean(log(two.loc[afs_val[3],"0.5:1","0","-1","or",1,]))

mean(log(two.loc.eu[afs_val[3],"0.5:1","0","-1","or",1,]))/mean(log(two.loc[afs_val[1],"0.5:1","0","-1","or",1,]))
mean(log(two.loc.eu[afs_val[3],"0.5:1","0","-1","or",1,]))/mean(log(two.loc[afs_val[1],"0.5:1","0.25","-1","or",1,]))
mean(log(two.loc.eu[afs_val[3],"0:1","0","-1","or",1,]))/mean(log(two.loc[afs_val[1],"0:1","0.25","-1","or",1,]))
## SA
loc1.SA=median(log(two.loc.b[afs_val[1],"0.5:1","0.25","-3","or",1,]))
loc2.SA=median(log(two.loc.b[afs_val[1],"0.5:1","0.25","-3","or",2,]))
## EU 
loc1.EU=median(log(two.loc.eu.b[afs_val[3],"0.5:1","0","-3","or",1,]))
loc2.EU=median(log(two.loc.eu.b[afs_val[3],"0.5:1","0","-3","or",2,]))

loc1.EU/loc2.EU
loc1.SA/loc2.SA

### plot two locus effects eu
quartz(width=12,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1","0:1","1:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val="0"
for (i in as.character(eff2[1:3])){
    d.s=0.2
    ylims=c(0,-4.5)
    xlims=c(1-d.s,length(H_levs)+d.s)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(log(t(two.loc.eu[afs_val[1],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-d.s,boxwex=d.s, col="red",ylim=ylims,xlim=xlims, yaxt="n",xaxt="n")
    boxplot(log(t(two.loc.eu[afs_val[2],H_levs,F.val,i,"or",1,])), at=1:length(H_levs),boxwex=d.s, col="green",ylim=ylims,xlim=xlims, yaxt="n",cex.axis=1.5,xaxt="n",add=T)
    #mtext("h",side=1,outer=F,line=2.5)
    mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)
    boxplot(log(t(two.loc.eu[afs_val[3],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)+d.s,boxwex=d.s, col="deepskyblue",ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    if (i=="-5"){
        axis(2,cex.axis=1.5)
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[1]~ "="~.(afs2[x,1])~","~loc[2]~"="~.(afs2[x,2]))))
legend("bottomleft",legend=leg_expr,title=c("avg. AFs"),fill=c("red","green","deepskyblue"),cex=1.5)
for (i in as.character(eff2[1:3])){
    d.s=0.2
    ylims=c(0,8)
    xlims=c(1-d.s,length(H_levs)+d.s)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(abs(log(t(two.loc.eu[afs_val[1],H_levs,F.val,i,"or",2,]))), at=1:length(H_levs)-d.s,boxwex=d.s, col="red",ylim=ylims,xlim=xlims, yaxt="n",xaxt="n")
    axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    boxplot(abs(log(t(two.loc.eu[afs_val[2],H_levs,F.val,i,"or",2,]))), at=1:length(H_levs),boxwex=d.s, col="green",ylim=ylims,xlim=xlims, yaxt="n",cex.axis=1.5,xaxt="n",add=T)
    mtext("h",side=1,outer=F,line=2.5)
    #mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)
    boxplot(abs(log(t(two.loc.eu[afs_val[3],H_levs,F.val,i,"or",2,]))), at=1:length(H_levs)+d.s,boxwex=d.s, col="deepskyblue",ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    if (i=="-5"){
        axis(2,cex.axis=1.5)
    }
}
mtext(expression(log(OR)),side=2,outer=T,line=2.0)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0)
dev.copy2pdf(file="logOR_100x1500_3r_2loc_h_afs_loc1_0.8_neg_eff_h2_0.3.pdf")

### plot two locus effects eu (b values)
quartz(width=12,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1","0:1","1:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val="0"
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
d.s=0.2
for (i in as.character(eff2b[1:3])){
    ylims=c(0,-3.75)
    xlims=c(1-d.s,length(H_levs)+d.s)
    for(j in 1:3){
        boxplot(log(t(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
    }
    mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)

    if (i=="-3"){
        axis(2,cex.axis=1.5)
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[1]~ "="~.(afs2b[x,1])~","~loc[2]~"="~.(afs2b[x,2]))))
legend("bottomleft",legend=leg_expr,title=c("avg. AFs"),fill=colors,cex=1.5)
for (i in as.character(eff2b[1:3])){
    ylims=c(0,6.5)
    xlims=c(1-d.s,length(H_levs)+d.s)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    for(j in 1:3){
        boxplot(log(t(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    mtext("h",side=1,outer=F,line=2.5)
    if (i=="-3"){
        axis(2,cex.axis=1.5)
    }
}
mtext(expression(log(OR)),side=2,outer=T,line=2.0)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0)
dev.copy2pdf(file="logOR_100x1500_3r_2loc_h_afs_loc1_0.8_neg_eff_0.5_h2_0.2.pdf")

### plot two locus effects eu eps (b values)
quartz(width=12,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1","0:1","1:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val="0"
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
d.s=0.2
for (i in as.character(eps.int)){
    ylims=c(0,-3.75)
    xlims=c(1-d.s,length(H_levs)+d.s)
    for(j in 1:3){
        boxplot(log(t(two.loc.eu.eps.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
    }
    mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)

    if (i==eps.int[1]){
        axis(2,cex.axis=1.5)
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[1]~ "="~.(afs2b[x,1])~","~loc[2]~"="~.(afs2b[x,2]))))
legend("bottomleft",legend=leg_expr,title=c("avg. AFs"),fill=colors,cex=1.5)
for (i in as.character(eps.int)){
    ylims=c(0,6.5)
    xlims=c(1-d.s,length(H_levs)+d.s)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    for(j in 1:3){
        boxplot(log(t(two.loc.eu.eps.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    mtext("h",side=1,outer=F,line=2.5)
    if (i==eps.int[1]){
        axis(2,cex.axis=1.5)
    }
}
mtext(expression(log(OR)),side=2,outer=T,line=2.0)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0)

### plot two locus PV eu (b values)
quartz(width=12,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1","0:1","1:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val="0"
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
d.s=0.2
for (i in as.character(eff2b[1:3])){
    ylims=c(0,90)
    xlims=c(1-d.s,length(H_levs)+d.s)
    for(j in 1:3){
        boxplot(-log10(t(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
    }
    mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)

    if (i=="-3"){
        axis(2,cex.axis=1.5)
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[1]~ "="~.(afs2b[x,1])~","~loc[2]~"="~.(afs2b[x,2]))))
legend("bottomleft",legend=leg_expr,title=c("avg. AFs"),fill=colors,cex=1.5)
for (i in as.character(eff2b[1:3])){
    ylims=c(0,90)
    xlims=c(1-d.s,length(H_levs)+d.s)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    for(j in 1:3){
        boxplot(-log10(t(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    mtext("h",side=1,outer=F,line=2.5)
    if (i=="-3"){
        axis(2,cex.axis=1.5)
    }
}
mtext(expression(-log[10](P)),side=2,outer=T,line=2.0)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0)

### plot two locus frequencies sa (b values)
quartz(width=12,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1","0:1","1:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val="0.25"
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
d.s=0.2
for (i in as.character(eff2b[1:3])){
    ylims=c(0,1.0)
    xlims=c(1-d.s,length(H_levs)+d.s)
    for(j in 1:3){
        boxplot(t(two.loc.b[afs_val[j],H_levs,F.val,i,"af.m",1,]), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(two.loc.b[afs_val[j],H_levs,F.val,i,"af.l",1,]), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(two.loc.b[afs_val[j],H_levs,F.val,i,"af.d",1,]), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
    mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)

    if (i=="-3"){
        axis(2,cex.axis=1.5)
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[1]~ "="~.(afs2b[x,1])~","~loc[2]~"="~.(afs2b[x,2]))))
legend("bottomleft",legend=leg_expr,title=c("avg. AFs"),fill=colors,cex=1.5)
for (i in as.character(eff2b[1:3])){
    ylims=c(0,1)
    xlims=c(1-d.s,length(H_levs)+d.s)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
   for(j in 1:3){
        boxplot(t(two.loc.b[afs_val[j],H_levs,F.val,i,"af.m",2,]), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(two.loc.b[afs_val[j],H_levs,F.val,i,"af.l",2,]), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(two.loc.b[afs_val[j],H_levs,F.val,i,"af.d",2,]), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    mtext("h",side=1,outer=F,line=2.5)
    if (i=="-3"){
        axis(2,cex.axis=1.5)
    }
}
mtext(expression(AFs),side=2,outer=T,line=2.0)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0)


### plot two locus frequencies eu (b values)
quartz(width=12,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1","0:1","1:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val="0"
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
d.s=0.2
for (i in as.character(eff2b[1:3])){
    ylims=c(0,1.0)
    xlims=c(1-d.s,length(H_levs)+d.s)
    for(j in 1:3){
        boxplot(t(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"af.m",1,]), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"af.l",1,]), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"af.d",1,]), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
    mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)

    if (i=="-3"){
        axis(2,cex.axis=1.5)
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[1]~ "="~.(afs2b[x,1])~","~loc[2]~"="~.(afs2b[x,2]))))
legend("bottomleft",legend=leg_expr,title=c("avg. AFs"),fill=colors,cex=1.5)
for (i in as.character(eff2b[1:3])){
    ylims=c(0,1)
    xlims=c(1-d.s,length(H_levs)+d.s)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
   for(j in 1:3){
        boxplot(t(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"af.m",2,]), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"af.l",2,]), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"af.d",2,]), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    mtext("h",side=1,outer=F,line=2.5)
    if (i=="-3"){
        axis(2,cex.axis=1.5)
    }
}
mtext(expression(AFs),side=2,outer=T,line=2.0)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0)
dev.copy2pdf(file="AFs_100x1500_3r_2loc_h_afs_loc1_0.8_neg_eff_0.5_h2_0.2.pdf")

### plot two locus lOR eu dominance (b values)
quartz(width=7,height=7)
layout(matrix(1:2, ncol = 1,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1","0:1","1:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val="0.25"
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
d.s=0.2
for (i in as.character(eff2b[3])){
    ylims=c(0,-3.0)
    xlims=c(1-d.s,length(H_levs)+d.s)
    for(j in 1:3){
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
    #mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)

    #if (i=="-3"){
        axis(2,cex.axis=1.5)
    #}
}
for (i in as.character(eff2b[3])){
    ylims=c(0,3.5)
    xlims=c(1-d.s,length(H_levs)+d.s)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
   for(j in 1:3){
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    mtext("h",side=1,outer=F,line=2.5,cex=1.5)
    #if (i=="-3"){
        axis(2,cex.axis=1.5)
    #}
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc2~"="~.(afs2b[x,2]))))
legend("bottomleft",legend=c(as.expression(bquote(loc1~"="~.(afs2b[1,1]))),leg_expr),fill=c(NA,colors),cex=1)
mtext(expression(log(OR)),side=2,outer=T,line=2.25,cex=1.5)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0,,cex=1.5)
dev.copy2pdf(file="lOR_dominance_sa_65x700_3r_2loc_h_afs_loc1_0.8_neg_eff_1_h2_0.15_F_0.25.pdf")
## pv sa dominance
quartz(width=7,height=7)
layout(matrix(1:2, ncol = 1,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
F.val="0.25"
for (i in as.character(eff2b[3])){
    ylims=c(0,69)
    xlims=c(1-d.s,length(H_levs)+d.s)
    for(j in 1:3){
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
    #mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)

    #if (i=="-3"){
        axis(2,cex.axis=1.5)
    #}
}
## leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc2~"="~.(afs2b[x,2]))))
## legend("topleft",legend=c(as.expression(bquote(loc1~"="~.(afs2b[1,1]))),leg_expr),fill=c(NA,colors),cex=1)

for (i in as.character(eff2b[3])){
    ylims=c(0,69)
    xlims=c(1-d.s,length(H_levs)+d.s)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
   for(j in 1:3){
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    mtext("h",side=1,outer=F,line=2.5,cex=1.5)
    #if (i=="-3"){
        axis(2,cex.axis=1.5)
    #}
}
mtext(expression(-log[10](P)),side=2,outer=T,line=2.25,cex=1.5)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0,,cex=1.5)
dev.copy2pdf(file="lpv_dominance_sa_65x700_3r_2loc_h_afs_loc1_0.8_neg_eff_1_h2_0.15_F_0.25.pdf")

### plot two locus lOR sa inbreeding (b values)
quartz(width=10,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.vals=c("0","0.25","0.5") 
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
d.s=0.15
for (i in F.vals){
    ylims=c(0,-3.5)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)
    for(j in 1:3){
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,i,3,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,i,3,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,i,3,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
    mtext(bquote(F[is]~"="~.(i)), side=3, line=1.0,cex=1.25)

    if (i==F.vals[1]){
        axis(2,cex.axis=1.5,at=c(0,-1,-2,-3))
    }
}
for (i in F.vals){
    ylims=c(0,3.5)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
   for(j in 1:3){
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,i,3,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,i,3,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,i,3,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    if (i==F.vals[1]){
        axis(2,cex.axis=1.5,at=c(0,1,2,3))
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc2~"="~.(afs2b[x,2]))))
legend("bottomleft",legend=c(as.expression(bquote(loc1~"="~.(afs2b[1,1]))),leg_expr),fill=c(NA,colors),cex=1)
mtext("h",side=1,outer=T,line=2.5,cex=1.25)
mtext(expression(log(OR)),side=2,outer=T,line=2.25,cex=1.25)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0,,cex=1.25)
dev.copy2pdf(file="lOR_F_sa_65x700_3r_2loc_h_afs_loc1_0.8_neg_eff_1_h2_0.15.pdf")
## pv sa inbreeding
quartz(width=10,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.vals=c("0","0.25","0.5") 
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
d.s=0.15
for (i in F.vals){
    ylims=c(0,69)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)
    for(j in 1:3){
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,i,3,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,i,3,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,i,3,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
          mtext(bquote(F[is]~"="~.(i)), side=3, line=1.0,cex=1.25)

    if (i==F.vals[1]){
        axis(2,cex.axis=1.5,at=seq(0,60,by=20))
    }
}
for (i in F.vals){
    ylims=c(0,69)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
   for(j in 1:3){
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,i,3,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,i,3,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,i,3,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    if (i==F.vals[1]){
        axis(2,cex.axis=1.5,at=seq(0,60,by=20))
    }
}
mtext("h",side=1,outer=T,line=2.5,cex=1.25)
## leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc2~"="~.(afs2b[x,2]))))
## legend("bottomleft",legend=c(as.expression(bquote(loc1~"="~.(afs2b[1,1]))),leg_expr),fill=c(NA,colors),cex=1)
mtext(expression(-log[10](P)),side=2,outer=T,line=2.25,cex=1.25)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0,,cex=1.25)
dev.copy2pdf(file="lpv_F_sa_65x700_3r_2loc_h_afs_loc1_0.8_neg_eff_1_h2_0.15.pdf")

### plot two locus lOR sa effect (b values)
quartz(width=10,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val=c("0.25") 
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
d.s=0.15
for (i in as.character(eff2b)){
    ylims=c(0,-3.5)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)
    for(j in 1:3){
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
    mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0,cex=1.25)

    if (i==as.character(eff2b[1])){
        axis(2,cex.axis=1.5,at=c(0,-1,-2,-3))
    }
}
for (i in as.character(eff2b)){
    ylims=c(0,5.5)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
   for(j in 1:3){
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(log(two.loc.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    if (i==as.character(eff2b[1])){
        axis(2,cex.axis=1.5,at=c(0,1,2,3,4,5))
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc2~"="~.(afs2b[x,2]))))
legend("topleft",legend=c(as.expression(bquote(loc1~"="~.(afs2b[1,1]))),leg_expr),fill=c(NA,colors),cex=1)
mtext("h",side=1,outer=T,line=2.5,cex=1.25)
mtext(expression(log(OR)),side=2,outer=T,line=2.25,cex=1.25)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0,,cex=1.25)
dev.copy2pdf(file="lOR_eff_sa_65x700_3r_2loc_h_afs_loc1_0.8_neg_eff_1_h2_0.15.pdf")
## pv sa eff
quartz(width=10,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val="0.25" 
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
d.s=0.15
for (i in as.character(eff2b)){
    ylims=c(0,50)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)
    for(j in 1:3){
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
          mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0,cex=1.25)

    if (i==as.character(eff2b[1])){
        axis(2,cex.axis=1.5,at=seq(0,60,by=20))
    }
}
for (i in as.character(eff2b)){
    ylims=c(0,90)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
   for(j in 1:3){
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(-log10(two.loc.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    if (i==as.character(eff2b[1])){
        axis(2,cex.axis=1.5,at=seq(0,80,by=20))
    }
}
mtext("h",side=1,outer=T,line=2.5,cex=1.25)
## leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc2~"="~.(afs2b[x,2]))))
## legend("bottomleft",legend=c(as.expression(bquote(loc1~"="~.(afs2b[1,1]))),leg_expr),fill=c(NA,colors),cex=1)
mtext(expression(-log[10](P)),side=2,outer=T,line=2.25,cex=1.25)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0,,cex=1.25)
dev.copy2pdf(file="lpv_eff_sa_65x700_3r_2loc_h_afs_loc1_0.8_neg_eff_1_h2_0.15_F_0.25.pdf")

### plot two locus lOR eu effect (b values)
quartz(width=10,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val=c("0") 
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
d.s=0.15
for (i in as.character(eff2b)){
    ylims=c(0,-2.75)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)
    for(j in 1:3){
        boxplot(t(log(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(log(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(log(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
    mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0,cex=1.25)

    if (i==as.character(eff2b[1])){
        axis(2,cex.axis=1.5,at=c(0,-1,-2,-3))
    }
}
for (i in as.character(eff2b)){
    ylims=c(0,4.25)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
   for(j in 1:3){
        boxplot(t(log(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(log(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(log(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    if (i==as.character(eff2b[1])){
        axis(2,cex.axis=1.5,at=c(0,1,2,3,4,5))
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc2~"="~.(afs2b[x,2]))))
legend("topright",legend=c(as.expression(bquote(loc1~"="~.(afs2b[1,1]))),leg_expr),fill=c(NA,colors),cex=1)
mtext("h",side=1,outer=T,line=2.5,cex=1.25)
mtext(expression(log(OR)),side=2,outer=T,line=2.25,cex=1.25)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0,,cex=1.25)
dev.copy2pdf(file="lOR_eff_eu_100x1500_6r_2loc_h_afs_loc1_0.8_neg_eff_1_h2_0.15.pdf")
## pv eu eff
quartz(width=10,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val="0" 
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
d.s=0.15
for (i in as.character(eff2b)){
    ylims=c(0,100)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)
    for(j in 1:3){
        boxplot(t(-log10(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(-log10(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(-log10(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
          mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0,cex=1.25)

    if (i==as.character(eff2b[1])){
        axis(2,cex.axis=1.5,at=seq(0,100,by=20))
    }
}
for (i in as.character(eff2b)){
    ylims=c(0,230)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
   for(j in 1:3){
        boxplot(t(-log10(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(-log10(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(-log10(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    if (i==as.character(eff2b[1])){
        axis(2,cex.axis=1.5,at=seq(0,200,by=50))
    }
}
mtext("h",side=1,outer=T,line=2.5,cex=1.25)
## leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc2~"="~.(afs2b[x,2]))))
## legend("bottomleft",legend=c(as.expression(bquote(loc1~"="~.(afs2b[1,1]))),leg_expr),fill=c(NA,colors),cex=1)
mtext(expression(-log[10](P)),side=2,outer=T,line=2.25,cex=1.25)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0,,cex=1.25)
dev.copy2pdf(file="lpv_eff_eu_100x1500_6r_2loc_h_afs_loc1_0.8_neg_eff_1_h2_0.15.pdf")
### plot two locus effects sa eps (b values)
quartz(width=10,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val="0.25"
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)

d.s=0.2
for (i in as.character(eps.int)){
    ylims=c(0,-4.5)
    for(j in 1:3){
        boxplot(log(t(two.loc.eps.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
    }
    mtext(bquote(ei~"="~.(i)), side=3, line=1.0,cex=1.25)

    if (i==eps.int[1]){
        axis(2,cex.axis=1.5,at=c(0,-1,-2,-3,-4))
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[2]~"="~.(afs2b[x,2]))))
legend("bottomright",legend=c(as.expression(bquote(loc1~"="~.(afs2b[1,1]))),leg_expr),fill=c(NA,colors),cex=1)
for (i in as.character(eps.int)){
    ylims=c(0,4.5)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    for(j in 1:3){
        boxplot(log(t(two.loc.eps.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    if (i==eps.int[1]){
        axis(2,cex.axis=1.5,at=c(0,1,2,3,4,5,6))
    }
}
mtext("h",side=1,outer=T,line=2.5,cex=1.25)
mtext(expression(log(OR)),side=2,outer=T,line=2.25,cex=1.25)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0,cex=1.25)
dev.copy2pdf(file="lOR_eps_sa_65x700_3r_2loc_h_afs_loc1_0.8_neg_eff_1_h2_0.15_F_0.25.pdf")
quartz(width=10,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val="0.25"
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)

d.s=0.2
for (i in as.character(eps.int)){
    ylims=c(0,70)
    for(j in 1:3){
        boxplot(-log10(t(two.loc.eps.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
    }
    mtext(bquote(ei~"="~.(i)), side=3, line=1.0,cex=1.25)

    if (i==eps.int[1]){
        axis(2,cex.axis=1.5,at=seq(0,80,by=20))
    }
}
## leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[1]~ "="~.(afs2b[x,1])~","~loc[2]~"="~.(afs2b[x,2]))))
## legend("bottomright",legend=c(as.expression(bquote(loc1~"="~.(afs2b[1,1]))),leg_expr),fill=c(NA,colors),cex=1)
for (i in as.character(eps.int)){
    ylims=c(0,90)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    for(j in 1:3){
        boxplot(-log10(t(two.loc.eps.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    if (i==eps.int[1]){
        axis(2,cex.axis=1.5,at=seq(0,80,by=20))
    }
}
mtext("h",side=1,outer=T,line=2.5,cex=1.25)
mtext(expression(-log[10](P)),side=2,outer=T,line=2.25,cex=1.25)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0,cex=1.25)
dev.copy2pdf(file="lpv_eps_sa_65x700_3r_2loc_h_afs_loc1_0.8_neg_eff_1_h2_0.15_F_0.25.pdf")

### plot two locus effects eu eps (b values)
quartz(width=10,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val="0"
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)

d.s=0.2
for (i in as.character(eps.int)){
    ylims=c(0,-3.75)
    for(j in 1:3){
        boxplot(log(t(two.loc.eu.eps.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
    }
    mtext(bquote(ei~"="~.(i)), side=3, line=1.0,cex=1.25)

    if (i==eps.int[1]){
        axis(2,cex.axis=1.5,at=c(0,-1,-2,-3))
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[1]~ "="~.(afs2b[x,1])~","~loc[2]~"="~.(afs2b[x,2]))))
legend("bottomright",legend=c(as.expression(bquote(loc1~"="~.(afs2b[1,1]))),leg_expr),fill=c(NA,colors),cex=1)
for (i in as.character(eps.int)){
    ylims=c(0,4.5)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    for(j in 1:3){
        boxplot(log(t(two.loc.eu.eps.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    if (i==eps.int[1]){
        axis(2,cex.axis=1.5,at=c(0,1,2,3,4,5,6))
    }
}
mtext("h",side=1,outer=T,line=2.5,cex=1.25)
mtext(expression(log(OR)),side=2,outer=T,line=2.25,cex=1.25)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0,cex=1.25)
dev.copy2pdf(file="lOR_eps_eu_100x1500_6r_2loc_h_afs_loc1_0.8_neg_eff_1_h2_0.15.pdf")
quartz(width=10,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val="0"
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
    xlims=c(1-d.s*1.5,length(H_levs)+d.s*1.5)

d.s=0.2
for (i in as.character(eps.int)){
    ylims=c(0,200)
    for(j in 1:3){
        boxplot(-log10(t(two.loc.eu.eps.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
    }
    mtext(bquote(ei~"="~.(i)), side=3, line=1.0,cex=1.25)

    if (i==eps.int[1]){
        axis(2,cex.axis=1.5,at=seq(0,200,by=50))
    }
}
## leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[1]~ "="~.(afs2b[x,1])~","~loc[2]~"="~.(afs2b[x,2]))))
## legend("bottomright",legend=c(as.expression(bquote(loc1~"="~.(afs2b[1,1]))),leg_expr),fill=c(NA,colors),cex=1)
for (i in as.character(eps.int)){
    ylims=c(0,200)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    for(j in 1:3){
        boxplot(-log10(t(two.loc.eu.eps.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    if (i==eps.int[1]){
        axis(2,cex.axis=1.5,at=seq(0,200,by=50))
    }
}
mtext("h",side=1,outer=T,line=2.5,cex=1.25)
mtext(expression(-log[10](P)),side=2,outer=T,line=2.25,cex=1.25)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0,cex=1.25)
dev.copy2pdf(file="lpv_eps_eu_100x1500_6r_2loc_h_afs_loc1_0.8_neg_eff_1_h2_0.15.pdf")

### plot two locus lOR eu dominance (b values)
quartz(width=7,height=7)
layout(matrix(1:2, ncol = 1,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.8:0.5","0.8:0.75","0.8:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1","0:1","1:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
F.val="0"
colors=c("red","green","deepskyblue")
bool.add=c(F,T,T)
d.s=0.2
for (i in as.character(eff2b[3])){
    ylims=c(0,-3.0)
    xlims=c(1-d.s,length(H_levs)+d.s)
    for(j in 1:3){
        boxplot(t(log(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(log(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(log(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"or",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
    #mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)

    #if (i=="-3"){
        axis(2,cex.axis=1.5)
    #}
}
for (i in as.character(eff2b[3])){
    ylims=c(0,3.5)
    xlims=c(1-d.s,length(H_levs)+d.s)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
   for(j in 1:3){
        boxplot(t(log(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(log(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(log(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"or",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    mtext("h",side=1,outer=F,line=2.5,cex=1.5)
    #if (i=="-3"){
        axis(2,cex.axis=1.5)
    #}
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc2~"="~.(afs2b[x,2]))))
legend("bottomleft",legend=c(as.expression(bquote(loc1~"="~.(afs2b[1,1]))),leg_expr),fill=c(NA,colors),cex=1)
mtext(expression(log(OR)),side=2,outer=T,line=2.25,cex=1.5)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0,,cex=1.5)
dev.copy2pdf(file="lOR_dominance_eu_100x1500_6r_2loc_h_afs_loc1_0.8_neg_eff_1_h2_0.15.pdf")
## pv eu dominance
quartz(width=7,height=7)
layout(matrix(1:2, ncol = 1,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
F.val="0"
colors=c("red","green","deepskyblue")
for (i in as.character(eff2b[3])){
    ylims=c(0,120)
    xlims=c(1-d.s,length(H_levs)+d.s)
    for(j in 1:3){
        boxplot(t(-log10(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(-log10(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(-log10(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"pv",1,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
    #mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)

    #if (i=="-3"){
        axis(2,cex.axis=1.5)
    #}
}
## leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc2~"="~.(afs2b[x,2]))))
## legend("bottom",legend=c(as.expression(bquote(loc1~"="~.(afs2b[1,1]))),leg_expr),fill=c(NA,colors),cex=1)
for (i in as.character(eff2b[3])){
    ylims=c(0,175)
    xlims=c(1-d.s,length(H_levs)+d.s)
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
   for(j in 1:3){
        boxplot(t(-log10(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=bool.add[j])
        boxplot(t(-log10(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
        boxplot(t(-log10(two.loc.eu.b[afs_val[j],H_levs,F.val,i,"pv",2,])), at=1:length(H_levs)-2*d.s+j*d.s,boxwex=d.s, col=colors[j],ylim=ylims,xlim=xlims, yaxt="n",xaxt="n",add=T)
    }
        axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    mtext("h",side=1,outer=F,line=2.5,cex=1.5)
    #if (i=="-3"){
        axis(2,cex.axis=1.5)
    #}
}
mtext(expression(-log[10](P)),side=2,outer=T,line=2.25,cex=1.5)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0,,cex=1.5)
dev.copy2pdf(file="pV_dominance_eu_100x1500_6r_2loc_h_afs_loc1_0.8_neg_eff_1_h2_0.15.pdf")

### plot two locus effects for F
quartz(width=12,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.75:0.5","0.75:0.75","0.75:0.9")
H_levs=c("0:0","0.5:0","0.5:0.5","0.5:1","0:1","1:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
eff_val="-1"
for (i in 1:(length(Fs2)-1)){
    d.s=0.2
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(log(t(two.loc[afs_val[1],H_levs,i,eff_val,"or",1,])), at=1:length(H_levs)-d.s,boxwex=d.s, col="red",ylim=c(0,-6),xlim=c(1-d.s,length(H_levs)+d.s), yaxt="n",xaxt="n")
    boxplot(log(t(two.loc[afs_val[2],H_levs,i,eff_val,"or",1,])), at=1:length(H_levs),boxwex=d.s, col="green",ylim=c(0,-6),xlim=c(1-d.s,length(H_levs)+d.s), yaxt="n",cex.axis=1.5,xaxt="n",add=T)
    #mtext("h",side=1,outer=F,line=2.5)
    mtext(bquote(F[is]~"="~.(Fs2[i])), side=3, line=1.0)
    boxplot(log(t(two.loc[afs_val[3],H_levs,i,eff_val,"or",1,])), at=1:length(H_levs)+d.s,boxwex=d.s, col="deepskyblue",ylim=c(0,-6),xlim=c(1-d.s,length(H_levs)+d.s), yaxt="n",xaxt="n",add=T)
    if (i==1){
        axis(2,cex.axis=1.5)
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[1]~ "="~.(afs2[x,1])~","~loc[2]~"="~.(afs2[x,2]))))
legend("bottomright",legend=leg_expr,title=c("avg. AFs"),fill=c("red","green","deepskyblue"),cex=1.5)
for (i in  1:(length(Fs2)-1)){
    d.s=0.2
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(abs(log(t(two.loc[afs_val[1],H_levs,i,eff_val,"or",2,]))), at=1:length(H_levs)-d.s,boxwex=d.s, col="red",ylim=c(0,8),xlim=c(1-d.s,length(H_levs)+d.s), yaxt="n",xaxt="n")
    axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    boxplot(abs(log(t(two.loc[afs_val[2],H_levs,i,eff_val,"or",2,]))), at=1:length(H_levs),boxwex=d.s, col="green",ylim=c(0,8),xlim=c(1-d.s,length(H_levs)+d.s), yaxt="n",cex.axis=1.5,xaxt="n",add=T)
    mtext("h",side=1,outer=F,line=2.5)
    #mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)
    boxplot(abs(log(t(two.loc[afs_val[3],H_levs,i,eff_val,"or",2,]))), at=1:length(H_levs)+d.s,boxwex=d.s, col="deepskyblue",ylim=c(0,8),xlim=c(1-d.s,length(H_levs)+d.s), yaxt="n",xaxt="n",add=T)
    if (i==1){
        axis(2,cex.axis=1.5)
    }
}
mtext(expression(log(OR)),side=2,outer=T,line=2.0)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0)
dev.copy2pdf(file="logOR_65x700_3r_2loc_h_afs_loc1_0.75_F_eff_1_h2_0.3.pdf")




### plot two locus effects
quartz(width=12,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.9:0.5","0.5:0.5","0.5:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))

for (i in as.character(eff2[2:4])){
    d.s=0.04
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(log(t(two.loc.b[afs_val[1],H_levs,"0",i,"or",1,])), at=hs2[H_levs,2]-d.s,boxwex=0.04, col="red",ylim=c(0,-3),xlim=c(-0.2,1.2), yaxt="n",xaxt="n")
    boxplot(log(t(two.loc.b[afs_val[2],H_levs,"0",i,"or",1,])), at=hs2[H_levs,2],boxwex=0.04, col="green",ylim=c(0,-3),xlim=c(-0.2,1.2), yaxt="n",cex.axis=1.5,xaxt="n",add=T)
    #mtext("h",side=1,outer=F,line=2.5)
    mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)
    boxplot(log(t(two.loc.b[afs_val[3],H_levs,"0",i,"or",1,])), at=hs2[H_levs,2]+d.s,boxwex=0.04, col="deepskyblue",ylim=c(0,-3),xlim=c(-0.2,1.2), yaxt="n",xaxt="n",add=T)
    if (i=="1"){
        axis(2,cex.axis=1.5)
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[1]~ "="~.(afs2[x,1])~","~loc[2]~"="~.(afs2[x,2]))))
legend("topright",legend=leg_expr,title=c("avg. AFs"),fill=c("red","green","deepskyblue"),cex=1.5)
for (i in as.character(eff2[2:4])){
    d.s=0.04
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(log(t(two.loc.b[afs_val[1],H_levs,"0",i,"or",2,])), at=hs2[H_levs,2]-d.s,boxwex=0.04, col="red",ylim=c(0,-6),xlim=c(-0.2,1.2), yaxt="n",xaxt="n")
    axis(1,at=hs2[H_levs,2],cex.axis=1.5,labels=H_labs)
    boxplot(log(t(two.loc.b[afs_val[2],H_levs,"0",i,"or",2,])), at=hs2[H_levs,2],boxwex=0.04, col="green",ylim=c(0,-6),xlim=c(-0.2,1.2), yaxt="n",cex.axis=1.5,xaxt="n",add=T)
    mtext("h",side=1,outer=F,line=2.5)
    #mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)
    boxplot(log(t(two.loc.b[afs_val[3],H_levs,"0",i,"or",2,])), at=hs2[H_levs,2]+d.s,boxwex=0.04, col="deepskyblue",ylim=c(0,-6),xlim=c(-0.2,1.2), yaxt="n",xaxt="n",add=T)
    if (i=="1"){
        axis(2,cex.axis=1.5)
    }
}
mtext(expression(log(OR)),side=2,outer=T,line=2.0)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0)
dev.copy2pdf(file="logOR_100x700_3r_2loc_h_afs_loc1_0.5.pdf")


quartz(width=12,height=7)
layout(matrix(1:6, ncol = 3,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=rownames(afs2[afs2$Var1==0.75,])
H_levs=c("0:0","0:0.5","0.5:0.5","0:1","1:0","0.5:1","1:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))

for (i in as.character(eff2[2:4])){
    d.s=0.2
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(log(t(two.loc.b[afs_val[1],H_levs,"0",i,"or",1,])), at=1:length(H_levs)-d.s,boxwex=d.s, col="red",ylim=c(0,-3),xlim=c(1-d.s,length(H_levs)+d.s), yaxt="n",xaxt="n")
    boxplot(log(t(two.loc.b[afs_val[2],H_levs,"0",i,"or",1,])), at=1:length(H_levs),boxwex=d.s, col="green",ylim=c(0,-3),xlim=c(1-d.s,length(H_levs)+d.s), yaxt="n",cex.axis=1.5,xaxt="n",add=T)
    #mtext("h",side=1,outer=F,line=2.5)
    mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)
    boxplot(log(t(two.loc.b[afs_val[3],H_levs,"0",i,"or",1,])), at=1:length(H_levs)+d.s,boxwex=d.s, col="deepskyblue",ylim=c(0,-3),xlim=c(1-d.s,length(H_levs)+d.s), yaxt="n",xaxt="n",add=T)
    if (i=="1"){
        axis(2,cex.axis=1.5)
    }
}
leg_expr=sapply(afs_val,function (x) as.expression(bquote(loc[1]~ "="~.(afs2[x,1])~","~loc[2]~"="~.(afs2[x,2]))))
legend("topright",legend=leg_expr,title=c("avg. AFs"),fill=c("red","green","deepskyblue"),cex=1.5)
for (i in as.character(eff2[2:4])){
    d.s=0.2
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(log(t(two.loc.b[afs_val[1],H_levs,"0",i,"or",2,])), at=1:length(H_levs)-d.s,boxwex=d.s, col="red",ylim=c(0,-6),xlim=c(1-d.s,length(H_levs)+d.s), yaxt="n",xaxt="n")
    axis(1,at=1:length(H_levs),cex.axis=1.5,labels=H_labs)
    boxplot(log(t(two.loc.b[afs_val[2],H_levs,"0",i,"or",2,])), at=1:length(H_levs),boxwex=d.s, col="green",ylim=c(0,-6),xlim=c(1-d.s,length(H_levs)+d.s), yaxt="n",cex.axis=1.5,xaxt="n",add=T)
    mtext("h",side=1,outer=F,line=2.5)
    #mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)
    boxplot(log(t(two.loc.b[afs_val[3],H_levs,"0",i,"or",2,])), at=1:length(H_levs)+d.s,boxwex=d.s, col="deepskyblue",ylim=c(0,-6),xlim=c(1-d.s,length(H_levs)+d.s), yaxt="n",xaxt="n",add=T)
    if (i=="1"){
        axis(2,cex.axis=1.5)
    }
}
mtext(expression(log(OR)),side=2,outer=T,line=2.0)
dev.copy2pdf(file="logOR_100x700_3r_2loc_h_afs_loc1_0.75.pdf")




c("0:0.5","0.5:0.5","1:0.5")
#0:0.5 0.5:0.5 1:0.5
quartz(width=12,height=7)
layout(matrix(1:4, ncol = 4), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.9:0.5","0.5:0.5","0.5:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1")
local_eff=attr(two.loc,"dimnames")$rel.eff
for (i in 1:length(local_eff)){
     d.s=0.03
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(abs(log(t(two.loc[afs_val[1],H_levs,"0",i,"or",1,]))), at=hs2[H_levs,2]-d.s,boxwex=0.03, col="red",ylim=c(0,6),xlim=c(-0.2,1.2), yaxt="n",cex.axis=1.5)
    mtext("h",side=1,outer=F,line=2.5)
    mtext(bquote(eff[2]~"="~.(local_eff[i])~eff[1]), side=3, line=1.0)
    boxplot(abs(log(t(two.loc[afs_val[2],H_levs,"0",i,"or",1,]))), at=hs2[H_levs,2],boxwex=0.03, col="green",ylim=c(0,6),xlim=c(-0.2,1.2), yaxt="n",xaxt="n",add=T)
    boxplot(abs(log(t(two.loc[afs_val[3],H_levs,"0",i,"or",1,]))), at=hs2[H_levs,2]+d.s,boxwex=0.03, col="deepskyblue",ylim=c(0,6),xlim=c(-0.2,1.2), yaxt="n",xaxt="n",add=T)
    if (i==1){
        axis(2,cex.axis=1.5)
    }
}
leg_expr=sapply(afs_val,function (x) paste0("loc1:",afs2[x,1],",loc2:",afs2[x,2]))
legend("bottomright",c("0.5","0.75","0.9"),title=c("avg. AFs"),fill=c("red","green","deepskyblue"),cex=1.5)
mtext(expression(log(OR)),side=2,outer=T,line=2.0)

quartz(width=12,height=7)
layout(matrix(1:4, ncol = 4), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
afs_val=c("0.9:0.5","0.5:0.5","0.5:0.9")
H_levs=c("0.5:0","0.5:0.5","0.5:1")
local_eff=attr(two.loc.b,"dimnames")$rel.eff
for (i in 1:length(local_eff)){
    
    d.s=0.03
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(abs(log(t(two.loc.b[afs_val[1],H_levs,"0",i,"or",1,]))), at=hs2[H_levs,2]-d.s,boxwex=0.03, col="red",ylim=c(0,6),xlim=c(-0.2,1.2), yaxt="n",cex.axis=1.5)
    mtext("h",side=1,outer=F,line=2.5)
    mtext(bquote(eff[2]~"="~.(local_eff[i])~eff[1]), side=3, line=1.0)
    boxplot(abs(log(t(two.loc.b[afs_val[2],H_levs,"0",i,"or",1,]))), at=hs2[H_levs,2],boxwex=0.03, col="green",ylim=c(0,6),xlim=c(-0.2,1.2), yaxt="n",xaxt="n",add=T)
    boxplot(abs(log(t(two.loc.b[afs_val[3],H_levs,"0",i,"or",1,]))), at=hs2[H_levs,2]+d.s,boxwex=0.03, col="lightblue",ylim=c(0,6),xlim=c(-0.2,1.2), yaxt="n",xaxt="n",add=T)
    if (i==1){
        axis(2,cex.axis=1.5)
    }
}
leg_expr=sapply(afs_val,function (x) paste0("loc1:",afs2[x,1],",loc2:",afs2[x,2]))
legend("bottomright",leg_expr,title=c("avg. AFs"),fill=c("red","green","lightblue"),cex=1.5)
mtext(expression(log(OR)),side=2,outer=T,line=2.0)










	###### logistic regression based on individual haplotypes
	
	
	resp <- rep(c(0,1),each=n.pools*n.replicates)
	
	geno.li = pos[,alllight]
	geno.dark = pos[,alldark]
	geno.lidark = cbind(geno.li,geno.dark )
	
	summary(glm(resp~t(geno.lidark)))
	
	
	logistic.odds.ratio_ajb[i,]  <- exp( glm(resp~t(geno.lidark))$coef[2:(1+n.loci)]  )
	

	############ is there correlation ? #########################
	#cors[i,1]<-fisher.test(as.factor(geno.dark[1,]),as.factor(geno.dark[8,]))$p
	#cors[i,2]<-cor(geno.dark[8,],geno.dark[1,])  # phi-coefficient
	#cors[i,3]<-fisher.test(as.factor(geno.lidark[8,]),as.factor(geno.lidark[1,]))$p
	#cors[i,4]<-cor(geno.lidark[8,],geno.lidark[1,])  # phi-coefficient
	#################################### not within dark pool but between light+dark pools ########################
	
	

	
	
	
#}

##############  Summarize our finding

o.r <- apply(odds.ratio_ajb,2,mean)
l.o.r <- apply(logistic.odds.ratio_ajb,2,mean)

o.eff<-order(eff)
o.o.r = order(o.r)
o.l.o.r = order(l.o.r)


gt.a=c(0.0,0.5,1.0)
gt.b=c(0.0,0.5,1.0)
gt.ab=expand.grid(gt.a,gt.b)
eff.a=1.0
eff.b=1.0
colors=c("red","blue","green")
quartz()
h.eps=c(0.5,0.5)
gt.eff=apply(gt.ab*2,1,function(x){(loc.eff(h.eps,x) - 0.5)})

plot(gt.a,(eff.a*gt.eff[1,1:3]+eff.b*gt.eff[2,1:3]),type="n",col="red",ylim=c(-1,1), main="Independent Loci",xaxt="n",xlab="GT",ylab="P")
for (i in 1:3){
    idxs=(3*(i-1)+1):(3*i)
    lines(gt.a,(eff.a*gt.eff[1,idxs]+eff.b*gt.eff[2,idxs]),type="o",col=colors[i])
}
axis(1,at=c(0,0.5,1),labels=c("aa","aA","AA"))
legend("bottomright",legend=c("bb","bB","BB"),col=colors,lty=1,pch=1)
dev.copy2pdf(file="independent_loci_reac_norm.pdf")

quartz()
plot(gt.a,(eff.a*gt.a+eff.b*gt.b[1]),type="n",ylim=c(-1.5,1.5), main="Positively Interacting Loci",xaxt="n",xlab="GT",ylab="P")
for (i in 1:3){
    idxs=(3*(i-1)+1):(3*i)
    lines(gt.a,(eff.a*gt.eff[1,idxs]+eff.b*gt.eff[2,idxs] + eff.b*gt.eff[1,idxs]*gt.eff[2,idxs]),type="o",col=colors[i])
}
axis(1,at=c(0,0.5,1),labels=c("aa","aA","AA"))
legend("bottomright",legend=c("bb","bB","BB"),col=colors,lty=1,pch=1)
dev.copy2pdf(file="pos_loci_reac_norm.pdf")

quartz()
plot(gt.a,(eff.a*gt.a+eff.b*gt.b[1]),type="n",ylim=c(-1.5,1.5), main="Negatively Interacting Loci",xaxt="n",xlab="GT",ylab="P")
for (i in 1:3){
    idxs=(3*(i-1)+1):(3*i)
lines(gt.a,(eff.a*gt.eff[1,idxs]+eff.b*gt.eff[2,idxs] - eff.b*gt.eff[1,idxs]*gt.eff[2,idxs]),type="o",col=colors[i])
}
axis(1,at=c(0,0.5,1),labels=c("aa","aA","AA"))
legend("bottomright",legend=c("bb","bB","BB"),col=colors,lty=1,pch=1)
dev.copy2pdf(file="neg_loci_reac_norm.pdf")


h.eps=c(0.5,0.5)
gt.eff=apply(gt.ab*2,1,function(x){(loc.eff(h.eps,x) - 0.5)})
quartz()
plot(gt.a,(eff.a*gt.eff[1,1:3]+eff.b*gt.eff[2,1:3]),type="o",col="red",ylim=c(-1,1), main="Independent Loci",xaxt="n",xlab="GT",ylab="P")
for (i in 1:3){
idxs=(3*(i-1)+1):(3*i)
h.eps=c(0.5,0.5)
gt.eff=apply(gt.ab*2,1,function(x){(loc.eff(h.eps,x) - 0.5)})
h.eps=c(0.5,0)
lines(gt.a,(eff.a*gt.eff[1,idxs]+eff.b*gt.eff[2,idxs]),type="o",col=colors[i],lty=1,pch=1)
gt.eff=apply(gt.ab*2,1,function(x){(loc.eff(h.eps,x) - 0.5)})
lines(gt.a,(eff.a*gt.eff[1,idxs]+eff.b*gt.eff[2,idxs]),type="o",col=colors[i],lty=2,pch=2)
h.eps=c(0.5,1.0)
gt.eff=apply(gt.ab*2,1,function(x){(loc.eff(h.eps,x) - 0.5)})
lines(gt.a,(eff.a*gt.eff[1,idxs]+eff.b*gt.eff[2,idxs]),type="o",col=colors[i],lty=3,pch=3)
}
dev.copy2pdf(file="independent_loci_reac_norm.pdf")

quartz()
plot(gt.a,(eff.a*gt.eff[1,1:3]+eff.b*gt.eff[2,1:3]),type="n",col="red",ylim=c(-1.5,1.5), main="Positively Interacting Loci",xaxt="n",xlab="GT",ylab="P")
for (i in 1:3){
idxs=(3*(i-1)+1):(3*i)
h.eps=c(0.5,0.5)
gt.eff=apply(gt.ab*2,1,function(x){(loc.eff(h.eps,x) - 0.5)})
h.eps=c(0.5,0)
lines(gt.a,(eff.a*gt.eff[1,idxs]+eff.b*gt.eff[2,idxs]+ eff.a*gt.eff[1,idxs]*gt.eff[2,idxs]),type="o",col=colors[i],lty=3,pch=1)
gt.eff=apply(gt.ab*2,1,function(x){(loc.eff(h.eps,x) - 0.5)})
lines(gt.a,(eff.a*gt.eff[1,idxs]+eff.b*gt.eff[2,idxs]+ eff.a*gt.eff[1,idxs]*gt.eff[2,idxs]),type="o",col=colors[i],lty=2,pch=2)
h.eps=c(0.5,1.0)
gt.eff=apply(gt.ab*2,1,function(x){(loc.eff(h.eps,x) - 0.5)})
lines(gt.a,(eff.a*gt.eff[1,idxs]+eff.b*gt.eff[2,idxs]+ eff.a*gt.eff[1,idxs]*gt.eff[2,idxs]),type="o",col=colors[i],lty=2,pch=3)
}
axis(1,at=c(0,0.5,1),labels=c("aa","aA","AA"))
legend("bottomright",legend=c("bb","bB","BB","B rec.","B add.", "B dom."),col=c(colors,"black","black","black"),lty=c(1,1,1,2,3,2),pch=c(NA,NA,NA,2,1,3))

### Runs for assessing parameter combinations
### Fixed afs
afs.eu=c(0.2,0.83)
afs.sa=c(0.17,0.47)
### Fis
###F.sa=0.25
F.sa=0.35
F.eu=0.0
### base effect
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
n.runs=1000

# parameters to go over
#eff2=c(3,2,1,0.5)
eff2=c(3,2,1) # for dominance
#hs2=c(0.0,0.25,0.5,0.75,1.0)
hs2=c(0.0,0.5,0.75,1.0,1.25)
hs2=expand.grid(hs2,hs2)
rownames(hs2)=apply(hs2,1,function (x) paste(x,collapse=":"))
eps.int=c(-1,-0.5,0,0.5,1)
eps.int=c(-1,-0.5,0,0.5,1) # for dominance
res.fields2=c("pv","or","ol","od","af.m","af.l","af.d")
two.loc.sa=array(0,dim=c(length(hs2[,1]),length(eff2),length(res.fields2),2,n.runs),dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"resfields"=res.fields2,"locus"=c(1,2),"run"=1:n.runs))
two.loc.eu=array(0,dim=c(length(hs2[,1]),length(eff2),length(res.fields2),2,n.runs),dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"resfields"=res.fields2,"locus"=c(1,2),"run"=1:n.runs))
two.loc.eps.sa=array(0,dim=c(length(hs2[,1]),length(eff2),length(eps.int),length(res.fields2),2,n.runs),dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"eps.int"=eps.int,"resfields"=res.fields2,"locus"=c(1,2),"run"=1:n.runs))
two.loc.eps.eu=array(0,dim=c(length(hs2[,1]),length(eff2),length(eps.int),length(res.fields2),2,n.runs),dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"eps.int"=eps.int,"resfields"=res.fields2,"locus"=c(1,2),"run"=1:n.runs))
two.loc.eps.xs=array(0,dim=c(length(hs2[,1]),length(eff2),length(eps.int),length(res.fields2),2,n.runs),dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"eps.int"=eps.int,"resfields"=res.fields2,"locus"=c(1,2),"run"=1:n.runs))

for(j in 1:length(hs2[,1])){
   for(l in 1:length(eff2)){
       eff=c(eff1/(1+abs(eff2[l])),eff1*eff2b[l]/(1+abs(eff2[l])))
       h.loci=unlist(hs2[j,])
       v.E=get.V.E(afs.cosm,eff,h.loci,h2)
       F.inbreed=F.sa
       n.freqs=afs.sa
       run.result=run_simulations(n.runs,v.E,total.ind=total.ind.sa,n.replicates=n.rep.sa,n.pools=n.pools.sa)
       two.loc.sa[j,l,,1,]=t(sapply(res.fields2, function (x) run.result[[x]][,1]))
       two.loc.sa[j,l,,2,]=t(sapply(res.fields2, function (x) run.result[[x]][,2]))
       F.inbreed=F.eu
       n.freqs=afs.eu
       run.result=run_simulations(n.runs,v.E,total.ind=total.ind.eu,n.replicates=n.rep.eu,n.pools=n.pools.eu)
       two.loc.eu[j,l,,1,]=t(sapply(res.fields2, function (x) run.result[[x]][,1]))
       two.loc.eu[j,l,,2,]=t(sapply(res.fields2, function (x) run.result[[x]][,2]))
   }
}

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
            run.result=run_simulations.eps(n.runs,v.E,total.ind=total.ind.sa,n.replicates=n.rep.sa,n.pools=n.pools.xs,eps.int=eps.int[l])
            two.loc.eps.xs[j,k,l,,1,]=t(sapply(res.fields2, function (x) run.result[[x]][,1]))
            two.loc.eps.xs[j,k,l,,2,]=t(sapply(res.fields2, function (x) run.result[[x]][,2]))

            F.inbreed=F.eu
            n.freqs=afs.eu
            run.result=run_simulations.eps(n.runs,v.E,total.ind=total.ind.eu,n.replicates=n.rep.eu,n.pools=n.pools.eu,eps.int=eps.int[l])
            two.loc.eps.eu[j,k,l,,1,]=t(sapply(res.fields2, function (x) run.result[[x]][,1]))
            two.loc.eps.eu[j,k,l,,2,]=t(sapply(res.fields2, function (x) run.result[[x]][,2]))
        }
    }
}

save(two.loc.eps.eu,two.loc.eps.sa,two.loc.eps.xs,file="two.loc.eu.sa_c.RData")

# just xs
for(j in 1:length(hs2[,1])){
    for(k in 1:length(eff2)){
        for(l in 1:length(eps.int)){
            eff=c(eff1/(1+abs(eff2[k])),eff1*eff2[k]/(1+abs(eff2[k])))
            h.loci=unlist(hs2[j,])
            v.E=get.V.E(afs.cosm,eff,h.loci,h2,eps.int=eps.int[l])
            F.inbreed=F.sa
            n.freqs=afs.sa
            run.result=run_simulations.eps(n.runs,v.E,total.ind=total.ind.sa,n.replicates=n.rep.sa,n.pools=n.pools.xs,eps.int=eps.int[l])
            two.loc.eps.xs[j,k,l,,1,]=t(sapply(res.fields2, function (x) run.result[[x]][,1]))
            two.loc.eps.xs[j,k,l,,2,]=t(sapply(res.fields2, function (x) run.result[[x]][,2]))
        }
    }
}

save(two.loc.eps.eu,two.loc.eps.sa,two.loc.eps.xs,file="two.loc.eu.sa.RData")
# get fract. of log(OR) for each combination lying in the right bin
bin.eu=c(0.25,0.5,0.75,1.0)
bin.sa=c(0.1,0.25,0.5,1.0,2.0)
tar.eu=-2.912427
tar.sa=-0.600489
tar.eu.sa=c(1.72,0.35)

perc.int <- function(target,intervals,values,values2=F){        
        if (values2[1] != F) {
            ratios=values[1,]/values2[1,]
            ratios2=values[2,]/values2[2,]
            return(sapply(intervals, function(x) { length(ratios[ abs(ratios - target[1]) <= abs(target[1]*x/2) &  abs(ratios2 - target[2]) <= abs(target[2]*x/2) ] ) })/length(ratios))
        }
        else {
            ratios=values[1,]/values[2,]
            return(sapply(intervals, function(x) { length(ratios[ abs(ratios - target) <= abs(target*x/2) ]) })/length(ratios))
        }
}
sqd <- function(target,values,values2=F){        
        if (values2[1] != F) {
            ratios=values[1,]/values2[1,]
            ratios2=values[2,]/values2[2,]
            return(mean((ratios-target[1])^2)+mean((ratios2-target[2])^2))
        }
        else {
            ratios=values[1,]/values[2,]
            return(mean((ratios-target)^2))
        }
}

p.eu = array(0,dim=c(length(hs2[,1]),length(eff2),length(bin.eu) ), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"bins"=bin.eu)) 
p.sa = array(0,dim=c(length(hs2[,1]),length(eff2),length(bin.eu) ), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"bins"=bin.eu))
p.eu.sa = array(0,dim=c(length(hs2[,1]),length(eff2),length(bin.eu) ), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"bins"=bin.eu)) 

m.eu = array(0,dim=c(length(hs2[,1]),length(eff2),2), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"loc"=c(1,2))) 
m.sa = array(0,dim=c(length(hs2[,1]),length(eff2),2), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"loc"=c(1,2)))


for(j in 1:length(hs2[,1])){
   for(l in 1:length(eff2)){
       p.eu[j,l,]=perc.int(tar.eu,bin.eu,log(two.loc.eu[j,l,"or",,]))
       p.sa[j,l,]=perc.int(tar.sa,bin.eu,log(two.loc.sa[j,l,"or",,]))
       p.eu.sa[j,l,]=perc.int(tar.eu.sa,bin.eu,log(two.loc.eu[j,l,"or",,]),log(two.loc.sa[j,l,"or",,]))
       m.eu[j,l,]=apply(log(two.loc.eu[j,l,"or",,]),1,median)
       m.sa[j,l,]=apply(log(two.loc.sa[j,l,"or",,]),1,median)
   }
}

p.eu.eps = array(0,dim=c(length(hs2[,1]),length(eff2),length(eps.int),length(bin.eu) ), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"eps.int"=eps.int,"bins"=bin.eu)) 
p.sa.eps =  array(0,dim=c(length(hs2[,1]),length(eff2),length(eps.int),length(bin.eu) ), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"eps.int"=eps.int,"bins"=bin.eu))
p.eu.sa.eps = array(0,dim=c(length(hs2[,1]),length(eff2),length(eps.int),length(bin.eu) ), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"eps.int"=eps.int,"bins"=bin.eu)) 
m.eu.eps = array(0,dim=c(length(hs2[,1]),length(eff2),length(eps.int),2), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"eps.int"=eps.int,"loc"=c(1,2)))
m.sa.eps = array(0,dim=c(length(hs2[,1]),length(eff2),length(eps.int),2), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"eps.int"=eps.int,"loc"=c(1,2)))
sqd.eu.eps = array(0,dim=c(length(hs2[,1]),length(eff2),length(eps.int)), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"eps.int"=eps.int)) 
sqd.sa.eps =  array(0,dim=c(length(hs2[,1]),length(eff2),length(eps.int) ), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"eps.int"=eps.int))
sqd.eu.sa.eps = array(0,dim=c(length(hs2[,1]),length(eff2),length(eps.int) ), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2,"eps.int"=eps.int)) 


for(j in 1:length(hs2[,1])){
   for(k in 1:length(eff2)){
       for(l in 1:length(eps.int)){
           p.eu.eps[j,k,l,]=perc.int(tar.eu,bin.eu,log(two.loc.eps.eu[j,k,l,"or",,]))
           p.sa.eps[j,k,l,]=perc.int(tar.sa,bin.eu,log(two.loc.eps.sa[j,k,l,"or",,]))
           p.eu.sa.eps[j,k,l,]=perc.int(tar.eu.sa,bin.eu,log(two.loc.eps.eu[j,k,l,"or",,]),log(two.loc.eps.sa[j,k,l,"or",,]))
           sqd.eu.eps[j,k,l]=sqd(tar.eu,log(two.loc.eps.eu[j,k,l,"or",,]))
           sqd.sa.eps[j,k,l]=sqd(tar.sa,log(two.loc.eps.sa[j,k,l,"or",,]))
           sqd.eu.sa.eps[j,k,l]=sqd(tar.eu.sa,log(two.loc.eps.eu[j,k,l,"or",,]),log(two.loc.eps.sa[j,k,l,"or",,]))
           m.eu.eps[j,k,l,]=apply(log(two.loc.eps.eu[j,k,l,"or",,]),1,median)
           m.sa.eps[j,k,l,]=apply(log(two.loc.eps.sa[j,k,l,"or",,]),1,median)
       }
   }
}
criteria=c(-2.51,-0.57,1.89,0.44)
criteria.xs=c(-2.51,-0.68,2.28,0.62)
criteria.xs=c(2.51,0.68,2.28,0.62)
criteria=c(2.51,0.57,1.89,0.44)
names(criteria)=c("eu","sa","tan_eusa","bab_eusa")
#var.crit.eps =  array(0,dim=c(length(hs2[,1]),length(eff2)-1,length(eps.int),length(criteria)), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2[1:(length(eff2)-1)],"eps.int"=eps.int,"criteria"=names(criteria)))
#var.crit.eps.xs =  array(0,dim=c(length(hs2[,1]),length(eff2)-1,length(eps.int),length(criteria)), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2[1:(length(eff2)-1)],"eps.int"=eps.int,"criteria"=names(criteria)))
load("two.loc.eu.sa.RData")
# use upper limit for OR
min.or=1.0
max.or=1.0
a=two.loc.eps.eu[,,,"or",,]
min.or=min(1.0,min(a[which(a != 0)]))
max.or=max(1.0,max(a[which(a != 0)]))
a=two.loc.eps.sa[,,,"or",,]
min.or=min(1.0,min(a[which(a != 0)]))
max.or=max(1.0,max(a[which(a != 0)]))
a=two.loc.eps.xs[,,,"or",,]
min.or=min(1.0,min(a[which(a != 0)]))
max.or=max(1.0,max(a[which(a != 0)]))
min.or=min.or/2
max.or=max.or*2
a=two.loc.eps.eu[,,,"or",,]
a[a <= min.or] = min.or; a[a >= max.or] = max.or
two.loc.eps.eu[,,,"or",,]=a
a=two.loc.eps.sa[,,,"or",,]
a[a <= min.or] = min.or; a[a >= max.or] = max.or
two.loc.eps.sa[,,,"or",,]=a
a=two.loc.eps.xs[,,,"or",,]
a[a <= min.or] = min.or; a[a >= max.or] = max.or
two.loc.eps.xs[,,,"or",,]=a

var.crit.eps =get.var(two.loc.eps.eu,two.loc.eps.sa,citeria,drop.eff=0)
var.crit.eps.xs =get.var(two.loc.eps.eu,two.loc.eps.xs,citeria.xs)

get.var <- function (two.loc.eps.1,two.loc.eps.2,citeria,drop.eff=1)
{
    var.crit.eps =  array(0,dim=c(length(hs2[,1]),length(eff2)-drop.eff,length(eps.int),length(criteria)), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2[1:(length(eff2)-drop.eff)],"eps.int"=eps.int,"criteria"=names(criteria)))
    for(j in 1:length(hs2[,1])){
        for(k in 1:(length(eff2)-drop.eff)){
            for(l in 1:length(eps.int)){
                a=log(two.loc.eps.1[j,k,l,"or",1,])/log(two.loc.eps.1[j,k,l,"or",2,])
                #b=quantile(a,c(0.005,0.995),na.rm=T); a=a[a >= b[1] & a <= b[2] ]
                                        #a[a>=10]=10; a[a<=-10]=-10; 
                var.crit.eps[j,k,l,1]=var((a[is.finite(a)]-criteria[1]))
                a=log(two.loc.eps.2[j,k,l,"or",1,])/log(two.loc.eps.2[j,k,l,"or",2,])
                #b=quantile(a,c(0.005,0.995),na.rm=T); a=a[a >= b[1] & a <= b[2] ]
                                        #a[a>=10]=10; a[a<=-10]=-10
                var.crit.eps[j,k,l,2]=var((a[is.finite(a)]-criteria[2]))
                a=log(two.loc.eps.1[j,k,l,"or",1,])/log(two.loc.eps.2[j,k,l,"or",1,])
                #b=quantile(a,c(0.005,0.995),na.rm=T); a=a[a >= b[1] & a <= b[2] ]
                                        #a[a>=10]=10; a[a<=-10]=-10
                var.crit.eps[j,k,l,3]=var((a[is.finite(a)]-criteria[3]))
                a=log(two.loc.eps.1[j,k,l,"or",2,])/log(two.loc.eps.2[j,k,l,"or",2,])
                #b=quantile(a,c(0.005,0.995),na.rm=T); a=a[a >= b[1] & a <= b[2] ]
                                        #a[a>=10]=10; a[a<=-10]=-10
                var.crit.eps[j,k,l,4]=var((a[is.finite(a)]-criteria[4]))
            }
        }
    }
    return(var.crit.eps)
}

library(reshape2)
## test.array=melt(var.crit.eps[,,,1],id=c("rel.eff"))
## boxplot(value~hs,data=melt(test.array))
med.var.crit=sapply(1:4,function(x) median(var.crit.eps[,,,x]))
mean.var.crit=sapply(1:4,function(x) mean(var.crit.eps[,,,x]))
med.var.crit.xs=sapply(1:4,function(x) median(var.crit.eps.xs[,,,x]))
mean.var.crit.xs=sapply(1:4,function(x) mean(var.crit.eps.xs[,,,x]))

#lh.two.loc.eps=array(0,dim=c(length(hs2[,1]),length(eff2)-1,length(eps.int),length(criteria)), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2[1:(length(eff2)-1)],"eps.int"=eps.int,"criteria"=names(criteria)))
lh.two.loc.eps=get.lh(two.loc.eps.eu,two.loc.eps.sa,med.var.crit,criteria)
lh.two.loc.eps.xs=get.lh(two.loc.eps.eu,two.loc.eps.xs,med.var.crit.xs,criteria.xs)
## lh.two.loc.eps=get.lh(two.loc.eps.eu,two.loc.eps.sa,var.crit.eps,criteria,bw=1)
## lh.two.loc.eps.xs=get.lh(two.loc.eps.eu,two.loc.eps.xs,var.crit.eps.xs,criteria.xs,bw=1)


# for overdominance
var.crit.eps =get.var(two.loc.eps.eu,two.loc.eps.sa,citeria,drop.eff=0)
med.var.crit=sapply(1:4,function(x) median(var.crit.eps[,,,x]))
mean.var.crit=sapply(1:4,function(x) mean(var.crit.eps[,,,x]))
lh.two.loc.eps=get.lh(two.loc.eps.eu,two.loc.eps.sa,med.var.crit,criteria,drop.eff=0)
var.crit.eps.xs =get.var(two.loc.eps.eu,two.loc.eps.xs,citeria,drop.eff=0)
med.var.crit.xs=sapply(1:4,function(x) median(var.crit.eps.xs[,,,x]))
mean.var.crit.xs=sapply(1:4,function(x) mean(var.crit.eps.xs[,,,x]))
lh.two.loc.eps.xs=get.lh(two.loc.eps.eu,two.loc.eps.xs,med.var.crit,criteria,drop.eff=0)

get.lh <- function (two.loc.eps.1,two.loc.eps.2,med.var.crit.int,criteria,drop.eff=1)
{
    lh.two.loc.eps=array(0,dim=c(length(hs2[,1]),length(eff2)-drop.eff,length(eps.int),length(criteria)), dimnames=list("hs"=rownames(hs2),"rel.eff"=eff2[1:(length(eff2)-drop.eff)],"eps.int"=eps.int,"criteria"=names(criteria)))
    for(j in 1:length(hs2[,1])){
        for(k in 1:(length(eff2)-drop.eff)){
            for(l in 1:length(eps.int)){
                a=log(two.loc.eps.1[j,k,l,"or",1,])/log(two.loc.eps.1[j,k,l,"or",2,])
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

kernel.dens <- function(a,criteria,crit,med.var.crit.int){
    n=length(a)
    bw=bw.nrd0((a-criteria[crit])/med.var.crit.int[crit])
    #return(mean((a-criteria[crit])^2))
    return(1/(2*pi*n*bw)*sum(exp(-(a-criteria[crit])^2/(2*bw*med.var.crit.int[crit]))))
    #return(density((a-criteria[crit])/med.var.crit.int[crit],bw=bw.nrd0((a-criteria[crit])/med.var.crit.int[crit]),kernel="gaussian"))
}

# get tan recessive
tanrec= which(hs2$Var1<0.5)
tandom= which(hs2$Var1>0.5)
tanadd=which(hs2$Var1 == 0.5)
lh.two.loc.eps.c12=lh.two.loc.eps[,,,1]*lh.two.loc.eps[,,,2]
lh.two.loc.eps.c34=lh.two.loc.eps[,,,3]*lh.two.loc.eps[,,,4]
lh.two.loc.eps.alls=lh.two.loc.eps[,,,1]*lh.two.loc.eps[,,,2]*lh.two.loc.eps[,,,3]*lh.two.loc.eps[,,,4]

sum(lh.two.loc.eps[tanrec,,,1])
sum(lh.two.loc.eps[tandom,,,1])
sum(lh.two.loc.eps[tanadd,,,1])
sum(lh.two.loc.eps[tanrec,,,2])
sum(lh.two.loc.eps[tandom,,,2])
sum(lh.two.loc.eps[tanadd,,,2])
sum(lh.two.loc.eps.c12[tandom,,])/sum(lh.two.loc.eps.c12[c(tanrec,tanadd),,])
sum(lh.two.loc.eps.c34[tandom,,])/sum(lh.two.loc.eps.c34[c(tanrec,tanadd),,])
sum(lh.two.loc.eps.alls[tandom,,])/sum(lh.two.loc.eps.alls[c(tanrec,tanadd),,])



x11()
a=log(two.loc.eps.eu[j,k,l,"or",1,])/log(two.loc.eps.eu[j,k,l,"or",2,])
b=quantile(a,c(0.01,0.99))
shapiro.test((a[a >= b[1] & a <= b[2]]))
hist((a[a >= b[1] & a <= b[2]]))
var(a[a >= b[1] & a <= b[2]])
qqplot(-rgamma(10000,shape=77,scale=1),a)
hist(-rgamma(10000,shape=20,scale=1))
qqnorm(a)
qqline(a)
shapiro.test(a)
shapiro.test(rnorm(100, mean = 5, sd = 3))
# gamma distribution
# shape = mean^2/var
# scale= var/mean
b=quantile(a,c(0.01,0.99))
a.mean=mean(a[a >= b[1] & a <= b[2]])
a.var=var(a[a >= b[1] & a <= b[2]])
a.shape=a.mean^2/a.var
a.scale=a.var/abs(a.mean)
qqplot(sign(a.mean)*rgamma(10000,shape=a.shape,scale=a.scale),a)

sum.lh.two = lh.two.loc.eps[,,,1]+lh.two.loc.eps[,,,2]+lh.two.loc.eps[,,,3]+lh.two.loc.eps[,,,4]
prod.lh.two = lh.two.loc.eps[,,,1]*lh.two.loc.eps[,,,2]*lh.two.loc.eps[,,,3]*lh.two.loc.eps[,,,4]
prod.lh.two.xs = lh.two.loc.eps.xs[,,,1]*lh.two.loc.eps.xs[,,,2]*lh.two.loc.eps.xs[,,,3]*lh.two.loc.eps.xs[,,,4]

# test kd estimation comparison
j=24; k=2; l=4 
a=log(two.loc.eps.eu[j,k,l,"or",1,])/log(two.loc.eps.eu[j,k,l,"or",2,])
a=(a-criteria[1])/med.var.crit[1]
b=log(two.loc.eps.sa[j,k,l,"or",1,])/log(two.loc.eps.sa[j,k,l,"or",2,])
b=(b-criteria[2])/med.var.crit[2]

lh.two.loc.eps[j,k,l,]
abc=bkde2D(cbind(a,b),c(bw.nrd0(a),bw.nrd0(b)),gridsize=c(250,250),range.x=list(c(0,1.0),c(0,10)))
abc=bkde(a, kernel = "normal", canonical = FALSE,bandwidth=bw.nrd0(a),range.x=c(-0.1,0.1))


lh.two.loc.eps=get.lh(two.loc.eps.eu,two.loc.eps.sa,med.var.crit,criteria,drop.eff=0)
lh.two.loc.eps.xs=get.lh(two.loc.eps.eu,two.loc.eps.xs,med.var.crit.xs,criteria.xs)
library(stats)
x11(width=10,height=9)
quartz(width=10,height=9)
library(reshape2,ggplot2)
source("/Volumes/Temp/Lukas/Tools/Scripts/R/image.scale.R")
effxeps=expand.grid(eff2[1:3],eps.int)
effxeps=apply(effxeps,1,function (x) paste("Eff:EI",x,collapse=":"))
plot.heatmap.lh(lh.two.loc.eps[,,,1]*lh.two.loc.eps[,,,2]*lh.two.loc.eps[,,,3]*lh.two.loc.eps[,,,4],doms=sort(unique(unlist(hs2[1]))))
dev.copy2pdf(file="lh_kerneldensity_parameters.pdf")
plot.heatmap.lh(lh.two.loc.eps[,,,1])
dev.copy2pdf(file="lh_kerneldensity_crit_eu_parameters.pdf")
plot.heatmap.lh(lh.two.loc.eps[,,,2])
dev.copy2pdf(file="lh_kerneldensity_crit_sa_parameters.pdf")
plot.heatmap.lh(lh.two.loc.eps[,,,3])
dev.copy2pdf(file="lh_kerneldensity_crit_eusa1_parameters.pdf")
plot.heatmap.lh(lh.two.loc.eps[,,,4])
dev.copy2pdf(file="lh_kerneldensity_crit_eusa2_parameters.pdf")
plot.heatmap.lh(lh.two.loc.eps[,,,3]*lh.two.loc.eps[,,,4])
dev.copy2pdf(file="lh_kerneldensity_crit_eusa_parameters.pdf")
plot.heatmap.lh(lh.two.loc.eps[,,,1]*lh.two.loc.eps[,,,2])
dev.copy2pdf(file="/Volumes/Temp/Lukas/Data/A7_pigmentation/Simulation/lh_kerneldensity_crit12_eu_sa_parameters.pdf")

plot.heatmap.lh(lh.two.loc.eps.xs[,,,1]*lh.two.loc.eps.xs[,,,2]*lh.two.loc.eps.xs[,,,3]*lh.two.loc.eps.xs[,,,4])
dev.copy2pdf(file="lh_xs_kerneldensity_parameters.pdf")
plot.heatmap.lh(lh.two.loc.eps.xs[,,,1])
dev.copy2pdf(file="lh_xs_kerneldensity_crit_eu_parameters.pdf")
plot.heatmap.lh(lh.two.loc.eps.xs[,,,2])
dev.copy2pdf(file="lh_xs_kerneldensity_crit_sa_parameters.pdf")
plot.heatmap.lh(lh.two.loc.eps.xs[,,,3])
dev.copy2pdf(file="lh_xs_kerneldensity_crit_eusa1_parameters.pdf")
plot.heatmap.lh(lh.two.loc.eps.xs[,,,4])
dev.copy2pdf(file="lh_xs_kerneldensity_crit_eusa2_parameters.pdf")
plot.heatmap.lh(lh.two.loc.eps.xs[,,,3]*lh.two.loc.eps.xs[,,,4])
dev.copy2pdf(file="lh_xs_kerneldensity_crit_eusa_parameters.pdf")
x11(width=10,height=9)
plot.heatmap.lh(lh.two.loc.eps.xs[,,,1]*lh.two.loc.eps.xs[,,,2])
dev.copy2pdf(file="/Volumes/Temp/Lukas/Data/A7_pigmentation/Simulation/lh_xs_kerneldensity_crit12_eu_sa_parameters.pdf")

x11(width=10,height=9)
doms=sort(unique(unlist(hs2[1])))
effxeps=expand.grid(eff2[1:3],eps.int)
effxeps=apply(effxeps,1,function (x) paste("Eff:EI",x,collapse=":"))
plot.heatmap.lh(lh.two.loc.eps[,,,1]*lh.two.loc.eps[,,,2]*lh.two.loc.eps[,,,3]*lh.two.loc.eps[,,,4],doms=sort(unique(unlist(hs2[1]))))
dev.copy2pdf(file="lh_kerneldensity_parameters_od.pdf")
plot.heatmap.lh(lh.two.loc.eps[,,,1],doms=doms)
dev.copy2pdf(file="lh_kerneldensity_crit_eu_parameters_od.pdf")
plot.heatmap.lh(lh.two.loc.eps[,,,2],doms=doms)
dev.copy2pdf(file="lh_kerneldensity_crit_sa_parameters_od.pdf")
plot.heatmap.lh(lh.two.loc.eps[,,,3],doms=doms)
dev.copy2pdf(file="lh_kerneldensity_crit_eusa1_parameters_od.pdf")
plot.heatmap.lh(lh.two.loc.eps[,,,4],doms=doms)
dev.copy2pdf(file="lh_kerneldensity_crit_eusa2_parameters_od.pdf")
plot.heatmap.lh(lh.two.loc.eps[,,,3]*lh.two.loc.eps[,,,4],doms=doms)
dev.copy2pdf(file="lh_kerneldensity_crit_eusa_parameters_od.pdf")
plot.heatmap.lh(lh.two.loc.eps[,,,1]*lh.two.loc.eps[,,,2],doms=doms)
dev.copy2pdf(file="/Volumes/Temp/Lukas/Data/A7_pigmentation/Simulation/lh_kerneldensity_crit12_eu_sa_parameters_od.pdf")



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
    abline(v=seq(1,length(effxeps)/3-1)*3+0.5,lw=3)
    abline(h=(1:length(hs2[,1])+0.5))
    abline(h=seq(1,length(hs2[,1])/length(doms)-1)*5+0.5,lw=3)
    axis(1,at=1:length(effxeps),labels=rep(abs(eff2[3:1]),length(eps.int)))
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
    mtext(eps.int,side=3,at=seq(0,length(eps.int)-1)*3+2.5,las=1,line=0.5,cex=1.25)
    mtext("epistatic interaction strength",side=3,at=length(effxeps)/2+0.5,line=2.0,cex=1.5)
                                        #axis(2, at=seq(0,1,,dim(db)[2]), labels=colnames(db))
    plot.rank=array(length(plot.dat)-rank(plot.dat),dim=dim(plot.dat))
    for(i in 1:length(plot.dat[,1])){
        for(j in 1:length(plot.dat[1,])){
            if (plot.rank[i,j] < 5){
                text(j,i,prettyNum(plot.dat[i,j],digits=2,drop0trailing=T),col="red",cex=0.6,font=2)
            }
            else if (plot.rank[i,j] < 10){  #(plot.dat[i,j] >= 1e-9){
                #if plot.dat[i,j] >= plot.dat
                text(j,i,prettyNum(plot.dat[i,j],digits=2,drop0trailing=T),col="red",cex=0.6,font=1)
            }
            else if (plot.rank[i,j] < 50 & plot.dat[i,j] >= 1e-9){
                text(j,i,prettyNum(plot.dat[i,j],digits=2,drop0trailing=T),col="black",cex=0.6,font=1)
            }
        }
    }
    par(mar=c(5,0,5,5.5))
    image.scale(c(1.01,0.01), col=cols[(length(cols)-print.ends):length(cols)], breaks=sort(mybreaks[1:(print.ends+2)]),horiz=FALSE, yaxt="n",xaxt="n",xlab="")
    #print.breaks=c(0,-1,seq(-3,-24,by=-3))
    axis(4,at=log(10^print.breaks),labels=(format(10^print.breaks,digits=2)), las=2)
    mtext("Likelihood", side=4, line=3.5,cex=1.5)
}

epsxeff=expand.grid(eps.int,eff2[1:3])
epsxeff=apply(epsxeff,1,function (x) paste("Eff:EI",x,collapse=":"))
lh.two.loc.eps.rear=aperm(lh.two.loc.eps,c(1,3,2,4))
lh.two.loc.eps.rear=lh.two.loc.eps.rear[,,3:1,]
plot.heatmap.lh.rear(lh.two.loc.eps.rear[,,,1]*lh.two.loc.eps.rear[,,,2])


lh.two.loc.eps.xs.rear=aperm(lh.two.loc.eps.xs,c(1,3,2,4))
lh.two.loc.eps.xs.rear=lh.two.loc.eps.xs.rear[,,3:1,]
plot.heatmap.lh.rear(lh.two.loc.eps.xs.rear[,,,1]*lh.two.loc.eps.xs.rear[,,,2])

plot.heatmap.lh.rear <- function(lh.array)
{
    ## takes rearanged lh array (hs, eps,int, eff)
    plot.dat=array(lh.array,dim=c(length(hs2[,1]),length(epsxeff)), dimnames=list("hs"=rownames(hs2),"epsxeff"=epsxeff))
    cols=c(colorRampPalette(c("beige","bisque","khaki","yellow","greenyellow","green1","green3","lightblue","darkblue"))(30))
    doms=c(0,0.25,0.5,0.75,1)
    ## cols=rev(heat.colors(18))
    ## cols=rev(terrain.colors(50))
    ## cols=rev(topo.colors(50))
    ## cols=colorRampPalette(brewer.pal(9,"Greens"))(50)
    mybreaks=c(log(10^seq(0,-15,length.out=25)),log(1e-20),log(1e-30),log(1e-50),log(1e-99),log(1e-999))
    print.breaks=seq(0,-15,length.out=15)
    layout(matrix(1:2, ncol = 2,byrow=T), widths = rep(7,1), respect = FALSE)
    par(mar=c(5,7,5,1))
    image(1:length(epsxeff),1:length(rownames(hs2)),log(t(plot.dat)),col=cols, breaks=sort(mybreaks),xaxt='n',yaxt='n',xlab="",ylab="")
    abline(v=(1:length(epsxeff)+0.5))
    abline(v=seq(1,length(epsxeff)/length(eps.int)-1)*length(eps.int)+0.5,lw=3)
    abline(h=(1:length(hs2[,1])+0.5))
    abline(h=seq(1,length(hs2[,1])/5-1)*5+0.5,lw=3)
    axis(1,at=1:length(epsxeff),labels=rep(eps.int,length(eff2[3:1])))
    mtext("epistatic interaction strength", side=1, line=2.5,cex=1.5)
    mtext("dominance",side=2, line=5.0,cex=1.25)
    mtext("locus", at=-0.75,side=3, line=1.0,cex=1.25)
    mtext("tan : bab1", at=-0.75,side=3, line=0.0,cex=1.25)
                                        #mtext(c(0,0.25,0.5,0.75,1),side=2,at=seq(0,length(hs2[,1])/5-1)*5+3,las=1,line=4)
    mtext(rep(doms,5),side=2,at=1:length(plot.dat[,1]),line=4,las=1,adj=c(0.5,0.5))
    mtext(":",side=2,at=1:length(plot.dat[,1]),las=1,line=2.75)
                                        #axis(2,at=1:length(hs2[,1]),labels=rep(c(0,0.25,0.5,0.75,1),5),las=1)
    axis(2,at=1:length(hs2[,1]),labels=NA,las=1)
    mtext(sapply(doms,function(x) rep(x,5)),side=2,at=1:length(hs2[,1]),las=1,line=1.5,adj=c(0.5,0.5))
#    mtext(sapply(c(1,0.75,0.5,0.25,0),function(x) rep(x,5)),side=2,at=1:length(hs2[,1]),las=1,line=1.5,adj=c(0.5,0.5))
    mtext(eff2[3:1],side=3,at=seq(0,length(eff2[3:1])-1)*length(eps.int)+3,las=1,line=0.5,cex=1.25)
    mtext("effect of bab1 relative to tan",side=3,at=length(epsxeff)/2+0.5,line=2.0,cex=1.5)
                                        #axis(2, at=seq(0,1,,dim(db)[2]), labels=colnames(db))
    plot.rank=array(length(plot.dat)-rank(plot.dat),dim=dim(plot.dat))
    for(i in 1:length(plot.dat[,1])){
        for(j in 1:length(plot.dat[1,])){
            if (plot.rank[i,j] < 5){
                text(j,i,prettyNum(plot.dat[i,j],digits=2,drop0trailing=T),col="red",cex=0.6,font=2)
            }
            else if (plot.rank[i,j] < 10){  #(plot.dat[i,j] >= 1e-9){
                #if plot.dat[i,j] >= plot.dat
                text(j,i,prettyNum(plot.dat[i,j],digits=2,drop0trailing=T),col="red",cex=0.6,font=1)
            }
            else if (plot.rank[i,j] < 100 & plot.dat[i,j] >= 1e-9){
                text(j,i,prettyNum(plot.dat[i,j],digits=2,drop0trailing=T),col="black",cex=0.6,font=1)
            }
        }
    }
    par(mar=c(5,0,5,5.5))
    image.scale(c(1.01,0.01), col=cols[(length(cols)-18):length(cols)], breaks=sort(mybreaks[1:20]),horiz=FALSE, yaxt="n",xaxt="n",xlab="")
    print.breaks=c(0,-1,seq(-3,-24,by=-3))
    axis(4,at=log(10^print.breaks),labels=(format(10^print.breaks,digits=2)), las=2)
    mtext("Likelihood", side=4, line=3.5,cex=1.5)
}







#0.077476894 0.001882348 0.071866123 0.007402966



sqd.eu.eps*sqd.sa.eps*sqd.eu.sa.eps
sqd.eu.eps+sqd.sa.eps+sqd.eu.sa.eps
comb.p.eu.sa=p.eu.eps[,,,]*p.sa.eps[,,,]*p.eu.sa.eps[,,,]
rank.p=p.eu.eps[,,,]
for (i in 1:length(bin.eu)) {
rank.p[,,,i]=rank(p.eu.eps[,,,i]*p.sa.eps[,,,i]*p.eu.sa.eps[,,,i])
which(comb.p.eu.sa == max(comb.p.eu.sa),arr.ind=T)
}
p.eu[,,"1"]*p.sa[,,"0.5"]
save(two.loc.eps.eu,two.loc.eps.sa,file="two.loc.eu.sa.RData")
setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Simulation")
load("two.loc.eu.sa.RData")
library("RColorBrewer")
display.brewer.all()

## Heatmap
## Do one bin at a time
## 5 heatmaps side by side (eps interactions)
bin.idx="1"
## rearrange into one big matrix
effxeps=expand.grid(eff2,eps.int)
effxeps=apply(effxeps,1,function (x) paste("Eff:EI",x,collapse=":"))
plot.comb.p.eu.sa=array(comb.p.eu.sa,dim=c(length(hs2[,1]),length(effxeps),length(bin.eu) ), dimnames=list("hs"=rownames(hs2),"effxeps"=effxeps,"bins"=bin.eu))
heatmap(plot.comb.p.eu.sa[,,bin.idx],Rowv=NA,Colv=NA)
image(plot.comb.p.eu.sa[,,bin.idx],col=heat.colors(25))
library(reshape2,ggplot2)
source("/Volumes/Temp/Lukas/Tools/Scripts/R/image.scale.R")
plot.dat=plot.comb.p.eu.sa[,,bin.idx]
plot.dat[which(plot.dat==0)]=1e-6
cols=c(colorRampPalette(c("khaki","yellow","greenyellow","green1","green4","blue"))(20))
cols=rev(heat.colors(18))
cols=rev(terrain.colors(50))
cols=rev(topo.colors(50))
cols=colorRampPalette(brewer.pal(9,"Greens"))(50)
mybreaks=c(seq(log(1),log(0.75),length.out=10),seq(log(0.65),log(0.25),length.out=6),seq(log(0.1),log(1e-2),length.out=3),log(1e-4),log(1e-6))
quartz(width=10,height=9)

x11(width=10,height=9)
plot.image.heat(plot.dat,cols,mybreaks)
dev.copy2pdf(file="combined_lh_parameters.pdf")
quartz(width=10,height=9)
plot.image.heat(array(p.eu.eps[,,,bin.idx],dim=c(length(hs2[,1]),length(effxeps) ), dimnames=list("hs"=rownames(hs2),"effxeps"=effxeps)),cols,mybreaks)
dev.copy2pdf(file="crit_eu_lh_parameters.pdf")
quartz(width=10,height=9)
plot.image.heat(array(p.sa.eps[,,,bin.idx],dim=c(length(hs2[,1]),length(effxeps) ), dimnames=list("hs"=rownames(hs2),"effxeps"=effxeps)),cols,mybreaks)
dev.copy2pdf(file="crit_sa_lh_parameters.pdf")
quartz(width=10,height=9)
plot.image.heat(array(p.eu.sa.eps[,,,bin.idx],dim=c(length(hs2[,1]),length(effxeps) ), dimnames=list("hs"=rownames(hs2),"effxeps"=effxeps)),cols,mybreaks)
dev.copy2pdf(file="crit_eu_sa_lh_parameters.pdf")



#function 
plot.image.heat <- function(plot.dat,cols,print.breaks=c(1.0,0.9,0.75,0.5,0.25,0.1,0.01,0.001)){
length(mybreaks)
exp(mybreaks)
length(cols)
plot.dat[which(plot.dat==0)]=1e-6
layout(matrix(1:2, ncol = 2,byrow=T), widths = rep(7,1), respect = FALSE)
par(mar=c(5,7,5,1))
image(1:length(effxeps),1:length(rownames(hs2)),log(t(plot.dat)),col=cols, breaks=sort(mybreaks),xaxt='n',yaxt='n',xlab="",ylab="")
#image(1:length(effxeps),1:length(rownames(hs2)),t(plot.dat),col=cols, breaks=sort(mybreaks),xaxt='n',yaxt='n',xlab="",ylab="")
abline(v=(1:length(effxeps)+0.5))
abline(v=seq(1,length(effxeps)/4-1)*4+0.5,lw=3)
abline(h=(1:length(hs2[,1])+0.5))
abline(h=seq(1,length(hs2[,1])/5-1)*5+0.5,lw=3)
axis(1,at=1:length(effxeps),labels=rep(eff2,length(eps.int)))
mtext("effect of locus 2 relative to locus 1", side=1, line=2.5,cex=1.5)
mtext("dominance",side=2, line=5.0,cex=1.25)
mtext("locus", at=-1,side=3, line=1.0,cex=1.25)
mtext("1   :   2", at=-1,side=3, line=0.0,cex=1.25)
#mtext(c(0,0.25,0.5,0.75,1),side=2,at=seq(0,length(hs2[,1])/5-1)*5+3,las=1,line=4)
mtext(rep(c(0,0.25,0.5,0.75,1),5),side=2,at=1:length(plot.dat[,1]),line=4,las=1,adj=c(0.5,0.5))
mtext(":",side=2,at=1:length(plot.dat[,1]),las=1,line=2.75)
#axis(2,at=1:length(hs2[,1]),labels=rep(c(0,0.25,0.5,0.75,1),5),las=1)
axis(2,at=1:length(hs2[,1]),labels=NA,las=1)
mtext(sapply(c(0,0.25,0.5,0.75,1),function(x) rep(x,5)),side=2,at=1:length(hs2[,1]),las=1,line=1.5,adj=c(0.5,0.5))
mtext(eps.int,side=3,at=seq(0,length(eps.int)-1)*4+2.5,las=1,line=0.5,cex=1.25)
mtext("epistatic interaction strength",side=3,at=length(effxeps)/2+0.5,line=2.0,cex=1.5)
#axis(2, at=seq(0,1,,dim(db)[2]), labels=colnames(db))
for(i in 1:length(plot.dat[,1])){
for(j in 1:length(plot.dat[1,])){
    if (plot.dat[i,j] >= 0.8){
text(j,i,round(plot.dat[i,j],2),col="white",cex=0.8)
     }
    else if (plot.dat[i,j] >= 0.01){
text(j,i,round(plot.dat[i,j],2),col="red",cex=0.8)
     }
  }
}
par(mar=c(5,0,5,5.5))
image.scale(c(1,0), col=cols[(length(cols)-15):length(cols)], breaks=sort(mybreaks)[(length(cols)-16):length(cols)],horiz=FALSE, yaxt="n",xaxt="n",xlab="")
print.breaks=c(1.0,0.9,0.75,0.5,0.25,0.1,0.01,0.001)
axis(4,at=log(sort(print.breaks)),labels=(sort(print.breaks)), las=2)
mtext("Likelihood", side=4, line=3.5,cex=1.5)
}

dev.copy2pdf(file="combined_lh_parameters.pdf")







plot.dat["0:0.25",]
plot.df=melt(plot.comb.p.eu.sa)
plot.df$s.value=plot.df$value/max(plot.df$value)
p=NA
ggplot(plot.df, aes( effxeps,hs)) + geom_tile(aes(fill = value), colour = "black")  + scale_fill_gradientn(colours=cols,values=c(1e-10,sort(mybreaks),1),trans="log",na.value="white",guide="legend",breaks= mybreaks, labels=round(mybreaks,2),name="joined P") + scale_x_discrete(labels=rep(eff2,length(eps.int))) + xlab("relative effect of locus 2") + ylab("dominance locus1:locus2")+ geom_text(map=aes(x=effxeps,y=hs,label=format(value,digits=3)),color="red")
#+ theme(plot.margin=unit(c(4,1,1,1),"lines")) + geom_text(map=aes(x = seq(2.5,2.5+4*(length(eps.int)-1),by=4), y = rep(Inf,5),label = paste("epist. int.:",eps.int)))
#+ geom_text(map=aes(x=effxeps,y=hs),label=round(plot.df$value,3),size=3,color="red")
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

# + geom_text(aes(x=effxeps,y=hs,label = round(value,3)),size=3,color="red")
#par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
for  (i in as.character(eps.int)){
    image(plot.comb.p.eu.sa[,,bin.idx],col=heat.colors(16))
    ## if (i==eps.int[1]){
    ##     axis(2,cex.axis=1.5,lables=rownames(hs2))
    ## }

}


### plot two locus effects
quartz(width=12,height=7)
layout(matrix(1:8, ncol = 4,byrow=T), widths = rep(1,3), respect = FALSE)
#par(mar = c(4.1, 0.0, 0.0, 2.1))
par(oma=c(5.1,5.1,4.1,3.1),mar=c(0,0,0,0))
H_levs=c("0.5:0","0.5:0.5","0.5:1")
H_labs=gsub("1","d",gsub("0","r",gsub("0.5","a",H_levs)))
ylimits=c(0.25,-4)
for (i in as.character(eff2)){
    d.s=0.04
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(log(t(two.loc.eps.eu[H_levs,i,"0","or",1,])), at=hs2[H_levs,2]-d.s,boxwex=0.04, col="red",ylim=ylimits,xlim=c(-0.2,1.2), yaxt="n",xaxt="n")
    boxplot(log(t(two.loc.eps.xs[H_levs,i,"0","or",1,])), at=hs2[H_levs,2],boxwex=0.04, col="green",ylim=ylimits,xlim=c(-0.2,1.2), yaxt="n",cex.axis=1.5,xaxt="n",add=T)
    #mtext("h",side=1,outer=F,line=2.5)
    mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)
    boxplot(log(t(two.loc.eps.sa[H_levs,i,"0","or",1,])), at=hs2[H_levs,2]+d.s,boxwex=0.04, col="deepskyblue",ylim=ylimits,xlim=c(-0.2,1.2), yaxt="n",xaxt="n",add=T)
    if (i==eff2[1]){
        axis(2,cex.axis=1.5)
    }
}
legend("bottomright",legend=c("EU","SA (ext)","SA"),fill=c("red","green","deepskyblue"),cex=1.5)
ylimits=c(0,5)
for (i in as.character(eff2)){
    d.s=0.04
    #par(mar = c(4.1, 0.0 + (i==1)*4.1, 4.1,0.0 + (i==length(Fs)) * 4.1))
    boxplot(log(t(two.loc.eps.eu[H_levs,i,"0","or",2,])), at=hs2[H_levs,2]-d.s,boxwex=0.04, col="red",ylim=ylimits,xlim=c(-0.2,1.2), yaxt="n",xaxt="n",add=F)
    axis(1,at=hs2[H_levs,2],cex.axis=1.5,labels=H_labs)
    boxplot(log(t(two.loc.eps.xs[H_levs,i,"0","or",2,])), at=hs2[H_levs,2],boxwex=0.04, col="green",ylim=ylimits,xlim=c(-0.2,1.2), yaxt="n",cex.axis=1.5,xaxt="n",add=T)
    mtext("h",side=1,outer=F,line=2.5)
    #mtext(bquote(eff[2]~"="~.(i) %*% eff[1]), side=3, line=1.0)
    boxplot(log(t(two.loc.eps.sa[H_levs,i,"0","or",2,])), at=hs2[H_levs,2]+d.s,boxwex=0.04, col="deepskyblue",ylim=ylimits,xlim=c(-0.2,1.2), yaxt="n",xaxt="n",add=T)
    if (i==eff2[1]){
        axis(2,cex.axis=1.5)
    }
}
mtext(expression(log(OR)),side=2,outer=T,line=2.0)
mtext(c("locus 2","locus 1"),at=c(0.25,0.75),side=4,outer=T,line=1.0)
dev.copy2pdf(file="logOR_comp_eu_sa_xs.pdf")


run.simulations.rnd= function(n.runs=50,h2=0.25,total.ind =2100,n.replicates=3,n.pools=100,afs.vals=c(0.5),h.vals=c(0.5),F.vals=c(0.0), random.cov=FALSE) {
#### Runs n simulations assigning afs, h and Fit values randomly for each run according to given arrays
#### ir random.cov is set to a number, a poisson distr random coverage around that umber is chosen and the allelel counts are sampled from a binomial
#######  Initialize
     d.replicates=total.ind/n.replicates # ind. per replicate
    ## odds.ratio_ld <- matrix(0,nrow=n.runs,ncol=n.loci)
    ## odds.ratio_l <- matrix(0,nrow=n.runs,ncol=n.loci)
    ## odds.ratio_d <- matrix(0,nrow=n.runs,ncol=n.loci)
                                    #log.odds.ratio_ajb <- matrix(0,nrow=n.runs,ncol=n.loci)
                                        #cors <- matrix(0,nrow=n.runs,ncol=4)
    ## af_avg <- matrix(0,nrow=n.runs,ncol=n.loci)
    ## af_light <- matrix(0,nrow=n.runs,ncol=n.loci)
    ## af_dark <- matrix(0,nrow=n.runs,ncol=n.loci)    
    ## setup local afs, h and values
    h.loci.loc <- vector("list",length=n.runs)
    n.freqs.loc <- vector("list",length=n.runs)
    F.inbreed.loc <- vector("list",length=n.runs)
    v.E.loc <- vector("list",length=n.runs)
    results = foreach (i=1:n.runs,.combine="rbind") %dopar% {
        ## get random afs, h and F values
        h.loci.loc[[i]] <- sample(h.vals, n.loci, replace = TRUE)
        n.freqs.loc[[i]] <- sample(afs.vals, n.loci, replace = TRUE)
        F.inbreed.loc[[i]] <- sample(F.vals, 1)
        v.E.loc[[i]] <- get.V.E(n.freqs.loc[[i]],eff,h.loci.loc[[i]],h2,eps.int=0.0)
        ## generate genotypes randomly (independent loci)
### alleles A a
        gt.freqs.loc = sapply(n.freqs.loc[[i]],function (p) c(p^2 + (p-p^2)*F.inbreed.loc[[i]],2*(p-p^2)*(1-F.inbreed.loc[[i]]),(1-p)^2 + (p-p^2)*F.inbreed.loc[[i]])) ### genotype frequencies for each locus, aa, aA, AA
        mean_eff=sum(gt.vals * gt.freqs.loc %*% eff)

        pos <- c()
                                        #pos <- replicate(total.ind,get.gt(gt.freqs,gt.vals))
                                        #pos <- matrix(pos,nrow=n.loci,ncol=total.ind,byrow=T)
        pos <- apply(gt.freqs.loc,2,function (x) get.gt.locus(x,gt.vals,total.ind))
        pos <- matrix(pos,nrow=n.loci,ncol=total.ind,byrow=T)  
        ## pigmentation phenotypes resulting from genotypes
        ## additive effects (eff centered around 0)
        pheno <- eff %*% (apply(pos,2,function(x) loc.eff(h.loci.loc[[i]],x)) - 0.5)
        ## environmental effect
        pheno <- pheno + rnorm(total.ind,0,sqrt(v.E.loc[[i]]))
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
            o.pheno <- sample.indices[order(sampled.pheno)]
            light   <- o.pheno[1:n.pools]   # take "n.pools" lightest colored individuals for light pool
            dark    <- o.pheno[((d.replicates-n.pools)+1):d.replicates] # take "n.pools" darkest individuals for dark pool
            alllight=c(alllight, light)
            alldark=c(alldark, dark)
                                        # par(mfrow=c(1,1)); b.points=seq(min(pheno[sample.indices]),max(pheno[sample.indices]),by=0.5); hist(pheno[sample.indices],breaks=b.points); hist(pheno[light], col="yellow",add=T,breaks=b.points); hist(pheno[dark], col="blue",add=T,breaks=b.points)
           
            ## genotypes for the pools
            geno.light  <-  pos[,light]
            geno.dark   <-  pos[,dark]
            ## allele counts
            tot.counts=n.pools*2
            tot.count.l=tot.counts
            tot.count.d=tot.counts
            ac.light=rowSums(geno.light)
            ac.dark =rowSums(geno.dark)
            ac.base =rowSums(pos[,sample.indices])
            ac.rep=rowSums(pos[,sample.indices])
            if (random.cov != FALSE){
                tot.counts.l=rpois(1,random.cov)
                tot.counts.d=rpois(1,random.cov)
                ac.light=sapply(ac.light,function (x) rbinom(1,tot.counts.l,x/tot.counts))
                ac.dark=sapply(ac.dark,function (x) rbinom(1,tot.counts.d,x/tot.counts))
            }
            dat[,,j]=array(c(ac.light,ac.dark, tot.count.l-ac.light, tot.count.d-ac.dark), dim=c( n.loci,4))
            dat_l[,,j]=array(c(ac.light,ac.base, tot.count.l-ac.light,d.replicates*2-ac.base), dim=c( n.loci,4))
            dat_d[,,j]=array(c(ac.dark,ac.base, tot.count.d-ac.dark,d.replicates*2-ac.base), dim=c( n.loci,4))
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
        return(c(p.val,odds.ratio_ld)) #,odds.ratio_l,odds.ratio_d,af_avg,af_light,af_dark))
    }
    
    ## return all tables in list structure
    l.names=c("pv","or") #,"ol","od","af.m","af.l","af.d")
    ret.list=list()
    for (i in 1:length(l.names)){
        ret.list[[l.names[i]]]=results[,((i-1)*n.loci+1):(i*n.loci)]
    }
    return(ret.list)    
}

eff=seq(-1,1,by=0.125)
n.loci=length(eff)
afs.vals=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
h.vals=c(0,0.25,0.5,0.75,1.0)
F.vals=c(0.0,0.25,0.5)
n.runs=5000
h2=0.25
res.afs=run.simulations.rnd(n.runs=n.runs,h2=h2,total.ind =2100,n.replicates=3,n.pools=100,afs.vals=afs.vals,h.vals=c(0.5),F.vals=c(0.0))

res.afs$lor <- log(res.afs$or)
res.afs$lpv <- -log10(res.afs$pv)
# get spearman correlations
fields=c("sp.lor","p.lor","sp.pv","p.pv")

res.afs$sp.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs$lor[x,],eff,method="spearman")))
res.afs$p.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs$lor[x,],eff,method="pearson")))
res.afs$sp.pv=sapply(1:n.runs,function (x) cor(res.afs$lpv[x,],abs(eff),method="spearman"))
res.afs$p.pv=sapply(1:n.runs,function (x) cor(res.afs$lpv[x,],abs(eff),method="pearson"))
boxplot(res.afs[fields])
n.runs=5000
res.afs.rc=run.simulations.rnd(n.runs=n.runs,h2=h2,total.ind =2100,n.replicates=3,n.pools=100,afs.vals=afs.vals,h.vals=c(0.5),F.vals=c(0.0),random.cov=50)
res.afs.rc$lor <- log(res.afs.rc$or)
res.afs.rc$lpv <- -log10(res.afs.rc$pv)
res.afs.rc$sp.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.rc$lor[x,],eff,method="spearman")))
res.afs.rc$p.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.rc$lor[x,],eff,method="pearson")))
res.afs.rc$sp.pv=sapply(1:n.runs,function (x) cor(res.afs.rc$lpv[x,],abs(eff),method="spearman"))
res.afs.rc$p.pv=sapply(1:n.runs,function (x) cor(res.afs.rc$lpv[x,],abs(eff),method="pearson"))
boxplot(res.afs.rc[fields])


res.afs.h=run.simulations.rnd(n.runs=n.runs,h2=h2,total.ind =2100,n.replicates=3,n.pools=100,afs.vals=afs.vals,h.vals=h.vals,F.vals=c(0.0))
res.afs.h$lor <- log(res.afs.h$or)
res.afs.h$lpv <- -log10(res.afs.h$pv)
res.afs.h$sp.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.h$lor[x,],eff,method="spearman")))
res.afs.h$p.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.h$lor[x,],eff,method="pearson")))
res.afs.h$sp.pv=sapply(1:n.runs,function (x) cor(res.afs.h$lpv[x,],abs(eff),method="spearman"))
res.afs.h$p.pv=sapply(1:n.runs,function (x) cor(res.afs.h$lpv[x,],abs(eff),method="pearson"))

res.afs.h.rc=run.simulations.rnd(n.runs=n.runs,h2=h2,total.ind =2100,n.replicates=3,n.pools=100,afs.vals=afs.vals,h.vals=h.vals,F.vals=c(0.0),random.cov=50)
res.afs.h.rc$lor <- log(res.afs.h.rc$or)
res.afs.h.rc$lpv <- -log10(res.afs.h.rc$pv)
res.afs.h.rc$sp.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.h.rc$lor[x,],eff,method="spearman")))
res.afs.h.rc$p.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.h.rc$lor[x,],eff,method="pearson")))
res.afs.h.rc$sp.pv=sapply(1:n.runs,function (x) cor(res.afs.h.rc$lpv[x,],abs(eff),method="spearman"))
res.afs.h.rc$p.pv=sapply(1:n.runs,function (x) cor(res.afs.h.rc$lpv[x,],abs(eff),method="pearson"))

## res.afs.F=run.simulations.rnd(n.runs=n.runs,h2=h2,total.ind =2100,n.replicates=3,n.pools=100,afs.vals=afs.vals,h.vals=c(0.0),F.vals=F.vals)
## res.afs.F$lor <- log(res.afs.F$or)
## res.afs.F$lpv <- -log10(res.afs.F$pv)
## res.afs.F$sp.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.F$lor[x,],eff,method="spearman")))
## res.afs.F$p.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.F$lor[x,],eff,method="pearson")))
## res.afs.F$sp.pv=sapply(1:n.runs,function (x) cor(res.afs.F$lpv[x,],abs(eff),method="spearman"))
## res.afs.F$p.pv=sapply(1:n.runs,function (x) cor(res.afs.F$lpv[x,],abs(eff),method="pearson"))

res.afs.hF25=run.simulations.rnd(n.runs=n.runs,h2=h2,total.ind =2100,n.replicates=3,n.pools=100,afs.vals=afs.vals,h.vals=h.vals,F.vals=c(0.25))
res.afs.hF25$lor <- log(res.afs.hF25$or)
res.afs.hF25$lpv <- -log10(res.afs.hF25$pv)
res.afs.hF25$sp.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.hF25$lor[x,],eff,method="spearman")))
res.afs.hF25$p.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.hF25$lor[x,],eff,method="pearson")))
res.afs.hF25$sp.pv=sapply(1:n.runs,function (x) cor(res.afs.hF25$lpv[x,],abs(eff),method="spearman"))
res.afs.hF25$p.pv=sapply(1:n.runs,function (x) cor(res.afs.hF25$lpv[x,],abs(eff),method="pearson"))

res.afs.hF25.rc=run.simulations.rnd(n.runs=n.runs,h2=h2,total.ind =2100,n.replicates=3,n.pools=100,afs.vals=afs.vals,h.vals=h.vals,F.vals=c(0.25),random.cov=50)
res.afs.hF25.rc$lor <- log(res.afs.hF25.rc$or)
res.afs.hF25.rc$lpv <- -log10(res.afs.hF25.rc$pv)
res.afs.hF25.rc$sp.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.hF25.rc$lor[x,],eff,method="spearman")))
res.afs.hF25.rc$p.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.hF25.rc$lor[x,],eff,method="pearson")))
res.afs.hF25.rc$sp.pv=sapply(1:n.runs,function (x) cor(res.afs.hF25.rc$lpv[x,],abs(eff),method="spearman"))
res.afs.hF25.rc$p.pv=sapply(1:n.runs,function (x) cor(res.afs.hF25.rc$lpv[x,],abs(eff),method="pearson"))

res.afs.hF5=run.simulations.rnd(n.runs=n.runs,h2=h2,total.ind =2100,n.replicates=3,n.pools=100,afs.vals=afs.vals,h.vals=h.vals,F.vals=c(0.5))
res.afs.hF5$lor <- log(res.afs.hF5$or)
res.afs.hF5$lpv <- -log10(res.afs.hF5$pv)
res.afs.hF5$sp.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.hF5$lor[x,],eff,method="spearman")))
res.afs.hF5$p.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.hF5$lor[x,],eff,method="pearson")))
res.afs.hF5$sp.pv=sapply(1:n.runs,function (x) cor(res.afs.hF5$lpv[x,],abs(eff),method="spearman"))
res.afs.hF5$p.pv=sapply(1:n.runs,function (x) cor(res.afs.hF5$lpv[x,],abs(eff),method="pearson"))

res.afs.hF5.rc=run.simulations.rnd(n.runs=n.runs,h2=h2,total.ind =2100,n.replicates=3,n.pools=100,afs.vals=afs.vals,h.vals=h.vals,F.vals=c(0.5),random.cov=50)
res.afs.hF5.rc$lor <- log(res.afs.hF5.rc$or)
res.afs.hF5.rc$lpv <- -log10(res.afs.hF5.rc$pv)
res.afs.hF5.rc$sp.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.hF5.rc$lor[x,],eff,method="spearman")))
res.afs.hF5.rc$p.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.hF5.rc$lor[x,],eff,method="pearson")))
res.afs.hF5.rc$sp.pv=sapply(1:n.runs,function (x) cor(res.afs.hF5.rc$lpv[x,],abs(eff),method="spearman"))
res.afs.hF5.rc$p.pv=sapply(1:n.runs,function (x) cor(res.afs.hF5.rc$lpv[x,],abs(eff),method="pearson"))

## res.afs.hF=run.simulations.rnd(n.runs=n.runs,h2=h2,total.ind =2100,n.replicates=3,n.pools=100,afs.vals=afs.vals,h.vals=h.vals,F.vals=F.vals)
## res.afs.hF$lor <- log(res.afs.hF$or)
## res.afs.hF$lpv <- -log10(res.afs.hF$pv)
## res.afs.hF$sp.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.hF$lor[x,],eff,method="spearman")))
## res.afs.hF$p.lor=sapply(1:n.runs,function (x) -1*(cor(res.afs.hF$lor[x,],eff,method="pearson")))
## res.afs.hF$sp.pv=sapply(1:n.runs,function (x) cor(res.afs.hF$lpv[x,],abs(eff),method="spearman"))
## res.afs.hF$p.pv=sapply(1:n.runs,function (x) cor(res.afs.hF$lpv[x,],abs(eff),method="pearson"))

wilcox.test(res.afs$p.lor,res.afs$p.pv)$p.value
t.test(res.afs$sp.lor,res.afs$sp.pv)
wilcox.test(res.afs.h$p.lor,res.afs.h$p.pv)
wilcox.test(res.afs.F$p.lor,res.afs.F$p.pv)
wilcox.test(res.afs.hF$p.lor,res.afs.hF$p.pv)$p.value

lor.cors=list(res.afs$sp.lor,res.afs$p.lor,res.afs.h$sp.lor,res.afs.h$p.lor,res.afs.F$sp.lor,res.afs.F$p.lor,res.afs.hF$sp.lor,res.afs.hF$p.lor)
lpv.cors=list(res.afs$sp.pv,res.afs$p.pv,res.afs.h$sp.pv,res.afs.h$p.pv,res.afs.F$sp.pv,res.afs.F$p.pv,res.afs.hF$sp.pv,res.afs.hF$p.pv)
lor.cors=list(res.afs$sp.lor,res.afs$p.lor,res.afs.h$sp.lor,res.afs.h$p.lor,res.afs.hF25$sp.lor,res.afs.hF25$p.lor,res.afs.hF5$sp.lor,res.afs.hF5$p.lor)
lpv.cors=list(res.afs$sp.pv,res.afs$p.pv,res.afs.h$sp.pv,res.afs.h$p.pv,res.afs.hF25$sp.pv,res.afs.hF25$p.pv,res.afs.hF5$sp.pv,res.afs.hF5$p.pv)
lor.cors.rc=list(res.afs.rc$sp.lor,res.afs.rc$p.lor,res.afs.h.rc$sp.lor,res.afs.h.rc$p.lor,res.afs.hF25.rc$sp.lor,res.afs.hF25.rc$p.lor,res.afs.hF5.rc$sp.lor,res.afs.hF5.rc$p.lor)
lpv.cors.rc=list(res.afs.rc$sp.pv,res.afs.rc$p.pv,res.afs.h.rc$sp.pv,res.afs.h.rc$p.pv,res.afs.hF25.rc$sp.pv,res.afs.hF25.rc$p.pv,res.afs.hF5.rc$sp.pv,res.afs.hF5.rc$p.pv)

x11()
posts=unlist(as.vector(sapply(1:4,function (x) c(x-0.2,x+0.2))))
bxwd=0.1
boxplot(lor.cors,boxwex=bxwd,at=posts-bxwd/2,xlim=c(0.5,4.5),ylim=c(0,1), xaxt='n',col=c("cyan"),outline=F,xlab='',ylab=expression(paste("Spearman/Pearson correlation ",rho)),add=F,notch=T,main="exact allele counts")
boxplot(lpv.cors,boxwex=bxwd,at=posts+bxwd/2,xlim=c(0,5),ylim=c(0,1), xaxt='n',col=c("red"),outline=F,xlab="",ylab="",add=T,notch=T)
legend("bottomleft",c(expression(paste(log(OR))),expression(paste(-log[10](P)))),fill=c("cyan","red"))
axis(1,at=posts,labels=rep(c("Pear","Spear"),4),cex.axis=0.75)
mtext(c("AF","AF,h",expression("AF,h,"),expression("AF,h,")),side=1,font=1,at=1:4,line=2.25)
mtext(c("","",expression(F[IT]~"=0.25"),expression(F[IT]~"=0.5")),side=1,font=1,at=1:4,line=3.25)
dev.copy2pdf(file="correlations_lorpval.pdf")
x11()
posts=unlist(as.vector(sapply(1:4,function (x) c(x-0.2,x+0.2))))
bxwd=0.1
boxplot(lor.cors.rc,boxwex=bxwd,at=posts-bxwd/2,xlim=c(0.5,4.5),ylim=c(0,1), xaxt='n',col=c("cyan"),outline=F,xlab='',ylab=expression(paste("Spearman/Pearson correlation ",rho)),add=F,notch=T,main="simulated pooled sequencing (mean coverage:50)")
boxplot(lpv.cors.rc,boxwex=bxwd,at=posts+bxwd/2,xlim=c(0,5),ylim=c(0,1), xaxt='n',col=c("red"),outline=F,xlab="",ylab="",add=T,notch=T)
legend("bottomright",c(expression(paste(log(OR))),expression(paste(-log[10](P)))),fill=c("cyan","red"))
axis(1,at=posts,labels=rep(c("Pear","Spear"),4),cex.axis=0.75)
mtext(c("AF","AF,h",expression("AF,h,"),expression("AF,h,")),side=1,font=1,at=1:4,line=2.25)
mtext(c("","",expression(F[IT]~"=0.25"),expression(F[IT]~"=0.5")),side=1,font=1,at=1:4,line=3.25)
dev.copy2pdf(file="correlations_lorpval_cov50.pdf")


x11(width=10,height=6)
layout(matrix(1:2, ncol = 2,byrow=T), widths = rep(1,2), respect = FALSE)
#par(oma=c(0,0,0,0),mar=c(0,0,0,0))
par(mar=c(5.1,3.1,4.1,0))
posts=unlist(as.vector(sapply(1:4,function (x) c(x-0.2,x+0.2))))
bxwd=0.1
boxplot(lor.cors,boxwex=bxwd,at=posts-bxwd/2,xlim=c(0.5,4.5),ylim=c(0,1), xaxt='n',col=c("cyan"),outline=F,xlab='',ylab="",add=F,notch=T,cex.axis=0.75)
boxplot(lpv.cors,boxwex=bxwd,at=posts+bxwd/2,xlim=c(0,5),ylim=c(0,1), xaxt='n', yaxt='n',col=c("red"),outline=F,xlab="",ylab="",add=T,notch=T)
legend("bottomleft",c(expression(paste(log(OR))),expression(paste(-log[10](P)))),fill=c("cyan","red"))
axis(1,at=posts,labels=rep(c("Pear","Spear"),4),cex.axis=0.75)
mtext(c("AF","AF,h",expression("AF,h,"),expression("AF,h,")),side=1,font=1,at=1:4,line=2.25)
mtext(expression(paste("Spearman/Pearson correlation (",rho,")")),side=2,line=1.75)
mtext(c("","",expression(F[IT]~"=0.25"),expression(F[IT]~"=0.5")),side=1,font=1,at=1:4,line=3.25)
mtext("exact allele counts",side=3,line=1)
par(mar=c(5.1,0,4.1,3.1))
boxplot(lor.cors.rc,boxwex=bxwd,at=posts-bxwd/2,xlim=c(0.5,4.5),ylim=c(0,1), yaxt="n",xaxt='n',col=c("cyan"),outline=F,xlab='',ylab="",add=F,notch=T)
boxplot(lpv.cors.rc,boxwex=bxwd,at=posts+bxwd/2,xlim=c(0,5),ylim=c(0,1), yaxt="n",xaxt='n',col=c("red"),outline=F,xlab="",ylab="",add=T,notch=T)
#legend("bottomright",c(expression(paste(log(OR))),expression(paste(-log[10](P)))),fill=c("cyan","red"))
axis(1,at=posts,labels=rep(c("Pear","Spear"),4),cex.axis=0.75)
axis(4,cex.axis=0.75)
mtext(c("AF","AF,h",expression("AF,h,"),expression("AF,h,")),side=1,font=1,at=1:4,line=2.25)
mtext(c("","",expression(F[IT]~"=0.25"),expression(F[IT]~"=0.5")),side=1,font=1,at=1:4,line=3.25)
mtext("simulated pooled sequencing",side=3,line=1.5)
mtext("mean coverage:50",side=3,line=0.5)
dev.copy2pdf(file="correlations_lorpval_all.pdf")




run_sims_extremes = function(n.runs=50,v.E=0.0,total.ind =700,n.pools=c(65,67,108,68)) {
  
  #######  Initialize
  gt.freqs = sapply(n.freqs,function (p) c(p^2 + (p-p^2)*F.inbreed,2*(p-p^2)*(1-F.inbreed),(1-p)^2 + (p-p^2)*F.inbreed)) 
  mean_eff=sum(gt.vals * gt.freqs %*% eff)
  results = foreach (i=1:n.runs,.combine="rbind") %do% {
    pos <- c()
    pos <- apply(gt.freqs,2,function (x) get.gt.locus(x,gt.vals,total.ind))
    pos <- matrix(pos,nrow=n.loci,ncol=total.ind,byrow=T)  
    pheno <- eff %*% (apply(pos,2,function(x) loc.eff(h.loci,x)) - 0.5)
    ## environmental effect
    pheno <- pheno + rnorm(total.ind,0,sqrt(v.E))
    ## order indices by pheno and take extreme ones
    o.pheno <- order(pheno,decreasing=T)
    vlight   <- o.pheno[1:n.pools[1]]   # take "n.pools" lightest colored individuals for light pool
    light   <- o.pheno[(n.pools[1]+1):(n.pools[1]+n.pools[2])]   # take  individuals for light pool
    dark    <- o.pheno[(total.ind-n.pools[4]-n.pools[3]):(total.ind-n.pools[4])] # takeindividuals for dark pool
    vdark    <- o.pheno[(total.ind-n.pools[4]+1):total.ind] # take darkest individuals for very dark pool
    ## genotypes for the pools
    geno.vlight  <-  pos[1,vlight]
    geno.vdark   <-  pos[1,vdark]
    geno.light  <-  pos[1,light]
    geno.dark   <-  pos[1,dark]
    ## allele counts       
    ac.light=sum(geno.light)
    ac.dark =sum(geno.dark)
    ac.vlight=sum(geno.vlight)
    ac.vdark =sum(geno.vdark)
    return(c(ac.vlight,ac.light,ac.dark,ac.vdark)/(2*n.pools))
  }
  
  return(c(colMeans(results),apply(results,2,var)))    
}

n.loci<-5
n.freqs <- rep(0.53,n.loci)
h2=0.15
h.loci=c(1.0,rep(0.5,n.loci-1))
eff=c(2.0,rep(0.1,n.loci-1))
v.E=get.V.E(n.freqs,eff,h.loci,h2)
F.inbreed=0.25
hs=c(0,0.5,1.0)
pch=c(0,1,2)
for(i in 1:3){
  h.loci[1]=hs[i]
  v.E=get.V.E(n.freqs,eff,h.loci,h2)
  results=run_sims_extremes(1000,v.E)
  points(rev(xs),results[1:4],pch=pch[i])
  lines(rev(xs),results[1:4],lty=3)
  arrows(rev(xs),unlist(results[1:4]-results[5:8]), rev(xs), unlist(results[1:4]+results[5:8]), length=0.025, angle=90, code=3, lty=3)
}
legend("topleft",c("0.0","0.5","1.0"),pch=pch,cex=0.75,lty=3,title=c("h light allele"))
dev.copy2pdf(file="bab_best4_dom_sims_afs.pdf")

# tan
n.freqs[1]=0.20
h2=0.15
h.loci=c(1.0,rep(0.5,n.loci-1))
eff=c(1.0,rep(0.1,n.loci-1))
v.E=get.V.E(n.freqs,eff,h.loci,h2)
F.inbreed=0.25
hs=c(0,0.5,1.0)
pch=c(0,1,2)
for(i in 1:3){
  h.loci[1]=hs[i]
  v.E=get.V.E(n.freqs,eff,h.loci,h2)
  results=run_sims_extremes(1000,v.E)
  points(rev(xs),results[1:4],pch=pch[i])
  lines(rev(xs),results[1:4],lty=3)
  arrows(rev(xs),unlist(results[1:4]-results[5:8]), rev(xs), unlist(results[1:4]+results[5:8]), length=0.025, angle=90, code=3, lty=3)
}
legend("bottomright",c("0.0","0.5","1.0"),pch=pch,cex=0.75,lty=3,title=c("h light allele"))
dev.copy2pdf(file="tan_best4_dom_sims_afs.pdf")
