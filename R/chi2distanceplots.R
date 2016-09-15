get_chisq_dist=function(obs, null,bins){
	l=length(bins)
	cat("l=", l, "\n")
	o=hist (obs, breaks=bins, plot=F)$counts
	cat("o=", length(o), "\n")
	bins=c(bins, max(null, na.rm=T))
	bins=unique(bins)
	e=hist(null, breaks=bins, plot=F)$counts
	e=e[1:(l-1)]
	cat("e=", length(e), "\n")
	diff=((o-e)^2)/e
	return(diff)
}

get_ab_ratio=function(c_pop){
	#calculates ratio of beta to alpha to obtain mean
	mean_cov=mean(c_pop)
	return((1-mean_cov)/mean_cov)
}

create_histdist=function(obs_pvals){
	#print(length(obs_pvals))
	bins = hist(obs_pvals, plot=F, breaks=1000)$breaks

	l=length(bins)-1
	histdist=data.frame(bins=bins[1:l])
	cnames=c("bin")
	allfiles = system("ls *nullP.out", intern=T)


	for (infile in allfiles)
		{	
#	parc=strsplit(infile,".0_nullP.out")[[1]][1]
#	parc=tail(strsplit(parc,"_")[[1]],1)
		parc=sub(".*out[_]?(\\d+)\\.0_nullP.out","\\1",infile,perl=TRUE)
		print (parc)
		cat("infile", infile, "parc", parc, "\n")
		alph = as.numeric(parc)
		ps=scan(infile)
#		cnames=append(cnames, paste("alph_", alph, sep=""))	
		cnames=append(cnames,alph)
		histdist = cbind(histdist,alph=get_chisq_dist(obs_pvals, ps, bins)[1:l])
		}
	names(histdist)=cnames
	return(histdist)
}

# sort columns by alpha value

setwd=("/Volumes/Temp/Lukas/Data/Vienna_2011/Polymorphisms_Filtered/CMH/FDR")
obs_pvals=scan("fem_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_masked_mc8_only_chroms.cmhout.obsPval")
histdist=create_histdist(obs_pvals)
histdist=histdist[c(names(histdist)[1],sort(as.numeric(names(histdist)[-1])))]
quartz()

x=1:(dim(histdist)[2]-1)
y=apply(histdist[,-1],2,sum,na.rm=T)
plot(x,y, xlab=expression(alpha),ylab=expression(paste(Chi^2," distance")), xaxt="n")
axis(1, at=x, labels=F) #names(histdist)[-1]) 
text(x,par("usr")[3],pos=1,offset=1.5,srt=45,labels=names(histdist)[-1],xpd=TRUE)
title(main=expression(paste("Vienna 2011 ",italic("D. sim")," filtered ")))
title(main=expression(paste(Chi^2," dist. of obs. to ", beta," distr. shuffled P values")))


dist_Pval=scan("fem_vie2011_pI25_pII25_pIII25_pIp_pIIp_pIIIp_filtered_q20_masked_mc8_only_chroms.cmhout60.0_nullP.out")
nl_oP = -1 * log10(obs_pvals)
nl_dP1 = -1 * log10(dist_Pval)
par(bg="white")
qqplot(nl_dP1,nl_oP,xlab=expression(-log[10](P[null])),ylab=expression(-log[10](P[obs])))
abline(0,1,col="red",lty=2)
title(main=expression(paste("Vienna 2011, filtered, ",alpha,"=60")))
dev.print(png,width=640,file="qqplot_vienna2010_cont_alpha60.png")

cov_ratio_pop=c(50/150,60/(270+60),50/(320+50))
get_ab_ratio(cov_ratio_pop)


