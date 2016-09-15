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
		parc=sub(".*_(\\d+)\\.0_nullP.out","\\1",infile,perl=TRUE)
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
