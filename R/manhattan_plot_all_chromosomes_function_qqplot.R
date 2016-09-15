##############################################################################
###Manhattan plot function 
###from http://gettinggeneticsdone.blogspot.co.at/2011/04/annotated-manhattan-plots-and-qq-plots.html
###Modified for Dmel chromosomes
##############################################################################

draw_line_manh <- function(coord,chr_offset,y=c(0,0))
{
    # should draw line at given chrom x coordinates
    
}

add_points_manhattan <- function(point_set, t, color="red",nplog=TRUE, ymax="max",ymin=0,charac="+")
    {
        #should draw points at given color onto a manhattan plot
        ps=point_set
        limitchromosomes = colnames(chrom_offset)
        if (nplog)
            {
                ymin<- 10^-(ymin)
                ps = subset(na.omit(ps), (P>0 & P<=ymin)) # remove na's, sort, and keep only 0<P<=10^ymin for nplog
                ps$P=-1*log10(ps$P)
                
            }
        else
            {
                if  (ymax == "max") ymax<-ceiling(max(ps$P))
                ps = subset(na.omit(ps), (P >=ymin & P <= ymax))
            }
        
        for(i in colnames(chrom_offset))
            {
                ps$BP[ps$CHR == i] = ps$BP[ps$CHR == i] + chrom_offset[i]
            }
        points(ps$BP,ps$P,col=color,pch=charac)
    }

get_offset_manh <- function(dataframe, ymax="max", limitchromosomes=c("X","2L","2R","3L","3R","4","2LHet","2RHet","3LHet","3RHet","XHet","YHet"), ymin=0, nplog=TRUE)
{
# returns named vector with offsets for plot 	
	d=dataframe
	if (!("CHR" %in% names(d) & "BP" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
    chromregion = NULL 
    if (length(limitchromosomes)>0) {
    	# check if region for chromosome given
    	if (length(limitchromosomes)==1 && grepl('\\w+\\:\\d+\\-\\d+',limitchromosomes,perl=TRUE)){
    		limitchromosomes = unlist(strsplit(limitchromosomes,"[:-]+"))
    		chromregion = as.numeric(limitchromosomes[2:3])
    		limitchromosomes <- limitchromosomes[1]
    	}
    	d=d[d$CHR %in% limitchromosomes, ]
    	}
        d$CHR <- factor(d$CHR,levels<-limitchromosomes)
        if (nplog)
            {
                ymin<- 10^-(ymin)
                d = subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=ymin)) # remove na's, sort, and keep only 0<P<=10^ymin for nplog
            }
        else
            {
                if  (ymax=="max") ymax<-ceiling(max(d$P))
                d = subset(na.omit(d[order(d$CHR, d$BP), ]), (P >=ymin & P <= ymax))
            }
    # if positions given extract only that region
		# order by chromosome vector
    if ( length(chromregion) > 0) d = subset(d, (d$BP >= chromregion[1]  & d$BP <= chromregion[2]))
	lastbase=0
	cc=0
	chr_offset=c()
	chrs=unique(d$CHR)
	for (i in chrs) {
		cc=cc+1
		if (cc==1) {
			chr_offset[1]=0
		} else {
		lastbase=lastbase+tail(subset(d,CHR==chrs[cc-1])$BP, 1)
		chr_offset[cc]=lastbase
		}
	}
	names(chr_offset)<-chrs
	return(chr_offset)
}

manhattan <- function(dataframe, colors=c("gray60","black"), ymax="max",ymin=0, ylab=expression(-log[10](italic(P))), limitchromosomes=c("X","2L","2R","3L","3R","4","2LHet","2RHet","3LHet","3RHet","XHet","YHet"), suggestiveline=NULL, genomewideline=NULL, annotate=NULL, cand_col="black",  pointscale="linear", ps_max = 1, ps_min = 0.1, p_ch=20, nplog=TRUE, replot=F,...) {
	# to only get a region give limitchromosomes=c("X:10000-200000")
	# to change the scaling of points give pointscale either as "linear", "none" to turn it off 
	# or to a threshold value underneath which cex=ps_min and above cex=ps_max
	# to change the printchar set p_ch to something else
	# to annotate a region or a point give in the following format
	# annotate=c("2L:10000-15000,2R:190000-20000","X:110000"), cand_col=c("red","green")
	# the cand_color vector must be of the same size as the annotate vector
	# if you do not want to have a logarithmic scale give nplog=FALSE and set y min to something low
	d=dataframe
	if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
        chromregion = NULL 
        if (length(limitchromosomes)>0) {
                                        # check if region for chromosome given
            if (length(limitchromosomes)==1 && grepl('\\w+\\:\\d+\\-\\d+',limitchromosomes,perl=TRUE)){
    		limitchromosomes = unlist(strsplit(limitchromosomes,"[:-]+"))
    		chromregion = as.numeric(limitchromosomes[2:3])
    		limitchromosomes <- limitchromosomes[1]
            }
            d=d[d$CHR %in% limitchromosomes, ]
    	   	
        }
                                         # order by chromosome vector
        d$CHR <- factor(d$CHR,levels<-limitchromosomes)
        if (nplog)
            {
                ymin<- 10^-(ymin)
                d = subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=ymin)) # remove na's, sort, and keep only 0<P<=10^ymin for nplog
                d$logp = -log10(d$P)
                dmax=max(d$logp)
            }
        else
            {
                if  (ymax=="max") ymax<-ceiling(max(d$P))
                d = subset(na.omit(d[order(d$CHR, d$BP), ]), (P >=ymin & P <= ymax))
                d$logp = d$P
                dmax=max(d$logp)
            }
    # if positions given extract only that region
        if ( length(chromregion) > 0) d = subset(d, (d$BP >= chromregion[1]  & d$BP <= chromregion[2]))
	d$pos=NA
	ticks=NULL
	lastbase=0
	 
	if (pointscale == "linear") { d$CEX= (ps_min + (ps_max-ps_min) * d$logp/dmax) }
	else if (pointscale == "none") {d$CEX = 1}
	else {
		threshold <- as.numeric(pointscale)
		d$CEX= ifelse(d$logp >= threshold, ps_max, ps_min)
	}
	numchroms=length(unique(d$CHR))
	colors <- rep(colors,numchroms)[1:numchroms]
	# setting the pointsizes in dependence of p Values
	if (ymax=="max") ymax<-ceiling(max(d$logp))
	#if (ymax<8) ymax<-8
	print("Setting up plot... ")
	if (numchroms==1) {
		d$pos=d$BP
		ticks=floor(length(d$pos))/2+1
		} 
	else {
		cc=0
		for (i in unique(d$CHR)) {
			cc=cc+1
			if (cc==1) {
				d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
			} else {
				lastbase=lastbase+tail(subset(d,CHR==unique(d$CHR)[cc-1])$BP, 1)
				d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
			}
			ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
		}
	}
	print("Plotting points... ")
	if (numchroms==1) {
		with(d, plot(pos, logp, ylim=c(0,ymax), ylab="", xlab=unique(d$CHR), col=colors, cex.lab=1.1, pch=p_ch, cex=d$CEX, ...))
		}	else {
      if (replot != T){
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=ylab, xlab="Chromosome", xaxt="n", type="n", pch=p_ch, cex=d$CEX, ...))
        axis(1, at=ticks, lab=unique(d$CHR), ...)
      }		
		icol=1
		for (i in unique(d$CHR)) {
			with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], pch=p_ch,cex=d[d$CHR==i, ]$CEX, ...))
			icol=icol+1
			}
		}

    if (!is.null(annotate)) {
    	nn=0
    	for(set in annotate){
    		nn=nn+1
    		print("Annotating candidates... ")
    		cand_regions = unlist(strsplit(set,"[,]+"))
    		for (i in cand_regions) {
    			region = unlist(strsplit(i,"[:-]+"))
    			if (length(region) == 2 ) {region[3] = region[2] }
    			chromosome = region[1]
    			chrom_region = as.numeric(region[2:3])
				d.annotate=d[which(d$CHR == chromosome & chrom_region[1] <= d$BP & d$BP <= chrom_region[2]),]
    			with(d.annotate, points(pos, logp, col=cand_col[nn], cex=CEX, pch=p_ch, ...)) 
    		

    	 	}
    	}
    }	
    if(!is.null(suggestiveline)) abline(h=suggestiveline, col="red",lty=3,lwd=1)
    if(!is.null(genomewideline)) abline(h=genomewideline, col="black",lty=1)
}

####Second function to deal with allele frequency changes rather than pvalues
manhattan_AF <- function(dataframe, colors="gray75", ymax="max", ymin="min", limitchromosomes=c("2L","2R","3L","3R","4","X","2LHet","2RHet","3LHet","3RHet","XHet","YHet"), suggestiveline=NULL, genomewideline=NULL, annotate=NULL, cand_col="black", cand_size=1, ...) {

    d=dataframe
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "AF" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and AF")
    chromregion = NULL 
    if (length(limitchromosomes)>0) {
    	# check if region for chromosome given
    	if (length(limitchromosomes)==1 && grepl('\\w+\\:\\d+\\-\\d+',limitchromosomes,perl=TRUE)){
    		limitchromosomes = unlist(strsplit(limitchromosomes,"[:-]+"))
    		chromregion = as.numeric(limitchromosomes[2:3])
    		limitchromosomes <- limitchromosomes[1]
    	}
    	d=d[d$CHR %in% limitchromosomes, ]
    	   	
    	}
    d=subset(na.omit(d[order(d$CHR, d$BP), ])) # remove na's, sort
    # if positions given extract only that region
    if ( length(chromregion) > 0) d = subset(d, d$BP >= chromregion[1]  & d$BP <= chromregion[2])
    d$pos=NA
    ticks=NULL
     
    lastbase=0
    numchroms=length(unique(d$CHR))
    colors <- rep(colors,numchroms)[1:numchroms]
    if (ymax=="max") ymax<-ceiling(max(d$AF))
    if (ymin=="min") ymin<-floor(min(d$AF))
    
    print("Setting up plot... ")
    if (numchroms==1) {
        d$pos=d$BP
        ticks=floor(length(d$pos))/2+1
    } else {
	  cc=0
        for (i in unique(d$CHR)) {
	    cc=cc+1
          if (cc==1) {
    			d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
    		} else {
    			lastbase=lastbase+tail(subset(d,CHR==unique(d$CHR)[cc-1])$BP, 1)
    			d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
    		}
    		ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
    	  }
    }
    
    print("Plotting points... ")
    if (numchroms==1) {
        with(d, plot(pos, d$AF, col=colors, ylim=c(ymin,ymax), ylab="", xlab=unique(d$CHR), cex.lab=1.1, ...))
    }	else {
        with(d, plot(pos, d$AF, ylim=c(ymin,ymax), ylab="Allele frequency change", xlab="Chromosome", xaxt="n", type="n", ...))
        axis(1, at=ticks, lab=unique(d$CHR), ...)
        icol=1
        for (i in unique(d$CHR)) {
            with(d[d$CHR==i, ],points(pos, d$AF, col=colors[icol], ...))
            icol=icol+1
    	  }
    }

    if (!is.null(annotate)) {
	  nn=0
	  for(set in annotate){
	  	nn=nn+1
	  	print("Annotating candidates... ")
        	d.annotate=d[which(d$SNP %in% set), ]
        	with(d.annotate, points(pos, AF, col=cand_col[nn], ...)) 
        }
    }
    
    if(!is.null(suggestiveline)) abline(h=suggestiveline, col="blue")
    if(!is.null(genomewideline)) abline(h=genomewideline, col="black")
}

##Make new column with minor allele freq for base population
GetMinorAllele <- function(j){
	if(j>0.5){
		return(1-j)
	} else {
		return(j)
	}
}

##Make new column based on minor allele freq in base for evolved populations
GetMinorAllele_evolved <- function(df){
	base <- df[1]
	evolved <- df[2]
	if(base>0.5){
		return(1-evolved)
	} else {
		return(evolved)
	}
}

##make column with rising allele frequency in base
GetRisingAllele <- function(df){
	base <- df[1]
	evolved <- df[2]
	if(base>evolved){
		return(1-base)
	} else {
		return(base)
	}
}
##############################################
# Taken and altered from Stephen Turner
# http://StephenTurner.us/
# http://GettingGeneticsDone.blogspot.com/
# See license at http://gettinggeneticsdone.blogspot.com/p/copyright.html 
# Make a pretty QQ plot of p-values
##############################################
qq = function(pvector, ...) {
	if (!is.numeric(pvector)) stop("D'oh! P value vector is not numeric.")
	pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]
	o = -log10(sort(pvector,decreasing=F))
	e = -log10( ppoints(length(pvector) ))
	plot(e,o,pch=19,cex=1, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))), xlim=c(0,max(e)), ylim=c(0,max(o)), ...)
	abline(0,1,col="red")
}



##############################################
#make dataframe for analyses
#############################################
# data <- with(All_SNPs_env, data.frame(Chromosome,Position,Base_All_freq,Cold_All_freq=F15Cold_All_freq,Hot_All_freq=F15Hot_All_freq,BaseCold_pval=Base.F15Cold_pval,BaseHot_pval=Base.F15Hot_pval,HotCold_pval=F15Hot.F15Cold_pval))
# data$ID <- paste(data$Chromosome,data$Position,sep="_")


##############################################
##Do manhattan plots
##############################################

########################
###pvalues
########################
# hotdf  <-  data[!is.na(data$BaseHot_pval),c("ID","Chromosome","Position","BaseHot_pval")]
# colddf  <-  data[!is.na(data$BaseCold_pval),c("ID","Chromosome","Position","BaseCold_pval")]
# hotcolddf  <-  data[!is.na(data$HotCold_pval),c("ID","Chromosome","Position","HotCold_pval")]
# colnames(hotdf) <- colnames(colddf) <-colnames(hotcolddf) <- c("SNP", "CHR", "BP", "P")

# ##get sets of SNPs according to overalp between hot and cold:
# hotcands <- hotdf[order(hotdf$P),][1:2000,]
# coldcands <- colddf[order(colddf$P),][1:2000,]

# #postscript format
# postscript("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Hot_chromosomes_NoChorion.eps",width=10,height=7.5)
# layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
# par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
# Ymax <- max(-log10(hotcands$P))
# for(chrom in c("X","2L","2R","3L","3R")){
	# hc <- hotcands$SNP[hotcands$CHR==chrom]
	# cc <- coldcands$SNP[coldcands$CHR==chrom]
	# hotonly <- setdiff(hc,cc)
	# coldonly <-  setdiff(cc,hc)
	# hotandcold <- intersect(hc,cc)
	# manhattan(hotdf, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col=c("red","lightblue","black"), ymax=Ymax)
# }
# mtext(c(expression(-log[10](italic(p))),"Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
# dev.off()

# postscript("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Cold_chromosomes_NoChorion.eps",width=10,height=7.5)
# layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
# par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
# Ymax <- max(-log10(hotcands$P))
# for(chrom in c("X","2L","2R","3L","3R")){
	# hc <- hotcands$SNP[hotcands$CHR==chrom]
	# cc <- coldcands$SNP[coldcands$CHR==chrom]
	# hotonly <- setdiff(hc,cc)
	# coldonly <-  setdiff(cc,hc)
	# hotandcold <- intersect(hc,cc)
	# manhattan(colddf, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col=c("red","lightblue","black"), ymax=Ymax)
# }
# mtext(c(expression(-log[10](italic(p))),"Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
# dev.off()


# #png format
# png("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Hot_chromosomes_NoChorion.png",width=1000,height=750)
# layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
# par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
# Ymax <- max(-log10(hotcands$P))
# for(chrom in c("X","2L","2R","3L","3R")){
	# hc <- hotcands$SNP[hotcands$CHR==chrom]
	# cc <- coldcands$SNP[coldcands$CHR==chrom]
	# hotonly <- setdiff(hc,cc)
	# coldonly <-  setdiff(cc,hc)
	# hotandcold <- intersect(hc,cc)
	# manhattan(hotdf, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col=c("red","lightblue","black"), ymax=Ymax)
# }
# mtext(c(expression(-log[10](italic(p))),"Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
# dev.off()

# png("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Cold_chromosomes_NoChorion.png",width=1000,height=750)
# layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
# par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
# Ymax <- max(-log10(hotcands$P))
# for(chrom in c("X","2L","2R","3L","3R")){
	# hc <- hotcands$SNP[hotcands$CHR==chrom]
	# cc <- coldcands$SNP[coldcands$CHR==chrom]
	# hotonly <- setdiff(hc,cc)
	# coldonly <-  setdiff(cc,hc)
	# hotandcold <- intersect(hc,cc)
	# manhattan(colddf, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col=c("red","lightblue","black"), ymax=Ymax)
# }
# mtext(c(expression(-log[10](italic(p))),"Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
# dev.off()
#expression(-log[10](italic(p)))
