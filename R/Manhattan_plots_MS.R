#################################################################################################################################
## Create database of all common SNPs -- uses the COVandFREQ file which is derived from the combined sync file with CMH values
#################################################################################################################################
#Vetgrid 13
setwd("/Volumes/Temp/Ray/Databases")
library(filehash)
dumpDF(read.table("/Volumes/Temp/Ray/CMH_python_output_NEW/Merged_CMH_files/Merged_CMH_Base-F15HotCold_NEW_COVandFREQ_NoChorion.txt", header=TRUE, sep="\t"), dbName="F0_F15_All_SNPs_COVandFREQ_no_chorion")
All_SNPs_env <- db2env(db="F0_F15_All_SNPs_COVandFREQ_no_chorion")
#check structure of environment
#ls.str(All_SNPs_env)


##############################################################################
###Manhattan plot function 
###from http://gettinggeneticsdone.blogspot.co.at/2011/04/annotated-manhattan-plots-and-qq-plots.html
###Modified for Dmel chromosomes
##############################################################################

manhattan <- function(dataframe, colors=c("gray60", "gray75"), ymax="max", limitchromosomes=c("2L","2LHet","2R","2RHet","3L","3LHet","3R","3RHet","4","X","XHet","YHet"), suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), annotate=NULL, cand_col="black", ...) {

    d=dataframe
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
    
    if (length(limitchromosomes)>0) d=d[d$CHR %in% limitchromosomes, ]
    d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
    d$logp = -log10(d$P)
    d$pos=NA
    ticks=NULL
    lastbase=0
    numchroms=length(unique(d$CHR))
    colors <- rep(colors,numchroms)[1:numchroms]
    if (ymax=="max") ymax<-ceiling(max(d$logp))
    if (ymax<8) ymax<-8
    
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
        with(d, plot(pos , logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
    }	else {
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
        axis(1, at=ticks, lab=unique(d$CHR), ...)
        icol=1
        for (i in unique(d$CHR)) {
            with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
            icol=icol+1
    	  }
    }

    if (!is.null(annotate)) {
	  nn=0
	  for(set in annotate){
	  	nn=nn+1
	  	print("Annotating candidates... ")
        	d.annotate=d[which(d$SNP %in% set), ]
        	with(d.annotate, points(pos, logp, col=cand_col[nn], ...)) 
        }
    }
    
    if (suggestiveline) abline(h=suggestiveline, col="blue")
    if (genomewideline) abline(h=genomewideline, col="red")
}

####Second function to deal with allele frequency changes rather than pvalues
manhattan <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", limitchromosomes=c("2L","2LHet","2R","2RHet","3L","3LHet","3R","3RHet","4","X","XHet","YHet"), suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), annotate=NULL, cand_col="black", cand_size=1, ......) {

    d=dataframe
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "AFC" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and AFC")
    
    if (length(limitchromosomes)>0) d=d[d$CHR %in% limitchromosomes, ]
    d=subset(na.omit(d[order(d$CHR, d$BP), ])) # remove na's, sort
    d$pos=NA
    ticks=NULL
    lastbase=0
    numchroms=length(unique(d$CHR))
    colors <- rep(colors,numchroms)[1:numchroms]
    if (ymax=="max") ymax<-ceiling(max(d$AFC))
    if (ymax<8) ymax<-8
    
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
        with(d, plot(pos, d$AFC, ylim=c(0,ymax), ylab="Allele frequency change", xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
    }	else {
        with(d, plot(pos, d$AFC, ylim=c(0,ymax), ylab="Allele frequency change", xlab="Chromosome", xaxt="n", type="n", ...))
        axis(1, at=ticks, lab=unique(d$CHR), ...)
        icol=1
        for (i in unique(d$CHR)) {
            with(d[d$CHR==i, ],points(pos, d$AFC, col=colors[icol], ...))
            icol=icol+1
    	  }
    }

    if (!is.null(annotate)) {
	  nn=0
	  for(set in annotate){
	  	nn=nn+1
	  	print("Annotating candidates... ")
        	d.annotate=d[which(d$SNP %in% set), ]
        	with(d.annotate, points(pos, d$AFC, col=cand_col[nn], ...)) 
        }
    }
    
    if (suggestiveline) abline(h=suggestiveline, col="blue")
    if (genomewideline) abline(h=genomewideline, col="red")
}

##############################################
#make dataframe for analyses
#############################################
data <- with(All_SNPs_env, data.frame(Chromosome,Position,Base_All_freq,Cold_All_freq=F15Cold_All_freq,Hot_All_freq=F15Hot_All_freq,BaseCold_pval=Base.F15Cold_pval,BaseHot_pval=Base.F15Hot_pval,HotCold_pval=F15Hot.F15Cold_pval))
data$ID <- paste(data$Chromosome,data$Position,sep="_")


##############################################
##Do manhattan plots
##############################################
hotdf  <-  data[!is.na(data$BaseHot_pval),c("ID","Chromosome","Position","BaseHot_pval")]
colddf  <-  data[!is.na(data$BaseCold_pval),c("ID","Chromosome","Position","BaseCold_pval")]
hotcolddf  <-  data[!is.na(data$HotCold_pval),c("ID","Chromosome","Position","HotCold_pval")]
colnames(hotdf) <- colnames(colddf) <-colnames(hotcolddf) <- c("SNP", "CHR", "BP", "P")

##get sets of SNPs according to overalp between hot and cold:
hotcands <- hotdf$SNP[order(hotdf$P)][1:2000]
coldcands <- colddf$SNP[order(colddf$P)][1:2000]
hotonly <- setdiff(hotcands,coldcands)
coldonly <- setdiff(coldcands,hotcands)
hotandcold <- intersect(hotcands,coldcands)

#postscript format
postscript("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/ALL_POPS_CMH_python_NEW_Base-F15Hot_1-9_2-4_3-5_CANDs2000_NoChorion.eps",width=10,height=7.5)
manhattan(hotdf, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), genomewideline=F, suggestiveline=F, main="", limitchromosomes=c("2L","2R","3L","3R","X"),cex=0.5, pch=20, cand_col=c("red","lightblue","black"))
dev.off()

postscript("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/ALL_POPS_CMH_python_NEW_Base-F15Cold_1-6_2-8_3-7_CANDs2000_NoChorion.eps",width=10,height=7.5)
manhattan(colddf, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), genomewideline=F, suggestiveline=F, main="", limitchromosomes=c("2L","2R","3L","3R","X"),cex=0.5, pch=20, cand_col=c("red","lightblue","black"))
dev.off()

postscript("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/ALL_POPS_CMH_python_NEW_Base-F15Hot_9-6_2-8_3-7_CANDs2000_NoChorion.eps",width=10,height=7.5)
manhattan(hotcolddf, annotate=list(hotcolddf$SNP[order(hotcolddf$P)][1:2000]), genomewideline=F, suggestiveline=F, main="", limitchromosomes=c("2L","2R","3L","3R","X"),cex=0.5, pch=20, cand_col="green")
dev.off()

#png format
png("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/ALL_POPS_CMH_python_NEW_Base-F15Hot_1-9_2-4_3-5_CANDs2000_NoChorion.png",width=600,height=400)
manhattan(hotdf, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), genomewideline=F, suggestiveline=F, main="", limitchromosomes=c("2L","2R","3L","3R","X"),cex=0.5, pch=20, cand_col=c("red","lightblue","black"))
dev.off()

png("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/ALL_POPS_CMH_python_NEW_Base-F15Cold_1-6_2-8_3-7_CANDs2000_NoChorion.png",width=600,height=400)
manhattan(colddf, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), genomewideline=F, suggestiveline=F, main="", limitchromosomes=c("2L","2R","3L","3R","X"),cex=0.5, pch=20, cand_col=c("red","lightblue","black"))
dev.off()

png("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/ALL_POPS_CMH_python_NEW_Base-F15Hot_9-6_2-8_3-7_CANDs2000_NoChorion.png",width=600,height=400)
manhattan(hotcolddf, annotate=list(hotcolddf$SNP[order(hotcolddf$P)][1:2000]), genomewideline=F, suggestiveline=F, main="", limitchromosomes=c("2L","2R","3L","3R","X"), pch=20, cex=0.5, cand_col="green")
dev.off()
