#################################################################################################################################
## Create database of all common SNPs -- uses the COVandFREQ file which is derived from the combined sync file with CMH values
#################################################################################################################################
#Vetgrid 13
setwd("/Volumes/Temp/Ray/Databases")
library(filehash)
#dumpDF(read.table("/Volumes/Temp/Ray/CMH_python_output_NEW/Merged_CMH_files/Merged_CMH_Base-F15HotCold_NEW_COVandFREQ_NoChorion.txt", header=TRUE, sep="\t"), dbName="F0_F15_All_SNPs_COVandFREQ_no_chorion")
All_SNPs_env <- db2env(db="F0_F15_All_SNPs_COVandFREQ_no_chorion")
#check structure of environment
#ls.str(All_SNPs_env)


##############################################################################
###Manhattan plot function 
###from http://gettinggeneticsdone.blogspot.co.at/2011/04/annotated-manhattan-plots-and-qq-plots.html
###Modified for Dmel chromosomes
##############################################################################

manhattan <- function(dataframe, colors="gray60", ymax="max", limitchromosomes=c("2L","2LHet","2R","2RHet","3L","3LHet","3R","3RHet","4","X","XHet","YHet"), suggestiveline=NULL, genomewideline=NULL, annotate=NULL, cand_col="black", ...) {
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
		with(d, plot(pos, logp, ylim=c(0,ymax), ylab="", xlab=unique(d$CHR), col=colors, cex.lab=1.1, ...))
		}	else {
		with(d, plot(pos, logp, ylim=c(0,ymax), ylab="expression(-log[10](italic(p)))", xlab="Chromosome", xaxt="n", type="n", ...))
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
    		
    if(!is.null(suggestiveline)) abline(h=suggestiveline, col="blue")
    if(!is.null(genomewideline)) abline(h=genomewideline, col="black")
}

####Second function to deal with allele frequency changes rather than pvalues
manhattan_AF <- function(dataframe, colors="gray75", ymax="max", ymin="min", limitchromosomes=c("2L","2LHet","2R","2RHet","3L","3LHet","3R","3RHet","4","X","XHet","YHet"), suggestiveline=NULL, genomewideline=NULL, annotate=NULL, cand_col="black", cand_size=1, ...) {

    d=dataframe
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "AF" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and AF")
    
    if (length(limitchromosomes)>0) d=d[d$CHR %in% limitchromosomes, ]
    d=subset(na.omit(d[order(d$CHR, d$BP), ])) # remove na's, sort
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
#make dataframe for analyses
#############################################
data <- with(All_SNPs_env, data.frame(Chromosome,Position,Base_All_freq,Cold_All_freq=F15Cold_All_freq,Hot_All_freq=F15Hot_All_freq,BaseCold_pval=Base.F15Cold_pval,BaseHot_pval=Base.F15Hot_pval,HotCold_pval=F15Hot.F15Cold_pval))
data$ID <- paste(data$Chromosome,data$Position,sep="_")


##############################################
##Do manhattan plots
##############################################

########################
###pvalues
########################
hotdf  <-  data[!is.na(data$BaseHot_pval),c("ID","Chromosome","Position","BaseHot_pval")]
colddf  <-  data[!is.na(data$BaseCold_pval),c("ID","Chromosome","Position","BaseCold_pval")]
hotcolddf  <-  data[!is.na(data$HotCold_pval),c("ID","Chromosome","Position","HotCold_pval")]
colnames(hotdf) <- colnames(colddf) <-colnames(hotcolddf) <- c("SNP", "CHR", "BP", "P")

##get sets of SNPs according to overalp between hot and cold:
hotcands <- hotdf[order(hotdf$P),][1:2000,]
coldcands <- colddf[order(colddf$P),][1:2000,]

#postscript format
postscript("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Hot_chromosomes_NoChorion.eps",width=10,height=7.5)
layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
Ymax <- max(-log10(hotcands$P))
for(chrom in c("X","2L","2R","3L","3R")){
	hc <- hotcands$SNP[hotcands$CHR==chrom]
	cc <- coldcands$SNP[coldcands$CHR==chrom]
	hotonly <- setdiff(hc,cc)
	coldonly <-  setdiff(cc,hc)
	hotandcold <- intersect(hc,cc)
	manhattan(hotdf, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col=c("red","lightblue","black"), ymax=Ymax)
}
mtext(c(expression(-log[10](italic(p))),"Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
dev.off()

postscript("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Cold_chromosomes_NoChorion.eps",width=10,height=7.5)
layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
Ymax <- max(-log10(hotcands$P))
for(chrom in c("X","2L","2R","3L","3R")){
	hc <- hotcands$SNP[hotcands$CHR==chrom]
	cc <- coldcands$SNP[coldcands$CHR==chrom]
	hotonly <- setdiff(hc,cc)
	coldonly <-  setdiff(cc,hc)
	hotandcold <- intersect(hc,cc)
	manhattan(colddf, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col=c("red","lightblue","black"), ymax=Ymax)
}
mtext(c(expression(-log[10](italic(p))),"Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
dev.off()


#png format
png("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Hot_chromosomes_NoChorion.png",width=1000,height=750)
layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
Ymax <- max(-log10(hotcands$P))
for(chrom in c("X","2L","2R","3L","3R")){
	hc <- hotcands$SNP[hotcands$CHR==chrom]
	cc <- coldcands$SNP[coldcands$CHR==chrom]
	hotonly <- setdiff(hc,cc)
	coldonly <-  setdiff(cc,hc)
	hotandcold <- intersect(hc,cc)
	manhattan(hotdf, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col=c("red","lightblue","black"), ymax=Ymax)
}
mtext(c(expression(-log[10](italic(p))),"Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
dev.off()

png("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Cold_chromosomes_NoChorion.png",width=1000,height=750)
layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
Ymax <- max(-log10(hotcands$P))
for(chrom in c("X","2L","2R","3L","3R")){
	hc <- hotcands$SNP[hotcands$CHR==chrom]
	cc <- coldcands$SNP[coldcands$CHR==chrom]
	hotonly <- setdiff(hc,cc)
	coldonly <-  setdiff(cc,hc)
	hotandcold <- intersect(hc,cc)
	manhattan(colddf, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col=c("red","lightblue","black"), ymax=Ymax)
}
mtext(c(expression(-log[10](italic(p))),"Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
dev.off()

###########################
###AFCs - minor allele freq
###########################
data$Base_All_minorfreq <- sapply(data$Base_All_freq,GetMinorAllele)
data$Hot_All_minorfreq <- apply(data.frame(data$Base_All_freq,data$Hot_All_freq),1,GetMinorAllele_evolved)
data$Cold_All_minorfreq <- apply(data.frame(data$Base_All_freq,data$Cold_All_freq),1,GetMinorAllele_evolved)
data$Hot_AFC <- data$Hot_All_minorfreq-data$Base_All_minorfreq
data$Cold_AFC <- data$Cold_All_minorfreq-data$Base_All_minorfreq

hotdf2  <-  data[!is.na(data$BaseHot_pval),c("ID","Chromosome","Position","Hot_AFC","BaseHot_pval")]
colddf2  <-  data[!is.na(data$BaseCold_pval),c("ID","Chromosome","Position","Cold_AFC","BaseCold_pval")]
colnames(hotdf2) <- colnames(colddf2) <- c("SNP", "CHR", "BP", "AF","P")

##get sets of SNPs according to overalp between hot and cold:
hotcands2 <- hotdf2[order(hotdf2$P),][1:2000,1:4]
coldcands2 <- colddf2[order(colddf2$P),][1:2000,1:4]

#postscript format
postscript("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Hot_chromosomes_NoChorion_AFC.eps",width=10,height=7.5)
layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
Ymax <- max(max(hotcands2$AF),max(coldcands2$AF))
Ymin <- min(min(hotcands2$AF),min(coldcands2$AF))
for(chrom in c("X","2L","2R","3L","3R")){
	hc <- hotcands2$SNP[hotcands2$CHR==chrom]
	cc <- coldcands2$SNP[coldcands2$CHR==chrom]
	hotonly <- setdiff(hc,cc)
	coldonly <-  setdiff(cc,hc)
	hotandcold <- intersect(hc,cc)
	manhattan_AF(hotdf2, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), genomewideline=0, main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col=c("red","lightblue","black"), ymax=Ymax, ymin=Ymin)
}
mtext(c("Allele frequency change","Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
dev.off()

postscript("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Cold_chromosomes_NoChorion_AFC.eps",width=10,height=7.5)
layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
Ymax <- max(max(hotcands2$AF),max(coldcands2$AF))
Ymin <- min(min(hotcands2$AF),min(coldcands2$AF))
for(chrom in c("X","2L","2R","3L","3R")){
	hc <- hotcands2$SNP[hotcands2$CHR==chrom]
	cc <- coldcands2$SNP[coldcands2$CHR==chrom]
	hotonly <- setdiff(hc,cc)
	coldonly <-  setdiff(cc,hc)
	hotandcold <- intersect(hc,cc)
	manhattan_AF(colddf2, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), genomewideline=0, main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col=c("red","lightblue","black"), ymax=Ymax, ymin=Ymin)
}
mtext(c("Allele frequency change","Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
dev.off()


#png format
png("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Hot_chromosomes_NoChorion_AFC.png",width=1000,height=750)
layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
Ymax <- max(max(hotcands2$AF),max(coldcands2$AF))
Ymin <- min(min(hotcands2$AF),min(coldcands2$AF))
for(chrom in c("X","2L","2R","3L","3R")){
	hc <- hotcands2$SNP[hotcands2$CHR==chrom]
	cc <- coldcands2$SNP[coldcands2$CHR==chrom]
	hotonly <- setdiff(hc,cc)
	coldonly <-  setdiff(cc,hc)
	hotandcold <- intersect(hc,cc)
	manhattan_AF(hotdf2, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), genomewideline=0, main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col=c("red","lightblue","black"), ymax=Ymax, ymin=Ymin)
}
mtext(c("Allele frequency change","Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
dev.off()

png("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Cold_chromosomes_NoChorion_AFC.png",width=1000,height=750)
layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
Ymax <- max(max(hotcands2$AF),max(coldcands2$AF))
Ymin <- min(min(hotcands2$AF),min(coldcands2$AF))
for(chrom in c("X","2L","2R","3L","3R")){
	hc <- hotcands2$SNP[hotcands2$CHR==chrom]
	cc <- coldcands2$SNP[coldcands2$CHR==chrom]
	hotonly <- setdiff(hc,cc)
	coldonly <-  setdiff(cc,hc)
	hotandcold <- intersect(hc,cc)
	manhattan_AF(colddf2, annotate=list(ho=hotonly,co=coldonly,hc=hotandcold), genomewideline=0, main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col=c("red","lightblue","black"), ymax=Ymax, ymin=Ymin)
}
mtext(c("Allele frequency change","Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
dev.off()

################################################
###Basal allele frequency - rising allele freq
#################################################
data$Base_All_risingfreqCold <- apply(data.frame(data$Base_All_freq,data$Cold_All_freq),1,GetRisingAllele)
data$Base_All_risingfreqHot <- apply(data.frame(data$Base_All_freq,data$Hot_All_freq),1,GetRisingAllele)

hotdf3  <-  data[!is.na(data$BaseHot_pval),c("ID","Chromosome","Position","Base_All_risingfreqHot","BaseHot_pval")]
colddf3  <-  data[!is.na(data$BaseCold_pval),c("ID","Chromosome","Position","Base_All_risingfreqCold","BaseCold_pval")]
colnames(hotdf3) <- colnames(colddf3) <- c("SNP", "CHR", "BP", "AF","P")

##get sets of SNPs according to overalp between hot and cold:
hotcands3 <- hotdf3[order(hotdf3$P),][1:2000,1:4]
coldcands3 <- colddf3[order(colddf3$P),][1:2000,1:4]

#postscript format
postscript("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Hot_chromosomes_NoChorion_BaseAF.eps",width=10,height=7.5)
layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
Ymax <- max(max(hotcands3$AF),max(coldcands3$AF))
Ymin <- min(min(hotcands3$AF),min(coldcands3$AF))
for(chrom in c("X","2L","2R","3L","3R")){
	hc <- hotcands3$SNP[hotcands3$CHR==chrom]
	manhattan_AF(hotdf3, annotate=list(hc), main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col="red", ymax=Ymax, ymin=Ymin)
}
mtext(c("Start Allele Frequency","Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
dev.off()

postscript("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Cold_chromosomes_NoChorion_BaseAF.eps",width=10,height=7.5)
layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
Ymax <- max(max(hotcands3$AF),max(coldcands3$AF))
Ymin <- min(min(hotcands3$AF),min(coldcands3$AF))
for(chrom in c("X","2L","2R","3L","3R")){
	cc <- coldcands3$SNP[coldcands3$CHR==chrom]
	manhattan_AF(colddf3, annotate=list(cc), main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col="lightblue", ymax=Ymax, ymin=Ymin)
}
mtext(c("Start Allele Frequency","Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
dev.off()


#png format
png("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Hot_chromosomes_NoChorion_BaseAF.png",width=1000,height=750)
layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
Ymax <- max(max(hotcands3$AF),max(coldcands3$AF))
Ymin <- min(min(hotcands3$AF),min(coldcands3$AF))
for(chrom in c("X","2L","2R","3L","3R")){
	hc <- hotcands3$SNP[hotcands3$CHR==chrom]
	manhattan_AF(hotdf3, annotate=list(hc), main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col="red", ymax=Ymax, ymin=Ymin)
}
mtext(c("Start Allele Frequency","Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
dev.off()

png("/Volumes/Temp/Ray/Plots/Plots_no_chorion_MS/Manhattan_BaseF15Cold_chromosomes_NoChorion_BaseAF.png",width=1000,height=750)
layout(matrix(c(0,1,1,0,2,2,3,3,4,4,5,5), 3, 4, byrow = TRUE))
par(mar=c(3.5,3.5,1,1),oma=c(3,3,1,1),mgp=c(2,1,0))
Ymax <- max(max(hotcands3$AF),max(coldcands3$AF))
Ymin <- min(min(hotcands3$AF),min(coldcands3$AF))
for(chrom in c("X","2L","2R","3L","3R")){
	cc <- coldcands3$SNP[coldcands3$CHR==chrom]
	manhattan_AF(colddf3, annotate=list(cc), main="", limitchromosomes=chrom, cex=0.75, pch=20, cand_col="lightblue", ymax=Ymax, ymin=Ymin)
}
mtext(c("Start Allele Frequency","Chromosome & position"),c(2,1),cex=1.2,line=1.3,outer=TRUE)
dev.off()



#expression(-log[10](italic(p)))