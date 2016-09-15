library(hexbin)
qqplot_hex <- function(nl_oP,nl_dP1,title,alpha=FALSE,beta=" by coverage"){
	s_nl_oP <- sort(nl_oP,decreasing=F)
	s_nl_dP <- sort(nl_dP1,decreasing=F)
	bin <- hexbin(s_nl_dP,s_nl_oP,xbins=400,xlab=expression(paste("-log(P",scriptstyle(null),")")),ylab=expression(paste("-log(P",scriptstyle(obs),")")))
	if (alpha){
		main_tit=bquote(paste(.(title)," ",alpha,"=",.(alpha)," ", beta,.(beta)))
		}
	else{
		main_tit=bquote(.(title))
	}
	pp<-plot(bin,legend=FALSE,style="constant.col",main=main_tit)
	
	hvp=hexViewpor(bin)
	hexVP.abline(hvp,0,1,col="red",lty=2)
}

source("/Volumes/Temp/Lukas/Tools/Scripts/R/manhattan_plot_all_chromosomes_function_qqplot.R")
setwd("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/")
data <- read.table("females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25_nlpl1.gwas",header=TRUE)
data$ID <- paste(data$CHR,data$BP,sep="_")
quartz(height=6,width=12)
par(bg="white")
chroms=c("X","2L","2R","3L","3R","4")#,"XHet","2LHet","2RHet","3LHet","3RHet")
manhattan(data,colors=c("black","slategrey"),limitchromosomes=chroms,main="Vienna (2010)", suggestiveline=c(-log10(4.134235e-07)))
dev.print(png,width=800,file="vie2010_ccr_manhattan.png")

source("/Volumes/Temp/Lukas/Tools/Scripts/R/manhattan_plot_all_chromosomes_function_qqplot.R")
setwd("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/")
data11 <- read.table("females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20_filt_mc8_chroms_only_mct8_mcv25_nlpl1.gwas",header=TRUE)
data11$ID <- paste(data$CHR,data$BP,sep="_")
quartz(height=6,width=12)
par(bg="white")
chroms=c("X","2L","2R","3L","3R","4")#,"XHet","2LHet","2RHet","3LHet","3RHet")
manhattan(data11,colors=c("black","slategrey"),limitchromosomes=chroms,main="Bolzano (2011)", suggestiveline=c(-log10(2.324404e-09)))
dev.print(png,width=800,file="boz2011_ccr_manhattan.png")
