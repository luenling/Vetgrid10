setwd("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/")
vie_snps=read.table("all_snps_samtools_I16_RPB.tab",na.strings = "NA")
colnames(vie_snps)=c("CHR","BPS","ALL","I16_0","I16_1","I16_2","I16_3","I16_12","I16_13","I16_14","I16_15","RPB")
vie_snps$n1=vie_snps$I16_0 + vie_snps$I16_1
vie_snps$n2=vie_snps$I16_2 + vie_snps$I16_3
vie_snps$SB=apply(vie_snps[,6:7],1,function(x) {min(x)/max(x)})
vie_snps$FS=apply(vie_snps[,4:7],1,function(x) {fisher.test(matrix(unlist(x),nrow=2,byrow=T))$p.value})
vie_snps$ATD=vie_snps$I16_14/vie_snps$n2
vie_snps$RTD=vie_snps$I16_12/vie_snps$n1
vie_snps$RV=(vie_snps$I16_13-vie_snps$n1*(vie_snps$RTD)^2 )/(vie_snps$n1-1)
vie_snps$AV=(vie_snps$I16_15-vie_snps$n2*(vie_snps$ATD)^2 )/(vie_snps$n2-1)
vie_snps$T=(vie_snps$RTD-vie_snps$ATD)/sqrt(vie_snps$RV/vie_snps$n1 + vie_snps$AV/vie_snps$n2 )
vie_snps$DF=(vie_snps$RV/vie_snps$n1 + vie_snps$AV/vie_snps$n2 )^2/(vie_snps$RV^2/(vie_snps$n1^2*(vie_snps$n1-1)) + vie_snps$AV^2/(vie_snps$n2^2*(vie_snps$n2-1)) )
vie_snps$TDB=1-pt(vie_snps$T,vie_snps$DF)
vie_snps$lTDB=-1*log10(vie_snps$TDB)
vie_snps$lTDB[is.infinite(vie_snps$lTDB)]=200
vie_snps$lFS=-1*log10(vie_snps$FS)
vie_snps$lFS[is.infinite(vie_snps$lFS)]=200
