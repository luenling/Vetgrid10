get_cov_ratio=function(coverage){
	coverage2=coverage[,2:(length(coverage[1,])-1)]
	indcs=seq(2,length(coverage2[1,]),by=2)
	return(coverage2[,indcs]/coverage2[,indcs-1])	}
coverage=read.table("/Volumes/Temp/Lukas/Data/SA_A7/Polymorphisms/sa_a7_VL1_VL4_VL5_L1_L4_L5_DVD1_DVD4_DVD5_VVD1_VVD4_VVD5_q20_filt_mc2_chroms_only.coverage", header=F, skip=1,nrows=14)
rownames(coverage)=coverage[,1]
get_cov_ratio(coverage)[c("X","2L","2R","3L","3R","4"),]
coverage2=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20.coverage", header=F, skip=1,nrows=20)
rownames(coverage2)=coverage2[,1]
apply(get_cov_ratio(coverage2)[c("X","2L","2R","3L","3R","4"),],1,median)
apply(get_cov_ratio(coverage)[c("X","2L","2R","3L","3R","4"),],1,median)
coverage3=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20.coverage", header=F, skip=1,nrows=19)
rownames(coverage3)=coverage3[,1]
apply(get_cov_ratio(coverage3)[c("X","2L","2R","3L","3R","4"),],1,median)
coverage3=read.table("/Volumes/Temp/Lukas/Data/Trident/Realigned/Polymorphisms/males_vie_pTD1_pTL2_pTL1_pTD2_q20_filt_mc1.5_chroms_only.coverage", header=F, skip=1,nrows=14)
rownames(coverage3)=coverage3[,1]
apply(get_cov_ratio(coverage3)[c("X","2L","2R","3L","3R","4"),],1,median)
coverage3=read.table("/Volumes/Temp/Lukas/Data/Trident/Italy/Polymorphisms/males_ita11_pTL1_pTL2_pTL3_pTD1_pTD2_pTD3_q20_filt_mc2_chroms_only.coverage", header=F, skip=2,nrows=19)
rownames(coverage3)=coverage3[,1]
apply(get_cov_ratio(coverage3)[c("X","2L","2R","3L","3R","4"),],1,median)
