setwd=("/Volumes/Temp/Lukas/Data/Vienna_2010_2011_joined/Comparison_Vie2010_2011")
popAB_joined=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/vie2010_cs25_el_p_ita2011_cs25_hl_pV_joined_averaged_af_nowolb.cmhout.af",header=F)
colnames(popAB_joined)=c("CHR","bps","A","Acs","Act","PA","Bcs","Bct","PB")
sigB=2.203e-08
sigA=2.997e-08
popAB_joined$Mct = (popAB_joined$Act + popAB_joined$Bct)/2
popA_sig = popAB_joined$PA < sigA
popB_sig = popAB_joined$PB < sigB
popAB_sig = popB_sig & popA_sig
num_A_sig = length(popAB_joined$PA[popA_sig])
# 436
num_B_sig = length(popAB_joined$PB[popB_sig])
#534
num_AB_sig = length(popAB_joined$PB[popB_sig & popA_sig])
#196
quartz()
hist(popAB_joined$Act[popAB_sig], breaks=seq(0,1,0.02), freq=F,col=rgb(0,0,1,1/4), xlim=c(0.2,1))
hist(popAB_joined$Bct[popAB_sig], breaks=seq(0,1,0.02), freq=F,col=rgb(0,1,0,1/4), xlim=c(0.2,1),add=T)
hist(popAB_joined$Bcs[popAB_sig], breaks=seq(0,1,0.02), freq=F,col=rgb(0,1,1,1/4), xlim=c(0.2,1),add=T)
#hist(popAB_joined$Act[! popA_sig], breaks=25, freq=F, col=rgb(0,1,0,1/4), xlim=c(0.2,1), add=T)
hist(popAB_joined$Mct,breaks=seq(0,1,0.02), freq=F, col=rgb(1,0,0,1/4), xlim=c(0.2,1), add=T)
hist(popAB_joined$Bct,breaks=seq(0,1,0.02), freq=F, col=rgb(1,0,0,1/4), xlim=c(0.2,1), add=T)
hist(popAB_joined$Act[popA_sig], breaks=seq(0,1,0.02), freq=F,col=rgb(0,0,1,1/4), xlim=c(0.2,1))
hist(popAB_joined$Bct[popB_sig], breaks=seq(0,1,0.02), freq=F,col=rgb(0,1,0,1/4), xlim=c(0.2,1),add=T)

hist_A=hist(popAB_joined$Act, breaks=seq(0.0,1,0.02),plot=F)
hist_B=hist(popAB_joined$Bct, breaks=seq(0.0,1,0.02),plot=F)
hist_M=hist(popAB_joined$Mct, breaks=seq(0.0,1,0.02),plot=F)
hist_A_sigA=hist(popAB_joined$Act[popA_sig], breaks=seq(0,1,0.02),plot=F)
hist_B_sigB=hist(popAB_joined$Bct[popB_sig], breaks=seq(0,1,0.02),plot=F)
hist_M_sigAB=hist(popAB_joined$Mct[popAB_sig], breaks=seq(0.6,1,0.02),plot=F)
hist_M_sigA=hist(popAB_joined$Mct[popA_sig], breaks=seq(0,1,0.02),plot=F)
hist_M_sigB=hist(popAB_joined$Mct[popB_sig], breaks=seq(0,1,0.02),plot=F)
hist_Mct_AB=hist(popAB_joined$Mct[popB_sig | popA_sig ], breaks=seq(0,1,0.02),plot=F)

popAB_Act=popAB_joined$Act
popAB_Act = data.frame(popAB_Act)
colnames(popAB_Act) = c("AF")
popAB_Act = popAB_Mct[order(popAB_Act$AF),]
popAB_Bct=popAB_joined$Bct
popAB_Bct = data.frame(popAB_Bct)
colnames(popAB_Bct) = c("AF")
popAB_Bct = popAB_Mct[order(popAB_Bct$AF),]
intA=seq(0.0,1,0.02)
popAB_Act$INT = findInterval(popAB_Act$AF,intA,rightmost.closed=T)
popAB_Bct$INT = findInterval(popAB_Bct$AF,intA,rightmost.closed=T)
hist_listA=list()
hist_listB=list()
for (i in seq(1,length(intA)-1)){
		hist_listA[[i]]=which(popAB_Act$INT == i)
		hist_listB[[i]]=which(popAB_Bct$INT == i)
}


sampling_frame = data.frame("zero"=numeric(0),"both"=numeric(0))
for (i in 1:100000){
	sampling_frame = rbind(sampling_frame,sample_fun(hist_A_sigA$counts,hist_B_sigB$counts,hist_listA,hist_listB))
}

sample_fun <- function(A,B,C,D){
	# A = hist with significant average allele freqs of A
	# B = hist with significant average allele freqs of B
	# C = list with indices of average allelefreq of snp in allefreq bins in A 
	# D = list with indices of average allelefreq of snp in allefreq bins in B   
	lenA=length(A)	
	results=c(0,0)
	trialA=c()
	trialB=c() 
	for (i in 1:lenA){
		if (A[i] == 0 & B[i] == 0){next}
		if (A[i] != 0){trialA = c(trialA, sample(C[[i]],A[i]))}
		if (B[i] != 0){trialB = c(trialB, sample(D[[i]],B[i]))}
	}
	results=c(length(union(trialA,trialB)),length(intersect(trialA,trialB)))
	return(results)
}

popB=read.table("/Volumes/Temp/Lukas/Data/Vienna_2011/Realigned/Poly_Real/CMH/females11_pI25_pII25_pIII25_pIhl_pIIhl_pIIIhl_realigned_q20_filt_mc8_chroms_only_mct8_mcv25.cmhout.af",header=F)
colnames(popB)=c("CHR","bps","A","csI","csII","csIII","bI","bII","bIII","P")
popA=read.table("/Volumes/Temp/Lukas/Data/Vienna_2010/Realigned/Real_poly/CMH/females_pI25_pII25_pIII25_pIel_pIIel_pIIIel_realigned_q20_filt_mc8_chroms_only_mct8_mcv25_nowolb.cmhout.af",header=F)
colnames(popA)=c("CHR","bps","A","csI","csII","csIII","bI","bII","bIII","P")
popA$cs=(popA$csI+popA$csII+popA$csIII)/3
popB$cs=(popB$csI+popB$csII+popB$csIII)/3
popA$b=(popA$bI+popA$bII+popA$bIII)/3
popB$b=(popB$bI+popB$bII+popB$bIII)/3
