setwd("/Volumes/Temp/Lukas/Sim_Mel_align/mauve_alignments")
mauve_algn=read.table("/Volumes/Temp/Lukas/Sim_Mel_align/mauve_alignments/Dsim_M252_draft_4_3R_align_parse.conserved",header=T)
mauve_ita_snps_algn=read.table("/Volumes/Temp/Lukas/Sim_Mel_align/mauve_alignments/Dmel_Dsim_M252_draft_4_3R_ita2011.conserved",header=T)

quartz()

plot(mauve_algn$win_center/10^6,mauve_algn$mismatch,col="blue",type="l",lty=1,xlab="bps [MB]",ylab="mismatched positions in D. sim./mel. aligment [percent]", main="3R", ylim=c(0,30))
lines(mauve_ita_snps_algn$win_center/10^6,mauve_ita_snps_algn$mismatch,col="red",type="l",lty=2)
legend("topleft",c("all positions","SNPs Italy 2011"),lty=c(1,2), col=c("blue","red"))
dev.copy2pdf(file="3R_mismatches_algn_itaSNPs.pdf")