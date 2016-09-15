library(optparse)
option_list <- list(
  make_option(c("-i", "--infile"),type="character", dest="infile",default=F,
              help="Input file to load"),
  make_option(c("-n", "--name"),type="character", dest="basename",default=F,
              help="Basename for figure title and output file")
  )
  
cl_options=parse_args(OptionParser(option_list = option_list))
if (cl_options$infile == F || cl_options$basename == F){
  print("infile and basename required")
  cl_options=parse_args(OptionParser(option_list = option_list),args=c("--help"))
}
# read table with mc bp
a=read.table(cl_options$infile)
# sort a
a=a[order(a[,1]),]
b = cbind((a[-1,1]+a[-length(a[,1]),1])/2,abs(diff(a[,2])/diff(a[,1])))
outfile=paste0(cl_options$basename,"_indelmasking.pdf")
pdf(file=outfile)
par(mar=c(5,4,4,5)+.1)
plot(a,xlab="indel threshold",ylab="# of masked bp",main=paste("InDel masking for",cl_options$basename),xlim=range(a[,1]))
par(new=T)
plot(b,axes=F,pch=4,col="red",xlab="",ylab="",xlim=range(a[,1]))
axis(4)
mtext("change of masked bps per minimal count change", side = 4,line = 3)
legend("topright",col=c("black","red"),pch=c(0,4),legend=c("bp masked","change of masked bp"))
dev.off()
