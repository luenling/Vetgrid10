# get read IDs of a input bam file 
samtools view $1 | awk '{print $1}' > `basename $1`.readIDs
gzip `basename $1`.readIDs
rm -f $1