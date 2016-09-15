#!/bin/bash
# turn on extended globbing behavior
# This script goes through a bam file and calculates the minimum, mean and maximum fragment length as needed for the simulans filtering script
# call with: bash get_min_mean_max_fraglength.sh infile.bam
# writes infile.fraglen with columns indicating the number of sequences read, the minimum, mean, and maximum fragmentlength
# if you only want to sample a limited amount of reads, eg. 100000, you can put a head -100000 | before the awk   
shopt -s extglob
INFILE=$1
OUTFILE=`basename $1 .bam`.fraglen
echo -e "number of\tminimal\tmean\tmaximal fragment length" > $OUTFILE
samtools view $INFILE | head -1000000 | awk 'BEGIN {sum=0.0; num=0; min=1000; max=0;}  $9 > 0 && $9 < 10000 && $7=="=" {sum+=$9; num++; if($9 < min){min=$9;}; if($9 > max) {max=$9;};} END {mean=sum/num; printf "%i\t%i\t%i\t%i\n",num,min,mean,max}' >> $OUTFILE
