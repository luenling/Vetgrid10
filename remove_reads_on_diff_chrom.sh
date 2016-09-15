#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
PICARD=/Volumes/Temp/Lukas/Tools/picard/picard-tools-1.79
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.18/samtools
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
BWA=/Volumes/Temp/Lukas/Tools/bwa-0.5.9/bwa
DIR=`pwd`
BASE_NAME=$1
#IN=`basename $1 .bam`
IN=${BASE_NAME}'.picnodupl.filtered.mq20'
[ -e ${IN}.bam ] || exit 1
echo Filtering out reads mapped to different chromosomes at `date` >> ${DIR}/${BASE_NAME}.log 
$SAMTOOLS view -h  $IN.bam |  awk ' (/^@/ || $7=="=") {print;}' | $SAMTOOLS  view -bSh - > ${IN}_chrfilt.bam
echo flagstat of ${IN}_chrfilt.bam at `date` >> ${DIR}/${BASE_NAME}.log 
$SAMTOOLS flagstat ${IN}_chrfilt.bam >> ${DIR}/${BASE_NAME}.log 

$SAMTOOLS index ${IN}_chrfilt.bam

echo idxstat of ${IN}_chrfilt.bam : >>  ${DIR}/${BASE_NAME}.log

samtools idxstats ${IN}_chrfilt.bam >> ${DIR}/${BASE_NAME}.log

echo filtering finished at `date` >> ${DIR}/${BASE_NAME}.log 
exit 0
