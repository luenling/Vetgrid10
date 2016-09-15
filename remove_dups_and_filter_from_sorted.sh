#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
PICARD=/Volumes/Temp/Lukas/Tools/picard/picard-tools-1.104
SNAPPY=/Volumes/Temp/Lukas/Tools/picard
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.19/samtools
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
BWA=/Volumes/Temp/Lukas/Tools/bwa-0.5.9/bwa
DIR=`pwd`
BASE_NAME=`basename $1 .bam`
### remove duplicates to get rid of PCR artefacts 
### prob. of having two reads starting at the same position is low 

echo Removing duplicates using Picard at `date` >> ${DIR}/${BASE_NAME}.log
echo nohup java -Xmx5g -Dsnappy.disable=true -jar ${PICARD}/MarkDuplicates.jar REMOVE_DUPLICATES=true I=${DIR}/${BASE_NAME}.bam O=${DIR}/${BASE_NAME}.picnodupl.bam M=${DIR}/${BASE_NAME}.picmetrics.txt VALIDATION_STRINGENCY=SILENT >> ${DIR}/${BASE_NAME}.log 

nohup java -Xmx5g -Dsnappy.disable=true -jar ${PICARD}/MarkDuplicates.jar REMOVE_DUPLICATES=true I=${DIR}/${BASE_NAME}.bam O=${DIR}/${BASE_NAME}.picnodupl.bam M=${DIR}/${BASE_NAME}.picmetrics.txt VALIDATION_STRINGENCY=SILENT &

wait $!

[ $? ] || { echo  MarkDuplicates exit with error $? at  `date`  >> ${DIR}/${BASE_NAME}.log   ; }

### filter bam file for correct read flags and min mapping quality

echo bam filtering using samtools at `date` >> ${DIR}/${BASE_NAME}.log
echo nohup $SAMTOOLS view  -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -bh  ${DIR}/${BASE_NAME}.picnodupl.bam \> ${DIR}/${BASE_NAME}.picnodupl.filtered.mq20.bam >> ${DIR}/${BASE_NAME}.log 

nohup $SAMTOOLS view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -bh ${DIR}/${BASE_NAME}.picnodupl.bam > ${DIR}/${BASE_NAME}.picnodupl.filtered.mq20.bam &

wait $!

[ $? ] || { echo  samtools exit with error $? at  `date`  >> ${DIR}/${BASE_NAME}.log ; }

rm -f ${DIR}/${BASE_NAME}.picnodupl.bam

echo flagstat of $DIR/${BASE_NAME}.picnodupl.filtered.mq20.bam : >>  ${DIR}/${BASE_NAME}.log
$SAMTOOLS flagstat  $DIR/${BASE_NAME}.picnodupl.filtered.mq20.bam >>  ${DIR}/${BASE_NAME}.log

$SAMTOOLS index $DIR/${BASE_NAME}.picnodupl.filtered.mq20.bam

echo idxstat of $DIR/${BASE_NAME}.picnodupl.filtered.mq20.bam : >>  ${DIR}/${BASE_NAME}.log

samtools idxstats $DIR/${BASE_NAME}.picnodupl.filtered.mq20.bam >>  ${DIR}/${BASE_NAME}.log

echo Sorting, removing duplicates and bam conversion finished at `date` >> ${DIR}/${BASE_NAME}.log

IN=${BASE_NAME}'.picnodupl.filtered.mq20'

[ -e ${IN}.bam ] || exit 1
echo Filtering out reads mapped to different chromosomes at `date` >> ${DIR}/${BASE_NAME}.log 
$SAMTOOLS view -h  $IN.bam |  awk ' (/^@/ || $7=="=") {print;}' | $SAMTOOLS  view -bSh - > ${IN}_chrfilt.bam
echo flagstat of ${IN}_chrfilt.bam at `date` >> ${DIR}/${BASE_NAME}.log 
$SAMTOOLS flagstat ${IN}_chrfilt.bam >> ${DIR}/${BASE_NAME}.log 

$SAMTOOLS index ${IN}_chrfilt.bam

echo idxstat of ${IN}_chrfilt.bam : >>  ${DIR}/${BASE_NAME}.log


samtools idxstats ${IN}_chrfilt.bam >> ${DIR}/${BASE_NAME}.log

rm -f ${IN}.bam

echo filtering finished at `date` >> ${DIR}/${BASE_NAME}.log 
exit 0
