#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
PICARD=/Volumes/Temp/Lukas/Tools/picard/picard-tools-1.104
SNAPPY=/Volumes/Temp/Lukas/Tools/picard
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.19/samtools
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
LOGFILE=contamination.log
INFILE=$1
OUT=`basename $INFILE .bam`
OUTDIR=''
if [ $2 ]; then 
    OUTDIR=$2'/'
    if [ ! -d $2 ]; then mkdir $2 ; fi
fi  


echo converting bam to fastq at `date` >> ${LOGFILE}
echo nohup java -Xmx4g -classpath $SNAPPY -jar ${PICARD}/SamToFastq.jar VALIDATION_STRINGENCY=SILENT INPUT=$INFILE FASTQ=${OUTDIR}${OUT}_1.fastq SECOND_END_FASTQ=${OUTDIR}${OUT}_2.fastq  >> ${LOGFILE}
nohup java -Xmx4g -classpath $SNAPPY -jar ${PICARD}/SamToFastq.jar VALIDATION_STRINGENCY=SILENT INPUT=$INFILE FASTQ=${OUTDIR}${OUT}_1.fastq SECOND_END_FASTQ=${OUTDIR}${OUT}_2.fastq &
