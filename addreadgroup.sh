#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
PICARD=/usr/local/share/java/picard.jar
#SNAPPY=/Volumes/Temp/Lukas/Tools/picard
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
FORMAT=ILLUMINA
LOGFILE=realign.log
INFILE=$1
LABEL=$2
if [ ! $LABEL ]; then LABEL="UNDEF"; fi
UNIT=$3
if [ ! $UNIT ]; then UNIT="UNDEF"; fi

IN_BASE=`basename $INFILE .bam`
SAMPLE=${IN_BASE#*tophat_}
SAMPLE=${SAMPLE%_accept*}

echo adding readgroup to $INFILE, label $LABEL, unit $UNIT at `date` >> $LOGFILE
echo java -jar ${PICARD} AddOrReplaceReadGroups INPUT=$INFILE OUTPUT=${IN_BASE}_RG.bam RGID="$SAMPLE" RGLB="$LABEL" RGPL="$FORMAT" RGSM="$SAMPLE" RGPU="$UNIT" VALIDATION_STRINGENCY=SILENT >> $LOGFILE
java -jar ${PICARD} AddOrReplaceReadGroups INPUT=$INFILE OUTPUT=${IN_BASE}_RG.bam RGLB="$LABEL" RGID="$SAMPLE" RGPL="$FORMAT" RGSM="$SAMPLE" RGPU="$UNIT" VALIDATION_STRINGENCY=SILENT 
echo indexing  ${IN_BASE}_RG.bam at `date` >> $LOGFILE
samtools index ${IN_BASE}_RG.bam
echo done addreadgroup $INFILE at `date` >> $LOGFILE


