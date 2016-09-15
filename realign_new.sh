#!/bin/bash
# turn on extended globbing behavior
# takes one required and two optional arguments: infile_name  ( logfile_name reference_genome_filename )
# by default writes infile_real.bam and realign.log 
shopt -s extglob
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
GATK=/Volumes/Temp/Lukas/Tools/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar
KNOWN=/Volumes/Temp/Lukas/Data/DGRP_vcf/freeze2.bins.indels.vcf.gz
INFILE=$1
IN_BASE=`basename $INFILE .bam`
LOGFILE=$2
if test $3 ; then REFGEN=$3; fi
if [ ! $LOGFILE ]; then LOGFILE="realign.log"; fi
if [ ! $REFGEN ]; then REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa; fi
echo at `date`: >> $LOGFILE
echo java -Xmx4g -jar $GATK  -dt NONE -R $REFGEN -nt 1 -T RealignerTargetCreator -I $INFILE -o ${IN_BASE}.intervals --maxIntervalSize 400 --minReadsAtLocus 10  -known $KNOWN >> $LOGFILE
java -Xmx4g -jar $GATK -R $REFGEN -nt 1 -T RealignerTargetCreator -I $INFILE -o ${IN_BASE}.intervals --maxIntervalSize 400 --minReadsAtLocus 10  -known $KNOWN 2>&1 >> ${IN_BASE}_realign.log
echo java -Xmx5g -jar $GATK -R $REFGEN -T IndelRealigner -dt NONE -I $INFILE -o ${IN_BASE}_real.bam -targetIntervals ${IN_BASE}.intervals  -entropy 0.05 -LOD 3 -maxConsensuses 50 -greedy 250 -model USE_READS -known $KNOWN >> $LOGFILE
java -Xmx5g -jar $GATK -R $REFGEN -T IndelRealigner -dt NONE -I $INFILE -o ${IN_BASE}_real.bam -targetIntervals ${IN_BASE}.intervals  -entropy 0.05 -LOD 3 -maxConsensuses 50 -greedy 250 -model USE_READS -known $KNOWN 2>&1 >> ${IN_BASE}_realign.log
echo done realigning $INFILE at `date` >> $LOGFILE
