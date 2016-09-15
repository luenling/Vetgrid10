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
BASE_NAME=females_pileup
VIENNA2010=/Volumes/Temp/Lukas/Data/Vienna_2010
POPICS25=`find $VIENNA2010/PopI/ColdShock25 -name "*chrfilt.bam"`
POPICS18=`find $VIENNA2010/PopI/ColdShock18 -name "*chrfilt_illu*.bam"`
POPIEL=`find $VIENNA2010/PopI/EarlyLate -name "*earlylate_merged.bam"`
POPIICS25=`find $VIENNA2010/PopII/ColdShock25 -name "*chrfilt.bam"`
POPIICS18=`find $VIENNA2010/PopII/ColdShock18/ -name "*chrfilt_illu*.bam"`
POPIIEL=`find $VIENNA2010/PopII/EarlyLate -name "*earlylate_merged.bam"`
POPIIICS25=`find $VIENNA2010/PopIII/ColdShock25 -name "*chrfilt.bam"`
POPIIIEL=`find $VIENNA2010/PopIII/EarlyLate -name "*earlylate_merged.bam"`
FEMFILE=females_pI25_pII25_pIII25_pI18_pII18_pIel_pIIel_pIIIel_unfiltered
### joining reads to create sam file
echo Mpileup of bam files for coldshock Vienna 2010  at `date` >> ${DIR}/${BASE_NAME}.log
echo $SAMTOOLS mpileup -6 -B -Q 0 -f $REFGEN $POPICS25 $POPIICS25 $POPIIICS25 $POPICS18 $POPIICS18 $POPIEL $POPIIEL $POPIIIEL \> $FEMFILE.mpileup >> ${DIR}/${BASE_NAME}.log
$SAMTOOLS mpileup -6 -B -Q 0 -f $REFGEN $POPICS25 $POPIICS25 $POPIIICS25 $POPICS18 $POPIICS18 $POPIEL $POPIIEL $POPIIIEL > $FEMFILE.mpileup
[ $? ] || { echo samtools pileup exited with error $? at `date`  >> ${DIR}/${BASE_NAME}.log; exit 100;  } 
echo Mpileup of bam files for coldshock Vienna 2010 finished at `date` >> ${DIR}/${BASE_NAME}.log 
echo java -Xmx12g -jar $POPO2/mpileup2sync.jar --input ${FEMFILE}.mpileup --output ${FEMFILE}_q20.sync --fastq-type sanger --min-qual 20 --threads 15  >> ${DIR}/${BASE_NAME}.log 
java -Xmx12g -jar $POPO2/mpileup2sync.jar --input ${FEMFILE}.mpileup --output ${FEMFILE}_q20.sync --fastq-type sanger --min-qual 20 --threads 15
[ $? ] || { echo mpileup2sync.jar exited with error $? at `date`  >> ${DIR}/${BASE_NAME}.log; exit 100;  } 
echo Syncfile generation for coldshock Vienna 2010 finished at `date` >> ${DIR}/${BASE_NAME}.log
exit 0