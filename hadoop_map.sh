#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.18/samtools
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
PICARD=/Volumes/Temp/Lukas/Tools/picard/picard-tools-1.79
BWA=/Volumes/Temp/Lukas/Tools/bwa-0.7.5a/bwa
#BWA_OPTIONS="\"-o 2 -n 0.01 -l200 -e 12 -d 12 \""
BWA_OPTIONS="\"-o 1 -n 0.01 -l 200 -e 12 -d 12\""
TRIMFILES="\""$1"\""
OUT_DIR=`pwd`
QUEUE=pg1
# pg1 pg2 pg3 
JOBDESC="\"Lukas mapping $2\""
BASE_NAME=$2

echo starting hadoop mapping on `date` >> ${BASE_NAME}.log
echo nohup /Volumes/cluster/DistMap_v1.1/distmap --reference-fasta $REFGEN --input $TRIMFILES --output $OUT_DIR --mapper bwa --mapper-path $BWA --picard-mergesamfiles-jar $PICARD/MergeSamFiles.jar --picard-sortsam-jar $PICARD/SortSam.jar --mapper-args $BWA_OPTIONS --output-format bam --queue-name $QUEUE --job-desc $JOBDESC  1\> distmap.out 2\> distmap.stderr.out \& >> ${BASE_NAME}.log


#nohup /Volumes/cluster/DistMap_v1.1/distmap --reference-fasta $REFGEN --input $TRIMFILES --output $OUT_DIR --mapper bwa --mapper-path $BWA --picard-mergesamfiles-jar $PICARD/MergeSamFiles.jar --picard-sortsam-jar $PICARD/SortSam.jar --mapper-args $BWA_OPTIONS --output-format bam --queue-name $QUEUE --job-desc $JOBDESC  1> distmap.out 2> distmap.stderr.out &

