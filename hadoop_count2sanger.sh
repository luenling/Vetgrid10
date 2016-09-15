#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.18/samtools
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
SCRIPTS=/Volumes/Temp/Lukas/Tools/Scripts
DIR=`pwd`
BASE_NAME=$1

echo Flagstat DistMap_output_Paired_end_reads.bam at `date` >> ${DIR}/${BASE_NAME}.log
nohup $SAMTOOLS flagstat DistMap_output_Paired_end_reads.bam  >> ${DIR}/${BASE_NAME}.log &
echo converting  DistMap_output_Paired_end_reads.bam  to sanger at `date` >> ${DIR}/${BASE_NAME}.log
python $SCRIPTS/bam_illumina2sanger.py --in DistMap_output_Paired_end_reads.bam --out ${BASE_NAME}_sanger.bam
echo Flagstat ${BASE_NAME}_sanger.bam at `date` >> ${DIR}/${BASE_NAME}.log
nohup $SAMTOOLS flagstat ${BASE_NAME}_sanger.bam  >> ${DIR}/${BASE_NAME}.log &
