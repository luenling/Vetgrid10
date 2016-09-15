#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
SCRIPTS=/Volumes/Temp/Lukas/Tools/Scripts
BASE_DIR=`pwd`
BASE_NAME=`basename $1 .bam`
#echo bamify at `date`>> ${BASE_NAME}.log
#bash $SCRIPTS/bamify_sort.sh $BASE_NAME

samtools flagstat ${BASE_NAME}.bam >> ${BASE_NAME}.log &
echo remove_dup_from_sorted.sh at `date`>> ${BASE_NAME}.log
bash $SCRIPTS/remove_dup_from_sorted.sh $BASE_NAME
echo remove_reads_on_diff_chrom.sh at `date`>> ${BASE_NAME}.log
bash $SCRIPTS/remove_reads_on_diff_chrom.sh $BASE_NAME
echo finished reomving reads on different chromosomes at `date` >> ${BASE_NAME}.log
FILE=${BASE_NAME}'.picnodupl.filtered.mq20_chrfilt.bam'
bash $SCRIPTS/add_readgroup_realign.sh $FILE ${BASE_NAME} ${BASE_NAME}.log
echo finished at `date` >> ${BASE_NAME}.log


