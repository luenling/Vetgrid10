#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
SCRIPTS=/Volumes/Temp/Lukas/Tools/Scripts
BASE_DIR=`pwd`
BASE_NAME=`basename $1 .bam`
#echo bamify at `date`>> ${BASE_NAME}.log
#bash $SCRIPTS/bamify_sort.sh $BASE_NAME

samtools flagstat ${BASE_NAME}.bam >> ${BASE_NAME}.log &
echo remove_dups_and_filter_from_sorted.sh at `date`>> ${BASE_NAME}.log
bash $SCRIPTS/remove_dups_and_filter_from_sorted.sh $BASE_NAME
FILE=${BASE_NAME}'.picnodupl.filtered.mq20_chrfilt.bam'
bash $SCRIPTS/realign.sh $FILE realign.log
echo finished at `date` >> ${BASE_NAME}.log


