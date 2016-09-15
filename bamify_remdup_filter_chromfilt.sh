#!/bin/bash
# turn on extended globbing behavior
# sorts, removes duplicates, filters and removes mates mapped to different chromosomes from a sam
shopt -s extglob
SCRIPTS=/Volumes/Temp/Lukas/Tools/Scripts
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
BASE_DIR=`pwd`
BASE_NAME=$1
echo bamify at `date`>> ${BASE_NAME}.log
bash $SCRIPTS/bamify_sort.sh $BASE_NAME
echo remove_dup_from_sorted.sh at `date`>> ${BASE_NAME}.log
bash $SCRIPTS/remove_dup_from_sorted.sh $BASE_NAME
echo remove_reads_on_diff_chrom.sh at `date`>> ${BASE_NAME}.log
bash $SCRIPTS/remove_reads_on_diff_chrom.sh $BASE_NAME
echo finished at `date` >> ${BASE_NAME}.log


