#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
SCRIPTS=/Volumes/Temp/Lukas/Tools/Scripts
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.18/samtools
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
BWA=/Volumes/Temp/Lukas/Tools/bwa-0.5.9/bwa
BASE_DIR=`pwd`
BASE_NAMES=$1

echo Mapping and Removing Duplicates at `date`: >> $BASE_DIR/map_run.log
for i in $BASE_NAMES; do
    echo Mapping $i at `date` >> ${BASE_DIR}/map_run.log
    cd $i
    DIR=`pwd`
    echo nohup bash $SCRIPTS/mapping_generic.sh $i \& >> ${BASE_DIR}/map_run.log
    nohup bash $SCRIPTS/mapping_generic.sh $i &
    wait $!
    [ $? ] || { echo Mapping $i exited with error $? at `date`  >> ${BASE_DIR}/map_run.log; exit 100;  }
    echo Removing Duplicates for $i at `date` >> ${BASE_DIR}/map_run.log
    echo nohup bash $SCRIPTS/remove_duplicates.sh $i \& >> ${BASE_DIR}/map_run.log 
    nohup  bash $SCRIPTS/remove_duplicates.sh $i &
    cd ${BASE_DIR}
done


