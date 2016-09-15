#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
SCRIPTS=/Volumes/Temp/Lukas/Tools/Scripts
FORMAT=SANGER
LOGFILE=realign.log
INFILE=$1
LABEL=`basename $INFILE .bam`
REF=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
BASE_NAME=`basename $INFILE .bam`
if  test $2  ; then REF=$2; fi


echo converting  $1 to sanger at `date` >> ${DIR}/${BASE_NAME}.log
echo python $SCRIPTS/bam_illumina2sanger.py --in $INFILE --out ${BASE_NAME}_sanger.bam
python $SCRIPTS/bam_illumina2sanger.py --in $INFILE --out ${BASE_NAME}_sanger.bam 
echo started adding readgroups at `date` >> $LOGFILE
echo bash $SCRIPTS/addreadgroup.sh ${BASE_NAME}_sanger.bam $LABEL >> ${DIR}/${BASE_NAME}.log
bash $SCRIPTS/addreadgroup.sh  ${BASE_NAME}_sanger.bam $LABEL
echo started at realigning at `date` >> ${DIR}/${BASE_NAME}.log
echo bash $SCRIPTS/realign.sh ${BASE_NAME}_sanger_RG.bam $LOGFILE $REF >> $LOGFILE
nohup bash $SCRIPTS/realign.sh ${BASE_NAME}_sanger_RG.bam $LOGFILE $REF &
