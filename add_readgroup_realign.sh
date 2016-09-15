#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
SCRIPTS=/Volumes/Temp/Lukas/Tools/Scripts
FORMAT=SANGER
LOGFILE=realign.log
INFILE=$1
LABEL=$2
if  test $3  ; then LOGFILE=$3; fi
echo started adding readgroups at `date` >> $LOGFILE
echo bash $SCRIPTS/addreadgroup.sh $INFILE $LABEL >> $LOGFILE
bash $SCRIPTS/addreadgroup.sh $INFILE $LABEL
IN_BASE=`basename $INFILE .bam`
echo started at realigning at `date` >> $LOGFILE
echo bash $SCRIPTS/realign.sh ${IN_BASE}_RG.bam $LOGFILE >> $LOGFILE
nohup bash $SCRIPTS/realign.sh ${IN_BASE}_RG.bam $LOGFILE &
