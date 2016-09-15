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
BASE_NAME=$1

# ### joining reads to create sam file
# echo Joining reads to sam file  at `date` >> ${DIR}/${BASE_NAME}.log
# echo $BWA sampe $REFGEN $DIR/${BASE_NAME}_trimmed_1.sai $DIR/${BASE_NAME}_trimmed_2.sai $DIR/${BASE_NAME}_trimmed_1 $DIR/${BASE_NAME}_trimmed_2 \> $DIR/${BASE_NAME}_trimmed_1and2.sam >> ${DIR}/${BASE_NAME}.log

# [ $? ] || { echo BWA sampe exited with error $? at `date`  >> ${DIR}/${BASE_NAME}.log; exit 100;  } 

# $BWA sampe $REFGEN $DIR/${BASE_NAME}_trimmed_1.sai $DIR/${BASE_NAME}_trimmed_2.sai $DIR/${BASE_NAME}_trimmed_1?(.gz) $DIR/${BASE_NAME}_trimmed_2?(.gz) > $DIR/${BASE_NAME}_trimmed_1and2.sam

# [ $? ] || { echo BWA sampe exited with error $? at `date`  >> ${DIR}/${BASE_NAME}.log; exit 100;  } 

# echo Finished joining at `date` >>  ${DIR}/${BASE_NAME}.log

# echo Checking number of mapped sequences: >>  ${DIR}/${BASE_NAME}.log
# #A=`zgrep -c '^@' $DIR/${BASE_NAME}_trimmed_1?(.gz)`
# A=`gunzip -c $DIR/${BASE_NAME}_trimmed_1?(.gz) | wc -l`
# #B=`zgrep -c '^@' $DIR/${BASE_NAME}_trimmed_2?(.gz)`
# B=`gunzip -c $DIR/${BASE_NAME}_trimmed_2?(.gz) | wc -l`
# C=`grep -cv '^@' $DIR/${BASE_NAME}_trimmed_1and2.sam`
# A=`echo $A/4 | bc`
# B=`echo $B/4 | bc`
# D=`expr $A + $B - $C`

# echo $DIR/${BASE_NAME}_trimmed_1 : $A  >>  ${DIR}/${BASE_NAME}.log
# echo $DIR/${BASE_NAME}_trimmed_2 : $B  >>  ${DIR}/${BASE_NAME}.log
# echo $DIR/${BASE_NAME}_trimmed_1and2.sam : $C  >>  ${DIR}/${BASE_NAME}.log
# echo Difference : $D >>  ${DIR}/${BASE_NAME}.log

# if [[ $D -eq 0 ]] ; then
#     echo removing $DIR/  $DIR/${BASE_NAME}_trimmed_1?(.gz) $DIR/${BASE_NAME}_trimmed_2?(.gz)  >>  ${DIR}/${BASE_NAME}.log
#     rm -f $DIR/${BASE_NAME}_trimmed_1?(.gz) $DIR/${BASE_NAME}_trimmed_2?(.gz) 
# else
#     echo not all sequences mapped >>  ${DIR}/${BASE_NAME}.log
#     exit 100
# fi

# echo convert sam to bam file and remove sam at `date` >> ${DIR}/${BASE_NAME}.log
# echo $SAMTOOLS view -Sbh $DIR/${BASE_NAME}_trimmed_1and2.sam \> $DIR/${BASE_NAME}_trimmed_1and2.bam  >>  ${DIR}/${BASE_NAME}.log
# $SAMTOOLS view -Sbh $DIR/${BASE_NAME}_trimmed_1and2.sam > $DIR/${BASE_NAME}_trimmed_1and2.bam
# echo $SAMTOOLS flagstats at `date`: >> ${DIR}/${BASE_NAME}.log
# $SAMTOOLS flagstat  $DIR/${BASE_NAME}_trimmed_1and2.bam >>  ${DIR}/${BASE_NAME}.log

# [ $? ] || { echo samtools exited with error $? at `date`  >> ${DIR}/${BASE_NAME}.log; exit 100;  }

# A=`$SAMTOOLS view $DIR/${BASE_NAME}_trimmed_1and2.bam | wc -l`
# echo lines in ${BASE_NAME}_trimmed_1and2.bam : $A >> ${DIR}/${BASE_NAME}.log

# if [[ $A -eq $C ]] ; then
#   echo removing  ${BASE_NAME}_trimmed_1and2.sam >> ${DIR}/${BASE_NAME}.log 
#   rm -f ${BASE_NAME}_trimmed_1and2.sam
# fi



### remove duplicates to get rid of PCR artefacts 
### prob. of having two reads starting at the same position is low 
# sort
echo Sort files using Picard at `date` >> ${DIR}/${BASE_NAME}.log
echo nohup java -Xmx4g -Dsnappy.disable=true -jar ${PICARD}/SortSam.jar I=$DIR/${BASE_NAME}.bam O=${DIR}/${BASE_NAME}.picsort.bam SO=coordinate VALIDATION_STRINGENCY=SILENT  >> ${DIR}/${BASE_NAME}.log 

nohup java -Xmx4g -Dsnappy.disable=true -jar ${PICARD}/SortSam.jar I=$DIR/${BASE_NAME}.bam O=${DIR}/${BASE_NAME}.picsort.bam SO=coordinate VALIDATION_STRINGENCY=SILENT &

wait $!
[ $? ] || { echo  SortSam exit with error $? at  `date`  >> ${DIR}/${BASE_NAME}.log   ; exit; }

echo Removing duplicates using Picard at `date` >> ${DIR}/${BASE_NAME}.log
echo nohup java -Xmx4g -Dsnappy.disable=true -jar ${PICARD}/MarkDuplicates.jar REMOVE_DUPLICATES=true I=${DIR}/${BASE_NAME}.picsort.bam O=${DIR}/${BASE_NAME}.picnodupl.bam M=${DIR}/${BASE_NAME}.picmetrics.txt VALIDATION_STRINGENCY=SILENT >> ${DIR}/${BASE_NAME}.log 

nohup java -Xmx4g -Dsnappy.disable=true -jar ${PICARD}/MarkDuplicates.jar REMOVE_DUPLICATES=true I=${DIR}/${BASE_NAME}.picsort.bam O=${DIR}/${BASE_NAME}.picnodupl.bam M=${DIR}/${BASE_NAME}.picmetrics.txt VALIDATION_STRINGENCY=SILENT &

wait $!

[ $? ] || { echo  MarkDuplicates exit with error $? at  `date`  >> ${DIR}/${BASE_NAME}.log   ; }

### convert into a bam file

echo bam filtering using samtools at `date` >> ${DIR}/${BASE_NAME}.log
echo nohup $SAMTOOLS view  -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -bh  ${DIR}/${BASE_NAME}.picnodupl.bam \> ${DIR}/${BASE_NAME}.picnodupl.filtered.mq20.bam >> ${DIR}/${BASE_NAME}.log 

nohup $SAMTOOLS view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -bh ${DIR}/${BASE_NAME}.picnodupl.bam > ${DIR}/${BASE_NAME}.picnodupl.filtered.mq20.bam &

wait $!

[ $? ] || { echo  samtools exit with error $? at  `date`  >> ${DIR}/${BASE_NAME}.log ; }


echo flagstat of $DIR/${BASE_NAME}.picnodupl.filtered.mq20.bam : >>  ${DIR}/${BASE_NAME}.log
$SAMTOOLS flagstat  $DIR/${BASE_NAME}.picnodupl.filtered.mq20.bam >>  ${DIR}/${BASE_NAME}.log

$SAMTOOLS index $DIR/${BASE_NAME}.picnodupl.filtered.mq20.bam

echo idxstat of $DIR/${BASE_NAME}.picnodupl.filtered.mq20.bam : >>  ${DIR}/${BASE_NAME}.log

samtools idxstats $DIR/${BASE_NAME}.picnodupl.filtered.mq20.bam >>  ${DIR}/${BASE_NAME}.log

echo Sorting, removing duplicates and bam conversion finished at `date` >> ${DIR}/${BASE_NAME}.log

exit 0