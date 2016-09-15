#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
PICARD=/Volumes/Temp/Lukas/Tools/picard/picard-tools-1.79
SNAPPY=/Volumes/Temp/Lukas/Tools/picard
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.18/samtools
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
BWA=/Volumes/Temp/Lukas/Tools/bwa-0.5.9/bwa
DIR=`pwd`
BASE_NAME=$1


echo Converting, filtering, sorting, and removing duplicates at `date` >>  ${DIR}/${BASE_NAME}.log

# echo Checking number of mapped sequences: >>  ${DIR}/${BASE_NAME}.log
# #A=`zgrep -c '^@' $DIR/${BASE_NAME}_trimmed_1?(.gz)`
# A=`gunzip -c $DIR/${BASE_NAME}_trimmed_1?(.gz) | wc -l`
# #B=`zgrep -c '^@' $DIR/${BASE_NAME}_trimmed_2?(.gz)`
# B=`gunzip -c $DIR/${BASE_NAME}_trimmed_2?(.gz) | wc -l`
# C=`grep -cv '^@' $DIR/${BASE_NAME}.sam`
# A=`echo $A/4 | bc`
# B=`echo $B/4 | bc`
# D=`expr $A + $B - $C`

# echo $DIR/${BASE_NAME}_trimmed_1 : $A  >>  ${DIR}/${BASE_NAME}.log
# echo $DIR/${BASE_NAME}_trimmed_2 : $B  >>  ${DIR}/${BASE_NAME}.log
# echo $DIR/${BASE_NAME}.sam : $C  >>  ${DIR}/${BASE_NAME}.log
# echo Difference : $D >>  ${DIR}/${BASE_NAME}.log

# if [[ $D -eq 0 ]] ; then
#     echo removing $DIR/  $DIR/${BASE_NAME}_trimmed_1?(.gz) $DIR/${BASE_NAME}_trimmed_2?(.gz)  >>  ${DIR}/${BASE_NAME}.log
#     rm -f $DIR/${BASE_NAME}_trimmed_1?(.gz) $DIR/${BASE_NAME}_trimmed_2?(.gz) 
# else
#     echo not all sequences mapped >>  ${DIR}/${BASE_NAME}.log
#     exit 100
# fi

echo joining and converting sam files to bams at `date` >> ${DIR}/${BASE_NAME}.log
A='' 
for i in Lukas_${BASE_NAME}*.sam; do 
    A=$A" I="$i
done
echo java -Xmx5g -classpath $SNAPPY -jar ${PICARD}/MergeSamFiles.jar $A O=${BASE_NAME}.bam USE_THREADING=TRUE VALIDATION_STRINGENCY=SILENT >> ${DIR}/${BASE_NAME}.log
nohup java -Xmx5g -classpath $SNAPPY -jar ${PICARD}/MergeSamFiles.jar $A O=${BASE_NAME}.bam USE_THREADING=TRUE VALIDATION_STRINGENCY=SILENT &
wait $!

[ $? ] || { echo Sorting/Merging exited with error $? at  `date`  >> ${DIR}/${BASE_NAME}.log   ; }

echo $SAMTOOLS flagstats at `date`: >> ${DIR}/${BASE_NAME}.log
$SAMTOOLS flagstat  $DIR/${BASE_NAME}.bam >>  ${DIR}/${BASE_NAME}.log

[ $? ] || { echo samtools exited with error $? at `date`  >> ${DIR}/${BASE_NAME}.log; exit 100;  }
C=0
for i in Lukas_${BASE_NAME}*.sam; do 
#    D=`echo $i | wc -m `
    D=`grep -cv '^@' $i`
    let "C += D"
done
A=`$SAMTOOLS view $DIR/${BASE_NAME}.bam | wc -l`
echo lines in ${BASE_NAME} sams : $C >> ${DIR}/${BASE_NAME}.log
echo lines in ${BASE_NAME}.bam : $A >> ${DIR}/${BASE_NAME}.log

if [[ $A -eq $C ]] ; then
  echo removing  Lukas_${BASE_NAME}\*.sam >> ${DIR}/${BASE_NAME}.log 
  rm -f Lukas_${BASE_NAME}*.sam
fi


exit 0