shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.19/samtools
PICARD=/Volumes/Temp/Lukas/Tools/picard/picard-tools-1.79
SNAPPY=/Volumes/Temp/Lukas/Tools/picard
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
DIR=`pwd`
BASE_NAME=$1

echo java -Xmx5g -classpath $SNAPPY -jar ${PICARD}/SortSam.jar I=${BASE_NAME}_trimmed_1and2.sam O=${BASE_NAME}_trimmed_1and2.bam SO=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=TRUE >> ${DIR}/${BASE_NAME}.log
java -Xmx5g -classpath $SNAPPY -jar ${PICARD}/SortSam.jar I=${BASE_NAME}_trimmed_1and2.sam O=${BASE_NAME}_trimmed_1and2.bam SO=coordinate VALIDATION_STRINGENCY=SILENT  CREATE_INDEX=TRUE
C=`grep -cv '^@' $DIR/${BASE_NAME}_trimmed_1and2.sam`
A=`$SAMTOOLS view $DIR/${BASE_NAME}_trimmed_1and2.bam | wc -l`

D=`expr $A - $C`
if [[ $D -eq 0 ]] ; then
    echo removing  $DIR/${BASE_NAME}_trimmed_1and2.sam >>  ${DIR}/${BASE_NAME}.log
    rm -f $DIR/${BASE_NAME}_trimmed_1and2.sam 
else
    echo not all sequences in bam >>  ${DIR}/${BASE_NAME}.log
    exit 100
fi

$SAMTOOLS index ${BASE_NAME}_trimmed_1and2.bam