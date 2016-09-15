#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.19/samtools
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
DIR=`pwd`
BASE_NAME=$1
TRIMFILES=$2
THREADS=10
BWA=/Volumes/Temp/Lukas/Tools/bwa-0.5.9/bwa

if [[ $1 == '-h' || $1 == '--help' ]]; then
    echo "usage: $0"
    echo "creates trimmed fq files,  and joins results to sam file."
    echo "everyhing is logged in the files BASE_NAME.log and BASE_NAME_1/2_bwa.log."
    echo "the script does hardly do any checking, so use with utmost care."
    exit
fi


# trimming reads
FILE=( ${TRIMFILES} ) 
#echo ${FILE[0]} ${FILE[1]}
echo time perl ${POPO}/basic-pipeline/trim-fastq_gz.pl --input1 ${DIR}/${FILE[0]} --input2 ${DIR}/${FILE[1]} --output ${BASE_NAME}_trimmed --quality-threshold 18 --fastq-type illumina --min-length 50 --no-5p-trim >>  ${DIR}/${BASE_NAME}.log

time perl ${POPO}/basic-pipeline/trim-fastq_gz.pl --input1 ${DIR}/${FILE[0]} --input2 ${DIR}/${FILE[1]} --output ${BASE_NAME}_trimmed --quality-threshold 18 --fastq-type illumina --min-length 50  --no-5p-trim 2>> ${DIR}/${BASE_NAME}.log

# ### joining reads to create sam file
echo Joining reads to sam file  at `date` >> ${DIR}/${BASE_NAME}.log
echo $BWA sampe $REFGEN $DIR/${BASE_NAME}_trimmed_1.sai $DIR/${BASE_NAME}_trimmed_2.sai $DIR/${BASE_NAME}_trimmed_1 $DIR/${BASE_NAME}_trimmed_2 \> $DIR/${BASE_NAME}_trimmed_1and2.sam >> ${DIR}/${BASE_NAME}.log 

$BWA sampe $REFGEN $DIR/${BASE_NAME}_trimmed_1.sai $DIR/${BASE_NAME}_trimmed_2.sai $DIR/${BASE_NAME}_trimmed_1 $DIR/${BASE_NAME}_trimmed_2 > $DIR/${BASE_NAME}_trimmed_1and2.sam

echo Finished joining at `date` >>  ${DIR}/${BASE_NAME}.log


echo Checking number of mapped sequences: >>  ${DIR}/${BASE_NAME}.log
A=`grep -c '@' $DIR/${BASE_NAME}_trimmed_1`
B=`grep -c '@' $DIR/${BASE_NAME}_trimmed_2`
C=`grep -vc '@' $DIR/${BASE_NAME}_trimmed_1and2.sam`
D=`expr $A + $B - $C`

echo $DIR/${BASE_NAME}_trimmed_1 : $A  >>  ${DIR}/${BASE_NAME}.log
echo $DIR/${BASE_NAME}_trimmed_2 : $B  >>  ${DIR}/${BASE_NAME}.log
echo $DIR/${BASE_NAME}_trimmed_1and2.sam : $C  >>  ${DIR}/${BASE_NAME}.log
echo Difference : $D >>  ${DIR}/${BASE_NAME}.log

if [[ $D -eq 0 ]] ; then
    echo removing $DIR/  $DIR/${BASE_NAME}_trimmed_1 $DIR/${BASE_NAME}_trimmed_2  >>  ${DIR}/${BASE_NAME}.log
    rm -f $DIR/${BASE_NAME}_trimmed_1 $DIR/${BASE_NAME}_trimmed_2 
else
    echo not all sequences mapped >>  ${DIR}/${BASE_NAME}.log
    exit 100
fi

bash /Volumes/Temp/Lukas/Tools/Scripts/sam_bam_sort_index.sh ${BASE_NAME}
# $SAMTOOLS view -Shb $DIR/${BASE_NAME}_trimmed_1and2.sam > $DIR/${BASE_NAME}_trimmed_1and2.bam

exit

