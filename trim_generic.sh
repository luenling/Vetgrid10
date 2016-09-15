#!/bin/bash
# turn on extended globbing behavior
# Time Stamp: <>
shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.18/samtools
SCRIPTS=/Volumes/Temp/Lukas/Tools/Scripts
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
BASE_DIR=`pwd`
BASE_NAMES=''
THREADS=10
BWA=/Volumes/Temp1/Lukas/Tools/bwa-0.5.9/bwa
#MIDDLE=_sequence_
ENDING="_sequence.txt|.fq|.txt"

if [[ $1 == '-h' || $1 == '--help' ]]; then
    echo "usage: $0 BASE_NAMES"
    echo "goes through subdirs and creates trimmed files for the 9* type read file names"
    echo "the script does hardly do any checking, so use with utmost care."
    exit
fi

BASE_NAMES=$1
for i in $BASE_NAMES; do
    echo trimming $i at `date` >> ${BASE_DIR}/map_run.log
    cd ./$i
    DIR=`pwd`
    FILE=( ${i}*_[12]+(${ENDING}) )
 #   BASE_NAME=${FILE[0]}
 #   BASE_NAME=${FILE[0]%_[12]+(${ENDING})}
    BASE_NAME=$i
    if [[ ! -e ${DIR}/${BASE_NAME}_trimmed_1 ]]; then
	FQ_TYPE=`perl $SCRIPTS/fastqFormatDetect.pl ${FILE[1]} 2>&1`
	if [[ $FQ_TYPE =~ sanger ]]; then
	    FQ_TYPE=sanger
	else
	    FQ_TYPE=illumina
	fi
	echo Trimming the files ${FILE[@]} at `date` >> ${DIR}/${BASE_NAME}.log
	echo nohup time perl ${POPO}/basic-pipeline/trim-fastq.pl --input1 ${DIR}/${FILE[0]} --input2 ${DIR}/${FILE[1]} --output ${BASE_NAME}_trimmed --quality-threshold 18 --fastq-type $FQ_TYPE --min-length 50 --no-5p-trim >>  ${DIR}/${BASE_NAME}.log 
	nohup time perl ${POPO}/basic-pipeline/trim-fastq.pl --input1 ${DIR}/${FILE[0]} --input2 ${DIR}/${FILE[1]} --output ${BASE_NAME}_trimmed --quality-threshold 18 --fastq-type $FQ_TYPE --min-length 50 --no-5p-trim  >> ${DIR}/${BASE_NAME}_trimstats.log 2>> ${DIR}/${BASE_NAME}.log & 
    fi
    cd $BASE_DIR
done
	# echo number of trimmed sequences: >> ${DIR}/${BASE_NAME}.log 
	# grep -c "@" ${BASE_NAME}_trimmed_* >> ${DIR}/${BASE_NAME}.log 
echo finished starting processes at `date` >> ${BASE_DIR}/map_run.log
exit 0
