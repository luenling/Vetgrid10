#!/bin/bash
# turn on extended globbing behavior
# Time Stamp: <>
shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
SCRIPTS=/Volumes/Temp/Lukas/Tools/Scripts
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
BASE_DIR=`pwd`
BASE_NAMES=''
MIDDLE=_sequence_
ENDING="_sequence.txt|.fq|.txt"

if [[ $1 == '-h' || $1 == '--help' ]]; then
    echo "usage: $0 \"BASE_NAMES\""
    echo "goes through subdirs and creates trimmed files for the files created by Run15-19"
    echo "the script does hardly any checking, so use with utmost care."
    echo "eg. bash trim_files_for_Runs.sh \"Run1 Run2 Run3\""
    exit
fi


BASE_NAMES=$1
for i in $BASE_NAMES; do
    echo trimming $i at `date` >> ${BASE_DIR}.log
    cd ./$i
    DIR=`pwd`
    FILE=( ${i:0:3}_[12]${MIDDLE}${i:4}?(${ENDING}) )
 #   BASE_NAME=${FILE[0]}
 #   BASE_NAME=${FILE[0]%_[12]+(${ENDING})}
    BASE_NAME=$i
    if [[ ! -e ${DIR}/${BASE_NAME}_trimmed_1 ]]; then
	# runs perl script to determine quality encoding
	# the script returns either "sanger" or "illumina"
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
echo finished starting processes at `date` >> ${BASE_DIR}.log
exit 0
