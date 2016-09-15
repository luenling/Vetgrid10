#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.18/samtools
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
DIR=`pwd`
BASE_NAME=$1
THREADS=10
BWA=/Volumes/Temp/Lukas/Tools/bwa-0.5.9/bwa

if [[ $1 == '-h' || $1 == '--help' ]]; then
    echo "usage: $0 BASE_NAME"
    echo "starts bwa alignment for BASE_NAME trimmed files, and joins results to sam file."
    echo "everyhing is logged in the files BASE_NAME.log and BASE_NAME_1/2_bwa.log."
    echo "the script does hardly do any checking, so use with utmost care."
    exit
fi


# trimming reads
if [[ ! -e ${DIR}/${BASE_NAME}_trimmed_1.gz && ! -e ${DIR}/${BASE_NAME}_trimmed_2.gz ]]; then
    echo "No trimmed files found for " ${BASE_NAME} >> ${DIR}/${BASE_NAME}.log
    exit 100
fi

### mapping trimmed reads to ref genome
echo mapping trimmed reads of $BASE_NAME with reference genome with  $THREADS  threads each at:  `date` >> ${DIR}/${BASE_NAME}.log
echo nohup nice time $BWA aln -l 200 -n 0.01 -o 2 -e 12 -d 12 -t $THREADS $REFGEN $DIR/${BASE_NAME}_trimmed_1.gz \> $DIR/${BASE_NAME}_trimmed_1.sai 2\> ${BASE_NAME}_1_bwa.log >> ${DIR}/${BASE_NAME}.log

nohup nice time $BWA aln -l 200 -n 0.01 -o 2 -e 12 -d 12 -t $THREADS $REFGEN $DIR/${BASE_NAME}_trimmed_1?(.gz) > $DIR/${BASE_NAME}_trimmed_1.sai 2> ${BASE_NAME}_1_bwa.log &

#nohup nice time sleep 10 >> ${DIR}/${BASE_NAME}.log &

echo Process ID $! at `date` >> ${DIR}/${BASE_NAME}.log
PIDS=( $! )

echo nohup nice time $BWA aln -l 200 -n 0.01 -o 2 -e 12 -d 12 -t $THREADS $REFGEN $DIR/${BASE_NAME}_trimmed_2.gz \> $DIR/${BASE_NAME}_trimmed_2.sai  2\> ${BASE_NAME}_2_bwa.log >> ${DIR}/${BASE_NAME}.log

nohup nice time $BWA aln -l 200 -n 0.01 -o 2 -e 12 -d 12 -t $THREADS $REFGEN $DIR/${BASE_NAME}_trimmed_2?(.gz) > $DIR/${BASE_NAME}_trimmed_2.sai 2> ${BASE_NAME}_2_bwa.log &

#nohup nice time sleep 10 >> ${DIR}/${BASE_NAME}.log &

echo Process ID $! at `date` >> ${DIR}/${BASE_NAME}.log

Pids=( ${PIDS[@]} $! )

### wait for the processes to finish
wait ${PIDS[0]}
[ $? ] || { echo BWA ${PIDS[0]} exited with error $? at `date`  >> ${DIR}/${BASE_NAME}.log; exit 100;  }
wait ${PIDS[1]}
[ $? ] || { echo BWA ${PIDS[1]} exited with error $? at `date`  >> ${DIR}/${BASE_NAME}.log; exit 100;  }

echo Mapping finished at `date` >> ${DIR}/${BASE_NAME}.log

exit

