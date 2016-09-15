#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
DIR=`pwd`
BASE_NAMES=$1
THREADS=10
BWA=/Volumes/Temp/Lukas/Tools/bwa-0.5.9/bwa

if [[ $1 == '-h' || $1 == '--help' ]]; then
    echo "usage: $0"
    echo "creates trimmed fq files, starts bwa alignment, and joins results to sam file."
    echo "everyhing is logged in the files BASE_NAME.log and BASE_NAME_1/2_bwa.log."
    echo "the script does hardly do any checking, so use with utmost care."
    exit
fi


# trimming reads
if [[ ! -e ${DIR}/${BASE_NAME}_trimmed_1 ]]; then
    echo "triming " ${BASE_NAME} >> ${DIR}/${BASE_NAME}.log
    echo time perl ${POPO}/basic-pipeline/trim-fastq.pl --input1 ${DIR}/${BASE_NAME}_1_sequence.fq --input2 ${DIR}/${BASE_NAME}_2_sequence.fq --output ${BASE_NAME}_trimmed --quality-threshold 18 --fastq-type illumina --min-length 50 >>  ${DIR}/${BASE_NAME}.log
    time perl ${POPO}/basic-pipeline/trim-fastq.pl --input1 ${DIR}/${BASE_NAME}_1_sequence.fq --input2 ${DIR}/${BASE_NAME}_2_sequence.fq --output ${BASE_NAME}_trimmed --quality-threshold 18 --fastq-type illumina --min-length 50 2>> ${DIR}/${BASE_NAME}.log
    echo removing ${DIR}/${BASE_NAME}_1_sequence.fq ${DIR}/${BASE_NAME}_2_sequence.fq >> ${DIR}/${BASE_NAME}.log
    rm -f ${DIR}/${BASE_NAME}_1_sequence.fq ${DIR}/${BASE_NAME}_2_sequence.fq
fi

### mapping trimmed reads to ref genome
echo "mapping trimmed reads of 9390__81KGWABXX_2_s_2 with reference genome with " $THREADS " threads each at: " `date` >> ${DIR}/${BASE_NAME}.log
echo nohup nice time $BWA aln -l 200 -n 0.01 -o 2 -e 12 -d 12 -t $THREADS $REFGEN $DIR/${BASE_NAME}_trimmed_1 \> $DIR/${BASE_NAME}_trimmed_1.sai 2\> ${BASE_NAME}_1_bwa.log >> ${DIR}/${BASE_NAME}.log

nohup nice time $BWA aln -l 200 -n 0.01 -o 2 -e 12 -d 12 -t $THREADS $REFGEN $DIR/${BASE_NAME}_trimmed_1 > $DIR/${BASE_NAME}_trimmed_1.sai 2> ${BASE_NAME}_1_bwa.log &

#nohup nice time sleep 10 >> ${DIR}/${BASE_NAME}.log &

echo Process ID $! at `date` >> ${DIR}/${BASE_NAME}.log
PIDS=( $! )

echo nohup nice time $BWA aln -l 200 -n 0.01 -o 2 -e 12 -d 12 -t $THREADS $REFGEN $DIR/${BASE_NAME}_trimmed_2 \> $DIR/${BASE_NAME}_trimmed_2.sai  2\> ${BASE_NAME}_2_bwa.log >> ${DIR}/${BASE_NAME}.log

nohup nice time $BWA aln -l 200 -n 0.01 -o 2 -e 12 -d 12 -t $THREADS $REFGEN $DIR/${BASE_NAME}_trimmed_2 > $DIR/${BASE_NAME}_trimmed_2.sai 2> ${BASE_NAME}_2_bwa.log &

#nohup nice time sleep 10 >> ${DIR}/${BASE_NAME}.log &

echo Process ID $! at `date` >> ${DIR}/${BASE_NAME}.log

Pids=( ${PIDS[@]} $! )

### wait for the processes to finish
wait ${PIDS[0]}
$? || { echo BWA ${PIDS[0]} exited with error $? at `date`  >> ${DIR}/${BASE_NAME}.log; exit 100;  }
wait ${PIDS[1]}
$? || { echo BWA ${PIDS[1]} exited with error $? at `date`  >> ${DIR}/${BASE_NAME}.log; exit 100;  }

echo Mapping finished at `date` >> ${DIR}/${BASE_NAME}.log

### joining reads to create sam file
echo Joining reads to sam file  at `date` >> ${DIR}/${BASE_NAME}.log
echo $BWA sampe $REFGEN $DIR/${BASE_NAME}_trimmed_1.sai $DIR/${BASE_NAME}_trimmed_2.sai $DIR/${BASE_NAME}_trimmed_1 $DIR/${BASE_NAME}_trimmed_2 \> $DIR/${BASE_NAME}_trimmed_1and2.sam >> ${DIR}/${BASE_NAME}.log 

$BWA sampe $REFGEN $DIR/${BASE_NAME}_trimmed_1.sai $DIR/${BASE_NAME}_trimmed_2.sai $DIR/${BASE_NAME}_trimmed_1 $DIR/${BASE_NAME}_trimmed_2 > $DIR/${BASE_NAME}_trimmed_1and2.sam

echo Finished joining at `date` >>  ${DIR}/${BASE_NAME}.log

echo Checking number of mapped sequences: >>  ${DIR}/${BASE_NAME}.log
A=`grep -c '@' $DIR/${BASE_NAME}_trimmed_1`
B=`grep -c '@' $DIR/${BASE_NAME}_trimmed_2`
C=`grep -cv '@' $DIR/${BASE_NAME}_trimmed_1and2.sam`
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

exit

