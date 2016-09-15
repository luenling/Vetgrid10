#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
DIR=`pwd`
perl ${POPO}/basic-pipeline/trim-fastq.pl --input1 ${DIR}/9390__81KGWABXX_2_s_2_1_sequence.fq --input2 ${DIR}/9390__81KGWABXX_2_s_2_2_sequence.fq --output 9390__81KGWABXX_2_s_2_trimmed --quality-threshold 18 --fastq-type illumina --min-length 50

exit