#!/bin/sh
# turn on extended globbing behavior
shopt -s extglob
DIR=/Volumes/Temp/Lukas/reference_genome
DIR_MEL=/Volumes/Temp/Lukas/reference_genome/mel
GENOMES=(dmel5.38-clean Acetobacter Wolbachia NC_001422)
# shorten IDs
for i in $GENOMES
do
   # awk '{print $1}'  ${DIR}/${i}.?(fa fasta fna) > ${DIR}/${i}_shortID.fa
    echo "${DIR}/${i}.?(fa fasta fna) ${DIR}/${i}_shortID.fa"
done