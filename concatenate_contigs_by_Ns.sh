#!/bin/bash
#----------
# author: Lukas Endler
# date: 08.10.2015 at 12:23
# takes a fasta file and tries to create a new one with all small contigs concatenated by Ns
# only works for genomes were the small ones are in a sequence (sorted by size) for a specific Dsim genome ;)
#--------------
FASTA=$1
FN=`basename $FASTA .fa`
FN=`basename $FN .fasta`
FN=`basename $FN .fna`

# create fasta index
if [[ ! -e ${FASTA}.fai ]]; then
    samtools faidx $1 
fi
# get small contigs
awk '$2 < 260000' ${FASTA}.fai | cut -f 1 | sed 's/BGI_25a_trimmed_with-5p-trim_1_PE-reads//g' | sed 's/\(.\)$/\1\$/g'  > small_contigs
# get lines of small contig headers and write sed commands to replace them with Ns
grep -n ">" $FASTA | grep -E -f small_contigs | cut -f 1 -d ":" | awk '{print $1"{r replacement"; print "d"; print"}"}' > sed.replacement
# replace them
sed -f sed.replacement < $FASTA > ${FN}_joined_by_N2.fasta
# normalize length
picard NormalizeFasta I=${FN}_joined_by_N2.fasta O=${FN}_joined_by_N.fasta_joined_by_N2_norm.fasta
