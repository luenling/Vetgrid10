#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
SNPEFF=/Volumes/Temp/Lukas/Tools/snpEff/snpEff.jar
SNPEFFCONF=/Volumes/Temp/Lukas/Tools/snpEff/snpEff.config
BASE_NAME=`basename $1 .vcf.gz`
gunzip -c $1 | java -Xmx4g -jar $SNPEFF -stats summary_snpeff_BDGP5.78 -v BDGP5.78 -c $SNPEFFCONF  > ${BASE_NAME}.BDPG5.78.snpeff.vcf 2> snpeff_bdpg5.78.stderr.log &
gunzip -c $1 | java -Xmx4g -jar $SNPEFF -stats summary_snpeff_dm5.48 -v dm5.48 -c $SNPEFFCONF  > ${BASE_NAME}.dm5.48.snpeff.vcf 2> snpeff_dm5.48.stderr.log &
