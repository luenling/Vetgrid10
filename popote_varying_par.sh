#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
POPOTE=/Volumes/Temp/Lukas/Tools/popoolationte-read-only/
PICARD=/Volumes/Temp/Lukas/Tools/picard/picard-tools-1.79
SNAPPY=/Volumes/Temp/Lukas/Tools/picard
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.19/samtools
REFGEN=/Volumes/Temp/Lukas/reference_genome/mel/masked_TEs/ref_genome_combined_TEs.TE_masked.fa
TEHIRE=/Volumes/Temp/Lukas/reference_genome/mel/masked_TEs/te_hirarchy_rk_paper.txt
BWA=/Volumes/Temp/Lukas/Tools/bwa-0.5.9/bwa

DIR=`pwd`
BASE_NAME=$1
# echo at `date`: >>  ${DIR}/${BASE_NAME}.log
# echo perl ${POPOTE}/identify-te-insertsites.pl --input ${BASE_NAME}_pe_sorted.bam --te-hierarchy-file $TEHIRE --te-hierarchy-level family --narrow-range 75 --min-count 3 --min-map-qual 15 --output ${BASE_NAME}_te_fwd_rev.txt >> ${DIR}/${BASE_NAME}.log

# perl ${POPOTE}/identify-te-insertsites.pl --input ${BASE_NAME}_pe_sorted.bam --te-hierarchy-file $TEHIRE --te-hierarchy-level family --narrow-range 75 --min-count 3 --min-map-qual 15 --output ${BASE_NAME}_te_fwd_rev.txt >>  ${DIR}/${BASE_NAME}_te.log

# echo at `date`: >>  ${DIR}/${BASE_NAME}.log
# echo perl ${POPOTE}/genomic-N-2gtf.pl --input $REFGEN > poly_n.gtf >> ${DIR}/${BASE_NAME}.log

# perl ${POPOTE}/genomic-N-2gtf.pl --input $REFGEN > poly_n.gtf

echo at `date`: >>  ${DIR}/${BASE_NAME}.log
echo perl ${POPOTE}/crosslink-te-sites.pl --directional-insertions ${BASE_NAME}_te_fwd_rev.txt --min-dist $2 --max-dist $3 --output ${BASE_NAME}_te_inserts.txt --single-site-shift 100 --poly-n poly_n.gtf --te-hierarchy $TEHIRE --te-hier-level order >>  ${DIR}/${BASE_NAME}.log

perl ${POPOTE}/crosslink-te-sites.pl --directional-insertions ${BASE_NAME}_te_fwd_rev.txt  --min-dist $2 --max-dist $3  --output ${BASE_NAME}_te_inserts_min_$2_max_$3.txt --single-site-shift 100 --poly-n poly_n.gtf --te-hierarchy $TEHIRE --te-hier-level order >> ${BASE_NAME}_te_min_$2_max_$3.log

# echo at `date`: >>  ${DIR}/${BASE_NAME}.log
# echo perl ${POPOTE}/update-teinserts-with-knowntes.pl --known known-te-insertions.txt --output te_insertions.updated.txt --te-hierarchy-file $TEHIRE --te-hierarchy-level family --max-dist 300 --te-insertions  ${BASE_NAME}_te_inserts.txt --single-site-shift 100 >>  ${DIR}/${BASE_NAME}.log

# perl ${POPOTE}/update-teinserts-with-knowntes.pl --known known-te-insertions.txt --output te_insertions.updated.txt --te-hierarchy-file $TEHIRE --te-hierarchy-level family --max-dist 300 --te-insertions  ${BASE_NAME}_te_inserts.txt --single-site-shift 100 >> ${BASE_NAME}_te.log


echo at `date`: >>  ${DIR}/${BASE_NAME}.log
echo perl ${POPOTE}/estimate-polymorphism.pl --sam-file  ${BASE_NAME}_pe_sorted.bam  --te-insert-file  ${BASE_NAME}_te_inserts_min_$2_max_$3.txt --te-hierarchy-file $TEHIRE --te-hierarchy-level family --min-map-qual 15 --output  ${BASE_NAME}_te_polymorphism_min_$2_max_$3.txt >>  ${DIR}/${BASE_NAME}.log
perl ${POPOTE}/estimate-polymorphism.pl --sam-file  ${BASE_NAME}_pe_sorted.bam  --te-insert-file  ${BASE_NAME}_te_inserts_min_$2_max_$3.txt --te-hierarchy-file $TEHIRE --te-hierarchy-level family --min-map-qual 15 --output  ${BASE_NAME}_te_polymorphism_min_$2_max_$3.txt >> ${BASE_NAME}_te_min_$2_max_$3.log

echo at `date`: >>  ${DIR}/${BASE_NAME}.log
echo perl ${POPOTE}/filter-teinserts.pl --te-insertions  ${BASE_NAME}_te_polymorphism_min_$2_max_$3.txt --output  ${BASE_NAME}_te_poly_filtered_min_$2_max_$3.txt --discard-overlapping --min-count 5 >>  ${DIR}/${BASE_NAME}.log
 perl ${POPOTE}/filter-teinserts.pl --te-insertions  ${BASE_NAME}_te_polymorphism_min_$2_max_$3.txt --output  ${BASE_NAME}_te_poly_filtered_min_$2_max_$3.txt --discard-overlapping --min-count 5 >> ${BASE_NAME}_te_min_$2_max_$3.log


grep -P '^(?:2L|2R|3L|3R|X|4)\s+' ${BASE_NAME}_te_poly_filtered_min_$2_max_$3.txt | sort -k1,1 -k2,2n > ${BASE_NAME}_te_poly_filtered_min_$2_max_$3_maj_chrom.txt

python -u /Volumes/Temp/Lukas/Tools/Scripts/create_bed_from_TE.py -i ${BASE_NAME}_te_poly_filtered_min_$2_max_$3_maj_chrom.txt > ${BASE_NAME}_te_poly_filtered_min_$2_max_$3_maj_chrom.bed