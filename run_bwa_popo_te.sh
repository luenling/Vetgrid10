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
echo $BWA bwasw -t 6 $REFGEN ../${BASE_NAME}_trimmed_1 \> ${BASE_NAME}_1.sam  >>  ${DIR}/${BASE_NAME}.log
nohup $BWA bwasw -t 6 $REFGEN ../${BASE_NAME}_trimmed_1 > ${BASE_NAME}_1.sam  2>>  ${DIR}/${BASE_NAME}_1_bwa.log &

echo Process ID $! at `date` >> ${DIR}/${BASE_NAME}.log
PIDS=( $! )

echo $BWA bwasw -t 6 $REFGEN ../${BASE_NAME}_trimmed_2 \> ${BASE_NAME}_2.sam >>  ${DIR}/${BASE_NAME}.log

nohup $BWA bwasw -t 6 $REFGEN ../${BASE_NAME}_trimmed_2 > ${BASE_NAME}_2.sam 2>>  ${DIR}/${BASE_NAME}_2_bwa.log &

echo Process ID $! at `date` >> ${DIR}/${BASE_NAME}.log

Pids=( ${PIDS[@]} $! )

### wait for the processes to finish
wait ${PIDS[0]}
[ $? ] || { echo BWA ${PIDS[0]} exited with error $? at `date`  >> ${DIR}/${BASE_NAME}.log; exit 100;  }
wait ${PIDS[1]}
[ $? ] || { echo BWA ${PIDS[1]} exited with error $? at `date`  >> ${DIR}/${BASE_NAME}.log; exit 100;  }

echo perl ${POPOTE}/samro.pl --sam1 ${BASE_NAME}_1.sam --sam2 ${BASE_NAME}_2.sam --fq1 ../${BASE_NAME}_trimmed_1 --fq2 ../${BASE_NAME}_trimmed_2 --output ${BASE_NAME}_pe.sam >>  ${DIR}/${BASE_NAME}.log


perl ${POPOTE}/samro.pl --sam1 ${BASE_NAME}_1.sam --sam2 ${BASE_NAME}_2.sam --fq1 ../${BASE_NAME}_trimmed_1 --fq2 ../${BASE_NAME}_trimmed_2 --output ${BASE_NAME}_pe.sam 

echo java -Xmx5g -classpath $SNAPPY -jar ${PICARD}/SortSam.jar I=${BASE_NAME}_pe.sam O=${BASE_NAME}_pe_sorted.bam SO=coordinate VALIDATION_STRINGENCY=SILENT >> ${DIR}/${BASE_NAME}.log
java -Xmx5g -classpath $SNAPPY -jar ${PICARD}/SortSam.jar I=${BASE_NAME}_pe.sam O=${BASE_NAME}_pe_sorted.bam SO=coordinate VALIDATION_STRINGENCY=SILENT

echo at `date`: >>  ${DIR}/${BASE_NAME}.log
echo perl ${POPOTE}/identify-te-insertsites.pl --input ${BASE_NAME}_pe_sorted.bam --te-hierarchy-file $TEHIRE --te-hierarchy-level family --narrow-range 75 --min-count 3 --min-map-qual 15 --output ${BASE_NAME}_te_fwd_rev.txt >> ${DIR}/${BASE_NAME}.log

perl ${POPOTE}/identify-te-insertsites.pl --input ${BASE_NAME}_pe_sorted.bam --te-hierarchy-file $TEHIRE --te-hierarchy-level family --narrow-range 75 --min-count 3 --min-map-qual 15 --output ${BASE_NAME}_te_fwd_rev.txt >>  ${DIR}/${BASE_NAME}_te.log

echo at `date`: >>  ${DIR}/${BASE_NAME}.log
echo perl ${POPOTE}/genomic-N-2gtf.pl --input $REFGEN > poly_n.gtf >> ${DIR}/${BASE_NAME}.log

perl ${POPOTE}/genomic-N-2gtf.pl --input $REFGEN > poly_n.gtf

echo at `date`: >>  ${DIR}/${BASE_NAME}.log
echo perl ${POPOTE}/crosslink-te-sites.pl --directional-insertions ${BASE_NAME}_te_fwd_rev.txt --min-dist 74 --max-dist 250 --output ${BASE_NAME}_te_inserts.txt --single-site-shift 100 --poly-n poly_n.gtf --te-hierarchy $TEHIRE --te-hier-level order >>  ${DIR}/${BASE_NAME}.log

perl ${POPOTE}/crosslink-te-sites.pl --directional-insertions ${BASE_NAME}_te_fwd_rev.txt --min-dist 74 --max-dist 250 --output ${BASE_NAME}_te_inserts.txt --single-site-shift 100 --poly-n poly_n.gtf --te-hierarchy $TEHIRE --te-hier-level order >> ${BASE_NAME}_te.log

# echo at `date`: >>  ${DIR}/${BASE_NAME}.log
# echo perl ${POPOTE}/update-teinserts-with-knowntes.pl --known known-te-insertions.txt --output te_insertions.updated.txt --te-hierarchy-file $TEHIRE --te-hierarchy-level family --max-dist 300 --te-insertions  ${BASE_NAME}_te_inserts.txt --single-site-shift 100 >>  ${DIR}/${BASE_NAME}.log

# perl ${POPOTE}/update-teinserts-with-knowntes.pl --known known-te-insertions.txt --output te_insertions.updated.txt --te-hierarchy-file $TEHIRE --te-hierarchy-level family --max-dist 300 --te-insertions  ${BASE_NAME}_te_inserts.txt --single-site-shift 100 >> ${BASE_NAME}_te.log


echo at `date`: >>  ${DIR}/${BASE_NAME}.log
echo perl ${POPOTE}/estimate-polymorphism.pl --sam-file  ${BASE_NAME}_pe_sorted.bam  --te-insert-file  ${BASE_NAME}_te_inserts.txt --te-hierarchy-file $TEHIRE --te-hierarchy-level family --min-map-qual 15 --output  ${BASE_NAME}_te_polymorphism.txt >>  ${DIR}/${BASE_NAME}.log
perl ${POPOTE}/estimate-polymorphism.pl --sam-file  ${BASE_NAME}_pe_sorted.bam  --te-insert-file  ${BASE_NAME}_te_inserts.txt --te-hierarchy-file $TEHIRE --te-hierarchy-level family --min-map-qual 15 --output  ${BASE_NAME}_te_polymorphism.txt >> ${BASE_NAME}_te.log

echo at `date`: >>  ${DIR}/${BASE_NAME}.log
echo perl ${POPOTE}/filter-teinserts.pl --te-insertions  ${BASE_NAME}_te_polymorphism.txt --output  ${BASE_NAME}_te_poly_filtered.txt --discard-overlapping --min-count 10 >>  ${DIR}/${BASE_NAME}.log
perl ${POPOTE}/filter-teinserts.pl --te-insertions  ${BASE_NAME}_te_polymorphism.txt --output  ${BASE_NAME}_te_poly_filtered.txt --discard-overlapping --min-count 10 >> ${BASE_NAME}_te.log

grep -P '^(?:2L|2R|3L|3R|X|4)\s+' ${BASE_NAME}_te_poly_filtered.txt | sort -k1,1 -k2,2n > ${BASE_NAME}_te_poly_filtered_maj_chrom.txt

python -u /Volumes/Temp/Lukas/Tools/Scripts/create_bed_from_TE.py -i ${BASE_NAME}_te_poly_filtered_maj_chrom.txt > ${BASE_NAME}_te_poly_filtered_maj_chrom.bed