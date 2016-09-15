REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
GATK=/Volumes/Temp/Lukas/Tools/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar
LOGFILE=annotator.log
BASE_NAME=$1
anno_files=`for i in ${BASE_NAME}_*.vcf; do echo -n  --variant $i" "; done`
echo java -Xmx4g -jar $GATK -R $REFGEN -T CombineVariants $anno_files --assumeIdenticalSamples --out ${BASE_NAME}_all.vcf >> $LOGFILE
java -Xmx4g -jar $GATK -R $REFGEN -T CombineVariants $anno_files --assumeIdenticalSamples --out ${BASE_NAME}_all.vcf 2>&1 >> $LOGFILE 

bash /Volumes/Temp/Lukas/Tools/Scripts/variantrecalibrator_indel.sh -o ${BASE_NAME}_recal -v ${BASE_NAME}_all.vcf -t /Volumes/Temp/Lukas/Data/DGRP_vcf/freeze2.bins.indels_noGT_v41.vcf.gz  2>&1 >> nohup.out

java -Xmx4g -jar $GATK -R $REFGEN -T ApplyRecalibration -input ${BASE_NAME}_all.vcf --ts_filter_level 95.0 -tranchesFile ${BASE_NAME}_recal.tranches -recalFile ${BASE_NAME}_recal.recal -mode INDEL -o ${BASE_NAME}_all_filter_95.vcf >> nohup.out 

bcftools view -M -f PASS -e '%MAX(MLPSAF)<0.05' -O z  ${BASE_NAME}_all_filter_95.vcf  > ${BASE_NAME}_all_filter_95_maxAF_gt_0.05.vcf.gz

bcftools view -M -f PASS -i '%MIN(DP)>=15' -O z ${BASE_NAME}_all_filter_95_maxAF_gt_0.05.vcf.gz > ${BASE_NAME}_all_filter_95_maxAF_gt_0.05_minDP_gt_15.vcf.gz

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' ${BASE_NAME}_all_filter_95_maxAF_gt_0.05_minDP_gt_15.vcf.gz | bgzip -c > ${BASE_NAME}_all_filter_95_maxAF_gt_0.05_minDP_gt_15.sync.gz

