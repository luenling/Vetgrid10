SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.19/samtools
BCFTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.19/bcftools/bcftools
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
SCRIPTS=/Volumes/Temp/Lukas/Tools/Scripts
CMH=$1
BAMS=$2
BASE_NAME=`basename $1 .cmhout`
echo "creating vcf file out of $1:"
echo python $SCRIPTS/sync2vcf_alt.py -s $CMH \| bgzip -c  \> ${BASE_NAME}.vcf.gz
python $SCRIPTS/sync2vcf_alt.py -s $CMH | bgzip -c  > ${BASE_NAME}.vcf.gz
echo tabix -p vcf ${BASE_NAME}.vcf.gz
tabix -p vcf ${BASE_NAME}.vcf.gz
echo $SAMTOOLS mpileup -g -S -I -R -d 1000 -B -Q 20 -f $REFGEN  $BAMS \| $BCFTOOLS view -l ${BASE_NAME}.vcf.gz - \| python $SCRIPTS/calc_biases_from_vcf.py \> ${BASE_NAME}_metrics.tab
$SAMTOOLS mpileup -g -S -I -d 1000 -B -Q 20 -f $REFGEN  $BAMS | $BCFTOOLS view -l ${BASE_NAME}.vcf.gz - | python $SCRIPTS/calc_biases_from_vcf.py > ${BASE_NAME}_metrics.tab
echo "awk \' (\$4 \>=30 \&\& \$5 <= 0.01 ) \|\| (\$8 \<= 8 \&\& $10 \>= 3.5) \'" ${BASE_NAME}_metrics.tab \> ${BASE_NAME}_bad_snps.tab
awk ' ($4 >=30 && $5 <= 0.01 ) || ($8 <= 8 && $10 >= 3.5) '  ${BASE_NAME}_metrics.tab > ${BASE_NAME}_bad_snps.tab
echo python $SCRIPTS/get_positions_from_sync.py -d -a ${BASE_NAME}_bad_snps.tab -b $CMH \> ${BASE_NAME}_filtered_SNPs.cmhout
python $SCRIPTS/get_positions_from_sync.py -d -a ${BASE_NAME}_bad_snps.tab -b $CMH > ${BASE_NAME}_filtered_SNPs.cmhout