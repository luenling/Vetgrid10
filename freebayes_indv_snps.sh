shopt -s extglob
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
FREEBAYES=/Volumes/Temp/Lukas/Tools/freebayes/bin/freebayes


echo $FREEBAYES -v `basename $1 .bam`_snps_freebayes.vcf -f $REFGEN -b $1 -i -n 3 -m 20 -q 20 -e 1 -F 0.2 -C 2 --min-coverage 10 -w -B 500 >> `basename $1 .bam`.log
 $FREEBAYES -v `basename $1 .bam`_snps_freebayes.vcf -f $REFGEN -b $1 -i -n 3 -m 20 -q 20 -e 1 -F 0.2 -C 2 --min-coverage 10 -w -B 500 