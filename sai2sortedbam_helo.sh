REFGEN=/Volumes/Temp/Helo/pigmentation/DATA/reference_genome/dmel_Phix_ref_genome.fasta
BASE_NAME=$1
bwa sampe $REFGEN ${BASE_NAME}_trimmed_1.sai ${BASE_NAME}_trimmed_2.sai ${BASE_NAME}_trimmed_1 ${BASE_NAME}_trimmed_2 >  ${BASE_NAME}.sam
samtools sort ${BASE_NAME}.sam | perl -F'\t' -lane 'if (exists($F[10])) {@q = split(//,$F[10]); $F[10] = join("",map {chr(ord($_)-31)} @q)}; print join("\t",@F);' | samtools -bhS - > ${BASE_NAME}.bam
rm ${BASE_NAME}.sam
