#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
GATK=/Volumes/Temp/Lukas/Tools/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar
LOGFILE=variant_recalibration.log

function USAGE ()
{
    echo "USAGE: $0 [-?oh] xx1.bam xx2.bam ..."
    echo "OPTIONS: -o  outprefix"
    echo "         -l  logfile (optional)"
    echo "         -t  truth & training set"
    echo "         -v  input vcf file"
    echo "         -g  gatk executable (default: $GATK)"
    echo "         -r  reference genome (optional, default: $REFGEN)"
    echo "         -h  this usage information"
    echo "appends to logfile"
    echo "EXAMPLE: $0 -o bla/out.vcf -l annotator.log -v input.vcf pop1.bam pop2.bam pop3.bam "
    echo ""
    exit $E_OPTERROR    # Exit and explain usage, if no argument(s) given.
}


#PROCESS ARGS
while getopts ":o:l:r:t:v:g:h" Option
do
    case $Option in
        o    ) OUT=$OPTARG;;
        l    ) LOGFILE=$OPTARG;;
        g    ) GATK=$OPTARG;;
        t    ) TRUTH=$OPTARG;;
        v    ) VCF=$OPTARG;;
        r    ) REFGEN=$OPTARG;;
        \?    ) USAGE; exit 0;;
        h    ) USAGE; exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               USAGE   # DEFAULT
	       exit
    esac
done

[[ -n $OUT && -n $TRUTH && -n $VCF ]] || { USAGE ; exit 1; }
shift $(($OPTIND - 1))

if [ ! $REFGEN ]; then REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa; fi

COMMAND="nohup java -Xmx6g -jar $GATK -nt 2\
   -T VariantRecalibrator \
   -R $REFGEN \
   -input $VCF \
   -resource:truth,known=false,training=true,truth=true,prior=12.0 $TRUTH \
  -an ReadPosRankSum -an FS -an MQRankSum -an DP\
   -mode INDEL \
   -recalFile $OUT.recal \
   -tranchesFile $OUT.tranches -tranche 97.5 -tranche 95.0 -tranche 92.5 -tranche 85.0 \
    -tranche 99.0 -tranche 99.9 -tranche 90.0 \
   -rscriptFile $OUT.plots.R" 
#  -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an DP -an BaseQRankSum \
#   --minNumBadVariants 1000 --percentBadVariants 0.005 --maxGaussians 4 \

echo $COMMAND 2\>\&1 \>\> $LOGFILE \& >> $LOGFILE


 $COMMAND 2>&1 >> $LOGFILE  
