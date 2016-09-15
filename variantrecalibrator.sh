#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
GATK=/Volumes/Temp/Lukas/Tools/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar
LOGFILE=variant_recalibration.log

function USAGE ()
{
    echo "USAGE: $0 [-?oh]"
    echo "OPTIONS: -o  outprefix"
    echo "         -l  logfile (optional)"
    echo "         -t  truth & training set"
    echo "         -f  training set (not truth=true)"
    echo "         -k  known set"
    echo "         -v  input vcf file"
    echo "         -g  gatk executable (default: $GATK)"
    echo "         -r  reference genome (optional, default: $REFGEN)"
    echo "         -h  this usage information"
    echo "appends to logfile"
    echo "EXAMPLE: $0 -o bla/out.vcf -l annotator.log -v input.vcf"
    echo ""
    exit $E_OPTERROR    # Exit and explain usage, if no argument(s) given.
}


#PROCESS ARGS
while getopts ":o:l:r:t:f:v:k:g:h" Option
do
    case $Option in
        o    ) OUT=$OPTARG;;
        l    ) LOGFILE=$OPTARG;;
        f    ) TRAIN=$OPTARG;;
        g    ) GATK=$OPTARG;;
        t    ) TRUTH=$OPTARG;;
        v    ) VCF=$OPTARG;;
        k    ) KNOWN=$OPTARG;;
        r    ) REFGEN=$OPTARG;;
        \?    ) USAGE; exit 0;;
        h    ) USAGE; exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               USAGE   # DEFAULT
	       exit
    esac
done

#echo  $OUT $TRUTH $VCF $KNOWN 
#exit

[[ -n $OUT && -n $TRUTH && -n $VCF && -n $KNOWN ]] || { USAGE ; exit 1; }
shift $(($OPTIND - 1))

if [ ! $REFGEN ]; then REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa; fi

COMMAND="nohup java -Xmx6g -jar $GATK -nt 3\
   -T VariantRecalibrator \
   -R $REFGEN \
   -input $VCF \
   -resource:truth,known=false,training=true,truth=true,prior=15.0 $TRUTH \
   -resource:train,known=false,training=true,truth=false,prior=10.0  $TRAIN \
   -resource:known,known=true,training=false,truth=false,prior=5.0 $KNOWN \
    -an MQRankSum -an ReadPosRankSum -an FS -an DP -an BaseQRankSum -an HaplotypeScore\
   -mode SNP \
   -recalFile $OUT.recal \
   -tranchesFile $OUT.tranches \
   -rscriptFile $OUT.plots.R" 
#  -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an DP -an BaseQRankSum \

echo $COMMAND 2\>\&1 \>\> $LOGFILE \& >> $LOGFILE


 $COMMAND 2>&1 >> $LOGFILE  
