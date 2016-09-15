#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
GATK=/Volumes/Temp/Lukas/Tools/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar
LOGFILE=variant_recalibration.log
FLEV=99
function USAGE ()
{
    echo "USAGE: $0 [-?oh]"
    echo "OPTIONS: -o  outprefix"
    echo "         -l  logfile (optional)"
    echo "         -t  tranches file"
    echo "         -f  Filterlevel (def: 99)"
    echo "         -c  Recalibrationfile"
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
while getopts ":o:l:r:t:f:v:c:g:h" Option
do
    case $Option in
        o    ) OUT=$OPTARG;;
        l    ) LOGFILE=$OPTARG;;
        g    ) GATK=$OPTARG;;
        t    ) TRFILE=$OPTARG;;
        v    ) VCF=$OPTARG;;
        r    ) REFGEN=$OPTARG;;
        c    ) RECAL=$OPTARG;;
        f    ) FLEV=$OPTARG;;
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

[[ -n $OUT && -n $VCF && -n $RECAL && -n $TRFILE ]] || { USAGE ; exit 1; }
shift $(($OPTIND - 1))

if [ ! $REFGEN ]; then REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa; fi

COMMAND="java -Xmx3g -jar $GATK \
   -T ApplyRecalibration \
   -R $REFGEN \
   -input $VCF \
   --ts_filter_level $FLEV \
   -tranchesFile $TRFILE \
   -recalFile $RECAL \
   -mode SNP \
   -o $OUT" 
#  -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an DP -an BaseQRankSum \

echo $COMMAND 2\>\&1 \>\> $LOGFILE \& >> $LOGFILE


 $COMMAND 2>&1 >> $LOGFILE  
