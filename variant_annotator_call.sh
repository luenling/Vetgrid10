#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
GATK=/Volumes/Temp/Lukas/Tools/GATK_2.3.9/GenomeAnalysisTK.jar
LOGFILE=annotator.log

function USAGE ()
{
    echo "USAGE: $0 [-?oh] xx1.bam xx2.bam ..."
    echo "OPTIONS: -o  outfile"
    echo "         -l  logfile (optional)"
    echo "         -v  vcf file"
    echo "         -g  gatk executable (default: $GATK)"
    echo "         -r  reference genome (optional, default: $REFGEN)"
    echo "         -h  this usage information"
    echo "appends to logfile"
    echo "EXAMPLE: $0 -o bla/out.vcf -l annotator.log -v input.vcf pop1.bam pop2.bam pop3.bam "
    echo ""
    exit $E_OPTERROR    # Exit and explain usage, if no argument(s) given.
}


#PROCESS ARGS
while getopts ":o:l:hr:v:g:" Option
do
    case $Option in
        o    ) OUT=$OPTARG;;
        l    ) LOGFILE=$OPTARG;;
        g    ) GATK=$OPTARG;;
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

[[ -n $OUT && -n $VCF ]] || { USAGE ; exit 1; }
shift $(($OPTIND - 1))
INPUTS=""
INFILES=""
for ARG in $*; do
    #[[ -f $ARG ]] || { echo file $ARG not found; exit 1; }
    INFILES=${INFILES}${ARG}" "
    INPUTS=$INPUTS" -I "$ARG
done

if [ ! $REFGEN ]; then REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa; fi

echo java -Xmx2g -jar $GATK \
   -R $REFGEN \
   -T VariantAnnotator \
   --group StandardAnnotation \
   -A ClippingRankSumTest -A BaseCounts -A GCContent\
   $INPUTS  \
   -o $OUT \
   --variant $VCF >> $LOGFILE

nohup java -Xmx2g -jar $GATK \
   -R $REFGEN \
   -T VariantAnnotator \
   --group StandardAnnotation \
   -A ClippingRankSumTest -A BaseCounts -A GCContent\
   $INPUTS  \
   -o $OUT \
   --variant $VCF >> $LOGFILE &