#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
GATK=/Volumes/Temp/Lukas/Tools/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar
LOGFILE=annotator.log

function USAGE ()
{
    echo "USAGE: $0 [-?oh] xx1.bam xx2.bam ..."
    echo "OPTIONS: -o  outfile"
    echo "         -l  logfile (optional)"
    echo "         -g  gatk executable (default: $GATK)"
    echo "         -r  reference genome (optional, default: $REFGEN)"
    echo "         -h  this usage information"
    echo "appends to logfile"
    echo "EXAMPLE: $0 -o bla/out.vcf -l annotator.log -v input.vcf pop1.bam pop2.bam pop3.bam "
    echo ""
    exit $E_OPTERROR    # Exit and explain usage, if no argument(s) given.
}


#PROCESS ARGS
while getopts ":o:l:hr:g:" Option
do
    case $Option in
        o    ) OUT=$OPTARG;;
        l    ) LOGFILE=$OPTARG;;
        g    ) GATK=$OPTARG;;
        r    ) REFGEN=$OPTARG;;
        \?    ) USAGE; exit 0;;
        h    ) USAGE; exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               USAGE   # DEFAULT
	       exit
    esac
done

CHROMS="2L 2R 3L 3R 4 X YHet"
[[ -n $OUT ]] || { USAGE ; exit 1; }
shift $(($OPTIND - 1))
INPUTS=""
INFILES=""
for ARG in $*; do
    #[[ -f $ARG ]] || { echo file $ARG not found; exit 1; }
    INFILES=${INFILES}${ARG}" "
    INPUTS=$INPUTS" -I "$ARG
done

if [ ! $REFGEN ]; then REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa; fi

for i in $CHROMS; do 
COMMAND="nohup java -Xmx3g -jar $GATK -dt NONE -L $i \
     -nct 2 -maxAltAlleles 3 -glm GENERALPLOIDYINDEL -deletions 1.1 \
    -mbq 20 -minIndelFrac 0.05 -minIndelCnt 20\
    -R $REFGEN -T UnifiedGenotyper \
    $INPUTS \
    -o ${OUT}_${i}.vcf \
    --output_mode EMIT_VARIANTS_ONLY \
   -stand_call_conf 10.0 \
   -stand_emit_conf 10.0 \
    --sample_ploidy 20" 
#  -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an DP -an BaseQRankSum \
echo $COMMAND >> $LOGFILE
nohup $COMMAND 2>&1 >> `basename $LOGFILE .log`_${i}.log &
done

# echo $COMMAND 2\>\&1 \>\> $LOGFILE \& >> $LOGFILE

# $COMMAND 2>&1 >> $LOGFILE &
