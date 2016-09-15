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
    echo "         -t  target regions file"
    echo "         -g  gatk executable (default: $GATK)"
    echo "         -r  reference genome (optional, default: $REFGEN)"
    echo "         -h  this usage information"
    echo "appends to logfile"
    echo "EXAMPLE: $0 -o bla/out.vcf -l annotator.log -v input.vcf pop1.bam pop2.bam pop3.bam "
    echo ""
    exit $E_OPTERROR    # Exit and explain usage, if no argument(s) given.
}


#PROCESS ARGS
while getopts ":o:l:hr:t:g:" Option
do
    case $Option in
        o    ) OUT=$OPTARG;;
        l    ) LOGFILE=$OPTARG;;
        g    ) GATK=$OPTARG;;
        t    ) TARGET=$OPTARG;;
        r    ) REFGEN=$OPTARG;;
        \?    ) USAGE; exit 0;;
        h    ) USAGE; exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               USAGE   # DEFAULT
	       exit
    esac
done

[[ -n $OUT && -n $TARGET ]] || { USAGE ; exit 1; }
shift $(($OPTIND - 1))
INPUTS=""
INFILES=""
for ARG in $*; do
    #[[ -f $ARG ]] || { echo file $ARG not found; exit 1; }
    INFILES=${INFILES}${ARG}" "
    INPUTS=$INPUTS" -I "$ARG
done

if [ ! $REFGEN ]; then REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa; fi

echo java -Xmx4g -jar $GATK -dt NONE\
     -nct 4 -nt 2 -L $TARGET -maxAltAlleles 3 -glm GENERALPLOIDYSNP\
    -R $REFGEN -T UnifiedGenotyper \
    $INPUTS \
    -o $OUT \
    --output_mode EMIT_VARIANTS_ONLY \
   -stand_call_conf 10.0 \
   -stand_emit_conf 10.0 \
    --sample_ploidy 25 \
    >> $LOGFILE

nohup java -Xmx4g -jar $GATK -dt NONE\
    -nct 4 -nt 2 -L $TARGET -maxAltAlleles 3 -glm GENERALPLOIDYSNP \
    -R $REFGEN -T UnifiedGenotyper \
    $INPUTS \
    -o $OUT \
    --output_mode EMIT_VARIANTS_ONLY \
   -stand_call_conf 10.0 \
   -stand_emit_conf 10.0 \
    --sample_ploidy 25 \
    2>&1 >> $LOGFILE &
