#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
PICARD=/Volumes/Temp/Lukas/Tools/picard/picard-tools-1.104
SNAPPY=/Volumes/Temp/Lukas/Tools/picard
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.19/samtools
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa

function USAGE ()
{
    echo "USAGE: $0 [-?oh] xx1.bam xx2.bam ..."
    echo "OPTIONS: -o  outfile"
    echo "         -l  logfile"
    echo "         -?  this usage information"
    echo "appends to logfile"
    echo "EXAMPLE: $0 -o bla/merged.bam -l merge.log "
    echo ""
    exit $E_OPTERROR    # Exit and explain usage, if no argument(s) given.
}


#PROCESS ARGS
while getopts ":o:l:?h" Option
do
    case $Option in
        o    ) OUT=$OPTARG;;
        l    ) LOGFILE=$OPTARG;;
        ?    ) USAGE; exit 0;;
        h    ) USAGE; exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               USAGE   # DEFAULT
	       exit
    esac
done

[[ -n $OUT && -n $LOGFILE ]] || { USAGE ; exit 1; }
shift $(($OPTIND - 1))
INPUTS=""
INFILES=""
for ARG in $*; do
    [[ -f $ARG ]] || { echo file $ARG not found; exit 1; }
    INFILES=${INFILES}${ARG}" "
    INPUTS=$INPUTS" I="$ARG
done

echo joining bam files at `date` >> ${LOGFILE}
echo nohup java -Xmx5g -classpath $SNAPPY -jar ${PICARD}/MergeSamFiles.jar $INPUTS O=${OUT} USE_THREADING=TRUE VALIDATION_STRINGENCY=SILENT >> ${LOGFILE}
java -Xmx5g -classpath $SNAPPY -jar ${PICARD}/MergeSamFiles.jar $INPUTS O=${OUT} USE_THREADING=TRUE VALIDATION_STRINGENCY=SILENT
$SAMTOOLS index ${OUT}
echo finished joining $INFILES to $OUT at `date` >> ${LOGFILE}
