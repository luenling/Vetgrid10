#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
POPO=/Volumes/Temp/Lukas/Tools/popoolation
POPO2=/Volumes/Temp/Lukas/Tools/popoolation2
PICARD=/Volumes/Temp/Lukas/Tools/picard/picard-tools-1.79
SNAPPY=/Volumes/Temp/Lukas/Tools/picard
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.18/samtools
REFGEN=/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa
BWA=/Volumes/Temp/Lukas/Tools/bwa-0.5.9/bwa

function USAGE ()
{
    echo "USAGE: $0 [-?oh] pI=xx1.bam pII=xx2.bam ..."
    echo "OPTIONS: -o  outdir"
    echo "         -p  prefix"
    echo "         -a  postfix"
    echo "         -l  logfile"
    echo "         -r  reference genome file"
    echo "         -?  this usage information"
    echo "appends to logfile, the mpileup file will be named prefix_pI_pII_postfix.mpile"
    echo "EXAMPLE: $0 -o ../blub -l merge.log -p females -a unfiltered pI=bla.bam pII=blub.bam"
    echo ""
    exit $E_OPTERROR    # Exit and explain usage, if no argument(s) given.
}


#PROCESS ARGS
while getopts ":o:l:p:a:r:?h" Option
do
    case $Option in
        o    ) OUTDIR=$OPTARG;;
        l    ) LOGFILE=$OPTARG;;
        p    ) PREF=$OPTARG;;
        a    ) POSTF=$OPTARG;;
        r    ) REFGEN=$OPTARG;;
        ?    ) USAGE; exit 0;;
        h    ) USAGE; exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               USAGE   # DEFAULT
	       exit
    esac
done

[[ -n $LOGFILE ]] || { USAGE ; exit 1; }
[[ -n $OUTDIR ]] || OUTDIR=`pwd`
[[ -d $OUTDIR ]] || { echo $OUTDIR is not a directory; exit 1; }

shift $(($OPTIND - 1))

POPS=""
INFILES=""
SHORT_INFS=""
for ARG in $*; do
    [[ $ARG =~ .+=.+ ]] || { echo infiles need to be of the form xx=filename; exit 1; }
    POP=${ARG//=*}
    INF=${ARG//*=}
    [[ $INF =~ .+\.bam$ ]] || { echo $INF not ending with \".bam\" ; exit 1; } 
    [[ -f $INF ]] || { echo file $ARG not found; exit 1; }
    INFILES=${INFILES}${INF}" "
    SHORT_INFS=${SHORT_INFS}`basename ${INF}`" "
    POPS=$POPS$POP"_"
done
MPILEUP=${PREF=females}"_"$POPS${POSTF=unfiltered}".mpileup"
SYNC=${PREF=females}"_"$POPS${POSTF=unfiltered}"_q20.sync"
echo creating pileup ${MPILEUP} out of $SHORT_INFS at `date` >> ${OUTDIR}/$LOGFILE 
echo $SAMTOOLS mpileup -6 -B -Q 0 -f $REFGEN $INFILES \>  ${MPILEUP}  >> ${OUTDIR}/${LOGFILE}
$SAMTOOLS mpileup -6 -B -Q 0 -f $REFGEN $INFILES > ${OUTDIR}/${MPILEUP}
[ $? ] || { echo Samtools pileup creation exited with error $? at `date` >>  ${OUTDIR}/${LOGFILE};  exit 1; }
echo finished pileup of $INFILES to $OUT at `date` >> ${OUTDIR}/${LOGFILE}
echo creating sync file at `date` >> ${OUTDIR}/${LOGFILE}
echo java -Xmx12g -jar $POPO2/mpileup2sync.jar --input ${OUTDIR}/${MPILEUP} --output ${OUTDIR}/$SYNC --fastq-type sanger --min-qual 20 --threads 10
java -Xmx12g -jar $POPO2/mpileup2sync.jar --input ${OUTDIR}/${MPILEUP} --output ${OUTDIR}/$SYNC --fastq-type sanger --min-qual 20 --threads 10
echo finished creation of $SYNC at `date` >> ${OUTDIR}/${LOGFILE}
[[ -d $OUTDIR/Indelfiltering ]] || mkdir ${OUTDIR}/Indelfiltering
echo creating gtf for Indel masking at `date` >> ${OUTDIR}/${LOGFILE}
INDELGTF=${OUTDIR}/Indelfiltering/`basename ${OUTDIR}/$MPILEUP .mpileup`"_INDELS"
for MC in "5 10 15"; do
    echo perl $POPO2/indel_filtering/identify-indel-regions.pl --input ${OUTDIR}/${MPILEUP} --output $INDELGTF"_mc"$MC".gtf" --indel-window 5 --min-count $MC >> ${OUTDIR}/Indelfiltering/indelfilt.log
perl $POPO2/indel_filtering/identify-indel-regions.pl --input ${OUTDIR}/${MPILEUP} --output $INDELGTF"_mc"$MC".gtf" --indel-window 5 --min-count $MC >> ${OUTDIR}/Indelfiltering/indelfilt.log
done
exit 0

