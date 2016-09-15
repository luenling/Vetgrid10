#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
SAMTOOLS=/Volumes/Temp/Lukas/Tools/samtools-0.1.18/samtools
SCRIPTS=/Volumes/Temp/Lukas/Tools/Scripts
REFMEL=/Volumes/Temp/Lukas/ref_genomes_contam/mel/dmel_all.fa
REFSIM=/Volumes/Temp/Lukas/ref_genomes_contam/simulans/dsim_all.fa
LOGFILE=sim_filtering.log
function USAGE ()
{
    echo "USAGE:  $0 [-?iomqxt]"
    echo "OPTIONS: -i  infile base name"
    echo "         -o output directory"
    echo "         -m mean fragment length"
    echo "         -x maximal fragment length"
    echo "         -t threads"
    echo "         -q quality format"
    echo "         -? help"
    echo "appends to logfile "
    echo ""
    exit $E_OPTERROR    # Exit and explain usage, if no argument(s) given.
}

OUTDIR='Filtering'

#PROCESS ARGS
while getopts ":i:o:m:x:t:q:?" Option
do
    case $Option in
        i    ) INFILE=$OPTARG;;
        o    ) OUTDIR=$OPTARG;;
        q    ) QUALITY=$OPTARG;;
        m    ) MEANFRAG=$OPTARG;;
        x    ) MAXFRAG=$OPTARG;;
        t    ) THREADS=$OPTARG;;
        ?    ) USAGE; exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               USAGE   # DEFAULT
    esac
done

if [ ! -d $OUTDIR ]; then mkdir $OUTDIR; fi

echo python $SCRIPTS/find-simulans-reads-v1.2_alt.py --mel-reference-fasta-file $REFMEL --sim-reference-fasta-file $REFSIM --read1 ${INFILE}1.fastq --read2 ${INFILE}2.fastq --fastq-type $QUALITY --output-dir $OUTDIR --thread $THREADS --mean-fragment-size $MEANFRAG --maximum-fragment-size $MAXFRAG >> ${LOGFILE}
nohup python $SCRIPTS/find-simulans-reads-v1.2_alt.py --mel-reference-fasta-file $REFMEL --sim-reference-fasta-file $REFSIM --read1 ${INFILE}1.fastq --read2 ${INFILE}2.fastq --fastq-type $QUALITY --output-dir $OUTDIR --thread $THREADS --mean-fragment-size $MEANFRAG --maximum-fragment-size $MAXFRAG &
wait $!