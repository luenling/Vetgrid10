#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob
function USAGE ()
{
    echo "USAGE:  ineff.sh [-?dh]"
    echo "OPTIONS: -d  date of logfiles (format: yyyymmdd)"
    echo "         -l  logfile"
    echo "         -?  this usage information"
    echo "appends to logfile "
    echo "EXAMPLE: ineff.sh -d 20060801 -h 2[12]"
    echo ""
    exit $E_OPTERROR    # Exit and explain usage, if no argument(s) given.
}

#PROCESS ARGS
while getopts ":d:h:?" Option
do
    case $Option in
        d    ) DT=$OPTARG;;
        l    ) HR=$OPTARG;;
        ?    ) USAGE; exit 0;;
        *    ) echo ""
               echo "Unimplemented option chosen."
               USAGE   # DEFAULT
    esac
done
echo $DT 
shift $(($OPTIND - 1))
echo $1
#  Decrements the argument pointer so it points to next argument.
#  $1 now references the first non option item supplied on the command line
#+ if one exists.
