#!/bin/bash
# usage: join_sync_files.sh infa.sync.gz infb.sync.gz > outfile.sync
# creates temp files with CHRBPS as key column, and joins on all common fields 
# can also keep all lines exclusively in the second file with the missing ones from the first replaced by 0:0:0:0:0:0
# and padding empty lines from  file 1 with 0s
# to do that add: -a $F -e "0:0:0:0:0:0" to the join line
# 
INFA=$1
INFB=$2
# use first or second file to join onto
F=1
gunzip -c $INFA | awk 'BEGIN {OFS = "\t"} {key=sprintf("%-25s%09d",tolower($1),$2); gsub(/ /,"-",key); print key,$0}'  > ${INFA}".temp"
gunzip -c $INFB | awk 'BEGIN {OFS = "\t"} {key=sprintf("%-25s%09d",tolower($1),$2); gsub(/ /,"-",key); print key,$0}'  > ${INFB}".temp"
wA=`head -1 ${INFA}.temp | wc -w`
wB=`head -1 ${INFB}.temp | wc -w`
#CMD='$a='$wA'; $b='$wB'; $c=""; for $i (2..$a) {$c .=" 1.".$i}; for $i (5..$b){$c.=" 2.".$i}; print $c;'
#CMD='$a='$wA'; $b='$wB'; $c=""; for $i (2..4){$c .=" '$F'.".$i} ; for $i (5..$a) {$c .=" 1.".$i}; for $i (5..$b){$c.=" 2.".$i}; print $c;'
#CMD='$a='$wA'; $b='$wB'; $c=""; for $i (2..3){$c .=" '$F'.".$i} ; for $i (4..$a) {$c .=" 1.".$i}; for $i (4..$b){$c.=" 2.".$i}; print $c;'
# create placeholder string for lines missing in file A
CMD='$a='$wA'; $b='$wB'; $c=""; for $i (2..3){$c .=" '$F'.".$i} ; for $i (4..$a) {$c .=" 1.".$i}; for $i (5..$b){$c.=" 2.".$i}; print $c;'
#echo $CMD
OUTFORM=$(perl -e "$CMD") 
#echo $OUTFORM
#join -t $'\t' -j 1 -a $F -e "-" -o "$OUTFORM"  ${INFA}".temp"  ${INFB}".temp"
# to not pad and only print common lines comment the last and uncomment the next line
join -t $'\t' -j 1 -o "$OUTFORM"  ${INFA}".temp"  ${INFB}".temp"
rm -f ${INFA}".temp"  ${INFB}".temp"

 
