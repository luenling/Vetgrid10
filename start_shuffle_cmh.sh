#!/bin/bash
# turn on extended globbing behavior
shopt -s extglob

SHUFFLES=( '2:3,3:2,5:4'  '2:1,3:2,5:4'  '2:3,5:6,6:5'  '1:2,2:1,3:2'  '2:1,3:2,6:5'  '1:3,4:5,4:6'  '2:1,4:6,5:6'  '2:1,4:5,5:6'  '2:3,4:6,5:6'  '1:3,4:6,5:6'  '2:3,3:1,6:5' )
LOGFILE=shuffle.log

for i in ${SHUFFLES[@]}; do

echo at `date`: >> $LOGFILE
echo nohup python /Volumes/Temp/Lukas/Tools/FDR/null_dist_shuffled_cmh.py -i $1 -c -t "${i}" \& >> $LOGFILE
nohup python /Volumes/Temp/Lukas/Tools/FDR/null_dist_shuffled_cmh.py -i $1 -c -t "${i}" &

done


