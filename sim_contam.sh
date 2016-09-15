#!/bin/bash
#----------
# author: Lukas Endler
# authored: 2016-04-21 13:56:58
# Time-stamp: <2016-05-02 12:01:03 lukasendler>
# takes a sync file and runs the simulans divergence checker
#--------------


COMP_DIV=/Volumes/Temp/Lukas/Tools/compare_divergence
BASE_NAME=`basename $1 .sync`
echo Comparing divergence at `date`: >>  ${BASE_NAME}.log
echo python ${COMP_DIV}/compare_divergence.py --input  $1 --div ${COMP_DIV}/europe_sim_cov100.div --pops 1,2,3,4,5,6,7,8 --mincount 1,1,1,1,1,1,1,1 --out ${FEMFILE}_q20.europe_sim_cov100.divergence.txt >> ${DIR}/${BASE_NAME}.log
python ${COMP_DIV}/compare_divergence.py --input  ${FEMFILE}_q20.sync --div ${COMP_DIV}/europe_sim_cov100.div --pops 1,2,3,4,5,6,7,8 --mincount 1,1,1,1,1,1,1,1 --out ${FEMFILE}_q20.europe_sim_cov100.divergence.txt
