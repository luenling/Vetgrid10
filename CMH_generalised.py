#!/usr/bin/python
import os
import re
import argparse
import numpy as np
from collections import defaultdict
import Sync_parser
from numpy import array, zeros
import rpy2.robjects as robjects
from rpy2.robjects import packages
import argparse


def cmh_test(three_way,tab_dims,ord_names_r,rep_names_r,alleles_r=robjects.Vector(["A1","A2"]),types_r=robjects.Vector(["ALL"]),score=False):
    #three_way_r=robjects.IntVector(three_way.flatten())
    three_way_r=robjects.Vector(three_way.flatten())
    three_way_array=robjects.r['array'](three_way_r,dim=tab_dims,dimnames=robjects.r['list'](Alleles=alleles_r, Classes=ord_names_r,Reps=rep_names))
    three_way_tab=robjects.r['as.table'](three_way_array)
    # R try with silent = true, returns nothing if error
    if score:
        result=robjects.r['try_CMHtest'](three_way_tab,types=types_r,overall=robjects.constants.TRUE,cscores=ord_names_r)
    else:
        result=robjects.r['try_CMHtest'](three_way_tab,types=types_r,overall=robjects.constants.TRUE,cscores="midrank")
    if len(result) > 1:
        return result.rx2('ALL').rx2('table').rx2('cor','Prob')[0]
    else:
        return "NA"
    # try:
    #     result=robjects.r['CMHtest'](three_way_tab,types=types_r,overall=robjects.constants.TRUE,cscores="midrank")
    #     return result.rx2('ALL').rx2('table').rx2('cor','Prob')[0]
    # except:
    #     return "NAN"
    # extracting the P value for the correlated test 
    


parser = argparse.ArgumentParser(description='calculate generalized CMH P values for correlation with an ordinal variable for rows of a sync file.') 
parser.add_argument("--reps", dest="replicates", help="Comma separated list of the ordinal value and replicate number seperated by \':\' of the samples in the CMH file. Samples with no : are omitted. The ordinal values with some sort of a sortable sequence (eg. A,B,C...) (def. \'a:1,a:2,a:3,b:1,b:2,b:3,c:1,c:2,c:3,d:1,d:2,d:3\')",default='a:1,a:2,a:3,b:1,b:2,b:3,c:1,c:2,c:3,d:1,d:2,d:3')
parser.add_argument("--score", dest="score",action="store_true", help="use the ordinal values as scores",default=False)
parser.add_argument("-i", dest="infile", help="Sync or cmhout file",required=True)
args = parser.parse_args()

replicates= vars(args)['replicates'].split(',')
reps=defaultdict(lambda: defaultdict(int) )
score=vars(args)['score']
for pos in range(0,len(replicates)) :
    try:
        (ord,rep)=replicates[pos].split(":")
    except:
        continue
    if score:
        ord = float(ord)
    reps[rep][ord]=pos
    
rep_names=sorted(reps.keys())
ord_names=sorted(reps[reps.keys()[0]].keys())
sort_reps=[ [ reps[x][y] for y in sorted(reps[x].keys()) ] for x in sorted(reps.keys())  ]
# setting up R env
robjects.packages.quiet_require("gnm")
robjects.packages.quiet_require("vcdExtra")
# r=robjects.r('library("gnm")')
# print "loaded gnm"
# r=robjects.r('library("vcdExtra")')
tab_dims=robjects.IntVector([2,len(ord_names),len(rep_names)])
ord_names_r=robjects.Vector(ord_names)
rep_names_r=robjects.Vector(rep_names)
alleles_r=robjects.Vector(["A1","A2"])
types_r=robjects.Vector(["ALL"])
robjects.r('''
        try_CMHtest <- function(three_way_tab,types,overall,cscores) {
        try( CMHtest(three_way_tab,types=types,overall=overall,cscores=cscores),
        silent=TRUE
        )
        }
        ''')

inf=open(vars(args)['infile'],'r')
for line in inf:
    if re.match("\s*\#",line):
        # comment lines
        continue
    line=line.rstrip()
    syncline=Sync_parser.SyncLineParser(line)
    syncline.get_reduced_seq()
    three_way=np.array([ syncline.reduced_seq[x] for x in sort_reps ],ndmin=3,dtype=float)
    three_way += 1e-6
    pV = cmh_test(three_way,tab_dims,ord_names_r,rep_names_r,alleles_r,types_r,score)
    if pV=="NA":
        print line+"\t"+pV
    else:
        print line+"\t{:.5g}".format(pV)


