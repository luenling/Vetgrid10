def combine_pops(pops,cov_p_ind,seqs):
    cov_p_ind=np.array([ cov_p_ind[x] for x in pops])
    seqs=[ seqs[x] for x in pops]
    ratio_cpi=cov_p_ind.min()/cov_p_ind
    scaled_counts=[ seqs[x]*ratio_cpi[x] for x in range(0,len(ratio_cpi)) ]
    return(np.array(sum(scaled_counts),dtype=int))

import sys, os, re
import numpy as np
import Sync_parser
from numpy import array, zeros
import math
import gzip
import argparse
#import rpy2.robjects as robjects
#from rpy2.robjects.packages import importr

parser = argparse.ArgumentParser(description="""
read sync file and combine two columns each, depending on the coverage per individual.
For each population combination it calculates the coverage per individual and the ratio of that to the minimal coverage per individual and then adds the nucleotide counts times that scaling factor.
It prints a sync file with combined popualitons.
""")

parser.add_argument("--pops","-p", dest="pops", help="the population pairs to combine \"|\" (eg. \"1+2,2+3,4+5\")",default=False)
parser.add_argument("--ind","-n", dest="inds", help="individuals per population (eg. \"65,75,100,140,110,202\")",default=False)

parser.add_argument("--in","-i", dest="infile", help="sync file to be read", required=True)
args = parser.parse_args()
infile = vars(args)['infile']
inds= np.array([ int(x) for x in vars(args)['inds'].split(",") ])
pops = [ [ int(y) -1 for y in x.split("+") ]  for x in vars(args)['pops'].split(",") ]

if re.search("\.b?gz",infile):
    inf = gzip.open(infile,'rb')
else:
    inf = open(infile,"r")

for line in inf:
    if re.match("\s*\#",line):
        # comment lines
        continue
    line=line.rstrip()
    syncline = Sync_parser.SyncLineParser(line)
    cov_p_ind=syncline.get_pop_coverages()/inds
    pop_counts = [ combine_pops(x,cov_p_ind,syncline.seqs)  for x in pops ]
    pop_countstr = "\t".join([  ":".join([ str(y) for y in x]) for x in pop_counts ])
    print syncline.chr+"\t"+str(syncline.pos)+"\t"+syncline.ref+"\t"+pop_countstr
inf.close()
        
                                                
