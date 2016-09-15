#!/usr/bin/env python
def get_sliding_averages(bps,pVals,winsize=5000):
    """
    gets 2 np arrays with positions and pValues and gives two arrays with averages over a sliding window and the window positions back. optional give window size 
    """
    pVals_NaN = np.isnan(pVals)
    pVals[pVals_NaN]=1.0
    pVals_ln = -1.0 * np.log(pVals)
    #    pVals_log = -1.0 * np.log10(pVals)
    # use for weighting
    p_rank=pVals_ln.argsort().argsort() + 1
    averages_w=[]
    pVals_w = pVals_ln * p_rank
    pVals_ln = np.square(pVals_ln)
    averages_ln=[]
    win_pos=[]
    start=0 # index for starting 
    numerator_w=0.0
    denominator_w=0.0
    numerator_ln=0.0
    # go over all positions
    for i in range(0,len(bps)):
        if (bps[i] < bps[start]):
            print "Warning: SNPs not at consecutive positions at position {0} and {1} in some godforsaken chromosome".format(bps[start],bps[i])
            break
        while ((bps[i]-bps[start]) > winsize):
            # numerator_log -= pVals_log[start]
            numerator_ln -= pVals_ln[start]
            if numerator_ln < 0:
                numerator_ln = 0
            numerator_w -=  pVals_w[start]
            denominator_w -= p_rank[start]
            # averages_log.append(float(numerator_log/(i-start)))
            averages_ln.append(float(np.math.sqrt(numerator_ln)/(i-start)))
            averages_w.append(float(numerator_w/denominator_w))
            win_pos.append(bps[start]+winsize/2) # only approximate
            start+=1
        if start==i and i > 0:
            # no more entries (snps more than ws apart)
            averages_ln.append(0)
            averages_w.append(0)
            win_pos.append((bps[start] + bps[start-1])/2)
        # numerator_log += pVals_log[i]
        numerator_ln += pVals_ln[i]
        numerator_w +=  pVals_w[i]
        denominator_w += p_rank[i]
        averages_ln.append(float(np.math.sqrt(numerator_ln)/(i-start+1)))
        averages_w.append(float(numerator_w/denominator_w))
        win_pos.append(bps[i]-winsize/2 if (bps[i] >= winsize) else bps[i]/2)
    return (win_pos,averages_ln,averages_w)
    

import sys,re
import collections
import numpy as np
from numpy import array, zeros
np.seterr(all=None, divide='ignore', over=None, under=None, invalid='ignore')
#from scipy.stats import chisquare
import argparse
parser = argparse.ArgumentParser(description='calculates a sliding window over sequence positions in a chm file and gives averages') 
parser.add_argument("--in", dest="infile", help="cmh file to read", required=True)
parser.add_argument("--out", dest="outfile", help="file to output data",required=True)
parser.add_argument("--ws", dest="windowsize", type=int, help="windowsize (default: 5000)",default=5000)
args = parser.parse_args()
infile = vars(args)['infile']
outfile = vars(args)['outfile']
winsize = vars(args)['windowsize']
chroms=[]
# intermediate lists of values
bps=[]
pVals=[] 
cmh_dict={}  #dictionary of chroms with positions and pVals of snps
inf = open(infile,"r")
#load dictionary
for line in inf:
    line.rstrip()
    fields=line.split()
    if not (fields[0] in chroms):
        chroms.append(fields[0])
        if (len(chroms) > 1):
            cmh_dict[chroms[-2]]={}
            cmh_dict[chroms[-2]]['bps']=np.array(bps)
            cmh_dict[chroms[-2]]['pVals']=np.array(pVals)
            bps=[]
            pVals=[]
    pVals.append(float(fields[-1]))
    bps.append(int(fields[1]))
# set last chromosome
cmh_dict[chroms[-1]]={}
cmh_dict[chroms[-1]]['bps']=np.array(bps)
cmh_dict[chroms[-1]]['pVals']=np.array(pVals)
inf.close()
outf = open(outfile,"w")
#go through chromosomes
for chrom in chroms:
    (win_pos,averages_ln,averages_w)=get_sliding_averages(cmh_dict[chrom]['bps'],cmh_dict[chrom]['pVals'],winsize)
    for i,pos in enumerate(win_pos):
        print >>outf,"{0}\t{1}\t{2}\t{3}".format(chrom,pos,averages_ln[i],averages_w[i])
outf.close
