#!/usr/bin/env python
def get_sliding_averages(bps,pVals,winsize=5000,quantile=0.25):
    """
    gets 2 np arrays with positions and pValues and gives two arrays with averages over a sliding window and the window positions back. optional give window size 
    """
    pVals_NaN = np.isnan(pVals)
    pVals[pVals_NaN]=1.0
    pVals_ln = -1.0 * np.log(pVals)
    averages_ln=[]
    # use for weighting
    p_rank=pVals_ln.argsort().argsort() + 1
    pVals_w = pVals_ln * p_rank
    averages_q=[]
    win_pos=[]
    start=0 # index for starting 
    numerator_ln=0.0
    # go over all positions
    for i in range(0,len(bps)):
        if (bps[i] < bps[start]):
            print "Warning: SNPs not at consecutive positions at position {0} and {1} in some godforsaken chromosome".format(bps[start],bps[i])
            break
        while ((bps[i]-bps[start]) > winsize):
            numerator_ln -= pVals_ln[start]
            start +=1
            if (bps[i]-bps[start-1]+1 == winsize):
                # rare case in which no window would fit inbetween (bps[start-1] and bps[i] just window size + 1 apart)
                break
            averages_ln.append(float(numerator_ln/(i-start)))
            win_pos.append(bps[start]+winsize/2) # only approximate
            #for calculating only a quantile
            # for weighted pVals
            win_values=pVals_w[start:i]
            win_weights=p_rank[start:i]
            sort_index = win_values.argsort()
            win_values=win_values[sort_index]
            win_weights=win_weights[sort_index]
            low_lim=int(len(win_values)*(1-quantile)) # index from which to start summation
            numerator = win_values[low_lim:].sum()
            denominator = win_weights[low_lim:].sum()
            #averages_q.append(float(numerator/(len(win_values)-low_lim )))
            averages_q.append(float(numerator/denominator))
        if start==i and i > 0:
            # no more entries (snps more than ws apart)
            # window value equals 0
            averages_ln.append(0)
            averages_q.append(0)
            win_pos.append((bps[start] + bps[start-1])/2)
        # numerator_log += pVals_log[i]
        numerator_ln += pVals_ln[i]
        averages_ln.append(float(numerator_ln/(i-start+1)))
        #for calcualting only a quantile
        win_values=pVals_ln[start:i+1]
        win_weights=p_rank[start:i+1]
        sort_index = win_values.argsort()
        win_values=win_values[sort_index]
        win_weights=win_weights[sort_index]
        low_lim=int(len(win_values)*(1-quantile)) # index from which to start summation
        numerator = win_values[low_lim:].sum()
        denominator = win_weights[low_lim:].sum()
        #averages_q.append(float(numerator/(len(win_values)-low_lim )))
        averages_q.append(float(numerator/denominator))
        win_pos.append(bps[i]-winsize/2 if (bps[i] >= winsize) else bps[i]/2)
    return (win_pos,averages_ln,averages_q)
    

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
parser.add_argument("--qu", dest="quantile", type=float, help="fraction of highest pValues taken for the average (default: 0.25)",default=0.25)
args = parser.parse_args()
infile = vars(args)['infile']
outfile = vars(args)['outfile']
winsize = vars(args)['windowsize']
quantile = vars(args)['quantile']
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
    (win_pos,averages_ln,averages_q)=get_sliding_averages(cmh_dict[chrom]['bps'],cmh_dict[chrom]['pVals'],winsize,quantile)
    for i,pos in enumerate(win_pos):
        print >>outf,"{0}\t{1}\t{2}\t{3}".format(chrom,pos,averages_ln[i],averages_q[i])
outf.close
