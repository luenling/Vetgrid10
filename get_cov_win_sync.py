import sys
import re
import gzip
import collections
import operator
import numpy as np
#from optparse import OptionParser, OptionGroup
import argparse
#########################################################   HELP   #########################################################################
#print
parser = argparse.ArgumentParser(description=
"""
    			D E S C R I P T I O N
Reads a sync file and writes a file with average coverage over windows of a given length. The windows are nonoverlapping and fixed to genomic coordinates starting from the same position at each chromosome.
most code stolen from Ray 
""") 

##########################################################  FUNCTIONS ######################################################################
def get_win_cov(win_cov, pop_num): # gets array with coverages and returns vector with mean
    win_cov=np.reshape(win_cov,(-1,pop_num))
    return np.mean(win_cov,0) # return mean as array

parser.add_argument('--sync','-s', dest="sync", help='Sync file.', required=True)
parser.add_argument('--popcols','-p', dest="popcols", help="Columns with allele count data to be outputted - note: column 1 = first column with count data, default=all columns")
parser.add_argument('--verbose','-v', dest="verb", action="store_true", help="print progress default: no output", default=False)
parser.add_argument("--output",'-o', dest="output", help="An output file (default=stdout)", required=False)
parser.add_argument("--ws",'-w', dest="ws", type=int, help="Window size (default=100)", default=100)
parser.add_argument("--start", dest="start", type=int, help="Starting position of first window on each Chromosome (default=5000)", default=5000)
parser.add_argument("--fst", dest="fst", action="store_true", help="File is fst file (from PoPoolation2) (default=False)", default=False)

args = parser.parse_args()
ws=vars(args)['ws']
start=vars(args)['start']
popcols=vars(args)['popcols']
fst=vars(args)['fst']

pop_num = 0
#columns with counts
if (not popcols): # if not set empty allpopcols array
	allpopcols = []
else:     # translate population columns to indeces
	allpopcols = popcols.split(",")
        if fst:
            allpopcols = [int(i)+2 for i in allpopcols ]
	else:
            allpopcols = [int(i)+4 for i in allpopcols ]
        pop_num = len(allpopcols)

verb=vars(args)['verb']

# create output file
output_file =  vars(args)['output']
if output_file:
    ofh = open(output_file,"w")
else:
    ofh=sys.stdout

# parse input file and generate output
infile = vars(args)['sync']
if re.search("\.b?gz",infile):
    inf = gzip.open(infile,'rb')
else:
    inf = open(infile,"r")


chrom=""
win=0
start_pos=start
end_pos=start_pos+ws
covs=[]
for r in inf:
    if re.match("\s*\#",r):
        # comment lines
        continue
    r=r.rstrip()
    r=r.rstrip()
    vals = r.split() # splits each row into columns separated by tab    
    if (not allpopcols):
        if fst:
            allpopcols=range(5,len(vals))
        else:
            allpopcols=range(3,len(vals))
        # if CMH file drop last 
        if not(re.search(":",vals[-1])):
            vals.pop()
            allpopcols.pop()
        pop_num = len(allpopcols)
    if chrom != vals[0]:
        if covs:
            if fst:
                mean_covs = "\t".join( [ "{:.5f}".format(x) for x in get_win_cov(covs, pop_num) ] )
            else:
                mean_covs = "\t".join( [ "{:.1f}".format(x) for x in get_win_cov(covs, pop_num) ] )
            win_pos = int((start_pos + end_pos)/2)
            print >> ofh, chrom+"\t"+str(win_pos)+"\t"+ mean_covs
            covs=[]
        start_pos=start
        end_pos=start_pos+ws
        win=0
        chrom=vals[0]
    curr_pos = int(vals[1])
    if curr_pos < start_pos: # before start position
            continue
    if curr_pos >= end_pos:
        # next window
        if covs:
            if fst:
                mean_covs = "\t".join( [ "{:.5f}".format(x) for x in get_win_cov(covs, pop_num) ] )
            else:
                mean_covs = "\t".join( [ "{:.1f}".format(x) for x in get_win_cov(covs, pop_num) ] )
            mean_covs = "\t".join( [ "{:.1f}".format(x) for x in get_win_cov(covs, pop_num) ] )
            win_pos = str(int((start_pos + end_pos)/2))
            print >> ofh, chrom+"\t"+win_pos+"\t"+ mean_covs
            covs=[]
        win = int((curr_pos-start)/ws)
        start_pos=start+win*ws
        end_pos=start_pos+ws                
        if fst:
            fsts=np.array([ vals[x].split("=")[1] for x in allpopcols ],dtype=float)
            covs.extend(fsts)
            print covs
        else:
            for i in allpopcols:
                alleles=np.array(vals[i].split(":")[0:3],dtype=int)
                covs.append(alleles.sum())
if covs:
    if fst:
        mean_covs = "\t".join( [ "{:.5f}".format(x) for x in get_win_cov(covs, pop_num) ] )
    else:
        mean_covs = "\t".join( [ "{:.1f}".format(x) for x in get_win_cov(covs, pop_num) ] )
    win_pos = str(int((start_pos + end_pos)/2))
    print >> ofh, chrom+"\t"+win_pos+"\t"+ mean_covs
    covs=[]

ofh.close()
f.close()


