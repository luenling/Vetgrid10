def read_snp_file(infile):
    """
    reads a sync like file with majmin alleles pV odds ratios and intervals and creates a chrom->bps->[alleles,[ [pV,oR,oRa,oRb] ...]]
    """
    snp_dict=defaultdict(defaultdict)  #dictionary of chroms with positions and pVals of snps
    if re.search("\.b?gz",infile):
	inf = gzip.open(infile,'rb')
    else:
        inf = open(infile,"r")
    #load dictionary
    pops=0
    for line in inf:
        if re.match("\#",line):
            continue
        line=line.rstrip()
        fields=line.split()
        a=np.array(fields[3:],dtype=float)
        a.shape=[len(a)/4,4]
        snp_dict[fields[0]][fields[1]]=[fields[2],a]
        if pops == 0:
            # gets first population count != 0
            pops=a.shape[0]
    inf.close()
    return (pops,snp_dict)

import sys, os, re, gzip
import numpy as np
#from scipy import stats
from collections import defaultdict
#np.seterr(all=None, divide='ignore', over=None, under=None, invalid='ignore')
import argparse
parser = argparse.ArgumentParser(description='reads pos from file A and only gets the lines with the same positions from file B and changes the maj minor alleles and ORs accordingly to file B. writes B + A to stdout or only the changed A if -oa ') 

parser.add_argument("--inA","-a", dest="inA", help="file with odds", required=True)
parser.add_argument("--incon","-i",action='store_true', dest="incon", help="print missing/insonsistant lines as nans ", default=False)
parser.add_argument("--oa", dest="oa",action='store_true', help="only change allele order and odds ratios in file a and print it", default=False)
parser.add_argument("--verbose", dest="verb",action='store_true', help="print crap on stderror", default=False)
parser.add_argument("--inB","-b", dest="inB", help="sync or cmh file to get positions from", required=True)

args = parser.parse_args()
inA = vars(args)['inA']
inB = vars(args)['inB']
oa = vars(args)['oa']
incon = vars(args)['incon']
verb = vars(args)['verb']

(pops,snp_dict)=read_snp_file(inA)
if re.search("\.b?gz",inB):
    inf = gzip.open(inB,'rb')
else:
    inf = open(inB,"r")

if incon:
    incostr=pops*4*"\tnan"
    #print incostr
    #print str(pops)

for line in inf:
    if re.match("\#",line):
        continue
    line=line.rstrip()
    fields=line.split()
    if ( fields[0] in snp_dict and fields[1] in snp_dict[ fields[0] ]):
        if snp_dict[ fields[0] ][ fields[1] ][0] == fields[2][::-1]:
            snp_dict[ fields[0] ][ fields[1] ][1][:,1:] = np.exp(-1*np.log(snp_dict[ fields[0] ][ fields[1] ][1][:,1:]))
        elif snp_dict[ fields[0] ][ fields[1] ][0] != fields[2]:
            if incon:
                #set to nan
                snp_dict[ fields[0] ][ fields[1] ][1] *= np.nan
            else:
                continue
        line_a= "\t".join([ str(x) for x in snp_dict[ fields[0] ][ fields[1] ][1].flatten()])
        if oa:
            print "\t".join(fields[0:3])+"\t"+line_a
        else:
            print "\t".join(fields[0:3])+"\t"+"\t".join(fields[3:])+"\t"+line_a
    elif incon:
        if oa:
            print "\t".join(fields[0:3])+incostr
        else:
            print "\t".join(fields[0:3])+"\t"+"\t".join(fields[3:])+incostr
inf.close()

