
def read_inversion_file(infile):
    """
    reads martins inversion file (format: CHR INV BPS ALLELE) and creates a dict with chrom->bps->Inv->Allele]
    """
    snp_dict=defaultdict(lambda: defaultdict(defaultdict))
    inf = open(infile,"r")
    #load dictionary
    for line in inf:
        line=line.rstrip()
        fields=line.split()
        snp_dict[fields[0]][fields[2]][fields[1]]=fields[3]
    inf.close()
    return snp_dict

import Sync_parser
import sys, os, re
#import numpy as np
#from scipy import stats
from collections import defaultdict
import argparse
import gzip
parser = argparse.ArgumentParser(description='reads martin\'s inversion file and a sync/cmh file and gives the fraction of the alleles fixed in the inversion.') 

parser.add_argument("--inv", dest="inv", help="file with alleles segregating with inversions", required=True)
parser.add_argument("--infile", dest="infile", help="sync or cmh file to get positions from, can be gzipped", required=True)
parser.add_argument("--header", dest="head",  action="store_true", help="write header (default=False)", default=False)
parser.add_argument("--count", dest="count",  action="store_true", help="add overall count column (default=False)", default=False)
args = parser.parse_args()
infile = vars(args)['infile']
inv=vars(args)['inv']
head=vars(args)['head']
count=vars(args)['count']
nucs={'A':0,'T':1,'C':2,'G':3}
snp_dict=read_inversion_file(inv)

if re.search("\.b?gz",infile):
    inf = gzip.open(infile,'rb')
else:
    inf = open(infile,"r")
for line in inf:
    line=line.rstrip()
    fields=line.split()
    if ( fields[0] in snp_dict and fields[1] in snp_dict[ fields[0] ]):
        #print line
        entry = Sync_parser.SyncLineParser(line)
        afs = entry.get_pop_allele_freqs(two=False)
        if head:            
            if count:
                pops = [ "pop"+str(x)+"\tcov"+str(x) for x in range(1,len(afs)+1) ]
            else:
                pops = [ "pop"+str(x) for x in range(1,len(afs)+1) ]
                
            # if entry.cmhp:
            #     pops.append("P")
            print "CHR\tINV\tBPS\tREF\tALT\t"+"\t".join(pops)
            head=False
        for x in snp_dict[ fields[0] ][ fields[1] ].keys():
            allele=snp_dict[ fields[0] ][ fields[1] ][x]
            afs_x=[str(i[nucs[allele]]) for i in afs]
            if count:
                cts_x = [str(i) for i in entry.get_pop_coverages() ]
                afs_x = ["\t".join((i,j)) for (i,j) in zip(afs_x,cts_x)]
            out_str=fields[0]+"\t"+x+"\t"+fields[1]+"\t"+entry.ref+"\t"+allele+"\t"+"\t".join(afs_x)
            print out_str
inf.close()




