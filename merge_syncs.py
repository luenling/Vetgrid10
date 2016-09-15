
def split_syncline(syncline):
	syncline = syncline.rstrip()
	entries=syncline.split()
	return entries

def get_next_line(inf):
	try:
		l=inf.next()
	except:
		l=False
	return l

def pad_sync(num):
	pentries=entries[0:3]
	pentries.extend(["0:0:0:0:0:0"]*num)
	return pentries


def catchup_chrom(inf,chrom):
	entries=[]
	while(True):
		l=get_next_line(inf)
		if (not l):
			entries=l
			break
		entries=split_syncline(l)		
		if (entries[0] == chrom):
			break
	return entries

def catchup_line_num(inf,line_num,chrom):
	entries=[]	
	while(True):
		l=get_next_line(inf)
		if (not l):
			entries=l
			break
		entries=split_syncline(l)
		if (int(entries[1]) >= line_num) or (entries[0] != chrom):
			break
	return entries

def read_snp_file(infile):
    """
    reads a cmh or sync like file and creates a dict with chrom->list(bps)
    """
    snp_dict=defaultdict(set)  #dictionary of chroms with positions and pVals of snps
    inf = open(infile,"r")
    #load dictionary
    for line in inf:
        line=line.rstrip()
        fields=line.split()
        snp_dict[fields[0]].add(int(fields[1]))
    inf.close()
    return snp_dict


import sys, os, re
import numpy as np
#from scipy import stats
from collections import defaultdict
#np.seterr(all=None, divide='ignore', over=None, under=None, invalid='ignore')
import argparse
parser = argparse.ArgumentParser(description='reads snps from file A and only gets the lines with the same positions from file B. writes to stdout') 

parser.add_argument("--inA", dest="inA", help="sync or cmh file A (best to take the smaller one)", required=True)
parser.add_argument("--inB", dest="inB", help="sync or cmh file B", required=True)
parser.add_argument("--out",  dest="out", help="outfile", required=True)

args = parser.parse_args()
inA = vars(args)['inA']
inB = vars(args)['inB']
outf = vars(args)['out']
# get dict of all snps
snp_dict=read_snp_file(inA)
# open files for reading again
inB=open(inB,'r')
inA=open(inA,'r')
out=open(outf,'w')
count=0
for line in inB:
    entriesB=split_syncline(line)
    # see if position in file A
    if ( (entriesB[0] in snp_dict) and (int(entriesB[1]) in snp_dict[entriesB[0]]) ):
        entriesA=split_syncline(inA.next())
        if (entriesA[0] != entriesB[0]):
            # one chrom lagging behind
            entriesA=catchup_chrom(inA,entriesB[0])
            # just for debugging, should not happen
            if  (entriesA[0] != entriesB[0]) or (not entriesA):
                print entriesA, entriesB
                break
        if (entriesA[1] != entriesB[1]):
            entriesA=catchup_line_num(inA,int(entriesB[1]),entriesB[0]) 
            if (not entriesA) or (entriesA[0] != entriesB[0]):
                print entriesA, entriesB
                break
        # line matches! woohoo!
        entries = entriesA
	entries.extend(entriesB[3:])
	out_line="\t".join(entries)
	count+=1
	print >>out, out_line

inA.close()
inB.close
out.close()
print "lines printed: "+str(count)

