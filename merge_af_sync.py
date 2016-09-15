import sys, os, re, gzip
import argparse
import numpy as np
from collections import defaultdict
def split_syncline(syncline):
	syncline = syncline.rstrip()
	entries=syncline.split()
	return entries

def pad_pileup(entries,num):
	pentries=entries[0:3]
	pentries.extend(["0","*"]*num)
	return pentries

def get_next_line(inf):
	try:
		l=inf.next()
	except:
		l=False
	return l

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
    reads a cmh or sync like file and creates a dict with chrom->set(bps)
    """
    snp_dict=defaultdict(set)  #dictionary of chroms with set of bps
    inf = open(infile,"r")
    #load dictionary
    for line in inf:
        line=line.rstrip()
        fields=line.split()
        snp_dict[fields[0]].add(int(fields[1]))
    inf.close()
    return snp_dict



#Author: Lukas Endler
parser = argparse.ArgumentParser(description='merge two allelefrequency sync files. put populations of one after the other. ommit positions existing in only one file. in both files chromosomes have to come in the same sequence, no chromosome can be missing, and the bps have to be in increasing order. Files can be gzipped') 
parser.add_argument("--infile1","-i1", dest="infile1", help="first sync file", required=True)
#parser.add_argument("--pV", dest="pV",  action="store_true", help="allele freq files contain pV column (default: False)", default=False)
parser.add_argument("--avg","-a", dest="avg", help="create averages of afs in populations in each files eg. 1-3,4-6 (default: None)", default=None)
#parser.add_argument("--type", dest="type", help="file type (mpileup/sync) default:sync", default="sync")
#parser.add_argument("--popnum", dest="pops", help="number of populations in pileup files n1:n2", default="0")
parser.add_argument("--infile2","-i2", dest="infile2", help="second sync file", required=True)
parser.add_argument("--out","-o", dest="outfile", help="output sync file, if \"stdout\" write to stdout", required=True)

args = parser.parse_args()
infile1 = vars(args)['infile1']
infile2 = vars(args)['infile2']
outfile = vars(args)['outfile']
#pV = vars(args)['pV']
avg= vars(args)['avg']
if avg:
    avg=[ x.split("-") for x in avg.split(",")]
    avg=[ [int(x[0])+1,int(x[1])+1]  for x in avg]

# file_type = vars(args)['type']
# if (file_type == "mpileup"):
# 	pops=vars(args)['pops'].split(:)
# 	if (len(pops) != 2): sys.exit("population numbers unclear")
# get dict of all snps
snp_dict=read_snp_file(infile1)
if re.search("\.b?gz",infile1):
    in1 = gzip.open(infile1,'rb')
else:
    in1 = open(infile1,"r")
if re.search("\.b?gz",infile2):
    in2 = gzip.open(infile2,'rb')
else:
    in2 = open(infile2,"r")
if outfile == "stdout" or outfile == "STDOUT":
    out=sys.stdout
else:
    out = open(outfile,"w")
count=0

for line in in2:
	entries2=split_syncline(line)
	# see if position in file A
	if not ( (entries2[0] in snp_dict) and (int(entries2[1]) in snp_dict[entries2[0]]) ):
		continue
	entries1=split_syncline(in1.next())
	if (entries1[0] != entries2[0]):
		# one chrom lagging behind
		entries1=catchup_chrom(in1,entries2[0])
		# just for debugging, should not happen
		if  (entries1[0] != entries2[0]) or (not entries1):
			print entries1, entries2
			break
	if (entries1[1] != entries2[1]):
		entries1=catchup_line_num(in1,int(entries2[1]),entries2[0]) 
		if (not entries1) or (entries1[0] != entries2[0]):
			print entries1, entries2
			break
        # line matches! woohoo!
	pV1=False
	pV2=False
        if (len(entries1[3:])%2):
            pV1 = entries1.pop()
	if (len(entries2[3:])%2):
            pV2 = entries2.pop()
        if avg:
            pop1=np.array(entries1[3:],dtype=float)
            pop2=np.array(entries2[3:],dtype=float)
            avg1 = [ pop1[pos[0]:pos[1]+1].mean() for pos in avg]
            avg2 = [ pop2[pos[0]:pos[1]+1].mean() for pos in avg]
            entries1 = entries1[0:3] + list(avg1)
            entries2 = entries2[0:3] + list(avg2)
            
        #check whether major and minor allele same or at least major are the same:
        if (entries1[2][0] != entries2[2][0]):
            if (entries1[2][0] == entries2[2][1]): # major of A is minor B (swap B)
		if (entries1[2][1] != entries2[2][0]): # min A not equal maj B
			entries1[2] = "|".join([ entries1[2],entries2[2][::-1] ]) # add reversed major minor of entrie 2
		entries2[3::2] = [ 1 - float(x) for x in entries2[3::2] ]    #1-maj af of second file
            elif (entries1[2][1] == entries2[2][0]): # major B is minor A
                entries1[2] = "|".join([ entries1[2][::-1], entries2[2] ]) # reverse major minor of entrie 1
		entries1[3::2] = [ 1 - float(x) for x in entries1[3::2] ]    #1-maj af of first file
            elif (entries1[2][1] == entries2[2][1]):
                entries1[2] = "|".join([ entries1[2][::-1], entries2[2][::-1] ]) # reverse both major minor
		entries1[3::2] = [ 1 - float(x) for x in entries1[3::2] ]    #1-maj af of first file
		entries2[3::2] = [ 1 - float(x) for x in entries2[3::2] ]    #1-maj af of first file
            else:
                print str(entries1)+"\n"+str(entries2)
                next
	elif (entries1[2][1] != entries2[2][1]):
            entries1[2] =entries1[2] + "|" +  entries2[2][1] # add second minor
        if pV1:
            entries1.append(pV1)
	if pV2:
            entries2.append(pV2)
	entries = entries1[0:]
        entries.extend(entries2[3:])
        entries = [ str(x) for x in entries ]
	out_line="\t".join(entries)
	count+=1
	print >>out, out_line
	#print out_line
out.close()
in1.close()
in2.close()
print "lines printed: "+str(count)
