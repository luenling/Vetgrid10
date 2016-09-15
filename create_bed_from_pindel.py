#import re,sys
#from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description="""
gets a pindel file and creates a bed file with Type:start:stop:support as the name.

""") 
parser.add_argument("-i", dest="infile", help="infile", required=True)

args = parser.parse_args()
infile = vars(args)['infile']

inf=open(infile,"r")

for line in inf:
    line = line.rstrip()
    fields= line.split()
    chrom = fields[7]
    start= fields[9]
    stop= fields[10]
    typ = fields[1]
    support = fields[16]
    length = fields[2]
    rest = "\t".join(fields[31:])
    name = ":".join([typ,chrom,start,stop,support])
    print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrom,int(start)-1, stop , name, typ, length, support, rest )

    

