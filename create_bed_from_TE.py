#import re,sys
#from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description="""
gets a TE frequency file from robert's popoo TE and creates a bed file with CHR:pos:strand:name:freq as the name.

""") 
parser.add_argument("-i", dest="infile", help="infile", required=True)

args = parser.parse_args()
infile = vars(args)['infile']

inf=open(infile,"r")

for line in inf:
    fields=line.split("\t")
    name = ":".join(fields[0:5])
    fam=fields[3]
    freq = fields[4]    
    in_ref = fields[6]
    ranges=[]
    if "F" in fields[2]:
        ranges.extend([ int(x) for x in fields[8:10] ])
    if "R" in fields[2]:
        ranges.extend([ int(x) for x in fields[15:17] ])
    print "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(fields[0],min(ranges)-1, max(ranges) , name, fam, freq, in_ref )

    

