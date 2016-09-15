######################### filtering a bam file by the ID's 

import pysam 
import sys
import os 
import collections
from optparse import OptionParser, OptionGroup
#########################################################   HELP   #########################################################################
usage="python %prog --input input.bam --out output.bam --dif reads_to_discard.bam"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

Note that the pysam package needs to be installed (type: sudo easy_install pysam) for this.
This script takes two bam files and produces discards all reads from the input file with ids identical to the ones in the second file.


""") 

parser.add_option("--input", dest="input", help="input BAM file")
parser.add_option("--out", dest="out", help="an ouput file")
parser.add_option("--dif", dest="dif", help="a bam file containing the reads to discard")

parser.add_option_group(group)
(options, args) = parser.parse_args()

## read ID's from the bam file
def get_ids(x):
	read_ids = set()
	samfile=pysam.Samfile(x,"rb")
	for l in samfile.fetch(until_eof=True):
		read_ids.add(l.qname)
	samfile.close()
	return read_ids

print "Reading IDs to discard from file "+options.dif+"..."
read_ids=get_ids(options.dif)
print str(len(read_ids))+" IDs read"

samfile=pysam.Samfile(options.input,"rb")
out=pysam.Samfile(options.out,"wb",template=samfile)
all_count=0
acc_count=0
for l in samfile.fetch(until_eof=True):
	all_count+=1
	if all_count%1000000 == 0 :
		print str(all_count)+" reads"
	if not(l.qname in read_ids):
		out.write(l)
		acc_count+=1
samfile.close() 
out.close()
print "total reads: "+str(all_count)
print "accepted reads: "+str(acc_count)
print "rejected reads: "+str(all_count-acc_count)
