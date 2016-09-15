######################### filtering a bam file by the ID's 

import pysam 
import sys
import os 
import collections
from optparse import OptionParser, OptionGroup
# edited by Lukas to support greater numbers of reads (> 200 Mio)
from Bio import trie
#########################################################   HELP   #########################################################################
usage="python %prog --input contaminated.bam --sim simulans_1.fastq --mel melanogaster_1.fastq"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

Note that the pysam package needs to be installed (type: sudo easy_install pysam) for this. This version also uses the trie data structure from the Biopython package, to support bigger numbers of reads. This script takes a set of reads from a text file (one read ID per line) of reads to be filtered out and  as the input and splits a BAM file accordingly into a filtered and not filtered BAM which are subsets of the original BAM.

""") 

parser.add_option("--input", dest="input", help="A BAM file")
parser.add_option("--filter", dest="filt", help="reads to filtered out")
parser.add_option_group(group)
(options, args) = parser.parse_args()

## read ID's from the fastq file
def get_ids(x):
    newtrie=trie.trie()
    filt_file= open(x,"r")
    for line in filt_file:
        line = line.rstrip()
        id = line.split()[0]
        newtrie[id]=len(id)
        if id=="":
            break
    return newtrie

print "reading "+ options.filt
filt_reads=get_ids(options.filt)
print "created trie for filtering"

## index BAM file if necessary
if not os.path.exists(options.input+".bai"):
	print "indexing "+options.input
	os.system("samtools index "+options.input)

samfile=pysam.Samfile(options.input,"rb")
inf = options.input
# remove .bam at end
if options.input[-4:] == ".bam":
    inf = inf[:-4]
filtout=pysam.Samfile(options.input+"_filt.bam","wb",template=samfile)
cleanout=pysam.Samfile(options.input+"_clean.bam","wb",template=samfile)
## split BAM file
for l in samfile.fetch(until_eof=True):

	if filt_reads.has_key(l.qname) > 0:
		filtout.write(l)
	else:
		cleanout.write(l)
filtout.close()
cleanout.close()
samfile.close() 
