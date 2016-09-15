######################### filtering a bam file by the ID's 

import pysam 
import sys
import os 
import glob 
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

Note that the pysam package needs to be installed (type: sudo easy_install pysam) for this. This version also uses the trie data structure from the Biopython package, to support bigger numbers of reads. This script uses the Fwd dataset of the melanogaster Fastq file (output from Ram's script; --mel ) and the Fwd dataset of the simulans dataset (--sim) as the input and splits a BAM file accordingly in a _sim and a _mel BAM which are subsets of the original BAM. Additionally, a _missed BAM will be created which may contain reads which were not in the two fastq input files (although it is supposed to be empty). This may for expample happen if some reads were not mapped in the gsnap approach and therefore were not in any of the fastq outputs.




""") 

parser.add_option("--input", dest="input", help="A BAM file or an expression globing to a list of BAM files")
parser.add_option("--mel", dest="mel", help="the melanogaster specific fwd FASTQ")
parser.add_option("--sim", dest="sim", help="the melanogaster specific fwd FASTQ")

parser.add_option_group(group)
(options, args) = parser.parse_args()

## read ID's from the fastq file
def get_ids(x):
	newtrie=trie.trie()
	fastq=open(x)
	while(True):
		id=fastq.readline()[1:-3]
		newtrie[id]=len(id)
		fastq.readline()
		fastq.readline()
		fastq.readline()
		if id=="":
			break
	return newtrie

print "creating trie out of "+options.mel
sys.stdout.flush()
mel=get_ids(options.mel)
print "creating trie out of "+options.sim
sys.stdout.flush()
sim=get_ids(options.sim)

for infile in glob.glob(options.input):
    ## index BAM file if necessary
    print "processing "+infile
    sys.stdout.flush()
    if not os.path.exists(infile+".bai"):
        print "indexing "+infile
        os.system("samtools index "+infile)
    
    samfile=pysam.Samfile(infile,"rb")

    melout=pysam.Samfile(infile+"_mel","wb",template=samfile)
    simout=pysam.Samfile(infile+"_sim","wb",template=samfile)
    missedout=pysam.Samfile(infile+"_missed","wb",template=samfile)
    ## split BAM file
    for l in samfile.fetch(until_eof=True):
        if mel.has_key(l.qname) > 0:
            melout.write(l)
        elif sim.has_key(l.qname) > 0:
            simout.write(l)
        else:
            missedout.write(l)
            
    melout.close()
    simout.close()
    missedout.close()
    samfile.close() 
