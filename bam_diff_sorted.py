######################### filtering a bam file by the ID's 

import pysam 
import sys
import os 
import collections
from optparse import OptionParser, OptionGroup
# edited by Lukas to support greater numbers of reads (> 200 Mio)
#from Bio import trie
#########################################################   HELP   #########################################################################
usage="python %prog --inA A.bam --inB B.bam"
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

Note that the pysam package needs to be installed (type: sudo easy_install pysam) for this.
Gives all reads in A that are not in B. Both files need to be sorted, need to have the same sequence of chromosomes.

""") 

parser.add_option("--inA", dest="inA", help="A BAM file")
parser.add_option("--inB", dest="inB", help="B BAM file")
parser.add_option_group(group)
(options, args) = parser.parse_args()

## read ID's from the fastq file
def get_ids(samB,lb,posA=0,tidA=-1):
    pos = lb.pos
    tid = lb.tid
    while lb.tid < tidA or ( lb.tid == tidA and lb.pos < posA):
        try:
            lb = samB.next()
        except:
            return ("EoF",False)
    ids=[lb.qname]
    for lbb in samB.fetch(until_eof=True):
        lb = lbb
        if lb.pos != pos:
            break
        ids.append(lb.qname)
    return ( lb, set(ids) )


## index BAM file if necessary
if not os.path.exists(options.inA+".bai"):
	print "indexing "+options.inA
	os.system("samtools index "+options.inA)

samA=pysam.Samfile(options.inA,"rb")
samB=pysam.Samfile(options.inB,"rb")
chroms={ samB.header['SQ'][x]['SN'] : x for x in range(0,len(samB.header['SQ'])) }
inf = options.inA
# remove .bam at end
if options.inA[-4:] == ".bam":
    inf = inf[:-4]
filtout=pysam.Samfile(inf+"_AnotB.bam","wb",template=samA)
#cleanout=pysam.Samfile(options.input+"_clean.bam","wb",template=samfile)
# read first read and pos of B
lb=samB.next()
posB=lb.pos
tidB=lb.tid
(lb,ids)=get_ids(samB,lb)
for la in samA.fetch(until_eof=True):
    posA=la.pos
    tidA=la.tid
    if (tidA > tidB or ( tidA == tidB and posA > posB)) and not (lb == "EoF"):
        # file A over current position in file B
        (lb,ids)=get_ids(samB,lb,posA,tidA)    
    if lb == "EoF":
        # eof of B reached
        filtout.write(la)
    elif tidA < tidB or ( tidA == tidB and posA < posB ) :
        filtout.write(la)
    elif tidA == tidB and posA == posB :
        if not ( la.qname in ids):
            filtout.write(la)
filtout.close()
samA.close()
samB.close() 
