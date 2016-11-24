######################### filtering a bam file by tags
from __future__ import print_function
import pysam 
import sys, re
# for fuzzy matching
import regex
import os 
import argparse
import collections
#########################################################   HELP   #########################################################################
parser = argparse.ArgumentParser(description='Filter a bam file by read barcodes giving bam files named after the base and the tag sequence. needs the oython packages pysam and regex to be installed. It allows to give a maximal number of substitutions, although that wont work for the first position of the barcode. the default value for permitted substitutions is 1.')

parser.add_argument("--in", dest="infile", help="bam file with tagged reads", required=True)
parser.add_argument("--tags", dest="tags", help="comma separated list of tags (eg: \"TGACCAAT,ACAGTGAT,GCCAATAT,CTTGTAAT\") or pairs of ids and tags (eg: \"light_RI:TGACCAAT,light_RII:ACAGTGAT,dark_RI:GCCAATAT,dark_RII:CTTGTAAT\") or name of tag file", required=True)
parser.add_argument("--subs", dest="subs", type=int, help="maximal number of substitutions in tag sequence (default=1)", default=1)
parser.add_argument("--enz1", dest="enz1",  help="read 1 restr. enzyme overlap (default=TGCAG)", default="TGCAG")
parser.add_argument("--enz2", dest="enz2", help="read 2 restr. enzyme overlap (default=CGG)", default="CGG")
parser.add_argument("--single", dest="single", help="single end reads (default=False)", default=False, action="store_true")

args = parser.parse_args()
infile = vars(args)['infile']
subs = vars(args)['subs']
single=vars(args)['single']
if  os.path.isfile(vars(args)['tags']):
    tags=list()
    with open(vars(args)['tags'], "r") as f:
#    with open(tagf, "r") as f:
        for line in f:
            line = line.rstrip()
            if (re.match("^\s*\#+",line) or re.match("^\s*$",line)):
                continue
            line=re.sub("\s+",":",line)
            tags.append(line)            
else:
    tags= vars(args)['tags'].split(",")
ids=list(tags)
enz2=vars(args)['enz2']
enz1=vars(args)['enz1']

for i,tag in enumerate(tags):
	if (re.search(":",tag)):
		(ids[i],tags[i])=tags[i].split(":")
# remove whitespace in tags, if there is any
tags=[ re.sub("\s","",x) for x in tags]
ids=[ re.sub("\s","",x) for x in ids]
outfiles=[]
## index BAM file if necessary
#if not os.path.exists(infile+".bai"):
#	print "indexing "+infile
#	os.system("samtools index "+infile)
base_name=os.path.basename(infile)
# get rid of trailing bam
base_name=re.sub("\.bam(?=$)","",base_name)
samfile=pysam.Samfile(infile,"rb",check_sq=False)
tag_files=collections.defaultdict()
tag_count=collections.defaultdict(int)
tag_pats=collections.defaultdict()
count=0
# open files for tag splitting
for i,tag in enumerate(tags):
	tag_files[tag]=pysam.Samfile(base_name+"_"+ids[i]+".bam","wb",template=samfile)
	tag_pats[tag]=regex.compile("(?:"+tag+enz1+"){s<="+str(subs)+"}")

enz2_pat=regex.compile("(?:"+enz2+"){s<="+str(subs)+"}")
## split BAM file
#import timeit
#start_time = timeit.default_timer()
# code you want to evaluate

for entry in samfile.fetch(until_eof=True):
    if entry.is_read2:
        continue
#    read_tag=entry.qname.split('#')[-1]
    count+=1
    if (count%100000 == 0):
        print(str(count)+' readpairs gone through', file=sys.stderr) 
    for tag in tags:
        if regex.match(tag_pats[tag],entry.seq):
            if single:
                tag_count[tag] += 1
                tag_files[tag].write(entry)
                break
            entry2 = samfile.next()
            if (entry2.qname != entry.qname):
                sys.exit("mates not in proper order!")
#            if (entry2.seq[0:len(enz2)] != enz2):
            if not regex.match(enz2_pat,entry2.seq):
                break
            tag_count[tag] += 1
            tag_files[tag].write(entry)
            tag_files[tag].write(entry2)            
            break
#elapsed = timeit.default_timer() - start_time

print("read pairs read:\t{}".format(count))
tot_tag=0
for i,tag in enumerate(tags):
    print( "tag {} id {}:\t{}".format(tag,ids[i],tag_count[tag]) )
    tot_tag += tag_count[tag]
    tag_files[tag].close
print("reads not taggged:\t{}".format(count-tot_tag))
samfile.close() 
