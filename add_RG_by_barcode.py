######################### filtering a bam file by tags
import pysam 
import sys, re
# for fuzzy matching
import regex
import os 
import argparse
import collections
#########################################################   HELP   #########################################################################
parser = argparse.ArgumentParser(description="""Go through a barcoded bam file and add read groups by barcodes, prints number of reads processed and grouped on stdOUT.
ID: run_name+barcode, or only run_name for reads with no matching barcode
PL: ILLUMINA
SM: either barcode or list of SAMPLE_IDs:barcode or NONE if no barcode
LB : SM _ barcode or NONE
needs the python packages pysam and regex to be installed. It allows to give a maximal number of substitutions, although that wont work for the first position of the barcode. the default value for permitted substitutions is 1.
By default leaves reads without matching barcode in the file. To obtain a file only containing barcoded reads, use flag \"--clean\". 
Can also translate Illumina base quality encoding to sanger, though it does no checking, so use with care!
""")

parser.add_argument("--in", dest="infile", help="bam file with tagged reads", required=True)
parser.add_argument("--tags", dest="tags", help="comma separated list of pairs of sample ids and tags (eg: \"light_RI:TGACCAAT,light_RII:ACAGTGAT,dark_RI:GCCAATAT,dark_RII:CTTGTAAT\") ", required=True)
parser.add_argument("--subs", dest="subs", type=int, help="maximal number of substitutions in tag sequence (default=1)", default=1)
parser.add_argument("--base", dest="base", help="basename of run", required=True)
parser.add_argument("--out", dest="out", help="output filename", default=None)
parser.add_argument("--platform", dest="platform", help="platform (def: \"ILLUMINA\")", default="ILLUMINA")
parser.add_argument("--clean", dest="clean", action="store_true", help="remove reads without matching barcode (default: FALSE)", default=False)
parser.add_argument("--sanger", dest="sanger", action="store_true", help="change Illumina encoding to Sanger (no checking) (default: FALSE)", default=False)
parser.add_argument("--index", dest="idx", action="store_true", help="index resulting bam file (default: FALSE)", default=False)

args = parser.parse_args()

infile = vars(args)['infile']
idx = vars(args)['idx']
sanger = vars(args)['sanger']
base_name = vars(args)['base']
out = vars(args)['out']
platform = vars(args)['platform']
clean = vars(args)['clean']
subs = vars(args)['subs']
tags= vars(args)['tags'].split(",")
# check output filename
if out == None:
    out=base_name+"_RG.bam"
if os.path.exists(out):
    sys.stderr.write(out+" exists, wont overwrite, fuck off")   
    sys.exit(1)
# get ids
ids=list(tags)
for i,tag in enumerate(tags):
	if (re.search(":",tag)):
		(ids[i],tags[i])=tags[i].split(":")
# remove whitespace in tags, if there is any
tags=[ re.sub("\s","",x) for x in tags]
ids=[ re.sub("\s","",x) for x in ids]
## index BAM file if necessary
#if not os.path.exists(infile+".bai"):
#	print "indexing "+infile
#	os.system("samtools index "+infile)
samfile=pysam.Samfile(infile,"rb")
tag_count=collections.defaultdict(int)
tag_pats=collections.defaultdict()
tag_id=collections.defaultdict()
count=0
# create translation table from illumina to sanger:
trtbl=''.join([ chr(i-31) if (i-31>0) else chr(0) for i in range(256)])
#create RG header:
new_header = samfile.header.copy()
if not 'RG' in new_header.keys():
	new_header['RG']=[]
for i,tag in enumerate(tags):
	tag_pats[tag]=regex.compile("(?:"+tag+"){s<="+str(subs)+"}")
	tag_id[tag]=base_name+"_"+ids[i]+"_"+tag
	RG_entry = { 'ID': tag_id[tag] , 'PL': platform , 'SM': ids[i], 'LB': ids[i]+"_"+tag }
	new_header['RG'].append(RG_entry)

if not clean:
	RG_entry = { 'ID': base_name , 'PL': platform , 'SM': "NONE", 'LB': "NONE" }
	tags.append('N')
	new_header['RG'].append(RG_entry)
	tag_id['N']=base_name

rg_file=pysam.Samfile(out,"wb",header = new_header)

## write RGs in BAM file

for entry in samfile.fetch(until_eof=True):
	read_tag=entry.qname.split('#')[-1]
	count+=1
	# if count > 10000:
	# 	break
	for i,tag in enumerate(tags):
		if tag == "N":
			entry.tags += [('RG',tag_id[tag])]
			tag_count[tag] += 1		
			if sanger:
				entry.qual=entry.qual.translate(trtbl)
			rg_file.write(entry)
		elif regex.match(tag_pats[tag],read_tag):
			entry.tags += [('RG',tag_id[tag])]
			tag_count[tag] += 1		
			if sanger:
				entry.qual=entry.qual.translate(trtbl)
			rg_file.write(entry)
			break
print "reads read:\t{}".format(count)
tot_tag=0
for tag in tags:
	print "tag {}:\t{}".format(tag,tag_count[tag])
	tot_tag += tag_count[tag]
print "reads not taggged:\t{}".format(count-tot_tag)
samfile.close()
rg_file.close()
# index new file
if idx:
	os.system("samtools index "+base_name+"_RG.bam")
