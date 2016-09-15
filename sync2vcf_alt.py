import sys
import gzip
import re
import collections
import math
import argparse
#import vcf

#Author: Martin Kapun, Lukas Endler
#version: 1.2

#########################################################   HELP   #########################################################################
#print
usage="python %prog --sync input.sync --snps subset.snps --vcf output.vcf"
parser = argparse.ArgumentParser(description='''		
H E L P:
_________

This script takes a file with polymorphic SNPs (e.g. the output of the CMH test) in sync file format as input and identifies the reference and an alternative allele (in case of more than two alleles, the two most common ones will be used) across all populations in the file. Then it creates an output in the VCF v.4.1 output format, which can be used as the input for SNPeff, etc. Note that only the columns REF, ALT and FILTER in the VCF file will be filled. all other columns will be left blank (with a \".\" ). If you only want to convert a subset of the SNP input to vcf you can use the parameter --snps and provide an input file, which has to have at least two columns (sync, vcf), specifying the Chromosome and the positions. SNPs in this file need to be a subset of the --sync input. Per default, this option is set to FALSE and does not need to be set. All provided files can be gzipped or bgzipped (need ending \".(b)gz\")\n

See output example:

##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO 
3R	4704124	.	T	C	.	PASS	.
3R	8566126	.	A	C	.	.	PASS
3R	1178449	.	C	T	.	.	PASS
3R	15220143	.	C	T	.	PASS	.
3R	13777597	.	G	A	.	PASS	.

	''', formatter_class=argparse.RawTextHelpFormatter)

#########################################################   CODE   #########################################################################

parser.add_argument("-s", "--sync", dest="sync", help="sync file, can be gzipped or bgzipped, can have non colon seperated entries such as pvalues (are not considered) ", required=True)
parser.add_argument("-p", "--snps", dest="snps", help="define a subset of the SNPs which should be converted to VCF, can be gzipped, can be vcf or sync like tab delimited format",default=False)
args = parser.parse_args()

snplist=set()
snps= vars(args)['snps']
sync= vars(args)['sync']
if snps:
	if re.search("\.b?gz",snps):
		inf = gzip.open(snps,'rb')
	else:
		inf = open(snps,"r")
	for l in inf:
		if (l.rstrip("\n") =="" or re.match("\s*\#",l) ):
			continue
		a=l.split()
		b=a[0]+"_"+a[1]
		snplist.add(b)
	inf.close()


if re.search("\.b?gz",sync):
	syncf = gzip.open(sync,'rb')
else:
	syncf = open(sync,"r")

print "##fileformat=VCFv4.1"
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
dict0={}
dict1={}
for l in syncf:
	if not (l.rstrip("\n") =="" or re.match("\s*\#",l) ):
		a=l.split()
		b=a[0]+"_"+a[1]
		if not (snps) or (b in snplist):
			for i in range(0,len(l.split())-3):
				#print i+3
				bases=["A","T","C","G"]
				comp=a[i+3]
				if comp != "-" and re.search(":",comp):
					dict0[str(i+3)]=dict(zip(bases,map(float,comp.split(":")[:4])))
			#print dict0
			if a[2]!="A":
				dict1["A"]=sum(value.get("A", 0) for value in dict0.values())
			if a[2]!="T":
				dict1["T"]=sum(value.get("T", 0) for value in dict0.values())
			if a[2]!="C":
				dict1["C"]=sum(value.get("C", 0) for value in dict0.values())
			if a[2]!="G":
				dict1["G"]=sum(value.get("G", 0) for value in dict0.values())
			print a[0]+"\t"+a[1]+"\t.\t"+a[2]+"\t"+max(dict1, key = lambda x: dict1.get(x))+"\t.\tPASS\t."
			dict0={}
			dict1={}
syncf.close()



