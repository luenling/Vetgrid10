from Bio import SeqIO
import sys, os, re
from collections import defaultdict
import argparse


def read_snp_file(infile,both=False):
    """
    reads a cmh or sync like file and creates a dict with chrom->set(bps) or chrom->bps->line
    """
    if (both):
        snp_dict=defaultdict(defaultdict)  #dictionary of chroms with positions and pVals of snps
    else:
        snp_dict=defaultdict(set)  #dictionary of chroms with positions and pVals of snps
    inf = open(infile,"r")
    #load dictionary
    for line in inf:
        if re.match("\#",line):
            continue
        line=line.rstrip()
        fields=line.split()
        if (len(fields) < 3): continue
        if both:
            snp_dict[fields[0]][fields[1]]="\t".join(fields[2:])
        else:
            snp_dict[fields[0]].add(int(fields[1]))        
    inf.close()
    return snp_dict


parser = argparse.ArgumentParser(description='reads snps from file and gets n seq above and below from fasta file') 

parser.add_argument("-s", dest="snpfile", help="file with snps", required=True)
parser.add_argument("-f", dest="fastafile", help="fastafile", default="/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa")
parser.add_argument("-n", dest="dist", help="distance", default=500)

args = parser.parse_args()
snpfile = vars(args)['snpfile']
fastafile = vars(args)['fastafile']
dist = int(vars(args)['dist'])

snp_dict=read_snp_file(snpfile)
#print snp_dict
# open fasta
fastafile="/Volumes/Temp/Lukas/reference_genome/ref_genome_combined.fa"
for seq in SeqIO.parse(open(fastafile,"ru"),"fasta"):
    #print seq.id 
    for pos in snp_dict[seq.id]:
        print ">"+seq.id+":"+str(pos)
        print seq[pos-dist-1:pos-1].seq.lower()+seq[pos-1]+seq[pos:pos+dist].seq.lower()

